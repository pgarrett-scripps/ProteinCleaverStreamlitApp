from collections import Counter

import requests
import streamlit as st
from peptacular.sequence import strip_modifications, calculate_sequence_length
from peptacular.digest import identify_cleavage_sites
from peptacular.mass import valid_mass_sequence
from peptacular.constants import AMINO_ACIDS
from peptacular.spans import calculate_span_coverage

from constants import *
from util import make_clickable, generate_peptide_df

# TODO: add color gradient to protein coverage to show the most covered regions

st.set_page_config(page_title="proteincleaver", page_icon=":knife:")

# CSS to inject contained in a string
hide_table_row_index = """
            <style>
            thead tr th:first-child {display:none}
            tbody th {display:none}
            </style>
            """
# Inject CSS with Markdown
st.markdown(hide_table_row_index, unsafe_allow_html=True)


@st.cache_data
def fetch_sequence_from_uniprot(accession_number):
    url = f"https://www.uniprot.org/uniprot/{accession_number}.fasta"
    response = requests.get(url)
    if response.status_code != 200:
        st.error(f"Error fetching sequence from UniProt: {response.status_code}")
        st.stop()
        return None
    return ''.join(response.text.split('\n')[1:])  # Remove the header line


with st.sidebar:
    st.title('Protein Cleaver 🔪')
    st.markdown("""
    **Protein Cleaver** is an app to compute protease-specific cleavage sites and peptides.
    """)
    st.markdown("""
    You can input your protein sequence in **FASTA format** (Don't worry if the sequence spans multiple lines), 
    or input the proteins accession number (e.g. P04406) / protein identifier (e.g. ALBU_HUMAN). 
    Get started and unveil the peptides!
    
    You can globally specify **static modifications** to be applied to all amino acids, or directly include them in the
    sequence. Modifications are specified by parenthesis, while terminal modifications use square brackets.
    
    Example: `[-13]MAS(1.2345)FRLFLLCLAGLVFVS[57.0]`
    """)

    protein_id = st.text_input(label='Protein accession number / identifier',
                               value='',
                               help='A protein accession number / identifier')

    raw_protein_sequence = None
    if protein_id:
        fetched_protein_sequence = fetch_sequence_from_uniprot(protein_id)
        if fetched_protein_sequence is not None:
            raw_protein_sequence = fetched_protein_sequence

    raw_protein_sequence = st.text_area(label="Protein sequence",
                                        value=raw_protein_sequence if raw_protein_sequence else DEFAULT_PROTEIN_SEQUENCE,
                                        help='An amino acid sequence to digest',
                                        max_chars=MAX_PROTEIN_INPUT_LENGTH,
                                        height=200)

    protein_sequence = raw_protein_sequence.replace(' ', '').replace('\n', '')
    stripped_protein_sequence = strip_modifications(protein_sequence)
    protein_length = calculate_sequence_length(protein_sequence)
    st.caption(f'Length: {protein_length}')

    if protein_length > MAX_PROTEIN_LEN:
        st.error(f'Protein sequence is too long. Please keep it under {MAX_PROTEIN_LEN} residues.')
        st.stop()

    # check if all residues are valid
    if not valid_mass_sequence(protein_sequence):
        st.error('Invalid amino acid sequence. Please check your input.')
        st.stop()

    c1, c2 = st.columns(2)
    proteases_selected = c1.multiselect(label="Proteases",
                                        options=list(VALID_PROTEASES.keys()),
                                        help='The proteases to use for digestion',
                                        default=DEFAULT_PROTEASES)
    custom_regex = c2.text_input(label='(Additional) Custom protease',
                                 value='',
                                 help='A custom regular expression to use for digestion. Will be used along with selected proteases')

    missed_cleavages = st.number_input(label='Max number of missed cleavages',
                                       min_value=MIN_MISSED_CLEAVAGES,
                                       max_value=MAX_MISSED_CLEAVAGES,
                                       value=DEFAULT_MISSED_CLEAVAGES,
                                       step=MISSED_CLEAVAGES_STEP,
                                       help='Number of missed cleavages to allow during digestion')

    enzyme_regexes = [VALID_PROTEASES[protease] for protease in proteases_selected]

    if custom_regex:
        enzyme_regexes.append(custom_regex)

    c1, c2 = st.columns(2)
    min_peptide_len = c1.number_input(label='Min peptide length',
                                      min_value=MIN_PEPTIDE_LEN,
                                      max_value=MAX_PEPTIDE_LEN,
                                      value=DEFAULT_MIN_PEPTIDE_LEN,
                                      step=PEPTIDE_LEN_STEP,
                                      help='Minimum peptide length (inclusive)')
    max_peptide_len = c2.number_input(label='Max peptide length',
                                      min_value=MIN_PEPTIDE_LEN,
                                      max_value=MAX_PEPTIDE_LEN,
                                      value=DEFAULT_MAX_PEPTIDE_LEN,
                                      step=PEPTIDE_LEN_STEP,
                                      help='Maximum peptide length (inclusive)')

    if min_peptide_len > max_peptide_len:
        st.error('Min length must be less than max length')
        st.stop()

    c3, c4 = st.columns(2)
    min_mass = c3.number_input(label='Min peptide mass',
                               min_value=MIN_PEPTIDE_MASS,
                               max_value=MAX_PEPTIDE_MASS,
                               value=DEFAULT_MIN_PEPTIDE_MASS,
                               step=PEPTIDE_MASS_STEP,
                               help='Minimum peptide mass (inclusive)')
    max_mass = c4.number_input(label='Max peptide mass',
                               min_value=MIN_PEPTIDE_MASS,
                               max_value=MAX_PEPTIDE_MASS,
                               value=DEFAULT_MAX_PEPTIDE_MASS,
                               step=PEPTIDE_MASS_STEP,
                               help='Maximum peptide mass (inclusive)')

    if min_mass > max_mass:
        st.error('Min mass must be less than max mass')
        st.stop()

    c1, c2 = st.columns(2)
    mass_type = c1.radio(label='Mass type',
                         options=MASS_TYPE_OPTIONS,
                         index=DEFAULT_MASS_TYPE_INDEX,
                         help='Mass type to use for calculations')
    is_mono = mass_type == 'monoisotopic'
    semi_enzymatic = c2.checkbox(label='Semi enzymatic?',
                                 help='Allow semi enzymatic peptides?')

    st.subheader('Static Modifications')

    # a selection for the user to specify the number of rows
    c1, c2 = st.columns(2)
    num_rows = st.number_input(label='Number of unique modifications',
                               min_value=MIN_STATIC_MODS,
                               max_value=MAX_STATIC_MODS,
                               value=DEFAULT_STATIC_MODS,
                               step=STATIC_MODS_STEP,
                               help='Add another modification row')

    # columns to lay out the inputs
    grid = st.columns([3, 2])


    def add_row(row):
        with grid[0]:
            st.multiselect(label='Amino acids',
                           key=f'residues{row}',
                           options=list(AMINO_ACIDS),
                           help='Select amino acids for which to apply the static modification',
                           default=['C'] if row == 0 else [])
        with grid[1]:
            st.number_input(label='Modification Mass (Da)',
                            step=0.00001, key=f'mass{row}',
                            help='The mass of the modification (in daltons)',
                            value=57.02146 if row == 0 else 0.0,
                            format='%.5f')


    # Loop to create rows of input widgets
    for r in range(num_rows):
        add_row(r)

    mods = {}
    for r in range(num_rows):
        for residue in st.session_state[f'residues{r}']:
            mods[residue] = "{:.5f}".format(st.session_state[f'mass{r}'])

sites = set()
for enzyme_regex in enzyme_regexes:
    sites.update(identify_cleavage_sites(protein_sequence, enzyme_regex))
sites = sorted(list(sites))

df = generate_peptide_df(protein_sequence, sites, missed_cleavages, min_peptide_len, max_peptide_len,
                         semi_enzymatic, mods, min_mass, max_mass, is_mono)

sequence_with_sites = stripped_protein_sequence
for site in sites[::-1]:
    sequence_with_sites = sequence_with_sites[:site] + '%' + sequence_with_sites[site:]

# color all | in red
sequence_with_sites = sequence_with_sites.replace('%', '<span style="color:red">%</span>')

spans = [(s, e, mc) for s, e, mc in df[['Start', 'End', 'MC']].values]
protein_cov_arr = calculate_span_coverage(spans, protein_length)

# color all covered amino acids in red
protein_coverage = ''
for i, aa in enumerate(stripped_protein_sequence):
    if protein_cov_arr[i] == 1:
        protein_coverage += '<span style="color:red">' + aa + '</span>'
    else:
        protein_coverage += aa

# calculate protein coverage at different MC
protein_cov_at_mcs = []
mcs = [mc for mc in range(0, missed_cleavages + 1)]
for mc in mcs:
    df_mc = df[df['MC'] <= mc]
    spans = [(s, e, mc) for s, e, mc in df_mc[['Start', 'End', 'MC']].values]
    cov = calculate_span_coverage(spans, protein_length)
    protein_cov_at_mcs.append(sum(cov) / len(cov) * 100)

# calculate protein coverage at different peptide lengths
protein_cov_at_lens = []
lens = [l for l in range(min_peptide_len, max_peptide_len + 1)]
for l in lens:
    df_len = df[df['AACount'] <= l]
    spans = [(s, e, mc) for s, e, mc in df_len[['Start', 'End', 'MC']].values]
    cov = calculate_span_coverage(spans, protein_length)
    protein_cov_at_lens.append(sum(cov) / len(cov) * 100)

# calculate protein coverage at different peptide Mass
protein_cov_at_mass = []
masses = [m for m in range(int(min_mass), int(max_mass) + 1, 100)]
for m in masses:
    df_mass = df[df['Mass'] <= m]
    spans = [(s, e, mc) for s, e, mc in df_mass[['Start', 'End', 'MC']].values]
    cov = calculate_span_coverage(spans, protein_length)
    protein_cov_at_mass.append(sum(cov) / len(cov) * 100)

protein_cov_perc = round(sum(protein_cov_arr) / len(protein_cov_arr) * 100, 2)


df.drop(columns=['StrippedPeptide'], inplace=True)
df.sort_values(by=['MC'], inplace=True)
df.drop_duplicates(subset=['PeptideSequence', 'IsSemi'], inplace=True)
df.sort_values(by=['Start', 'AACount'], inplace=True)
df = df[['PeptideSequence', 'Start', 'End', 'AACount', 'MC', 'IsSemi', 'Mass']]
df_download = df.to_csv(index=False)
df_clickable = df.copy(deep=True)
df_clickable['PeptideSequence'] = [make_clickable(peptide, mass_type) for peptide in df_clickable['PeptideSequence']]


t1, t2, t3, t4, t5 = st.tabs(['Digestion Metrics', 'Cleavage & Coverage', 'Motif Analysis', 'Wiki', 'Help'])

with t1:
    st.subheader('Digestion Overview')
    c1, c2, c3 = st.columns(3)
    c1.metric('Cleavage Sites', len(sites))
    c2.metric('Total Peptides', len(df_clickable))
    c3.metric('Unique Peptides', len(df_clickable.drop_duplicates(subset=['PeptideSequence'])))
    c1, c2, c3 = st.columns(3)
    c1.metric('Semi Peptides', len(df_clickable[df_clickable['IsSemi'] == True]))
    c2.metric('Enzymatic Peptides', len(df_clickable[df_clickable['IsSemi'] == False]))
    c3.metric('Protein Coverage', f'{protein_cov_perc}%')

    st.subheader('Peptide Results')
    df_html = df_clickable.to_html(escape=False)

    st.caption('Click on a sequence to see the fragment ions!')
    st.write(df_html, unsafe_allow_html=True, use_container_width=True)
    st.write(' ')
    st.download_button('Download CSV', df_download, 'digestion.csv', 'text/csv', use_container_width=True)


with t2:
    st.subheader('Cleavage Sites', help='Red % are cleavage sites')
    st.markdown(sites)
    st.markdown(sequence_with_sites, unsafe_allow_html=True)

    st.subheader('Protein Coverage', help='Red amino acids are covered by peptides')
    st.markdown(protein_coverage, unsafe_allow_html=True)

    st.subheader('Coverage Dynamics: Missed Cleavages')
    st.line_chart(data={'Missed Cleavages': mcs, 'Protein Coverage (%)': protein_cov_at_mcs},
                  x='Missed Cleavages', y='Protein Coverage (%)')

    st.caption('Coverage Dynamics: Peptide Lengths')
    st.line_chart(data={'Peptide Length': lens, 'Protein Coverage (%)': protein_cov_at_lens},
                  x='Peptide Length', y='Protein Coverage (%)')

    st.caption('"Coverage Dynamics: Peptide Masses')
    st.line_chart(data={'Peptide Mass': masses, 'Protein Coverage (%)': protein_cov_at_mass},
                  x='Peptide Mass', y='Protein Coverage (%)')

with t3:
    st.subheader('Peptide Pattern Identification')
    st.write('This page identifies and tallies occurrences of your specified pattern within the digested peptides.'
             ' It provides a count of how many times each pattern appears in every peptide.')

    site_regex = st.text_input('Enter Pattern Regex', '(K)')

    if site_regex:
        sites = identify_cleavage_sites(stripped_protein_sequence, site_regex)

        site_counts = []
        for row in df[['Start', 'End']].values:
            site_counts.append(sum([1 for site in sites if row[0] <= site < row[1]]))
        df['PatternMatch'] = site_counts

        st.subheader('Peptides')
        st.dataframe(df)

        counter = Counter(site_counts)

        st.subheader('Frequency Analysis', help='Number of peptides with N number of motif matches')
        for k, v in sorted(counter.items()):
            st.metric(f'Peptides with {k} motif matches', v)

        st.subheader('Coverage Analysis', help='Coverage of protein based on peptides with N number of motif matches')
        for k, v in sorted(counter.items()):
            df_tmp = df[df['PatternMatch'] == k]
            tmp_spans = [(s, e, mc) for s, e, mc in df_tmp[['Start', 'End', 'MC']].values]
            cov = calculate_span_coverage(tmp_spans, protein_length)
            st.metric(f'Protein Coverage with {k} motif matches', f'{round(sum(cov) / len(cov) * 100, 2)}%')
with t4:
    st.markdown(PROTEASE_WIKI)

with t5:
    st.markdown(HELP)

    st.subheader('Proteases:')
    st.write(VALID_PROTEASES)
