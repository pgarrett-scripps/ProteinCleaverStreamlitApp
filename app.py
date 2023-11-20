import pickle
from collections import Counter

import streamlit as st
import numpy as np
import requests
from peptacular.sequence import strip_modifications, calculate_sequence_length
from peptacular.digest import identify_cleavage_sites
from peptacular.mass import valid_mass_sequence
from peptacular.constants import AMINO_ACIDS
from peptacular.spans import calculate_span_coverage
import matplotlib as mpl

from constants import *
from wiki import *
from util import make_clickable, generate_peptide_df, bin_aa_counts, coverage_string, create_colorbar

# TODO: add color gradient to protein coverage to show the most covered regions

st.set_page_config(page_title="proteincleaver", page_icon=":knife:", layout="wide")

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

    protein_id = st.text_input(label='Protein accession number / identifier',
                               value='',
                               help='A protein accession number / identifier')

    raw_sequence = None
    if protein_id:
        fetched_protein_sequence = fetch_sequence_from_uniprot(protein_id)
        if fetched_protein_sequence is not None:
            raw_sequence = fetched_protein_sequence

    raw_sequence = st.text_area(label="Protein sequence",
                                value=raw_sequence if raw_sequence else DEFAULT_PROTEIN_SEQUENCE,
                                help='An amino acid sequence to digest',
                                max_chars=MAX_PROTEIN_INPUT_LENGTH,
                                height=200)

    protein_sequence = raw_sequence.replace(' ', '').replace('\n', '')
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
                                 help='A custom regular expression to use for digestion. Will be used along with '
                                      'selected proteases')

    c1, c2 = st.columns(2)
    missed_cleavages = c1.number_input(label='Max missed cleavages',
                                       min_value=MIN_MISSED_CLEAVAGES,
                                       max_value=MAX_MISSED_CLEAVAGES,
                                       value=DEFAULT_MISSED_CLEAVAGES,
                                       step=MISSED_CLEAVAGES_STEP,
                                       help='Number of missed cleavages to allow during digestion')

    mass_type = c2.selectbox(label='Mass type', options=MASS_TYPE_OPTIONS, index=DEFAULT_MASS_TYPE_INDEX,
                             help='Mass type to use for calculations')
    is_mono = mass_type == 'monoisotopic'

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

    semi_enzymatic = c1.checkbox(label='Semi enzymatic?',
                                 help='Allow semi enzymatic peptides?')
    infer_charge = c2.checkbox(label='Infer charge', value=True, help='Infer charge of peptides based on L and K count')

    min_charge, max_charge, min_mz, max_mz = None, None, None, None
    if infer_charge:
        c1, c2 = st.columns(2)
        min_charge = c1.number_input(label='Min charge',
                                     min_value=MIN_CHARGE,
                                     max_value=MAX_CHARGE,
                                     value=DEFAULT_MIN_CHARGE,
                                     step=CHARGE_STEP,
                                     help='Minimum peptide charge (inclusive)')

        max_charge = c2.number_input(label='Max charge',
                                     min_value=MIN_CHARGE,
                                     max_value=MAX_CHARGE,
                                     value=DEFAULT_MAX_CHARGE,
                                     step=CHARGE_STEP,
                                     help='Maximum peptide charge (inclusive)')

        c1, c2 = st.columns(2)
        min_mz = c1.number_input(label='Min m/z',
                                 min_value=MIN_MZ,
                                 max_value=MAX_MZ,
                                 value=DEFAULT_MIN_MZ,
                                 step=MZ_STEP,
                                 help='Minimum peptide m/z (inclusive)')

        max_mz = c2.number_input(label='Max m/z',
                                 min_value=MIN_MZ,
                                 max_value=MAX_MZ,
                                 value=DEFAULT_MAX_MZ,
                                 step=MZ_STEP,
                                 help='Maximum peptide m/z (inclusive)')

        if min_charge > max_charge:
            st.error('Min charge must be less than max charge')
            st.stop()

    remove_non_proteotypic = st.checkbox(label='Remove non-proteotypic peptides', value=False,
                                         help='Remove peptides that are not proteotypic')

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


    def add_row(r):
        with grid[0]:
            st.multiselect(label='Amino acids',
                           key=f'residues{r}',
                           options=list(AMINO_ACIDS),
                           help='Select amino acids for which to apply the static modification',
                           default=['C'] if r == 0 else [])
        with grid[1]:
            st.number_input(label='Modification Mass (Da)',
                            step=0.00001, key=f'mass{r}',
                            help='The mass of the modification (in daltons)',
                            value=57.02146 if r == 0 else 0.0,
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
                         semi_enzymatic, mods, min_mass, max_mass, is_mono, infer_charge, min_charge,
                         max_charge, min_mz, max_mz)

rt_model = pickle.load(open("rt_model.pkl", "rb"))
df['RT'] = rt_model.predict(np.array([bin_aa_counts(strip_modifications(seq)) for seq in df['Sequence']]))
df['RT'] = df['RT'].round(3)

im_model = pickle.load(open("im_model.pkl", "rb"))
df['IM'] = im_model.predict(
    np.array([bin_aa_counts(strip_modifications(seq), c) for seq, c in df[['Sequence', 'Charge']].values]))
df['IM'] = df['IM'].round(3)

proteotypic_model = pickle.load(open("proteotypic_model.pkl", "rb"))
df['Proteotypic'] = proteotypic_model.predict(
    np.array([bin_aa_counts(strip_modifications(seq)) for seq in df['Sequence']])).astype(bool)

if remove_non_proteotypic:
    df = df[df['Proteotypic']]

# Start the HTML string for the site indexes
site_indexes_html = '<span style="font-family: Courier New, monospace; font-size: 16px;">'
for index in sites:
    site_indexes_html += f'<span style="background-color:#f0f0f0; font-weight:900; color:red; padding:2px; ' \
                         f'margin:1px; border:1px solid #ffcc00; border-radius:3px;">{index}</span>'
    site_indexes_html += ' '
site_indexes_html += '</span>'

sequence_with_sites = '<span style="font-family: Courier New, monospace; font-size: 16px;">'
for i, aa in enumerate(stripped_protein_sequence):
    # Add the amino acid with its original index
    sequence_with_sites += f'<span title="Index: {i + 1}" style="background-color:#f0f0f0; font-weight:900; ' \
                           f'color:#333; padding:2px; margin:1px; border:1px solid #cccccc; ' \
                           f'border-radius:3px;">{aa}</span>'

    # Check if the next position is a cleavage site and insert '%' character
    if i + 1 in sites:
        # Highlight '%' character in blue
        sequence_with_sites += f'<span style="background-color:#e0e0ff; font-weight:900; color:red; font-weight:bold;' \
                               f' padding:2px; margin:1px; border:1px solid #a0a0ff; border-radius:3px;">%</span>'
sequence_with_sites += '</span>'

# Create a colormap
cmap = mpl.colormaps.get_cmap(CMAP)

spans = [(s, e, mc) for s, e, mc in df[['Start', 'End', 'MC']].values]
protein_cov_arr = calculate_span_coverage(spans, protein_length, accumulate=True)
protein_coverage = coverage_string(protein_cov_arr, stripped_protein_sequence, cmap)

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
    df_len = df[df['Len'] <= l]
    spans = [(s, e, mc) for s, e, mc in df_len[['Start', 'End', 'MC']].values]
    cov = calculate_span_coverage(spans, protein_length)
    protein_cov_at_lens.append(sum(cov) / len(cov) * 100)

# calculate protein coverage at different peptide Mass
protein_cov_at_mass = []
masses = [m for m in range(int(min_mass), int(max_mass) + 1, 100)]
for m in masses:
    df_mass = df[df['NeutralMass'] <= m]
    spans = [(s, e, mc) for s, e, mc in df_mass[['Start', 'End', 'MC']].values]
    cov = calculate_span_coverage(spans, protein_length)
    protein_cov_at_mass.append(sum(cov) / len(cov) * 100)

df.drop(columns=['StrippedPeptide'], inplace=True)
df.sort_values(by=['MC'], inplace=True)
df.drop_duplicates(subset=['Start', 'Sequence', 'Semi'], inplace=True)
df.sort_values(by=['Start', 'Len'], inplace=True)

t1, t2, t3, t4, t5, t6 = st.tabs(['Digestion Metrics', 'Cleavage & Coverage', 'Motif Analysis', 'Wiki', 'Help',
                                  'Machine Learning'])

with t1:
    c1, c2, c3, c4 = st.columns(4)
    c1.metric('Total Peptides', len(df))
    c2.metric('Semi Peptides', len(df[df['Semi']]))
    c3.metric('Enzymatic Peptides', len(df[~df['Semi']]))
    c4.metric('Unique Peptides', len(df['Sequence'].unique()))

    st.subheader('Peptides')
    clickable = st.checkbox('Peptide Fragmenter Links', value=True)

    if clickable:
        df_clickable = df.copy(deep=True)
        df_clickable['Sequence'] = [make_clickable(peptide, mass_type) for peptide in
                                    df_clickable['Sequence']]
        st.caption('Click on a sequence to see the fragment ions!')
        st.write(df_clickable.to_html(escape=False), unsafe_allow_html=True, use_container_width=True)
    else:
        st.dataframe(df, use_container_width=True)

with t2:
    c1, c2 = st.columns(2)
    c1.metric('Cleavage Sites', len(sites))

    protein_cov_arr_bin = calculate_span_coverage(spans, protein_length, accumulate=False)
    protein_cov_perc = round(sum(protein_cov_arr_bin) / len(protein_cov_arr_bin) * 100, 2)
    c2.metric('Protein Coverage', f'{protein_cov_perc}%')

    st.subheader('Site Indexes', help='Red % are cleavage sites')
    st.markdown(site_indexes_html, unsafe_allow_html=True)
    st.write("")
    st.subheader('Sites', help='Red amino acids are cleavage sites')
    st.markdown(sequence_with_sites, unsafe_allow_html=True)

    st.subheader('Sequence Coverage', help='Red amino acids are covered by peptides')
    st.markdown(protein_coverage, unsafe_allow_html=True)

    # Example usage in a Streamlit app
    f = create_colorbar(max(protein_cov_arr), cmap)
    st.pyplot(f)

    st.subheader('Coverage vs Missed Cleavages')
    st.line_chart(data={'Missed Cleavages': mcs, 'Protein Coverage (%)': protein_cov_at_mcs},
                  x='Missed Cleavages', y='Protein Coverage (%)')

    st.caption('Coverage vs Peptide Lengths')
    st.line_chart(data={'Peptide Length': lens, 'Protein Coverage (%)': protein_cov_at_lens},
                  x='Peptide Length', y='Protein Coverage (%)')

    st.caption('"Coverage vs Peptide Masses')
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

        cov_site_mat = [0] * len(stripped_protein_sequence)
        for row in df[['Start', 'End', 'PatternMatch']].values:

            if row[2] == 0:
                continue

            for i in range(row[0], row[1] + 1):
                if i in sites:
                    if row[2] == 1:
                        cov_site_mat[i - 1] = 1
                    else:
                        if cov_site_mat[i - 1] == 0:
                            cov_site_mat[i - 1] = row[2]
                        else:
                            cov_site_mat[i - 1] = min(row[2], cov_site_mat[i - 1])

        clickable2 = st.checkbox('Peptide Fragmenter Links', value=False, key=1)
        st.subheader('Peptides')
        if clickable2:
            df_clickable = df.copy(deep=True)
            df_clickable['Sequence'] = [make_clickable(peptide, mass_type) for peptide in df_clickable['Sequence']]
            st.caption('Click on a sequence to see the fragment ions!')
            st.write(df_clickable.to_html(escape=False), unsafe_allow_html=True, use_container_width=True)
        else:
            st.dataframe(df, use_container_width=True)

        counter = Counter(site_counts)

        st.subheader('Coverage Analysis', help='Coverage of protein based on peptides with N number of motif matches')
        sites = [i - 1 for i in sites if i != 0]
        protein_coverage = coverage_string(cov_site_mat, stripped_protein_sequence, cmap, sites)
        st.markdown(protein_coverage, unsafe_allow_html=True)
        f = create_colorbar(max(cov_site_mat), cmap, label='Min Number of Motif Matches')
        st.pyplot(f)
        with st.expander('Show Coverage Table'):
            for i, (k, v) in enumerate(sorted(counter.items())):
                df_tmp = df[df['PatternMatch'] == k]
                tmp_spans = [(s, e, mc) for s, e, mc in df_tmp[['Start', 'End', 'MC']].values]
                cov = calculate_span_coverage(tmp_spans, protein_length, accumulate=True)

                cov_bin = calculate_span_coverage(tmp_spans, protein_length, accumulate=False)

                c1, c2 = st.columns(2)
                c1.metric(f'Protein Coverage with {k} motif matches', f'{round(sum(cov_bin) / len(cov_bin) * 100, 2)}%')
                c2.metric(f'Peptides with {k} motif matches', v)

                protein_coverage = coverage_string(cov, stripped_protein_sequence, cmap)
                st.markdown(protein_coverage, unsafe_allow_html=True)

                f = create_colorbar(max(cov), cmap)
                st.pyplot(f)

                if i >= 5:
                    st.warning('Warning: High motif match counts may result in long runtimes. Stopping...')
                    break

with t4:
    st.markdown(PROTEASE_WIKI)

with t5:
    st.markdown(HELP)

    st.subheader('Proteases:')
    st.write(VALID_PROTEASES)

with t6:
    st.subheader('Download Models')


    def get_model_file_as_byte_stream(path):
        with open(path, 'rb') as file:
            byte_stream = file.read()
        return byte_stream


    # download models
    c1, c2, c3 = st.columns(3)
    c1.download_button(
        label='Download RT Model',
        data=get_model_file_as_byte_stream('rt_model.pkl'),
        file_name="rt_model.pkl",
        mime='application/octet-stream'
    )
    c2.download_button(
        label='Download IM Model',
        data=get_model_file_as_byte_stream('im_model.pkl'),
        file_name="im_model.pkl",
        mime='application/octet-stream'
    )
    c3.download_button(
        label='Download Proteotypic Model',
        data=get_model_file_as_byte_stream('proteotypic_model.pkl'),
        file_name="proteotypic_model.pkl",
        mime='application/octet-stream'
    )

    st.write('---')

    st.subheader('Example Usage')

    st.code(MODEL_CODE)

    st.write('---')

    st.subheader('About')

    with st.expander('Proteotypic Peptide Prediction Model', expanded=True):
        st.markdown(PROTEOTYPIC_MODEL_HELP)

    with st.expander('RT Prediction Model', expanded=False):
        st.markdown(RT_MODEL_HELP)

    with st.expander('IM Prediction Model', expanded=False):
        st.markdown(IM_MODEL_HELP)
