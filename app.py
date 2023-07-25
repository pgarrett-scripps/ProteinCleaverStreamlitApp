import numpy as np
import streamlit as st
from peptacular.protein import calculate_protein_coverage
from peptacular.sequence import identify_cleavage_sites, add_static_mods, strip_modifications
from peptacular.mass import calculate_mass
from peptacular.constants import AMINO_ACIDS

from constants import MAX_PROTEIN_LEN, VALID_PROTEASES, EXAMPLE_PROTEIN, MAX_PEPTIDE_MASS, MAX_PEPTIDE_LEN, \
    MAX_MISSED_CLEAVAGES, PROTEASE_WIKI, HELP
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

with st.sidebar:
    st.title('Protein Cleaver 🔪')
    st.markdown("""
    **Protein Cleaver** is an app to compute protease-specific cleavage sites and peptides.
    """)
    st.markdown("""
    You can input your protein sequence in the well-known **FASTA format**. Don't worry if the sequence spans multiple 
    lines - Protein Cleaver will handle it seamlessly. Get started and unveil the peptides!
    """)

    sequence = st.text_area("Protein sequence", value=EXAMPLE_PROTEIN, help='An amino acid sequence to digest',
                            max_chars=MAX_PROTEIN_LEN)
    sequence = sequence.replace(' ', '').replace('\n', '')
    st.caption(f'Length: {len(sequence)}')

    if len(sequence) > MAX_PROTEIN_LEN:
        st.error(f'Protein sequence is too long. Please keep it under {MAX_PROTEIN_LEN} characters.')
        st.stop()

    # check if all residues are valid
    for aa in sequence:
        if aa not in AMINO_ACIDS:
            st.error(f'Invalid amino acid: {aa}')
            st.stop()

    c1, c2 = st.columns(2)
    proteases_selected = c1.multiselect("Select Proteases", options=list(VALID_PROTEASES.keys()),
                                        help='The proteases to use for digestion', default=['trypsin'])
    enzyme_regexes = [VALID_PROTEASES[protease] for protease in proteases_selected]
    missed_cleavages = c2.number_input('Max number of missed cleavages', min_value=0,
                                       max_value=MAX_MISSED_CLEAVAGES, value=1, step=1,
                                       help='Maximum number of missed cleavages to allow')

    c1, c2 = st.columns(2)
    min_len = c1.number_input('Min peptide length', min_value=1, max_value=MAX_PEPTIDE_LEN, value=7, step=1,
                              help='Minimum peptide length (inclusive)')
    max_len = c2.number_input('Max peptide length', min_value=1, max_value=MAX_PEPTIDE_LEN, value=20, step=1,
                              help='Maximum peptide length (inclusive)')

    if min_len > max_len:
        st.error('Min length must be less than max length')
        st.stop()

    c3, c4 = st.columns(2)
    min_mass = c3.number_input('Min peptide mass', min_value=0, max_value=MAX_PEPTIDE_MASS, value=200, step=100,
                               help='Minimum peptide mass (inclusive)')
    max_mass = c4.number_input('Max peptide mass', min_value=0, max_value=MAX_PEPTIDE_MASS, value=3000, step=100,
                               help='Maximum peptide mass (inclusive)')

    if min_mass > max_mass:
        st.error('Min mass must be less than max mass')
        st.stop()

    c1, c2 = st.columns(2)
    mass_type = c1.radio('Mass type', options=['monoisotopic', 'average'], index=0,
                         help='Mass type to use for calculations')
    is_mono = mass_type == 'monoisotopic'
    semi_enzymatic = c2.checkbox('Semi enzymatic?', help='Allow semi enzymatic peptides?')

    st.subheader('Static Modifications')

    # a selection for the user to specify the number of rows
    c1, c2 = st.columns(2)
    num_rows = st.number_input('Number of unique modifications', min_value=0, value=1, step=1, help='Add another modification row')

    # columns to lay out the inputs
    grid = st.columns([3, 2])

    def add_row(row):
        if row  == 0:
            with grid[0]:
                st.multiselect('Amino acids', key=f'residues{row}', options=list(AMINO_ACIDS),
                               help='Select amino acids for which to apply the static modification', default=['C'])
            with grid[1]:
                st.number_input('Modification Mass (Da)', step=0.00001, key=f'mass{row}',
                                help='The mass of the modification (in daltons)', value=57.02146, format='%.5f')
        else:
            with grid[0]:
                st.multiselect('Amino acids', key=f'residues{row}', options=list(AMINO_ACIDS),
                               help='Select amino acids for which to apply the static modification')
            with grid[1]:
                st.number_input('Modification Mass (Da)', step=0.00001, key=f'mass{row}',
                                help='The mass of the modification (in daltons)', format='%.5f')


    # Loop to create rows of input widgets
    for r in range(num_rows):
        add_row(r)

    mods = {}
    for r in range(num_rows):
        for residue in st.session_state[f'residues{r}']:
            mods[residue] = "{:.5f}".format(st.session_state[f'mass{r}'])

sites = set()
for enzyme_regex in enzyme_regexes:
    sites.update(identify_cleavage_sites(sequence, enzyme_regex))
sites = sorted(list(sites))

df = generate_peptide_df(sequence, sites, missed_cleavages, min_len, max_len, semi_enzymatic)

sequence_with_sites = sequence
for site in sites[::-1]:
    sequence_with_sites = sequence_with_sites[:site] + '%' + sequence_with_sites[site:]

# color all | in red
sequence_with_sites = sequence_with_sites.replace('%', '<span style="color:red">%</span>')

# make bool
df['IsSemi'] = df['IsSemi'].apply(lambda x: bool(x))
df.drop_duplicates(inplace=True)
df['PeptideSequence'] = df['PeptideSequence'].apply(lambda x: add_static_mods(x, mods))

df['Mass'] = [round(calculate_mass(sequence, monoisotopic=is_mono), 5) for sequence in df['PeptideSequence']]
df = df[(df['Mass'] >= min_mass) & (df['Mass'] <= max_mass)]

protein_cov_arr = calculate_protein_coverage(sequence, [strip_modifications(seq) for seq in df['PeptideSequence']])

# given a array of 0 and 1 which represent amino acids that are either not covered or covered higjlight the protein sequence
protein_coverage = ''
for i, aa in enumerate(sequence):
    if protein_cov_arr[i] == 1:
        protein_coverage += '<span style="color:red">' + aa + '</span>'
    else:
        protein_coverage += aa

protein_cov_perc = round(sum(protein_cov_arr) / len(protein_cov_arr) * 100, 2)

# keep first fuplicate
df.sort_values(by=['MC'], inplace=True)
df.drop_duplicates(subset=['PeptideSequence', 'IsSemi'], inplace=True)
df.sort_values(by=['Start', 'AACount'], inplace=True)
df = df[['PeptideSequence', 'Start', 'End', 'AACount', 'MC', 'IsSemi', 'Mass']]

df_download = df.to_csv(index=False)

# link is the column with hyperlinks
df['PeptideSequence'] = [make_clickable(peptide, mass_type) for peptide in df['PeptideSequence']]


st.warning('This is a work in progress. Please report any issues or suggestions to pgarrett@scripps.edu.')

t1, t2, t3, t4 = st.tabs(['Results', 'Protease Regexes', 'Wiki', 'Help'])


with t1:

    c1, c2, c3 = st.columns(3)
    c1.metric('Cleavage Sites', len(sites))
    c2.metric('Total Peptides', len(df))
    c3.metric('Unique Peptides', len(df.drop_duplicates(subset=['PeptideSequence'])))
    c1, c2, c3 = st.columns(3)
    c1.metric('Semi Peptides', len(df[df['IsSemi'] == True]))
    c2.metric('Enzymatic Peptides', len(df[df['IsSemi'] == False]))
    c3.metric('Protein Coverage', f'{protein_cov_perc}%')

    with st.expander('Protein Coverage'):
        st.markdown(protein_coverage, unsafe_allow_html=True)
        st.caption('Red amino acids are covered by peptides, black are not')

    with st.expander('Cleavage Sites'):
        st.markdown(sites)
        st.markdown(sequence_with_sites, unsafe_allow_html=True)
        st.caption('Cleavage sites are marked with a red %')

    st.subheader('Peptides')
    df = df.to_html(escape=False)

    st.caption('Click on a sequence to see the fragment ions!')
    st.write(df, unsafe_allow_html=True, use_container_width=True)
    st.write(' ')
    st.download_button('Download CSV', df_download, 'digestion.csv', 'text/csv', use_container_width=True)

with t2:
    st.write(VALID_PROTEASES)

with t3:
    st.markdown(PROTEASE_WIKI)

with t4:
    st.markdown(HELP)