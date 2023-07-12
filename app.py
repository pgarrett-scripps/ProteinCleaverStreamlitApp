import copy

import pandas as pd
import streamlit as st
from peptacular.sequence import digest_sequence, identify_cleavage_sites, calculate_mass, add_static_mods
from peptacular.constants import PROTEASES, AMINO_ACIDS

st.set_page_config(page_title="proteindigestor", page_icon=":knife:")

# CSS to inject contained in a string
hide_table_row_index = """
            <style>
            thead tr th:first-child {display:none}
            tbody th {display:none}
            </style>
            """

# Inject CSS with Markdown
st.markdown(hide_table_row_index, unsafe_allow_html=True)
VALID_PROTEASES = copy.deepcopy(PROTEASES)
VALID_PROTEASES.pop('non-specific', None)
with st.expander('Protease Regexes'):
    st.write(VALID_PROTEASES)

example_protein = "HRNGEMYACEQEHDKEPHMKIMPHGSGGFFPLVQFGRHFGQLKNLKRPAVHVDTEVLYWCNTRCEFLMWAFDCRIDPRDWGMDHMHCRESRCYASFRG" \
                  "TRGFDNLFYPAKHLEMHGTMISIMQWFQANGDKTLHSTYKFMSPCSGEKRMYQSWKWGEKPRCYSTQHVYCAVDKRMSVWSKCFSQGKALGTKESLNN" \
                  "VDDHHDLKQCVMISSWSTPYCKIPNCAAEQWMETMTMPWDWPPMFIKIVIASDRCVVLHPQLGLHAHGMTRWATTVRKGKIGFYDPGPNMCYWQQQWL" \
                  "FTVAGS"

LINK = 'https://peptidefragmenter.streamlit.app/'

with st.sidebar:
    st.title('Protein Digestion')
    st.markdown("""Digests a protein into peptides according to a protease of your choice.""")

    sequence = st.text_input("sequence to be digested", value=example_protein)
    sequence = sequence.replace(' ', '').replace('\n', '')

    #check if all residues are valid
    for aa in sequence:
        if aa not in AMINO_ACIDS:
            st.error(f'Invalid amino acid: {aa}')
            st.stop()

    protease_selected = st.selectbox("Select the Protease", options=list(VALID_PROTEASES.keys()))
    enzyme_regexes = VALID_PROTEASES[protease_selected]

    c1, c2 = st.columns(2)
    min_len = c1.number_input('Min peptide length', min_value=1, value=7, step=1)
    max_len = c2.number_input('Max peptide length', min_value=1, value=25, step=1)

    if min_len > max_len:
        st.error('Min length must be less than max length')
        st.stop()

    c3, c4 = st.columns(2)
    min_mass = c3.number_input('Min peptide mass', min_value=0, value=200, step=100)
    max_mass = c4.number_input('Max peptide mass', min_value=0, value=3000, step=100)

    if min_mass > max_mass:
        st.error('Min mass must be less than max mass')
        st.stop()

    missed_cleavages = st.number_input('Max number of missed cleavages', min_value=0, value=1, step=1)
    semi_enzymatic = st.checkbox('Semi Enzymatic?')

    st.subheader('Static Modifications')

    # a selection for the user to specify the number of rows
    c1, c2 = st.columns(2)
    num_rows = st.number_input('Add Modification', min_value=1, value=1, step=1)

    # columns to lay out the inputs
    grid = st.columns([3,2])

    # Function to create a row of widgets (with row number input to assure unique keys)
    def add_row(row):
        with grid[0]:
            st.multiselect('Residues', key=f'residues{row}', options=list('ACDEFGHIKLMNPQRSTVWY'))
        with grid[1]:
            st.number_input('Mass', step=0.01, key=f'mass{row}')


    # Loop to create rows of input widgets
    for r in range(num_rows):
        add_row(r)

    mods = {}
    for r in range(num_rows):
        for residue in st.session_state[f'residues{r}']:
            mods[residue] = st.session_state[f'mass{r}']

sites = identify_cleavage_sites(sequence, enzyme_regexes)
sites = sorted(sites, reverse=True)
sequence_with_sites = sequence
for site in sites:
    sequence_with_sites = sequence_with_sites[:site] + '%' + sequence_with_sites[site:]

# color all | in red
sequence_with_sites = sequence_with_sites.replace('%', '<span style="color:red">%</span>')
st.markdown(sequence_with_sites, unsafe_allow_html=True)

digested_sequences, mc_status, semi_status = digest_sequence(sequence, enzyme_regexes, missed_cleavages, min_len,
                                                             max_len, semi_enzymatic, info=True)
df = pd.DataFrame()
df['Sequence'] = list(digested_sequences)
df['MC'] = list(mc_status)
df['Semi'] = list(semi_status)
#make bool
df['Semi'] = df['Semi'].apply(lambda x: bool(x))
df['NeutralMass'] = [calculate_mass(sequence) for sequence in df['Sequence']]
df = df[(df['NeutralMass'] >= min_mass) & (df['NeutralMass'] <= max_mass)]
df['Len'] = df['Sequence'].apply(len)
df.drop_duplicates(inplace=True)
#df['Start'] = df['Sequence'].apply(lambda x: sequence.find(x))
df.sort_values(by=['NeutralMass'], inplace=True)

df['Sequence'] = df['Sequence'].apply(lambda x: add_static_mods(x, mods))

df_download = df.to_csv(index=False)


def make_clickable(sequence):
    # target _blank to open new window
    # extract clickable text to display for your link
    link = LINK + f'?sequence={sequence}'
    print(link)
    return f'<a target="_blank" href="{link}">{sequence}</a>'

# link is the column with hyperlinks
df['Sequence'] = df['Sequence'].apply(make_clickable)
#df['Link'] = df['Sequence'].apply(lambda x: LINK + f'?sequence={x}')
df = df.to_html(escape=False)

st.write(df, unsafe_allow_html=True, use_container_width=True)
st.write(' ')
st.download_button('Download CSV', df_download, 'digestion.csv', 'text/csv', use_container_width=True)



