import copy

import pandas as pd
import streamlit as st
from peptacular.sequence import identify_cleavage_sites, calculate_mass, add_static_mods, \
    _digest_sequence, get_left_semi_sequences, get_right_semi_sequences
from peptacular.constants import PROTEASES, AMINO_ACIDS

st.set_page_config(page_title="proteindigestor", page_icon=":knife:")

st.warning('This is a work in progress. Please report any issues or suggestions to pgarrett@scripps.edu.')

def digest_sequence(sequence: str, enzyme_regex: str, missed_cleavages: int, min_len: int,
                    max_len: int, semi_enzymatic: bool):
    """
    A function that digests a given amino acid sequence based on the provided positive and negative regular expressions and
    the number of missed cleavages. It returns a list of tuples containing the digested_sequences and the number of missed
    cleavages. The returned digested_sequences will be filtered based on the provided minimum and maximum length.
    """


    digested_sequences, mc_status, semi_status, start_indexes = [], [], [], []

    cleavage_sites = identify_cleavage_sites(sequence, enzyme_regex)
    sequences, spans = _digest_sequence(sequence, cleavage_sites, missed_cleavages, min_len, max_len)

    digested_sequences.extend(sequences)
    mc_status.extend([span[2] for span in spans])
    semi_status.extend([0 for _ in range(len(sequences))])
    start_indexes.extend([span[0] for span in spans])

    if semi_enzymatic is True:
        semi_sequences, semi_mc_status, semi_semi_status, semi_starts = [], [], [], []
        for sequence, mc, start in zip(sequences, mc_status, start_indexes):

            semi_left_sequences, semi_left_offsets = get_left_semi_sequences(sequence, min_len, max_len, info=True)
            semi_sequences.extend(semi_left_sequences)
            semi_mc_status.extend([mc for _ in range(len(semi_left_sequences))])
            semi_semi_status.extend([1 for _ in range(len(semi_left_sequences))])
            semi_starts.extend([start for offset in semi_left_offsets])

            semi_right_sequences, semi_right_offsets = get_right_semi_sequences(sequence, min_len, max_len, info=True)
            semi_sequences.extend(semi_right_sequences)
            semi_mc_status.extend([mc for _ in range(len(semi_right_sequences))])
            semi_semi_status.extend([1 for _ in range(len(semi_right_sequences))])
            semi_starts.extend([start + len(sequence) + offset for offset in semi_right_offsets])

        digested_sequences.extend(semi_sequences)
        mc_status.extend(semi_mc_status)
        semi_status.extend(semi_semi_status)
        start_indexes.extend(semi_starts)

    # filter based on min / max len
    flags = [max_len >= len(sequence) >= min_len for sequence in digested_sequences]
    digested_sequences = [sequence for sequence, flag in zip(digested_sequences, flags) if flag is True ]
    mc_status = [mc for mc, flag in zip(mc_status, flags) if flag is True ]
    semi_status = [semi for semi, flag in zip(semi_status, flags) if flag is True ]
    start_indexes = [semi for semi, flag in zip(start_indexes, flags) if flag is True ]

    return digested_sequences, mc_status, semi_status, start_indexes

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

    st.markdown("""Protein sequence can be in fasta format, where the sequence is split into multiple lines.""")

    sequence = st.text_input("sequence to be digested", value=example_protein, help='Enter a protein sequence')
    sequence = sequence.replace(' ', '').replace('\n', '')

    #check if all residues are valid
    for aa in sequence:
        if aa not in AMINO_ACIDS:
            st.error(f'Invalid amino acid: {aa}')
            st.stop()

    protease_selected = st.selectbox("Select the Protease", options=list(VALID_PROTEASES.keys()), help='Select a protease')
    enzyme_regexes = VALID_PROTEASES[protease_selected]

    c1, c2 = st.columns(2)
    min_len = c1.number_input('Min peptide length', min_value=1, value=7, step=1, help='Minimum peptide length')
    max_len = c2.number_input('Max peptide length', min_value=1, value=25, step=1, help='Maximum peptide length')

    if min_len > max_len:
        st.error('Min length must be less than max length')
        st.stop()

    c3, c4 = st.columns(2)
    min_mass = c3.number_input('Min peptide mass', min_value=0, value=200, step=100, help='Minimum peptide mass')
    max_mass = c4.number_input('Max peptide mass', min_value=0, value=3000, step=100, help = 'Maximum peptide mass')

    if min_mass > max_mass:
        st.error('Min mass must be less than max mass')
        st.stop()

    missed_cleavages = st.number_input('Max number of missed cleavages', min_value=0, value=1, step=1, help='Maximum number of missed cleavages')
    semi_enzymatic = st.checkbox('Semi Enzymatic?', help='Allow semi enzymatic peptides?')

    st.subheader('Static Modifications')

    # a selection for the user to specify the number of rows
    c1, c2 = st.columns(2)
    num_rows = st.number_input('Add Modification', min_value=1, value=1, step=1, help='Add another modification row')

    # columns to lay out the inputs
    grid = st.columns([3,2])

    # Function to create a row of widgets (with row number input to assure unique keys)
    def add_row(row):
        with grid[0]:
            st.multiselect('Residues', key=f'residues{row}', options=list('ACDEFGHIKLMNPQRSTVWY'), help='Select residues to modify')
        with grid[1]:
            st.number_input('Mass', step=0.01, key=f'mass{row}', help='Modification mass')


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

st.subheader('Sequence with cleavage sites')
st.markdown(sequence_with_sites, unsafe_allow_html=True)

digested_sequences, mc_status, semi_status, start_indexes = digest_sequence(sequence, enzyme_regexes, missed_cleavages, min_len,
                                                             max_len, semi_enzymatic)
df = pd.DataFrame()
df['Sequence'] = list(digested_sequences)
df['MC'] = list(mc_status)
df['Semi'] = list(semi_status)
df['Start'] = start_indexes

#make bool
df['Semi'] = df['Semi'].apply(lambda x: bool(x))
df.drop_duplicates(inplace=True)
df['Sequence'] = df['Sequence'].apply(lambda x: add_static_mods(x, mods))

df['NeutralMass'] = [calculate_mass(sequence) for sequence in df['Sequence']]
df = df[(df['NeutralMass'] >= min_mass) & (df['NeutralMass'] <= max_mass)]
df['Len'] = df['Sequence'].apply(len)

# keep first fuplicate
df.sort_values(by=['MC'], inplace=True)
df.drop_duplicates(subset=['Sequence', 'Semi'], inplace=True)
df.sort_values(by=['Start','Len'], inplace=True)

df_download = df.to_csv(index=False)


def make_clickable(sequence):
    # target _blank to open new window
    # extract clickable text to display for your link
    link = LINK + f'?sequence={sequence}'
    return f'<a target="_blank" href="{link}">{sequence}</a>'

# link is the column with hyperlinks
df['Sequence'] = df['Sequence'].apply(make_clickable)
df = df.to_html(escape=False)

st.subheader('Peptide Results')
st.caption('Click on a sequence to see the fragment ions!')
st.write(df, unsafe_allow_html=True, use_container_width=True)
st.write(' ')
st.download_button('Download CSV', df_download, 'digestion.csv', 'text/csv', use_container_width=True)



