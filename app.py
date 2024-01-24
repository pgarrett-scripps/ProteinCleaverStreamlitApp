import uuid
from collections import Counter

import pandas as pd
import streamlit as st
from peptacular.sequence import strip_modifications, calculate_sequence_length
from peptacular.digest import identify_cleavage_sites
from peptacular.mass import valid_mass_sequence
from peptacular.constants import AMINO_ACIDS
from peptacular.spans import calculate_span_coverage
import matplotlib as mpl
import regex as reg

from constants import *
from wiki import *
from util import generate_peptide_df, coverage_string, create_colorbar, generate_app_url, fetch_sequence_from_uniprot, \
    make_clickable

st.set_page_config(page_title="proteincleaver", page_icon=":knife:", layout="wide")

# Parse query parameters
params = st.query_params
query_peptide_sequence = params.get('protein_sequence', DEFAULT_PROTEIN_SEQUENCE)
query_proteases = params.get('proteases', ';'.join(DEFAULT_PROTEASES)).split(',')
query_custom_regex = params.get('custom_regex', '')
query_missed_cleavages = int(params.get('missed_cleavages', DEFAULT_MISSED_CLEAVAGES))
query_mass_type = params.get('mass_type', DEFAULT_MASS_TYPE)
query_min_peptide_len = int(params.get('min_peptide_len', DEFAULT_MIN_PEPTIDE_LEN))
query_max_peptide_len = int(params.get('max_peptide_len', DEFAULT_MAX_PEPTIDE_LEN))
query_min_mass = float(params.get('min_mass', DEFAULT_MIN_PEPTIDE_MASS))
query_max_mass = float(params.get('max_mass', DEFAULT_MAX_PEPTIDE_MASS))
query_semi_enzymatic = params.get('semi_enzymatic', 'False').lower() == 'true'
query_infer_charge = params.get('infer_charge', 'False').lower() == 'true'
query_min_charge = int(params.get('min_charge', DEFAULT_MIN_CHARGE))
query_max_charge = int(params.get('max_charge', DEFAULT_MAX_CHARGE))
query_min_mz = float(params.get('min_mz', DEFAULT_MIN_MZ))
query_max_mz = float(params.get('max_mz', DEFAULT_MAX_MZ))
query_remove_non_proteotypic = params.get('remove_non_proteotypic', 'False').lower() == 'true'
query_n_term_static_mod = float(params.get('n_term_static_mod', 0.0))
query_c_term_static_mod = float(params.get('c_term_static_mod', 0.0))
query_num_static_mods = int(params.get('num_static_mods', DEFAULT_STATIC_MODS))
query_n_term_var_mod = float(params.get('n_term_var_mod', 0.0))
query_c_term_var_mod = float(params.get('c_term_var_mod', 0.0))
query_max_var_mods = int(params.get('max_var_mods', DEFAULT_MAX_VAR_MODS))
query_num_variable_mods = int(params.get('num_variable_mods', DEFAULT_VAR_MODS))
query_static_mods_str = params.get('static_mods', 'C:57.02146')
query_static_mods = [(s.split(':')[0], float(s.split(':')[1])) for s in query_static_mods_str.split(';') if s]
query_variable_mods_str = params.get('variable_mods', '')
query_variable_mods = [(s.split(':')[0], float(s.split(':')[1])) for s in query_variable_mods_str.split(';') if s]

# CSS to inject contained in a string
hide_table_row_index = """
            <style>
            thead tr th:first-child {display:none}
            tbody th {display:none}
            </style>
            """
# Inject CSS with Markdown
st.markdown(hide_table_row_index, unsafe_allow_html=True)

if 'protein_sequence_key' not in st.session_state:
    st.session_state['protein_sequence_key'] = 'protein_sequence'

with st.sidebar:
    st.title('Protein Cleaver 🔪')

    c1, c2 = st.columns([3, 1])
    protein_id = c1.text_input(label='Protein accession number / identifier',
                               value='',
                               help='A protein accession number / identifier')

    reset = c2.button('Reset Sequence')

    if reset:
        st.session_state['protein_sequence_key'] = 'protein_sequence' + str(uuid.uuid4())

    raw_sequence = query_peptide_sequence
    fetched_sequence = None
    if protein_id:
        response = fetch_sequence_from_uniprot(protein_id)
        if response.status_code != 200:
            st.error(f"Error fetching sequence from UniProt: {response.status_code}")
            st.stop()
        fetched_sequence = ''.join(response.text.split('\n')[1:])
        raw_sequence = fetched_sequence

    raw_sequence = st.text_area(label="Protein sequence",
                                value=raw_sequence,
                                help='An amino acid sequence to digest',
                                max_chars=MAX_PROTEIN_INPUT_LENGTH,
                                height=200,
                                key=st.session_state['protein_sequence_key'])

    if fetched_sequence and raw_sequence != fetched_sequence:
        st.warning('Input protein sequence does not match fetched protein sequence from UniProt. '
                   'Ensure that any changes made are deliberate')

    protein_sequence = raw_sequence.replace(' ', '').replace('\n', '')

    if not protein_sequence:
        st.error('Please enter a protein sequence')
        st.stop()

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
                                        default=query_proteases)

    custom_regex = c2.text_input(label='(Additional) Custom protease',
                                 value=query_custom_regex,
                                 help='A custom regular expression to use for digestion. Will be used along with '
                                      'selected proteases. For example a regex expression for trypsin would look like: '
                                      '([KR])')

    c1, c2 = st.columns(2)
    missed_cleavages = c1.number_input(label='Max missed cleavages',
                                       min_value=MIN_MISSED_CLEAVAGES,
                                       max_value=MAX_MISSED_CLEAVAGES,
                                       value=query_missed_cleavages,
                                       step=MISSED_CLEAVAGES_STEP,
                                       help='Number of missed cleavages to allow during digestion')

    mass_type = c2.selectbox(label='Mass type', options=MASS_TYPE_OPTIONS,
                             index=MASS_TYPE_OPTIONS.index(query_mass_type),
                             help='Mass type to use for calculations')
    is_mono = mass_type == 'monoisotopic'

    enzyme_regexes = [VALID_PROTEASES[protease] for protease in proteases_selected]

    if custom_regex:
        enzyme_regexes.append(custom_regex)

    c1, c2 = st.columns(2)
    min_peptide_len = c1.number_input(label='Min length',
                                      min_value=MIN_PEPTIDE_LEN,
                                      max_value=MAX_PEPTIDE_LEN,
                                      value=query_min_peptide_len,
                                      step=PEPTIDE_LEN_STEP,
                                      help='Minimum peptide length (inclusive)')

    max_peptide_len = c2.number_input(label='Max length',
                                      min_value=MIN_PEPTIDE_LEN,
                                      max_value=MAX_PEPTIDE_LEN,
                                      value=query_max_peptide_len,
                                      step=PEPTIDE_LEN_STEP,
                                      help='Maximum peptide length (inclusive)')

    if min_peptide_len > max_peptide_len:
        st.error('Min length must be less than max length')
        st.stop()

    c3, c4 = st.columns(2)
    min_mass = c3.number_input(label='Min neutral mass',
                               min_value=MIN_PEPTIDE_MASS,
                               max_value=MAX_PEPTIDE_MASS,
                               value=query_min_mass,
                               step=PEPTIDE_MASS_STEP,
                               help='Minimum peptide neutral mass (inclusive)')

    max_mass = c4.number_input(label='Max neutral mass',
                               min_value=MIN_PEPTIDE_MASS,
                               max_value=MAX_PEPTIDE_MASS,
                               value=query_max_mass,
                               step=PEPTIDE_MASS_STEP,
                               help='Maximum peptide neutral mass (inclusive)')

    if min_mass > max_mass:
        st.error('Min mass must be less than max mass')
        st.stop()

    c1, c2 = st.columns(2)

    semi_enzymatic = c1.checkbox(label='Semi enzymatic?',
                                 help='Allow semi enzymatic peptides?',
                                 value=query_semi_enzymatic)
    infer_charge = c2.checkbox(label='Infer charge',
                               value=query_infer_charge,
                               help='Infer charge of peptides based on L and K count')

    min_charge, max_charge, min_mz, max_mz = None, None, None, None
    if infer_charge:
        c1, c2 = st.columns(2)
        min_charge = c1.number_input(label='Min charge',
                                     min_value=MIN_CHARGE,
                                     max_value=MAX_CHARGE,
                                     value=query_min_charge,
                                     step=CHARGE_STEP,
                                     help='Minimum peptide charge (inclusive)')

        max_charge = c2.number_input(label='Max charge',
                                     min_value=MIN_CHARGE,
                                     max_value=MAX_CHARGE,
                                     value=query_max_charge,
                                     step=CHARGE_STEP,
                                     help='Maximum peptide charge (inclusive)')

        c1, c2 = st.columns(2)
        min_mz = c1.number_input(label='Min m/z',
                                 min_value=MIN_MZ,
                                 max_value=MAX_MZ,
                                 value=query_min_mz,
                                 step=MZ_STEP,
                                 help='Minimum peptide m/z (inclusive)')

        max_mz = c2.number_input(label='Max m/z',
                                 min_value=MIN_MZ,
                                 max_value=MAX_MZ,
                                 value=query_max_mz,
                                 step=MZ_STEP,
                                 help='Maximum peptide m/z (inclusive)')

        if min_charge > max_charge:
            st.error('Min charge must be less than max charge')
            st.stop()

    remove_non_proteotypic = st.checkbox(label='Remove non-proteotypic peptides',
                                         value=query_remove_non_proteotypic,
                                         help='Remove peptides that are not proteotypic')

    with st.expander('Static Modifications'):

        c1, c2 = st.columns(2)
        n_term_static_mod = c1.number_input(label='N-term mod',
                                            value=query_n_term_static_mod,
                                            help='Apply a static modification to the N-terminus')

        c_term_static_mod = c2.number_input(label='C-term mod',
                                            value=query_c_term_static_mod,
                                            help='Apply a static modification to the C-terminus')

        num_static_mods = st.number_input(label='Number of unique static modifications',
                                          min_value=MIN_STATIC_MODS,
                                          max_value=MAX_STATIC_MODS,
                                          value=query_num_static_mods,
                                          step=STATIC_MODS_STEP,
                                          help='Add another modification row')

        # columns to lay out the inputs
        grid = st.columns([3, 2])


        def add_static_modification(r):

            aas = list(query_static_mods[r][0]) if r < len(query_static_mods) else []
            mod = query_static_mods[r][1] if r < len(query_static_mods) else 0.0

            with grid[0]:
                st.multiselect(label='Amino acids',
                               key=f'static_mod_residue{r}',
                               options=list(AMINO_ACIDS),
                               help='Select amino acids for which to apply the static modification',
                               default=aas)
            with grid[1]:
                st.number_input(label='Modification Mass (Da)',
                                step=0.00001,
                                key=f'static_mod_mass{r}',
                                help='The mass of the modification (in daltons)',
                                value=mod,
                                format='%.5f')


        # Loop to create rows of input widgets
        for r in range(num_static_mods):
            add_static_modification(r)

        static_mods = {}
        for r in range(num_static_mods):
            mod = "{:.5f}".format(st.session_state[f'static_mod_mass{r}'])
            for residue in st.session_state[f'static_mod_residue{r}']:
                static_mods[residue] = mod

    with st.expander('Variable Modifications'):

        c1, c2 = st.columns(2)
        n_term_var_mod = c1.number_input(label='N-term var mod',
                                         value=query_n_term_var_mod,
                                         help='Apply a variable modification to the N-terminus')
        c_term_var_mod = c2.number_input(label='C-term var mod',
                                         value=query_c_term_var_mod,
                                         help='Apply a variable modification to the C-terminus')

        max_var_mods = st.number_input(label='Max var mods',
                                       min_value=MIN_MAX_VAR_MODS,
                                       max_value=MAX_MAX_VAR_MODS,
                                       value=query_max_var_mods,
                                       step=MAX_VAR_MOD_STEP,
                                       help='Maximum number of variable modifications per peptide')

        num_variable_mods = st.number_input(label='Unique var mods',
                                            min_value=MIN_VAR_MODS,
                                            max_value=MAX_VAR_MODS,
                                            value=query_num_variable_mods,
                                            step=VAR_MOD_STEP,
                                            help='Add another modification row')

        # columns to lay out the inputs
        grid = st.columns([3, 2])

        def add_variable_modification(r):
            aas = list(query_variable_mods[r][0]) if r < len(query_variable_mods) else []
            mod = query_variable_mods[r][1] if r < len(query_variable_mods) else 0.0

            with grid[0]:
                st.multiselect(label='Amino acids',
                               key=f'var_mod_residue{r}',
                               options=list(AMINO_ACIDS),
                               help='Select amino acids for which to apply the variable modification',
                               default=aas)
            with grid[1]:
                st.number_input(label='Modification Mass (Da)',
                                step=0.00001, key=f'var_mod_mass{r}',
                                help='The mass of the modification (in daltons)',
                                value=mod,
                                format='%.5f')


        # Loop to create rows of input widgets
        for r in range(num_variable_mods):
            add_variable_modification(r)

        var_mods = {}
        for r in range(num_variable_mods):
            mod = "{:.5f}".format(st.session_state[f'var_mod_mass{r}'])
            for residue in st.session_state[f'var_mod_residue{r}']:
                var_mods[residue] = mod

    url = generate_app_url(protein_id, protein_sequence, proteases_selected, custom_regex, missed_cleavages, mass_type,
                           min_peptide_len, max_peptide_len, min_mass, max_mass, semi_enzymatic, infer_charge,
                           min_charge, max_charge, min_mz, max_mz, remove_non_proteotypic, n_term_static_mod,
                           c_term_static_mod, num_static_mods, n_term_var_mod, c_term_var_mod, max_var_mods,
                           num_variable_mods, static_mods, var_mods)

sites = set()
for enzyme_regex in enzyme_regexes:
    sites.update(identify_cleavage_sites(protein_sequence, enzyme_regex))
sites = sorted(list(sites))

with st.expander('Edit Sites'):
    sites = st.multiselect(label="Sites",
                            options=list(range(len(stripped_protein_sequence)+1)),
                            help='The proteases to use for digestion',
                            default=sites)


df = generate_peptide_df(protein_sequence, sites, missed_cleavages, min_peptide_len, max_peptide_len,
                         semi_enzymatic, static_mods, min_mass, max_mass, is_mono, infer_charge, min_charge,
                         max_charge, min_mz, max_mz, var_mods, max_var_mods, n_term_static_mod, c_term_static_mod,
                         n_term_var_mod, c_term_var_mod, remove_non_proteotypic)

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

st.write(f'##### [Analysis URL]({url}) (copy me and send to your friends!)')

t1, t2, t3, t5 = st.tabs(['Digestion Metrics', 'Cleavage & Coverage', 'Motif Analysis', 'Help'])

with t1:
    st.header('Digestion Metrics')
    c1, c2, c3 = st.columns(3)
    c1.metric('Total Peptides', len(df))
    c2.metric('Semi Peptides', len(df[df['Semi']]))
    c3.metric('Enzymatic Peptides', len(df[~df['Semi']]))

    st.subheader('Peptides')

    df['Link'] = [make_clickable(peptide, mass_type) for peptide in df['Sequence']]

    st.dataframe(
        df,
        column_config={
            "Link": st.column_config.LinkColumn(
                display_text="View Ions"),
        },
        hide_index=True,
    )



with t2:
    st.header('Cleavage & Coverage')

    c1, c2 = st.columns(2)
    c1.metric('Cleavage Sites', len(sites))

    protein_cov_arr_bin = calculate_span_coverage(spans, protein_length, accumulate=False)
    protein_cov_perc = round(sum(protein_cov_arr_bin) / len(protein_cov_arr_bin) * 100, 2)
    c2.metric('Protein Coverage', f'{protein_cov_perc}%')

    st.subheader('Site Indexes')
    st.markdown(site_indexes_html, unsafe_allow_html=True)
    st.write("")
    st.subheader('Sites')
    st.markdown(sequence_with_sites, unsafe_allow_html=True)

    st.subheader('Sequence Coverage')
    st.markdown(protein_coverage, unsafe_allow_html=True)

    # Example usage in a Streamlit app
    f = create_colorbar(max(protein_cov_arr), cmap)
    st.pyplot(f)

    st.caption('Coverage vs Missed Cleavages')
    st.line_chart(data={'Missed Cleavages': mcs, 'Protein Coverage (%)': protein_cov_at_mcs},
                  x='Missed Cleavages', y='Protein Coverage (%)')

    st.caption('Coverage vs Peptide Lengths')
    st.line_chart(data={'Peptide Length': lens, 'Protein Coverage (%)': protein_cov_at_lens},
                  x='Peptide Length', y='Protein Coverage (%)')

    st.caption('"Coverage vs Peptide Masses')
    st.line_chart(data={'Peptide Mass': masses, 'Protein Coverage (%)': protein_cov_at_mass},
                  x='Peptide Mass', y='Protein Coverage (%)')

with t3:
    st.header('Motif Analysis')
    c1, c2, c3 = st.columns(3)
    motif_regex = c1.text_input('Motifs Regex', '(K)')

    st.cache_data()
    def get_motif_sites(motif_regex, stripped_protein_sequence):
        motif_sites = list(reg.finditer(motif_regex, stripped_protein_sequence, overlapped=True))
        return motif_sites

    if motif_regex:
        motif_sites = get_motif_sites(motif_regex, stripped_protein_sequence)

        def count_motifs(row):
            return sum([1 for site in motif_sites if row['Start'] <= site.start() < row['End']])

        df['Motifs'] = df.apply(count_motifs, axis=1)
        motif_cov_indexes = {i for site in motif_sites for i in range(site.start(), site.end())}

        motif_cov_array = [0] * len(stripped_protein_sequence)
        for row in df[['Start', 'End', 'Motifs']].values:

            if row[2] == 0:
                continue

            for i in range(row[0], row[1] + 1):
                if i in motif_cov_indexes:
                    if row[2] == 1:
                        motif_cov_array[i] = 1
                    else:
                        if motif_cov_array[i] == 0:
                            motif_cov_array[i] = row[2]
                        else:
                            motif_cov_array[i] = min(row[2], motif_cov_array[i])


        min_moitifs = c2.number_input('Min Motifs', min_value=0, max_value=max(df['Motifs']), value=0)
        max_motifs = c3.number_input('Max Motifs', min_value=0, max_value=max(df['Motifs']), value=max(df['Motifs']))
        df = df[(df['Motifs'] >= min_moitifs) & (df['Motifs'] <= max_motifs)]

        st.subheader('Peptides')

        # Make the Link column the last column int he dataframe
        df = df[[c for c in df if c not in ['Link']] + ['Link']]

        st.dataframe(
            df,
            column_config={
                "Link": st.column_config.LinkColumn(
                    display_text="View Ions"),
            },
            hide_index=True,
        )

        counter = Counter(df['Motifs'])


        st.subheader('Motif Site Coverage', help='The color corresponds to the peptide with the fewest number of motif '
                                                 'matches (excluding 0 matches). Example: Lets assume that the first '
                                                 'site is covered by two peptides, the first with one match and the '
                                                 'second with two matches, then this site will be displayed as 1.')
        protein_coverage = coverage_string(motif_cov_array, stripped_protein_sequence, cmap, motif_cov_indexes)
        st.markdown(protein_coverage, unsafe_allow_html=True)
        f = create_colorbar(max(motif_cov_array), cmap, label='Min Number of Motif Matches')
        st.pyplot(f)
        with st.expander('Show Coverage Table'):
            for i, (k, v) in enumerate(sorted(counter.items())):
                df_tmp = df[df['Motifs'] == k]
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

# with t4:
#    st.markdown(PROTEASE_WIKI)

with t5:
    st.header('Help')

    st.subheader('General')

    with st.expander('Protein Cleaver Overview'):
        st.markdown(HELP)

    with st.expander('Column Descriptions'):
        st.markdown(COLUMN_DESCRIPTIONS)

    with st.expander('Protease Regexes'):
        st.subheader('Protease Regexes')
        data = [{'Name': k, 'Regex': v} for k, v in VALID_PROTEASES.items()]
        protease_df = pd.DataFrame(data)
        st.table(protease_df)

    with st.expander('Contact'):
        st.markdown(CONTACT)

    st.subheader('Models')

    with st.expander('IM Model'):
        st.markdown(IM_MODEL_HELP)

    with st.expander('RT Model'):
        st.markdown(RT_MODEL_HELP)

    with st.expander('Proteotypic Model'):
        st.markdown(PROTEOTYPIC_MODEL_HELP)

    with st.expander('How to use ML Models?'):
        def get_model_file_as_byte_stream(path):
            with open(path, 'rb') as file:
                byte_stream = file.read()
            return byte_stream


        st.subheader('Download Models')
        # download models
        c1, c2, c3 = st.columns(3)
        c1.download_button(
            label='RT Model',
            data=get_model_file_as_byte_stream('rt_model.pkl'),
            file_name="rt_model.pkl",
            mime='application/octet-stream'
        )
        c2.download_button(
            label='IM Model',
            data=get_model_file_as_byte_stream('im_model.pkl'),
            file_name="im_model.pkl",
            mime='application/octet-stream'
        )
        c3.download_button(
            label='Proteotypic Model',
            data=get_model_file_as_byte_stream('proteotypic_model.pkl'),
            file_name="proteotypic_model.pkl",
            mime='application/octet-stream'
        )

        st.subheader('Example Code')
        st.code(MODEL_CODE)
