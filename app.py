import pickle
from collections import Counter

import pandas as pd
import streamlit as st
import streamlit_permalink as stp
import numpy as np
import requests
import peptacular as pt
import matplotlib as mpl

from constants import *
from wiki import *
from util import make_clickable, generate_peptide_df, bin_aa_counts, coverage_string, create_colorbar


def fetch_sequence_from_uniprot(accession_number):
    url = f"https://www.uniprot.org/uniprot/{accession_number}.fasta"
    response = requests.get(url)
    if response.status_code != 200:
        st.error(f"Error fetching sequence from UniProt: {response.status_code}")
        return None
    return ''.join(response.text.split('\n')[1:])  # Remove the header line


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

with st.sidebar:


    # Inject custom CSS to set the width of the sidebar
    st.markdown(
        """
        <style>
            section[data-testid="stSidebar"] {
                width: 600px !important; # Set the width to your desired value
            }
        </style>
        """,
        unsafe_allow_html=True,
    )

    st.title('Protein Cleaver 🔪')

    c1, c2 = st.columns(2, vertical_alignment='bottom')
    with c1:
        protein_id = st.text_input(label='Protein accession number / identifier',
                                value='',
                                placeholder='P12345',
                                help='A protein accession number / identifier',
                                max_chars=20,
                                key='protein_id')

    with c2:
        # btn to fetch sequence from uniprot
        fetch_sequence = st.button(label='Fetch from UniProt',
                                help='Fetch protein sequence from UniProt using the accession number',
                                use_container_width=True,
                                key='fetch_sequence')
        
    if fetch_sequence:
        raw_sequence = None
        if protein_id:
            fetched_protein_sequence = fetch_sequence_from_uniprot(protein_id)
            if fetched_protein_sequence is not None:
                # update URL value
                st.query_params['protein_sequence'] = fetched_protein_sequence
                st.rerun()
        else:
            st.error('Please enter a valid protein accession number')

    raw_sequence = stp.text_area(label="Protein sequence",
                                value=DEFAULT_PROTEIN_SEQUENCE,
                                help='An amino acid sequence to digest',
                                key='protein_sequence',
                                max_chars=MAX_PROTEIN_INPUT_LENGTH,
                                height=200)

    protein_sequence = raw_sequence.replace(' ', '').replace('\n', '')
    stripped_protein_sequence = pt.strip_mods(protein_sequence)
    protein_length = pt.sequence_length(protein_sequence)
    st.caption(f'Length: {protein_length}')

    if protein_length > MAX_PROTEIN_LEN:
        st.error(f'Protein sequence is too long. Please keep it under {MAX_PROTEIN_LEN} residues.')
        st.stop()

    # check if all residues are valid
    try:
        mass = pt.mass(protein_sequence)
    except Exception as e:
        st.error('Invalid amino acid sequence. Please check your input.')
        st.stop()

    c1, c2 = st.columns(2)
    with c1:
        proteases_selected = stp.multiselect(label="Proteases",
                                            options=list(VALID_PROTEASES.keys()),
                                            help='The proteases to use for digestion',
                                            default=DEFAULT_PROTEASES,
                                            key='protease')

    with c2:
        custom_regex = stp.text_input(label='(Additional) Custom protease',
                                    value='',
                                    help='A custom regular expression to use for digestion. Will be used along with '
                                        'selected proteases',
                                        key='custom_regex')

    c1, c2 = st.columns(2)
    with c1:
        missed_cleavages = stp.number_input(label='Max missed cleavages',
                                        min_value=MIN_MISSED_CLEAVAGES,
                                        max_value=MAX_MISSED_CLEAVAGES,
                                        value=DEFAULT_MISSED_CLEAVAGES,
                                        step=MISSED_CLEAVAGES_STEP,
                                        help='Number of missed cleavages to allow during digestion',
                                        key='missed_cleavages')

    with c2:
        mass_type = stp.selectbox(label='Mass type', options=MASS_TYPE_OPTIONS, index=DEFAULT_MASS_TYPE_INDEX,
                                help='Mass type to use for calculations', key='mass_type')
        is_mono = mass_type == 'monoisotopic'

    enzyme_regexes = [VALID_PROTEASES[protease] for protease in proteases_selected]

    if custom_regex:
        enzyme_regexes.append(custom_regex)

    c1, c2 = st.columns(2)
    with c1:
        min_peptide_len = stp.number_input(label='Min length',
                                        min_value=MIN_PEPTIDE_LEN,
                                        max_value=MAX_PEPTIDE_LEN,
                                        value=DEFAULT_MIN_PEPTIDE_LEN,
                                        step=PEPTIDE_LEN_STEP,
                                        help='Minimum peptide length (inclusive)',
                                        key='min_peptide_len')

    with c2:
        max_peptide_len = stp.number_input(label='Max length',
                                        min_value=MIN_PEPTIDE_LEN,
                                        max_value=MAX_PEPTIDE_LEN,
                                        value=DEFAULT_MAX_PEPTIDE_LEN,
                                        step=PEPTIDE_LEN_STEP,
                                        help='Maximum peptide length (inclusive)',
                                        key='max_peptide_len')

    if min_peptide_len > max_peptide_len:
        st.error('Min length must be less than max length')
        st.stop()

    c3, c4 = st.columns(2)
    with c3:
        min_mass = stp.number_input(label='Min neutral mass',
                                min_value=MIN_PEPTIDE_MASS,
                                max_value=MAX_PEPTIDE_MASS,
                                value=DEFAULT_MIN_PEPTIDE_MASS,
                                step=PEPTIDE_MASS_STEP,
                                help='Minimum peptide neutral mass (inclusive)',
                                key='min_peptide_mass')

    with c4:
        max_mass = stp.number_input(label='Max neutral mass',
                                min_value=MIN_PEPTIDE_MASS,
                                max_value=MAX_PEPTIDE_MASS,
                                value=DEFAULT_MAX_PEPTIDE_MASS,
                                step=PEPTIDE_MASS_STEP,
                                help='Maximum peptide neutral mass (inclusive)',
                                key='max_peptide_mass')

    if min_mass > max_mass:
        st.error('Min mass must be less than max mass')
        st.stop()

    semi_enzymatic = stp.checkbox(label='Semi enzymatic?',
                            help='Allow semi enzymatic peptides?',
                            value=False,
                            key='semi_enzymatic')

    with st.expander('Advanced Options (Peptide Properties)'):
           
        infer_charge = stp.checkbox(label='Infer Charge', 
                                    value=False, 
                                    help='Infer charge of peptides based on +/- 1 of (#Lysines and #Arginines + 1)',
                                    key='infer_charge')

        min_charge, max_charge, min_mz, max_mz = None, None, None, None
        plus_minus_charge = 0
        if infer_charge:
            c1, c2 = st.columns(2)
            with c1:
                min_charge = stp.number_input(label='Min charge',
                                            min_value=MIN_CHARGE,
                                            max_value=MAX_CHARGE,
                                            value=DEFAULT_MIN_CHARGE,
                                            step=CHARGE_STEP,
                                            help='Minimum peptide charge (inclusive)',
                                            key='min_charge')

                min_mz = stp.number_input(label='Min m/z',
                                min_value=MIN_MZ,
                                max_value=MAX_MZ,
                                value=DEFAULT_MIN_MZ,
                                step=MZ_STEP,
                                help='Minimum peptide m/z (inclusive)',
                                key='min_mz')

            with c2:
                max_charge = stp.number_input(label='Max charge',
                                            min_value=MIN_CHARGE,
                                            max_value=MAX_CHARGE,
                                            value=DEFAULT_MAX_CHARGE,
                                            step=CHARGE_STEP,
                                            help='Maximum peptide charge (inclusive)',
                                            key='max_charge')

                max_mz = stp.number_input(label='Max m/z',
                                            min_value=MIN_MZ,
                                            max_value=MAX_MZ,
                                            value=DEFAULT_MAX_MZ,
                                            step=MZ_STEP,
                                            help='Maximum peptide m/z (inclusive)',
                                            key='max_mz')

            if min_charge > max_charge:
                st.error('Min charge must be less than max charge')
                st.stop()

            plus_minus_charge = stp.number_input(label='Charge +/-',
                                                min_value=-1,
                                                max_value=1,
                                                value=0,
                                                step=1,
                                                help='Charge +/- for m/z calculations',
                                                key='plus_minus_charge')

        infer_retention_time = stp.checkbox(label='Infer RT',
                                            value=True,
                                            help='Use retention time model to predict retention times for peptides',
                                            key='infer_retention_time')
        
        if infer_retention_time:
            c1, c2 = st.columns(2, vertical_alignment='bottom')
            with c1:
                retention_time = stp.number_input(label='RT (min)',
                                                min_value=0,
                                                value=60,
                                                help='Retention time of the peptide in minutes',
                                                key='retention_time')

            with c2:
                filter_invalid_rt = stp.checkbox(label='Filter invalid RT',
                                                value=False,
                                                help='Filter out peptides with invalid retention times',
                                                key='filter_invalid_rt')


        infer_ion_mobility = stp.checkbox(label='Infer Ion Mobility',
                                            value=True,
                                            help='Use ion mobility model to predict ion mobility for peptides',
                                            key='infer_ion_mobility')

        infer_proteotypic = stp.checkbox(label='Infer Proteotypic-ness',
                                            value=True,
                                            help='Use proteotypic model to predict proteotypic peptides',
                                            key='infer_proteotypic')
        
        if infer_proteotypic:
            c1, c2 = st.columns(2, vertical_alignment='bottom')
            with c1:
                remove_non_proteotypic = stp.checkbox(label='Filter Non Proteotypic Peptides', value=False,
                                                    help='Remove peptides that are not proteotypic')
                
            with c2:
                score_threshold = stp.number_input(label='Proteotypic Score Threshold',
                                                    min_value=0.0,
                                                    max_value=1.0,
                                                    value=0.5,
                                                    step=0.01,
                                                    help='Proteotypic score threshold for filtering peptides',
                                                    key='proteotypic_score_threshold')

    with st.expander('Static Modifications'):

        df = pd.DataFrame({'Amino Acids': ['C'], 'Mass (Da)': [57.02146]})

        df = stp.data_editor(df,
                            use_container_width=True,
                            hide_index=True,
                            key='static_mods',
                            num_rows = 'dynamic',
                            column_config={
                                'Amino Acids': stp.column_config.SelectboxColumn(
                                    options=list(pt.AMINO_ACIDS),
                                    help='Select amino acids for which to apply the static modification',
                                    default='C',
                                ),
                                'Mass (Da)': stp.column_config.NumberColumn(
                                    format='%.5f',
                                    help='The mass of the modification (in daltons)',
                                    default=57.02146,
                                    
                                )
                            },)

        # get mods from df
        mods = {}
        for i, row in df.iterrows():
            for residue in row['Amino Acids']:
                mods[residue] = "{:.5f}".format(row['Mass (Da)'])

sites = set()
for enzyme_regex in enzyme_regexes:
    sites.update(pt.get_cleavage_sites(protein_sequence, enzyme_regex))
sites = sorted(list(sites))

df = generate_peptide_df(protein_sequence, sites, missed_cleavages, min_peptide_len, max_peptide_len,
                         semi_enzymatic, mods, min_mass, max_mass, is_mono, infer_charge, min_charge,
                         max_charge, min_mz, max_mz, plus_minus_charge)


if len(df) == 0:
    st.error('No peptides found. Please adjust your parameters.')
    st.stop()

if infer_retention_time:
    rt_model = pickle.load(open("rt_model.pkl", "rb"))
    df['RT'] = rt_model.predict(np.array([bin_aa_counts(pt.strip_mods(seq)) for seq in df['Sequence']]))
    df['RT'] = df['RT'].round(3)

    #scale to retention time
    df['RT'] = df['RT'] * retention_time 

    if filter_invalid_rt:
        df = df[df['RT'] <= retention_time]
        # greater than 0
        df = df[df['RT'] > 0]

if 'Charge' in df.columns:
    im_model = pickle.load(open("im_model.pkl", "rb"))
    df['IM'] = im_model.predict(
        np.array([bin_aa_counts(pt.strip_mods(seq), c) for seq, c in df[['Sequence', 'Charge']].values]))
    df['IM'] = df['IM'].round(3)

if infer_proteotypic:
    proteotypic_model = pickle.load(open("proteotypic_model.pkl", "rb"))

    probs = proteotypic_model.predict_proba(
        np.array([bin_aa_counts(pt.strip_mods(seq)) for seq in df['Sequence']]))
    
    df['score'] = probs[:, 1]
    df['Proteotypic'] = probs[:, 1] >= score_threshold

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
protein_cov_arr = pt.calculate_span_coverage(spans, protein_length, accumulate=True)
protein_coverage = coverage_string(protein_cov_arr, stripped_protein_sequence, cmap)

# calculate protein coverage at different MC
protein_cov_at_mcs = []
mcs = [mc for mc in range(0, missed_cleavages + 1)]
for mc in mcs:
    df_mc = df[df['MC'] <= mc]
    spans = [(s, e, mc) for s, e, mc in df_mc[['Start', 'End', 'MC']].values]
    cov = pt.calculate_span_coverage(spans, protein_length)
    protein_cov_at_mcs.append(sum(cov) / len(cov) * 100)

# calculate protein coverage at different peptide lengths
protein_cov_at_lens = []
lens = [l for l in range(min_peptide_len, max_peptide_len + 1)]
for l in lens:
    df_len = df[df['Len'] <= l]
    spans = [(s, e, mc) for s, e, mc in df_len[['Start', 'End', 'MC']].values]
    cov = pt.calculate_span_coverage(spans, protein_length)
    protein_cov_at_lens.append(sum(cov) / len(cov) * 100)

# calculate protein coverage at different peptide Mass
protein_cov_at_mass = []
masses = [m for m in range(int(min_mass), int(max_mass) + 1, 100)]
for m in masses:
    df_mass = df[df['NeutralMass'] <= m]
    spans = [(s, e, mc) for s, e, mc in df_mass[['Start', 'End', 'MC']].values]
    cov = pt.calculate_span_coverage(spans, protein_length)
    protein_cov_at_mass.append(sum(cov) / len(cov) * 100)

df.drop(columns=['StrippedPeptide'], inplace=True)
df.sort_values(by=['MC'], inplace=True)
if 'Charge' in df.columns:
    df.drop_duplicates(subset=['Start', 'Sequence', 'Semi', 'Charge'], inplace=True)
else:
    df.drop_duplicates(subset=['Start', 'Sequence', 'Semi'], inplace=True)
df.sort_values(by=['Start', 'Len'], inplace=True)




t1, t2, t3, t5 = st.tabs(['Peptides', 'Cleavage & Coverage', 'Motif Analysis', 'Help'])


with t1:

    c1, c2, c3, c4 = st.columns(4)
    c1.metric('Total Peptides', len(df))
    c2.metric('Semi Peptides', len(df[df['Semi']]))
    c3.metric('Enzymatic Peptides', len(df[~df['Semi']]))
    c4.metric('Unique Peptides', len(df['Sequence'].unique()))

    df['Frag_Ions'] = [make_clickable(peptide, mass_type) for peptide in df['Sequence']]
                                
    st.dataframe(df, use_container_width=True, hide_index=True,
                    column_config={
                        'Start': stp.column_config.NumberColumn(
                            help='Start index of the peptide in the protein sequence',
                            format='%d'
                        ),
                        'End': stp.column_config.NumberColumn(
                            help='End index of the peptide in the protein sequence',
                            format='%d'
                        ),
                        'MC': stp.column_config.NumberColumn(
                            help='Number of missed cleavages for this peptide',
                            format='%d'
                        ),
                        'Sequence': stp.column_config.TextColumn(
                            help='Peptide sequence',
                        ),
                        'Semi': stp.column_config.CheckboxColumn(
                            help='Is this peptide semi-enzymatic?',
                        ),
                        'Len': stp.column_config.NumberColumn(
                            help='Length of the peptide',
                            format='%d'
                        ),
                        'NeutralMass': stp.column_config.NumberColumn(
                            help='Neutral mass of the peptide',
                            format='%.5f'
                        ),
                        'Charge': stp.column_config.NumberColumn(
                            help='Charge of the peptide',
                            format='%d'
                        ),
                        'Mz': stp.column_config.NumberColumn(
                            help='m/z of the peptide',
                            format='%.5f'
                        ),
                        'RT': stp.column_config.NumberColumn(
                            help='Retention time of the peptide',
                            format='%.3f'
                        ),
                        'IM': stp.column_config.NumberColumn(
                            help='Ion mobility of the peptide',
                            format='%.3f'
                        ),
                        'Proteotypic': stp.column_config.CheckboxColumn(
                            help='Is this peptide proteotypic?',
                        ),
                        'score': stp.column_config.NumberColumn(
                            help='Proteotypic score of the peptide',
                            format='%.3f'
                        ),
                        'Frag_Ions': stp.column_config.LinkColumn(
                            help='Link to Peptide Fragmenter',
                            display_text='🔗',
                            width ='small',
                            
                        )
                    })

with t2:
    c1, c2 = st.columns(2)
    c1.metric('Cleavage Sites', len(sites))

    protein_cov_arr_bin = pt.calculate_span_coverage(spans, protein_length, accumulate=False)
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
    site_regex = st.text_input('Motifs Regex', 'N[ST][^P]')

    if site_regex:
        site_ranges = list(pt.get_regex_match_range(stripped_protein_sequence, site_regex))
        site_counts = []
        for row in df[['Start', 'End']].values:
            peptide_start, peptide_end = row[0], row[1]
            motif_cnt = 0
            for site_start, site_end in site_ranges:
                # add motif count if the site is within the peptide sequence
                if peptide_start <= site_start < peptide_end or peptide_start < site_end <= peptide_end:
                    motif_cnt += 1
            site_counts.append(motif_cnt)
        df['Motifs'] = site_counts

        if len(site_counts) == 0:
            min_motifs, max_motifs = st.slider(label='Min/Max Motif Matches', min_value=0, max_value=max(site_counts), value=(0, max(site_counts)))
            df = df[(df['Motifs'] >= min_motifs) & (df['Motifs'] <= max_motifs)]
        tmp_spans = [(s, e, mc) for s, e, mc in df[['Start', 'End', 'MC']].values]
        cov_site_mat = pt.calculate_span_coverage(tmp_spans, protein_length, accumulate=True)
        
        indexes_to_keep = set()
        for s, e in site_ranges:
            for i in range(s, e):
                indexes_to_keep.add(i)

        # set the coverage matrix to 0 for all indexes not in the site ranges
        for i in range(len(cov_site_mat)):
            if i not in indexes_to_keep:
                cov_site_mat[i] = 0

        # drop Frag_Ions column
        df.drop(columns=['Frag_Ions'], inplace=True)
        counter = Counter(site_counts)

        st.subheader('Coverage Analysis', help='Coverage of protein based on peptides with N number of motif matches')
        #sites = [i - 1 for i in sites if i != 0]
        print(list(indexes_to_keep))
        protein_coverage = coverage_string(cov_site_mat, stripped_protein_sequence, cmap, sites=list(indexes_to_keep))
        st.markdown(protein_coverage, unsafe_allow_html=True)
        f = create_colorbar(max(cov_site_mat), cmap, label='Min Number of Motif Matches')
        st.pyplot(f)

        st.dataframe(df, use_container_width=True, hide_index=True)


#with t4:
#    st.markdown(PROTEASE_WIKI)

with t5:

    st.subheader('General')

    with st.expander('Protein Cleaver Overview'):
        st.markdown(HELP)

    with st.expander('Column Descriptions'):
        st.markdown(COLUMN_DESCRIPTIONS)

    with st.expander('Protease Regexes'):
        st.subheader('Protease Regexes')
        data = [{'Name': k, 'Regex': v} for k, v in PROTEASES.items()]
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


