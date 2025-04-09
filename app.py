import pickle
from collections import Counter

import pandas as pd
import streamlit as st
import streamlit_permalink as stp
import numpy as np
import peptacular as pt
import matplotlib as mpl
from streamlit_js_eval import get_page_location

from constants import *
from wiki import *
from util import (make_clickable, 
                  generate_peptide_df,
                    bin_aa_counts, 
                    coverage_string, 
                    create_colorbar, 
                    get_query_params_url,
                    shorten_url,
                    get_site_index_html,
                    get_sequence_site_html)
from input_options import get_input_options, validate_input_option

st.set_page_config(page_title='Protein Cleaver', page_icon='🔪', layout='wide')

def inject_header_ccs():
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

def display_title():
    st.markdown(
        f"""
        <div style='text-align: center; padding: 15px; top-margin: 0px'>
            <h3 style='margin: 0; font-size: 1.5em; color: #333;'>Protein Cleaver 🔪</h3>
            <p style='font-size: 1.0em; line-height: 1.6; color: #555;'>
                Analyze and digest protein sequences using specified protease(s). 
                Generate peptides, calculate thier properties, and visualize resulting sequence coverage. Protein must be 
                <a href="https://peptacular.readthedocs.io/en/latest/modules/getting_started.html#proforma-notation" 
                target="_blank" style='color: #007BFF; text-decoration: none;'>proforma2.0 compliant</a>.
            </p>
        </div>
    """,
        unsafe_allow_html=True,
    )


# side bar and header
with st.sidebar:
    inject_header_ccs()
    display_title()
    input_options = get_input_options()

is_valid = validate_input_option(input_options)
df = generate_peptide_df(input_options)

if df is None or len(df) == 0:
    st.error('No peptides found. Please check your input options.')
    is_valid = False

if is_valid:
    top_window, bottom_window = st.container(), st.container()

    with bottom_window:
        page_loc = get_page_location()

    with top_window:

        title_c, _, button_c = st.columns([2, 1, 1])
        help_msg = "This page's URL automatically updates with your input and can be shared with others. You can optionally use the Generate TinyURL button to create a shortened URL."
        title_c.header("Results", help=help_msg)

        # Create a colormap
        cmap = mpl.colormaps.get_cmap(CMAP)

        spans = [(s, e, mc) for s, e, mc in df[['Start', 'End', 'MC']].values]
        protein_cov_arr = pt.calculate_span_coverage(spans, input_options.sequence_length, accumulate=True)
        protein_coverage = coverage_string(protein_cov_arr, input_options.stripped_protein_sequence, cmap)

        # calculate protein coverage at different MC
        protein_cov_at_mcs = []
        mcs = [mc for mc in range(0, input_options.missed_cleavages + 1)]
        for mc in mcs:
            df_mc = df[df['MC'] <= mc]
            spans = [(s, e, mc) for s, e, mc in df_mc[['Start', 'End', 'MC']].values]
            cov = pt.calculate_span_coverage(spans, input_options.sequence_length)
            protein_cov_at_mcs.append(sum(cov) / len(cov) * 100)

        # calculate protein coverage at different peptide lengths
        protein_cov_at_lens = []
        lens = [l for l in range(input_options.min_len, input_options.max_len + 1)]
        for l in lens:
            df_len = df[df['Len'] <= l]
            spans = [(s, e, mc) for s, e, mc in df_len[['Start', 'End', 'MC']].values]
            cov = pt.calculate_span_coverage(spans, input_options.sequence_length)
            protein_cov_at_lens.append(sum(cov) / len(cov) * 100)

        # calculate protein coverage at different peptide Mass
        protein_cov_at_mass = []
        masses = [m for m in range(int(input_options.min_mass), int(input_options.max_mass) + 1, 100)]
        for m in masses:
            df_mass = df[df['NeutralMass'] <= m]
            spans = [(s, e, mc) for s, e, mc in df_mass[['Start', 'End', 'MC']].values]
            cov = pt.calculate_span_coverage(spans, input_options.sequence_length)
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

            df['Frag_Ions'] = [make_clickable(peptide, input_options.mass_type) for peptide in df['Sequence']]
                                        
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
            c1.metric('Cleavage Sites', len(input_options.cleavage_sites))

            protein_cov_arr_bin = pt.calculate_span_coverage(spans, input_options.sequence_length, accumulate=False)
            protein_cov_perc = round(sum(protein_cov_arr_bin) / len(protein_cov_arr_bin) * 100, 2)
            c2.metric('Protein Coverage', f'{protein_cov_perc}%')

            st.subheader('Site Indexes')
            site_indexes_html = get_site_index_html(input_options)
            st.markdown(site_indexes_html, unsafe_allow_html=True)
            st.write("")
            
            st.subheader('Sites')
            sequence_with_sites = get_sequence_site_html(input_options)
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
            site_regex = stp.text_input('Motifs Regex', 'N[ST][^P]')

            if site_regex:
                site_ranges = list(pt.get_regex_match_range(input_options.stripped_protein_sequence, site_regex))
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
                cov_site_mat = pt.calculate_span_coverage(tmp_spans, input_options.sequence_length, accumulate=True)
                
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
                protein_coverage = coverage_string(cov_site_mat, input_options.stripped_protein_sequence, cmap, sites=list(indexes_to_keep))
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

        if page_loc and "origin" in page_loc:
            url_origin = page_loc["origin"]
            if button_c.button("Generate TinyURL", key="generate_tinyurl", type="primary"):
                url_params = {k: st.query_params.get_all(k) for k in st.query_params.keys()}
                page_url = f"{url_origin}{get_query_params_url(url_params)}"
                short_url = shorten_url(page_url)

                @st.dialog(title="Share your results")
                def url_dialog(url):
                    st.write(f"Shortened URL: {url}")

                url_dialog(short_url)

