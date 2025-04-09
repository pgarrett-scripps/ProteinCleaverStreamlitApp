from collections import Counter
from typing import List

import pandas as pd
import peptacular as pt
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import matplotlib as mpl
from urllib.parse import quote_plus
import requests

from constants import LINK
from input_options import InputOptions
import pickle
import numpy as np

def make_clickable(sequence, mass_type):
    # target _blank to open new window
    # extract clickable text to display for your link
    link = LINK + f'?peptide={sequence}&mass_type={mass_type}'
    return link


def generate_peptide_df(options: InputOptions) -> pd.DataFrame:
    spans = list(pt.build_enzymatic_spans(max_index=options.sequence_length, 
                                          enzyme_sites=options.cleavage_sites, 
                                          missed_cleavages=options.missed_cleavages, 
                                          min_len=1, 
                                          max_len=None))
    
    df = pd.DataFrame(spans, columns=['Start', 'End', 'MC'])
    df['Sequence'] = [pt.span_to_sequence(options.protein_annotation, span) for span in spans]
    df['Semi'] = 0

    if options.semi_enzymatic is True:
        semi_spans = list(pt.build_semi_spans(spans, options.min_len, options.max_len))
        semi_df = pd.DataFrame(semi_spans, columns=['Start', 'End', 'MC'])
        semi_df['Sequence'] = [pt.span_to_sequence(options.protein_annotation, span) for span in semi_spans]
        semi_df['Semi'] = 1
        df = pd.concat([df, semi_df], ignore_index=True)

    df['Len'] = df['End'] - df['Start']

    df = df.sort_values(by=['Start', 'End'])
    # make bool
    df['Semi'] = df['Semi'].apply(lambda x: bool(x))

    df.drop_duplicates(inplace=True)

    df = df[(df['Len'] >= options.min_len) & (df['Len'] <= options.max_len)]

    df['Sequence'] = df['Sequence'].apply(lambda x: pt.apply_static_mods(x, options.static_modifications))

    df['NeutralMass'] = [round(pt.mass(sequence, monoisotopic=options.is_monoisotopic), 5) for sequence in df['Sequence']]
    df = df[(df['NeutralMass'] >= options.min_mass) & (df['NeutralMass'] <= options.max_mass)]
    df['StrippedPeptide'] = df['Sequence'].apply(pt.strip_mods)

    if options.infer_charge is True:
        # charge should be the Lysine and Arginine count + 1
        df['base_charge'] = df['StrippedPeptide'].apply(lambda x: x.count('K') + x.count('R') + 1)
        df['Charge'] = df['base_charge'].apply(lambda c: [i for i in range(c-options.plus_minus_charge, c+options.plus_minus_charge + 1)])
        df = df.explode('Charge')
        df['Charge'] = df['Charge'].astype(int)

        df['Mz'] = (df['NeutralMass'] + df['Charge'] * pt.PROTON_MASS) / df['Charge']
        df = df[(df['Charge'] >= options.min_charge) & (df['Charge'] <= options.max_charge)]
        df = df[(df['Mz'] >= options.min_mz) & (df['Mz'] <= options.max_mz)]

        #drop base_charge
        df.drop(columns=['base_charge'], inplace=True)


    if options.infer_retention_time:
        rt_model = pickle.load(open("rt_model.pkl", "rb"))
        df['RT'] = rt_model.predict(np.array([bin_aa_counts(pt.strip_mods(seq)) for seq in df['Sequence']]))
        df['RT'] = df['RT'].round(3)

        #scale to retention time
        df['RT'] = df['RT'] * options.retention_time 

        if options.filter_invalid_rt:
            df = df[df['RT'] <= options.retention_time]
            # greater than 0
            df = df[df['RT'] > 0]

    if 'Charge' in df.columns:
        im_model = pickle.load(open("im_model.pkl", "rb"))
        df['IM'] = im_model.predict(
            np.array([bin_aa_counts(pt.strip_mods(seq), c) for seq, c in df[['Sequence', 'Charge']].values]))
        df['IM'] = df['IM'].round(3)

    if options.infer_proteotypic:
        proteotypic_model = pickle.load(open("proteotypic_model.pkl", "rb"))

        probs = proteotypic_model.predict_proba(
            np.array([bin_aa_counts(pt.strip_mods(seq)) for seq in df['Sequence']]))
        
        df['score'] = probs[:, 1]
        df['Proteotypic'] = probs[:, 1] >= options.score_threshold

        if options.remove_non_proteotypic:
            df = df[df['Proteotypic']]


    return df


def bin_aa_counts(seq, charge=None):
    amino_acids = 'ACDEFGHIKLMNPQRSTVWY'
    aa_counts = Counter(seq)

    enc = [aa_counts.get(aa, 0) for aa in amino_acids]
    if charge is not None:
        enc.append(charge)
    return enc


def coverage_string(protein_cov_arr, stripped_protein_sequence, cmap, sites=None):
    # Find the maximum coverage value
    max_coverage = max(protein_cov_arr)

    # Color all covered amino acids based on coverage and show index on hover, using a monospace font
    protein_coverage = '<span style="font-family: Courier New, monospace; font-size: 16px;">'
    for i, aa in enumerate(stripped_protein_sequence):
        coverage = protein_cov_arr[i]
        if coverage > 0 and max_coverage > 0:
            # Normalize the coverage based on the maximum value
            normalized_coverage = coverage / max_coverage

            # Get color from colormap
            color = cmap(normalized_coverage)
            hex_color = mcolors.to_hex(color)

            protein_coverage += f'<span title="Index: {i + 1}; Coverage: {coverage}" style="background-color:#e0e0ff; color:{hex_color}; font-weight:900; padding:3px; margin:1px; border:1px solid #a0a0ff; border-radius:3px;">{aa}</span>'
        else:
            if sites and i in sites:
                protein_coverage += f'<span title="Index: {i + 1}" style="background-color:red; color:#333; font-weight:900; padding:3px; margin:1px; border:1px solid #cccccc; border-radius:3px;">{aa}</span>'
            else:
                # Style the non-covered amino acid for a polished look with a tooltip
                protein_coverage += f'<span title="Index: {i + 1}" style="background-color:#f0f0f0; color:#333; font-weight:900; padding:3px; margin:1px; border:1px solid #cccccc; border-radius:3px;">{aa}</span>'
    protein_coverage += '</span>'

    return protein_coverage

def create_colorbar(max_coverage, cmap, label='Coverage'):
    # Create a figure and a subplot for the colorbar
    fig, ax = plt.subplots(figsize=(15, 1))
    fig.subplots_adjust(bottom=0.5)

    # Create a normalization from 0 to max_coverage
    norm = mpl.colors.Normalize(vmin=0, vmax=max_coverage)

    # Create a colorbar in the subplot
    cb = mpl.colorbar.ColorbarBase(ax, cmap=cmap, norm=norm, orientation='horizontal')
    cb.set_label(label)

    # Display the colorbar in Streamlit
    return fig


def listify(o=None):
    if o is None:
        res = []
    elif isinstance(o, list):
        res = o
    elif isinstance(o, str):
        res = [o]
    else:
        res = [o]
    return res


def shorten_url(url: str) -> str:
    """Shorten a URL using TinyURL."""
    api_url = f"http://tinyurl.com/api-create.php?url={url}"

    try:
        response = requests.get(api_url)
        response.raise_for_status()
        return response.text
    except requests.RequestException as e:
        return f"Error: {e}"


def get_query_params_url(params_dict):
    """
    Create url params from alist of parameters and a dictionary with values.

    Args:
        params_list (str) :
            A list of parameters to get the value of from `params_dict`
        parmas_dict (dict) :
            A dict with values for the `parmas_list .
        **kwargs :
            Extra keyword args to add to the url
    """
    return "?" + "&".join(
        [
            f"{key}={quote_plus(str(value))}"
            for key, values in params_dict.items()
            for value in listify(values)
        ]
    )


def get_site_index_html(input_options):
    site_indexes_html = '<span style="font-family: Courier New, monospace; font-size: 16px;">'
    for cleavage_site in input_options.cleavage_sites:
        site_indexes_html += f'<span style="background-color:#f0f0f0; font-weight:900; color:red; padding:2px; ' \
                            f'margin:1px; border:1px solid #ffcc00; border-radius:3px;">{cleavage_site}</span>'
        site_indexes_html += ' '
    site_indexes_html += '</span>'
    return site_indexes_html

def get_sequence_site_html(input_options):

    sequence_with_sites = '<span style="font-family: Courier New, monospace; font-size: 16px;">'
    for i, aa in enumerate(input_options.stripped_protein_sequence):
        # Add the amino acid with its original index
        sequence_with_sites += f'<span title="Index: {i + 1}" style="background-color:#f0f0f0; font-weight:900; ' \
                            f'color:#333; padding:2px; margin:1px; border:1px solid #cccccc; ' \
                            f'border-radius:3px;">{aa}</span>'

        # Check if the next position is a cleavage site and insert '%' character
        if i + 1 in input_options.cleavage_sites:
            # Highlight '%' character in blue
            sequence_with_sites += f'<span style="background-color:#e0e0ff; font-weight:900; color:red; font-weight:bold;' \
                                f' padding:2px; margin:1px; border:1px solid #a0a0ff; border-radius:3px;">%</span>'
    sequence_with_sites += '</span>'

    return sequence_with_sites


