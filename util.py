from collections import Counter
from typing import List

import pandas as pd
from peptacular.mass import calculate_mass
from peptacular.sequence import span_to_sequence, calculate_sequence_length, apply_static_modifications, \
    strip_modifications
from peptacular.spans import build_enzymatic_spans, build_semi_spans
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import matplotlib as mpl

from constants import LINK, CMAP


def make_clickable(sequence, mass_type):
    # target _blank to open new window
    # extract clickable text to display for your link
    link = LINK + f'?sequence={sequence}&mass_type={mass_type}'
    return f'<a target="_blank" href="{link}">{sequence}</a>'


def generate_peptide_df(sequence: str, cleavage_sites: List, missed_cleavages: int, min_len: int,
                        max_len: int, semi_enzymatic: bool, static_mods: dict, min_mass: float, max_mass: float,
                        is_mono: bool, infer_charge: bool, min_charge: int, max_charge: int, min_mz: float,
                        max_mz: float):
    cleavage_sites = sorted(cleavage_sites)
    spans = build_enzymatic_spans(calculate_sequence_length(sequence), cleavage_sites, missed_cleavages, 1, None)
    df = pd.DataFrame(spans, columns=['Start', 'End', 'MC'])
    df['Sequence'] = [span_to_sequence(sequence, span) for span in spans]
    df['Semi'] = 0

    if semi_enzymatic is True:
        semi_spans = build_semi_spans(spans, min_len, max_len)
        semi_df = pd.DataFrame(semi_spans, columns=['Start', 'End', 'MC'])
        semi_df['Sequence'] = [span_to_sequence(sequence, span) for span in semi_spans]
        semi_df['Semi'] = 1
        df = pd.concat([df, semi_df], ignore_index=True)

    df['Len'] = df['End'] - df['Start']

    df = df.sort_values(by=['Start', 'End'])
    # make bool
    df['Semi'] = df['Semi'].apply(lambda x: bool(x))

    df.drop_duplicates(inplace=True)

    df = df[(df['Len'] >= min_len) & (df['Len'] <= max_len)]

    df['Sequence'] = df['Sequence'].apply(lambda x: apply_static_modifications(x, static_mods))

    df['NeutralMass'] = [round(calculate_mass(sequence, monoisotopic=is_mono), 5) for sequence in df['Sequence']]
    df = df[(df['NeutralMass'] >= min_mass) & (df['NeutralMass'] <= max_mass)]
    df['StrippedPeptide'] = df['Sequence'].apply(strip_modifications)

    if infer_charge is True:
        # charge should be the Lysine and Arginine count + 1
        df['Charge'] = df['StrippedPeptide'].apply(lambda x: x.count('K') + x.count('R') + 1)
        df['Mz'] = df['NeutralMass'] / df['Charge'] + df['Charge'] * 1.00727646677
        df = df[(df['Charge'] >= min_charge) & (df['Charge'] <= max_charge)]
        df = df[(df['Mz'] >= min_mz) & (df['Mz'] <= max_mz)]

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