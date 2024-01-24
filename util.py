import pickle
from collections import Counter
from typing import List

import numpy as np
import pandas as pd
import requests
from peptacular.mass import calculate_mass
from peptacular.sequence import span_to_sequence, calculate_sequence_length, apply_static_modifications, \
    strip_modifications, apply_variable_modifications
from peptacular.spans import build_enzymatic_spans, build_semi_spans
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import matplotlib as mpl
from peptacular.term.modification import add_n_term_modification, add_c_term_modification

from constants import LINK, BASE_URL


def fetch_sequence_from_uniprot(accession_number):
    url = f"https://www.uniprot.org/uniprot/{accession_number}.fasta"
    response = requests.get(url)
    return response


def generate_peptide_df(sequence: str, cleavage_sites: List, missed_cleavages: int, min_len: int,
                        max_len: int, semi_enzymatic: bool, static_mods: dict, min_mass: float, max_mass: float,
                        is_mono: bool, infer_charge: bool, min_charge: int, max_charge: int, min_mz: float,
                        max_mz: float, var_mods: dict, max_var_mods: int, n_term_static_mod: float,
                        c_term_static_mod: float, n_term_var_mod: float, c_term_var_mod: float,
                        remove_non_proteotypic: bool):
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

    # Apply variable modifications to each sequence in the DataFrame
    def apply_var_mods(sequence: str) -> str:

        var_seqs = apply_variable_modifications(sequence, var_mods, max_var_mods)

        if n_term_var_mod and c_term_var_mod:
            n_term_seq = add_n_term_modification(sequence, n_term_var_mod)
            c_term_seq = add_c_term_modification(sequence, c_term_var_mod)
            n_c_term_seq = add_c_term_modification(n_term_seq, c_term_var_mod)

            n_term_seqs = apply_variable_modifications(n_term_seq, var_mods, max_var_mods)
            c_term_seqs = apply_variable_modifications(c_term_seq, var_mods, max_var_mods)
            n_c_term_seqs = apply_variable_modifications(n_c_term_seq, var_mods, max_var_mods)

            return ';'.join(list(set(var_seqs + n_term_seqs + c_term_seqs + n_c_term_seqs)))

        elif n_term_var_mod:
            n_term_seq = add_n_term_modification(sequence, n_term_var_mod)
            n_term_seqs = apply_variable_modifications(n_term_seq, var_mods, max_var_mods)
            return ';'.join(list(set(var_seqs + n_term_seqs)))
        elif c_term_var_mod:
            c_term_seq = add_c_term_modification(sequence, c_term_var_mod)
            c_term_seqs = apply_variable_modifications(c_term_seq, var_mods, max_var_mods)
            return ';'.join(list(set(var_seqs + c_term_seqs)))
        else:
            return ';'.join(var_seqs)

    # Apply the apply_var_mods function to the 'Sequence' column
    df['Sequence'] = df['Sequence'].apply(apply_var_mods)

    # expand the sequence column into multiple rows (sequences are separated by ';')
    df = df.assign(Sequence=df.Sequence.str.split(';')).explode('Sequence')

    def apply_static_mods(sequence: str) -> str:

        sequence = apply_static_modifications(sequence, static_mods)

        if n_term_static_mod:
            sequence = add_n_term_modification(sequence, n_term_static_mod)

        if c_term_static_mod:
            sequence = add_c_term_modification(sequence, c_term_static_mod)

        return sequence

    df['Sequence'] = df['Sequence'].apply(apply_static_mods)

    # drop duplicates

    df['NeutralMass'] = [round(calculate_mass(sequence, monoisotopic=is_mono), 5) for sequence in df['Sequence']]
    df = df[(df['NeutralMass'] >= min_mass) & (df['NeutralMass'] <= max_mass)]
    df['StrippedPeptide'] = df['Sequence'].apply(strip_modifications)

    df.sort_values(by=['MC'], inplace=True)
    df.drop_duplicates(subset=['Start', 'Sequence'], inplace=True)
    df.sort_values(by=['Start', 'Len'], inplace=True)

    if infer_charge is True:
        # charge should be the Lysine and Arginine count + 1
        df['Charge'] = df['StrippedPeptide'].apply(lambda x: x.count('K') + x.count('R') + 1)
        df['Mz'] = df['NeutralMass'] / df['Charge'] + df['Charge'] * 1.00727646677
        df = df[(df['Charge'] >= min_charge) & (df['Charge'] <= max_charge)]
        df = df[(df['Mz'] >= min_mz) & (df['Mz'] <= max_mz)]

    if len(df) == 0:
        return df

    rt_model = pickle.load(open("rt_model.pkl", "rb"))
    df['RT'] = rt_model.predict(np.array([bin_aa_counts(strip_modifications(seq)) for seq in df['Sequence']]))
    df['RT'] = df['RT'].round(3)

    if infer_charge:
        im_model = pickle.load(open("im_model.pkl", "rb"))
        df['IM'] = im_model.predict(
            np.array([bin_aa_counts(strip_modifications(seq), c) for seq, c in df[['Sequence', 'Charge']].values]))
        df['IM'] = df['IM'].round(3)

    proteotypic_model = pickle.load(open("proteotypic_model.pkl", "rb"))
    df['Proteotypic'] = proteotypic_model.predict(
        np.array([bin_aa_counts(strip_modifications(seq)) for seq in df['Sequence']])).astype(bool)

    if remove_non_proteotypic:
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


def generate_app_url(protein_id, protein_sequence, proteases, custom_regex, missed_cleavages, mass_type,
                     min_peptide_len, max_peptide_len, min_mass, max_mass, semi_enzymatic, infer_charge,
                     min_charge, max_charge, min_mz, max_mz, remove_non_proteotypic, n_term_static_mod,
                     c_term_static_mod, num_static_mods, n_term_var_mod, c_term_var_mod, max_var_mods,
                     num_variable_mods, static_mods, variable_mods):

    # flip the dictionary
    static_mods_rev = {}
    for aa, mod in static_mods.items():
        static_mods_rev.setdefault(mod, []).append(aa)
    static_mod_str = ';'.join([f'{"".join(aa)}:{mod}' for mod, aa in static_mods_rev.items()])

    variable_mods_rev = {}
    for aa, mod in variable_mods.items():
        variable_mods_rev.setdefault(mod, []).append(aa)
    variable_mod_str = ';'.join([f'{"".join(aa)}:{mod}' for mod, aa in variable_mods_rev.items()])

    params = {
        'protein_id': protein_id,
        'protein_sequence': protein_sequence,
        'proteases': ';'.join(proteases),
        'custom_regex': custom_regex,
        'missed_cleavages': missed_cleavages,
        'mass_type': mass_type,
        'min_peptide_len': min_peptide_len,
        'max_peptide_len': max_peptide_len,
        'min_mass': min_mass,
        'max_mass': max_mass,
        'semi_enzymatic': semi_enzymatic,
        'infer_charge': infer_charge,
        'min_charge': min_charge,
        'max_charge': max_charge,
        'min_mz': min_mz,
        'max_mz': max_mz,
        'remove_non_proteotypic': remove_non_proteotypic,
        'n_term_static_mod': n_term_static_mod,
        'c_term_static_mod': c_term_static_mod,
        'num_static_mods': num_static_mods,
        'n_term_var_mod': n_term_var_mod,
        'c_term_var_mod': c_term_var_mod,
        'max_var_mods': max_var_mods,
        'num_variable_mods': num_variable_mods,
        'static_mods': static_mod_str,
        'variable_mods': variable_mod_str
    }
    query_string = '&'.join([f'{key}={value}' for key, value in params.items() if value is not None])
    return f'{BASE_URL}?{query_string}'


def make_clickable(sequence, mass_type):
    # target _blank to open new window
    # extract clickable text to display for your link
    link = LINK + f'?sequence={sequence}&mass_type={mass_type}'
    return link
