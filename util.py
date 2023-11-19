from typing import List

import pandas as pd
from peptacular.mass import calculate_mass
from peptacular.sequence import span_to_sequence, calculate_sequence_length, apply_static_modifications, \
    strip_modifications
from peptacular.spans import build_enzymatic_spans, build_semi_spans

from constants import LINK


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
    df['PeptideSequence'] = [span_to_sequence(sequence, span) for span in spans]
    df['IsSemi'] = 0

    if semi_enzymatic is True:
        semi_spans = build_semi_spans(spans, min_len, max_len)
        semi_df = pd.DataFrame(semi_spans, columns=['Start', 'End', 'MC'])
        semi_df['PeptideSequence'] = [span_to_sequence(sequence, span) for span in semi_spans]
        semi_df['IsSemi'] = 1
        df = pd.concat([df, semi_df], ignore_index=True)

    df['AACount'] = df['End'] - df['Start']

    df = df.sort_values(by=['Start', 'End'])
    # make bool
    df['IsSemi'] = df['IsSemi'].apply(lambda x: bool(x))

    df.drop_duplicates(inplace=True)

    df = df[(df['AACount'] >= min_len) & (df['AACount'] <= max_len)]

    df['PeptideSequence'] = df['PeptideSequence'].apply(lambda x: apply_static_modifications(x, static_mods))

    df['Mass'] = [round(calculate_mass(sequence, monoisotopic=is_mono), 5) for sequence in df['PeptideSequence']]
    df = df[(df['Mass'] >= min_mass) & (df['Mass'] <= max_mass)]
    df['StrippedPeptide'] = df['PeptideSequence'].apply(strip_modifications)

    if infer_charge is True:
        # charge should be the Lysine and Arginine count + 1
        df['Charge'] = df['StrippedPeptide'].apply(lambda x: x.count('K') + x.count('R') + 1)
        df['Mz'] = df['Mass'] / df['Charge'] + df['Charge'] * 1.00727646677
        df = df[(df['Charge'] >= min_charge) & (df['Charge'] <= max_charge)]
        df = df[(df['Mz'] >= min_mz) & (df['Mz'] <= max_mz)]

    return df
