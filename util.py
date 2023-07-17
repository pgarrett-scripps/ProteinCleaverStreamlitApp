from typing import List

import pandas as pd
from peptacular.sequence import identify_cleavage_sites
from peptacular.spans import get_enzymatic_spans, get_semi_spans

from constants import LINK


def make_clickable(sequence):
    # target _blank to open new window
    # extract clickable text to display for your link
    link = LINK + f'?sequence={sequence}'
    return f'<a target="_blank" href="{link}">{sequence}</a>'


def generate_peptide_df(sequence: str, cleavage_sites: List, missed_cleavages: int, min_len: int,
                        max_len: int, semi_enzymatic: bool):
    cleavage_sites = sorted(cleavage_sites)
    spans = get_enzymatic_spans(len(sequence), cleavage_sites, missed_cleavages, None, None)
    df = pd.DataFrame(spans, columns=['Start', 'End', 'MC'])
    df['Sequence'] = df.apply(lambda x: sequence[x['Start']:x['End']], axis=1)
    df['Semi'] = 0
    df = df[(df['Sequence'].str.len() >= min_len) & (df['Sequence'].str.len() <= max_len)]

    if semi_enzymatic is True:
        semi_spans = get_semi_spans(spans, min_len, max_len)
        semi_df = pd.DataFrame(semi_spans, columns=['Start', 'End', 'MC'])
        semi_df['Sequence'] = semi_df.apply(lambda x: sequence[x['Start']:x['End']], axis=1)
        semi_df['Semi'] = 1
        df = pd.concat([df, semi_df], ignore_index=True)

    df = df.sort_values(by=['Start', 'End'])
    df['Len'] = df['Sequence'].str.len()

    return df
