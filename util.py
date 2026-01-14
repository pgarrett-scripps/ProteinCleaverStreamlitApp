from typing import Any
from collections import Counter

import pandas as pd
import peptacular as pt
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import matplotlib as mpl
from urllib.parse import quote_plus
import requests

from input_options import InputOptions
import pickle
import numpy as np


def generate_peptide_df(options: InputOptions) -> pd.DataFrame:
    """Generate peptide dataframe with unique peptides."""
    enzymatic_spans: list[pt.Span] = []
    for protease in options.proteases:
        peptide_spans = options.protein_annotation.digest(
            enzyme=protease,
            missed_cleavages=options.missed_cleavages,
            min_len=options.min_len,
            max_len=options.max_len,
            semi=False
        )
        enzymatic_spans.extend(list(peptide_spans))
    
    enzymatic_peptides = [options.protein_annotation[span] for span in enzymatic_spans]
    enzymatic_df = pd.DataFrame(enzymatic_spans, columns=["Start", "End", "MC"])
    enzymatic_df["Sequence"] = enzymatic_peptides
    enzymatic_df["Semi"] = False
    
    if options.semi_enzymatic:
        semi_enzymatic_spans: list[pt.Span] = []
        for protease in options.proteases:
            peptide_spans = options.protein_annotation.digest(
                enzyme=protease,
                missed_cleavages=options.missed_cleavages,
                min_len=options.min_len,
                max_len=options.max_len,
                semi=True
            )
            semi_enzymatic_spans.extend(list(peptide_spans))
        
        semi_enzymatic_peptides = [options.protein_annotation[span] for span in semi_enzymatic_spans]
        semi_enzymatic_df = pd.DataFrame(semi_enzymatic_spans, columns=["Start", "End", "MC"])
        semi_enzymatic_df["Sequence"] = semi_enzymatic_peptides
        semi_enzymatic_df["Semi"] = True
        
        # Concatenate and remove duplicates
        df = pd.concat([enzymatic_df, semi_enzymatic_df], ignore_index=True)
        df = df.drop_duplicates(subset=["Sequence"], keep="first")  # Keep enzymatic version
    else:
        # Remove duplicates from multiple proteases
        df = enzymatic_df.drop_duplicates(subset=["Sequence"], keep="first")
    

    df["Len"] = df["End"] - df["Start"]

    df = df.sort_values(by=["Start", "End"])
    df.drop_duplicates(inplace=True)
    df = df[(df["Len"] >= options.min_len) & (df["Len"] <= options.max_len)]

    internal_static_mods = {k: [v] for k, v in options.static_modifications.items()}
    nterm_static_mods = {k: [v] for k, v in options.n_term_modifications.items()}
    cterm_static_mods = {k: [v] for k, v in options.c_term_modifications.items()}

    df["Sequence"] = df["Sequence"].apply(
        lambda x: pt.build_mods(x, internal_static=internal_static_mods, nterm_static=nterm_static_mods, cterm_static=cterm_static_mods)[0]
    )

    df["NeutralMass"] = pt.mass( df["Sequence"].to_list(), monoisotopic=options.is_monoisotopic, method='sequential')
    df = df[
        (df["NeutralMass"] >= options.min_mass)
        & (df["NeutralMass"] <= options.max_mass)
    ]
    df["StrippedPeptide"] = df["Sequence"].apply(pt.strip_mods)

    if options.infer_charge is True:
        # charge should be the Lysine and Arginine count + 1
        df["base_charge"] = df["StrippedPeptide"].apply(
            lambda x: x.count("K") + x.count("R") + 1
        )
        df["Charge"] = df["base_charge"].apply(
            lambda c: [
                i
                for i in range(
                    c - options.plus_minus_charge, c + options.plus_minus_charge + 1
                )
            ]
        )
        df = df.explode("Charge")
        df["Charge"] = df["Charge"].astype(int)

        df["Mz"] = (df["NeutralMass"] + df["Charge"] * pt.PROTON_MASS) / df["Charge"]
        df = df[
            (df["Charge"] >= options.min_charge) & (df["Charge"] <= options.max_charge)
        ]
        df = df[(df["Mz"] >= options.min_mz) & (df["Mz"] <= options.max_mz)]

        # drop base_charge
        df.drop(columns=["base_charge"], inplace=True)

    if options.infer_retention_time:
        rt_model = pickle.load(open("models/rt_model.pkl", "rb"))
        df["RT"] = rt_model.predict(
            np.array([bin_aa_counts(pt.strip_mods(seq)) for seq in df["Sequence"]])
        )
        df["RT"] = df["RT"].round(3)

        # scale to retention time
        df["RT"] = df["RT"] * options.retention_time

        if options.filter_invalid_rt:
            df = df[df["RT"] <= options.retention_time]
            # greater than 0
            df = df[df["RT"] > 0]

    if "Charge" in df.columns:
        im_model = pickle.load(open("models/im_model.pkl", "rb"))
        df["IM"] = im_model.predict(
            np.array(
                [
                    bin_aa_counts(pt.strip_mods(seq), c)
                    for seq, c in df[["Sequence", "Charge"]].values
                ]
            )
        )
        df["IM"] = df["IM"].round(3)

    if options.infer_proteotypic:
        proteotypic_model = pickle.load(open("models/proteotypic_model.pkl", "rb"))

        probs = proteotypic_model.predict_proba(
            np.array([bin_aa_counts(pt.strip_mods(seq)) for seq in df["Sequence"]])
        )

        df["score"] = probs[:, 1]
        df["Proteotypic"] = probs[:, 1] >= options.score_threshold

        if options.remove_non_proteotypic:
            df = df[df["Proteotypic"]]

    return df


def bin_aa_counts(seq: str, charge: int | None = None):
    amino_acids = "ACDEFGHIKLMNPQRSTVWY"
    aa_counts = Counter(seq)

    enc = [aa_counts.get(aa, 0) for aa in amino_acids]
    if charge is not None:
        enc.append(charge)
    return enc


def coverage_string(protein_cov_arr: list[int], stripped_protein_sequence: str, cmap: Any, sites: list[int] | None = None):
    # Find the maximum coverage value
    max_coverage = max(protein_cov_arr)

    # Color all covered amino acids based on coverage and show index on hover, using a monospace font
    protein_coverage = (
        '<span style="font-family: Courier New, monospace; font-size: 16px;">'
    )
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
    protein_coverage += "</span>"

    return protein_coverage


def create_colorbar(max_coverage: int, cmap: Any, label: str ="Coverage"):
    # Create a figure and a subplot for the colorbar
    fig, ax = plt.subplots(figsize=(15, 1))
    fig.subplots_adjust(bottom=0.5)

    # Create a normalization from 0 to max_coverage
    norm = mpl.colors.Normalize(vmin=0, vmax=max_coverage)

    # Create a colorbar in the subplot
    cb = mpl.colorbar.ColorbarBase(ax, cmap=cmap, norm=norm, orientation="horizontal")
    cb.set_label(label)

    # Display the colorbar in Streamlit
    return fig


def listify(o: Any = None) -> list[Any]:
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


def get_query_params_url(params_dict: dict[str, Any], **kwargs: Any) -> str:
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


def get_site_index_html(input_options: InputOptions) -> str:
    site_indexes_html = (
        '<span style="font-family: Courier New, monospace; font-size: 16px;">'
    )
    for cleavage_site in input_options.cleavage_sites:
        site_indexes_html += (
            f'<span style="background-color:#f0f0f0; font-weight:900; color:red; padding:2px; '
            f'margin:1px; border:1px solid #ffcc00; border-radius:3px;">{cleavage_site}</span>'
        )
        site_indexes_html += " "
    site_indexes_html += "</span>"
    return site_indexes_html


def get_sequence_site_html(input_options: InputOptions):
    sequence_with_sites = (
        '<span style="font-family: Courier New, monospace; font-size: 16px;">'
    )
    for i, aa in enumerate(input_options.stripped_protein_sequence):
        # Add the amino acid with its original index
        sequence_with_sites += (
            f'<span title="Index: {i + 1}" style="background-color:#f0f0f0; font-weight:900; '
            f"color:#333; padding:2px; margin:1px; border:1px solid #cccccc; "
            f'border-radius:3px;">{aa}</span>'
        )

        # Check if the next position is a cleavage site and insert '%' character
        if i + 1 in input_options.cleavage_sites:
            # Highlight '%' character in blue
            sequence_with_sites += (
                f'<span style="background-color:#e0e0ff; font-weight:900; color:red; font-weight:bold;'
                f' padding:2px; margin:1px; border:1px solid #a0a0ff; border-radius:3px;">%</span>'
            )
    sequence_with_sites += "</span>"

    return sequence_with_sites
