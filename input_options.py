from dataclasses import dataclass
import streamlit as st
import streamlit_permalink as stp
import constants as c

import help_messages as h
from typing import Optional

import pandas as pd
import peptacular as pt
import requests
from functools import cached_property


def add_dumb_mobile_buffer():
    # blank space for mobile cause streamlit is weird
    st.markdown("<br>", unsafe_allow_html=True)
    st.markdown("<br>", unsafe_allow_html=True)
    st.markdown("<br>", unsafe_allow_html=True)
    st.markdown("<br>", unsafe_allow_html=True)
    st.markdown("<br>", unsafe_allow_html=True)
    st.markdown("<br>", unsafe_allow_html=True)


def fetch_sequence_from_uniprot(accession_number: str) -> Optional[str]:
    url = f"https://www.uniprot.org/uniprot/{accession_number}.fasta"
    response = requests.get(url)
    if response.status_code != 200:
        st.error(f"Error fetching sequence from UniProt: {response.status_code}")
        return None
    return "".join(response.text.split("\n")[1:])  # Remove the header line


@dataclass
class InputOptions:
    protein_id: str
    protein_sequence: str
    proteases_selected: list[str]
    custom_regex: str
    missed_cleavages: int
    mass_type: str
    min_len: int
    max_len: int
    min_mass: float
    max_mass: float
    semi_enzymatic: bool
    infer_charge: bool
    min_charge: Optional[int]
    max_charge: Optional[int]
    min_mz: Optional[float]
    max_mz: Optional[float]
    plus_minus_charge: Optional[int]
    infer_retention_time: bool
    retention_time: Optional[float]
    filter_invalid_rt: Optional[bool]
    infer_ion_mobility: bool
    infer_proteotypic: bool
    remove_non_proteotypic: Optional[bool]
    score_threshold: Optional[float]
    static_modifications: dict[str, float]
    n_term_modifications: dict[str, float]
    c_term_modifications: dict[str, float]

    @cached_property
    def proteases(self):
        proteases: list[str] = []
        for protease in self.proteases_selected:
            proteases.append(protease)
        if self.custom_regex:
            proteases.append(self.custom_regex)
        return proteases

    @property
    def is_monoisotopic(self):
        return self.mass_type == "monoisotopic"

    @property
    def is_average(self):
        return self.mass_type == "average"

    @cached_property
    def sequence_length(self):
        return pt.sequence_length(self.protein_sequence)

    @cached_property
    def stripped_protein_sequence(self):
        return pt.strip_mods(self.protein_sequence)

    @cached_property
    def neutral_mass(self):
        return pt.mass(self.protein_sequence, monoisotopic=self.is_monoisotopic)

    @cached_property
    def cleavage_sites(self):
        sites: set[int] = set()
        for enzyme_regex in self.proteases:
            sites.update(pt.cleavage_sites(self.protein_sequence, enzyme_regex))
        return sites

    @cached_property
    def protein_annotation(self):
        return pt.parse(self.protein_sequence)


def validate_input_option(options: InputOptions) -> bool:
    if options.min_len > options.max_len:
        st.error(
            f"Min length: {options.min_len} must be less than max length: {options.max_len}"
        )
        return False

    if (
        options.min_charge is not None
        and options.max_charge is not None
        and options.min_charge > options.max_charge
    ):
        st.error(
            f"Min charge: {options.min_charge} must be less than max charge: {options.max_charge}"
        )
        return False

    if (
        options.min_mz is not None
        and options.max_mz is not None
        and options.min_mz > options.max_mz
    ):
        st.error(
            f"Min m/z: {options.min_mz} must be less than max m/z: {options.max_mz}"
        )
        return False

    try:
        if options.sequence_length > c.MAX_PROTEIN_INPUT_LENGTH:
            st.error(
                f"Protein sequence length: {options.sequence_length} exceeds max length: {c.MAX_PROTEIN_INPUT_LENGTH}"
            )
            return False
    except Exception as e:
        st.error(f"Error parsing protein sequence: {e}")
        return False

    if options.min_mass > options.max_mass:
        st.error(
            f"Min mass: {options.min_mass} must be less than max mass: {options.max_mass}"
        )
        return False

    try:
        _ = options.neutral_mass
    except Exception as e:
        st.error(f"Error Calculating Protein Mass: {e}")
        return False

    return True


def get_input_options():
    st.subheader("Options", divider="grey")

    c1, c2 = st.columns(2, vertical_alignment="bottom")
    with c1:
        protein_id = st.text_input(
            label="Protein accession number / identifier",
            value=c.DEFAULT_PROTEIN_ID,
            placeholder=c.DEFAULT_PROTEIN_ID_PLACEHOLDER,
            help=h.PROTEIN_ID_HELP,
            max_chars=20,
            key="protein_id",
        )

    with c2:
        # btn to fetch sequence from uniprot
        fetch_sequence = st.button(
            label="Fetch from UniProt",
            help="Fetch protein sequence from UniProt using the accession number",
            use_container_width=True,
            key="fetch_sequence",
        )

    if fetch_sequence:
        if protein_id:
            fetched_protein_sequence = fetch_sequence_from_uniprot(protein_id)
            if fetched_protein_sequence is not None:
                st.query_params["protein_sequence"] = fetched_protein_sequence
                st.rerun()
        else:
            st.error("Please enter a protein accession number")

    protein_sequence = stp.text_area(
        label="Protein Sequence (Proforma 2.0)",
        value=c.DEFAULT_PROTEIN_SEQUENCE,
        help=h.PROTEIN_SEQUENCE_HELP,
        key="protein_sequence",
        max_chars=c.MAX_PROTEIN_INPUT_LENGTH,
        height=200,
    )

    # remove all new lines / spaces from the sequence
    protein_sequence = protein_sequence.replace("\n", "").replace(" ", "")

    c1, c2 = st.columns(2, vertical_alignment="bottom")

    with c1:
        proteases_selected = stp.multiselect(
            label="Proteases",
            options=c.VALID_PROTEASES.keys(),
            help=h.PROTEASES_HELP,
            default=c.DEFAULT_PROTEASES,
            key="protease",
        )

    with c2:
        custom_regex = stp.text_input(
            label="(Additional) Custom protease",
            value=c.DEFAULT_CUSTOM_REGEX,
            help=h.CUSTOM_REGEX_HELP,
            key="custom_regex",
        )

    c1, c2 = st.columns(2, vertical_alignment="bottom")

    with c1:
        missed_cleavages = stp.number_input(
            label="Max missed cleavages",
            min_value=c.MIN_MISSED_CLEAVAGES,
            max_value=c.MAX_MISSED_CLEAVAGES,
            value=c.DEFAULT_MISSED_CLEAVAGES,
            step=c.MISSED_CLEAVAGES_STEP,
            help=h.MISSED_CLEAVAGES_HELP,
            key="missed_cleavages",
        )

    with c2:
        mass_type = stp.selectbox(
            label="Mass type",
            options=c.MASS_TYPE_OPTIONS,
            index=c.DEFAULT_MASS_TYPE_INDEX,
            help=h.MASS_TYPE_HELP,
            key="mass_type",
        )

    c1, c2 = st.columns(2, vertical_alignment="bottom")

    with c1:
        min_len = stp.number_input(
            label="Min length",
            min_value=c.MIN_PEPTIDE_LEN,
            max_value=c.MAX_PEPTIDE_LEN,
            value=c.DEFAULT_MIN_PEPTIDE_LEN,
            step=c.PEPTIDE_LEN_STEP,
            help=h.MIN_PEPTIDE_LEN_HELP,
            key="min_peptide_len",
        )

    with c2:
        max_len = stp.number_input(
            label="Max length",
            min_value=c.MIN_PEPTIDE_LEN,
            max_value=c.MAX_PEPTIDE_LEN,
            value=c.DEFAULT_MAX_PEPTIDE_LEN,
            step=c.PEPTIDE_LEN_STEP,
            help=h.MAX_PEPTIDE_LEN_HELP,
            key="max_peptide_len",
        )

    c1, c2 = st.columns(2, vertical_alignment="bottom")

    with c1:
        min_mass = stp.number_input(
            label="Min neutral mass",
            min_value=c.MIN_PEPTIDE_MASS,
            max_value=c.MAX_PEPTIDE_MASS,
            value=c.DEFAULT_MIN_PEPTIDE_MASS,
            step=c.PEPTIDE_MASS_STEP,
            help=h.MIN_PEPTIDE_MASS_HELP,
            key="min_peptide_mass",
        )

    with c2:
        max_mass = stp.number_input(
            label="Max neutral mass",
            min_value=c.MIN_PEPTIDE_MASS,
            max_value=c.MAX_PEPTIDE_MASS,
            value=c.DEFAULT_MAX_PEPTIDE_MASS,
            step=c.PEPTIDE_MASS_STEP,
            help=h.MAX_PEPTIDE_MASS_HELP,
            key="max_peptide_mass",
        )

    semi_enzymatic = stp.checkbox(
        label="Semi enzymatic?",
        help=h.SEMI_ENZYMATIC_HELP,
        value=c.DEFAULT_SEMI_ENZYMATIC,
        key="semi_enzymatic",
    )

    with st.expander("Peptide Properties", expanded=False):
        st.subheader("Charge", divider="grey")
        infer_charge = stp.checkbox(
            label="Infer Charge",
            value=c.DEFAULT_INFER_CHARGE,
            help=h.INFER_CHARGE_HELP,
            key="infer_charge",
        )
        min_charge = max_charge = min_mz = max_mz = plus_minus_charge = None
        if infer_charge:
            c1, c2 = st.columns(2, vertical_alignment="bottom")
            with c1:
                min_charge = stp.number_input(
                    label="Min charge",
                    min_value=c.MIN_CHARGE,
                    max_value=c.MAX_CHARGE,
                    value=c.DEFAULT_MIN_CHARGE,
                    step=c.CHARGE_STEP,
                    help=h.MIN_CHARGE_HELP,
                    key="min_charge",
                )

                min_mz = stp.number_input(
                    label="Min m/z",
                    min_value=c.MIN_MZ,
                    max_value=c.MAX_MZ,
                    value=c.DEFAULT_MIN_MZ,
                    step=c.MZ_STEP,
                    help=h.MIN_MZ_HELP,
                    key="min_mz",
                )

            with c2:
                max_charge = stp.number_input(
                    label="Max charge",
                    min_value=c.MIN_CHARGE,
                    max_value=c.MAX_CHARGE,
                    value=c.DEFAULT_MAX_CHARGE,
                    step=c.CHARGE_STEP,
                    help=h.MAX_CHARGE_HELP,
                    key="max_charge",
                )

                max_mz = stp.number_input(
                    label="Max m/z",
                    min_value=c.MIN_MZ,
                    max_value=c.MAX_MZ,
                    value=c.DEFAULT_MAX_MZ,
                    step=c.MZ_STEP,
                    help=h.MAX_MZ_HELP,
                    key="max_mz",
                )
            plus_minus_charge = stp.number_input(
                label="Charge +/-",
                min_value=-1,
                max_value=1,
                value=c.DEFAULT_PLUS_MINUS_CHARGE,
                step=1,
                help=h.PLUS_MINUS_CHARGE_HELP,
                key="plus_minus_charge",
            )

        st.subheader("Retention Time", divider="grey")
        infer_retention_time = stp.checkbox(
            label="Infer RT",
            value=c.DEFAULT_INFER_RT,
            help=h.INFER_RT_HELP,
            key="infer_retention_time",
        )

        retention_time = filter_invalid_rt = None
        if infer_retention_time:
            c1, c2 = st.columns(2, vertical_alignment="bottom")
            with c1:
                retention_time = stp.number_input(
                    label="RT (min)",
                    min_value=0,
                    value=c.DEFAULT_RT,
                    help=h.RT_HELP,
                    key="retention_time",
                )
            with c2:
                filter_invalid_rt = stp.checkbox(
                    label="Filter invalid RT",
                    value=c.DEFAULT_FILTER_INVALID_RT,
                    help=h.FILTER_INVALID_RT_HELP,
                    key="filter_invalid_rt",
                )

        st.subheader("Ion Mobility", divider="grey")
        infer_ion_mobility = stp.checkbox(
            label="Infer Ion Mobility",
            value=c.DEFAULT_INFER_ION_MOBILITY,
            help=h.INFER_ION_MOBILITY_HELP,
            key="infer_ion_mobility",
        )

        st.subheader("Proteotypic-ness", divider="grey")
        infer_proteotypic = stp.checkbox(
            label="Infer Proteotypic-ness",
            value=c.DEFAULT_INFER_PROTEOTYPIC,
            help=h.INFER_PROTEOTYPIC_HELP,
            key="infer_proteotypic",
        )
        remove_non_proteotypic = score_threshold = None
        if infer_proteotypic:
            c1, c2 = st.columns(2, vertical_alignment="bottom")
            with c1:
                remove_non_proteotypic = stp.checkbox(
                    label="Filter Non Proteotypic Peptides",
                    value=c.DEFAULT_REMOVE_NON_PROTEOTYPIC,
                    help=h.REMOVE_NON_PROTEOTYPIC_HELP,
                )
            with c2:
                score_threshold = stp.number_input(
                    label="Proteotypic Score Threshold",
                    min_value=0.0,
                    max_value=1.0,
                    value=c.DEFAULT_PROTEOTYPIC_THRESHOLD,
                    step=0.01,
                    help=h.PROTEOTYPIC_THRESHOLD_HELP,
                    key="proteotypic_score_threshold",
                )

    with st.expander("Static Modifications", expanded=False):
        st.subheader("Internal Modifications", divider="grey")
        df = stp.data_editor(
            pd.DataFrame(c.DEFAULT_STATIC_MODS_DATA),
            use_container_width=True,
            hide_index=True,
            key="static_mods",
            num_rows="dynamic",
            column_config={
                "Amino Acids": stp.column_config.SelectboxColumn(
                    options=list("ACDEFGHIKLMNPQRSTVWY"),
                    help=h.AMINO_ACIDS_HELP,
                    default="C",
                ),
                "Mass (Da)": stp.column_config.NumberColumn(
                    format="%.5f",
                    help=h.MOD_MASS_HELP,
                    default=57.02146,
                ),
            },
        )
        static_modifications = {
            row["Amino Acids"]: row["Mass (Da)"] for _, row in df.iterrows()
        }

        all_amino_acids = list("ACDEFGHIKLMNPQRSTVWY")

        st.subheader("N-Terminal Modifications", divider="grey")
        nc1, nc2 = st.columns(2, vertical_alignment="bottom")
        with nc1:
            n_term_mod_aas = st.multiselect(
                label="Target Amino Acids",
                options=all_amino_acids,
                help="Select amino acids to apply N-terminal modification. Leave empty to apply to any amino acid.",
                key="n_term_mod_aas",
            )
        with nc2:
            n_term_mod_mass = stp.number_input(
                label="Mass (Da)",
                value=0.0,
                step=1.0,
                format="%.5f",
                help="Mass to add to N-terminus",
                key="n_term_mod_mass",
            )
        
        n_term_modifications = {}
        if n_term_mod_mass != 0.0:
            target_aas = n_term_mod_aas if n_term_mod_aas else all_amino_acids
            n_term_modifications = {aa: n_term_mod_mass for aa in target_aas}

        st.subheader("C-Terminal Modifications", divider="grey")
        cc1, cc2 = st.columns(2, vertical_alignment="bottom")
        with cc1:
            c_term_mod_aas = st.multiselect(
                label="Target Amino Acids",
                options=all_amino_acids,
                help="Select amino acids to apply C-terminal modification. Leave empty to apply to any amino acid.",
                key="c_term_mod_aas",
            )
        with cc2:
            c_term_mod_mass = stp.number_input(
                label="Mass (Da)",
                value=0.0,
                step=1.0,
                format="%.5f",
                help="Mass to add to C-terminus",
                key="c_term_mod_mass",
            )

        c_term_modifications = {}
        if c_term_mod_mass != 0.0:
            target_aas = c_term_mod_aas if c_term_mod_aas else all_amino_acids
            c_term_modifications = {aa: c_term_mod_mass for aa in target_aas}

    return InputOptions(
        protein_id=protein_id,
        protein_sequence=protein_sequence,
        proteases_selected=proteases_selected,
        custom_regex=custom_regex,
        missed_cleavages=missed_cleavages,
        mass_type=mass_type,
        min_len=min_len,
        max_len=max_len,
        min_mass=min_mass,
        max_mass=max_mass,
        semi_enzymatic=semi_enzymatic,
        infer_charge=infer_charge,
        min_charge=min_charge,
        max_charge=max_charge,
        min_mz=min_mz,
        max_mz=max_mz,
        plus_minus_charge=plus_minus_charge,
        infer_retention_time=infer_retention_time,
        retention_time=retention_time,
        filter_invalid_rt=filter_invalid_rt,
        infer_ion_mobility=infer_ion_mobility,
        infer_proteotypic=infer_proteotypic,
        remove_non_proteotypic=remove_non_proteotypic,
        score_threshold=score_threshold,
        static_modifications=static_modifications,
        n_term_modifications=n_term_modifications,
        c_term_modifications=c_term_modifications,
    )
