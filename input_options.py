from dataclasses import dataclass
import streamlit as st
import streamlit_permalink as stp
from constants import *
from help_messages import *
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

def fetch_sequence_from_uniprot(accession_number):
    url = f"https://www.uniprot.org/uniprot/{accession_number}.fasta"
    response = requests.get(url)
    if response.status_code != 200:
        st.error(f"Error fetching sequence from UniProt: {response.status_code}")
        return None
    return ''.join(response.text.split('\n')[1:])  # Remove the header line


@dataclass
class InputOptions:
    protein_id: str
    protein_sequence: str
    proteases_selected: list
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
    static_modifications: dict

    @cached_property
    def proteases(self):
        proteases = []
        for protease in self.proteases_selected:
            proteases.append(VALID_PROTEASES[protease])
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
        sites = set()
        for enzyme_regex in self.proteases:
            sites.update(pt.get_cleavage_sites(self.protein_sequence, enzyme_regex))
        return sites
    
    @cached_property
    def protein_annotation(self):
        return pt.parse(self.protein_sequence)
        



def validate_input_option(options: InputOptions) -> bool:
    if options.min_len > options.max_len:
        st.error(f'Min length: {options.min_len} must be less than max length: {options.max_len}')
        return False

    if options.min_charge is not None and options.max_charge is not None and options.min_charge > options.max_charge:
        st.error(f'Min charge: {options.min_charge} must be less than max charge: {options.max_charge}')
        return False
    
    if options.min_mz is not None and options.max_mz is not None and options.min_mz > options.max_mz:
        st.error(f'Min m/z: {options.min_mz} must be less than max m/z: {options.max_mz}')
        return False
    
    try:
        if options.sequence_length > MAX_PROTEIN_INPUT_LENGTH:
            st.error(f'Protein sequence length: {options.sequence_length} exceeds max length: {MAX_PROTEIN_INPUT_LENGTH}')
            return False
    except Exception as e:
        st.error(f'Error parsing protein sequence: {e}')
        return False
    
    if options.min_mass > options.max_mass:
        st.error(f'Min mass: {options.min_mass} must be less than max mass: {options.max_mass}')
        return False
    
    try:
        _ = options.neutral_mass
    except Exception as e:
        st.error(f'Error Calculating Protein Mass: {e}')
        return False
    
    return True


    
    

def get_input_options():
    st.subheader("Options", divider="grey")

    c1, c2 = st.columns(2, vertical_alignment='bottom')
    with c1:
        protein_id = st.text_input(label='Protein accession number / identifier',
                                value=DEFAULT_PROTEIN_ID,
                                placeholder=DEFAULT_PROTEIN_ID_PLACEHOLDER,
                                help=PROTEIN_ID_HELP,
                                max_chars=20,
                                key='protein_id')

    with c2:
        # btn to fetch sequence from uniprot
        fetch_sequence = st.button(label='Fetch from UniProt',
                                help='Fetch protein sequence from UniProt using the accession number',
                                use_container_width=True,
                                key='fetch_sequence')
        
    if fetch_sequence:
        if protein_id:
            fetched_protein_sequence = fetch_sequence_from_uniprot(protein_id)
            if fetched_protein_sequence is not None:
                st.query_params['protein_sequence'] = fetched_protein_sequence
                st.rerun()
        else:
            st.error('Please enter a protein accession number')

    protein_sequence = stp.text_area(label="Protein Sequence (Proforma 2.0)",
                                value=DEFAULT_PROTEIN_SEQUENCE,
                                help=PROTEIN_SEQUENCE_HELP,
                                key='protein_sequence',
                                max_chars=MAX_PROTEIN_INPUT_LENGTH,
                                height=200)

    c1, c2 = st.columns(2, vertical_alignment="bottom")

    with c1:
        proteases_selected = stp.multiselect(
            label="Proteases",
            options=list(VALID_PROTEASES.keys()),
            help=PROTEASES_HELP,
            default=DEFAULT_PROTEASES,
            key="protease",
        )

    with c2:
        custom_regex = stp.text_input(
            label="(Additional) Custom protease",
            value=DEFAULT_CUSTOM_REGEX,
            help=CUSTOM_REGEX_HELP,
            key="custom_regex",
        )


    c1, c2 = st.columns(2, vertical_alignment="bottom")

    with c1:
        missed_cleavages = stp.number_input(
            label="Max missed cleavages",
            min_value=MIN_MISSED_CLEAVAGES,
            max_value=MAX_MISSED_CLEAVAGES,
            value=DEFAULT_MISSED_CLEAVAGES,
            step=MISSED_CLEAVAGES_STEP,
            help=MISSED_CLEAVAGES_HELP,
            key="missed_cleavages",
        )

    with c2:
        mass_type = stp.selectbox(
            label="Mass type",
            options=MASS_TYPE_OPTIONS,
            index=DEFAULT_MASS_TYPE_INDEX,
            help=MASS_TYPE_HELP,
            key="mass_type",
        )


    c1, c2 = st.columns(2, vertical_alignment="bottom")

    with c1:
        min_len = stp.number_input(
            label="Min length",
            min_value=MIN_PEPTIDE_LEN,
            max_value=MAX_PEPTIDE_LEN,
            value=DEFAULT_MIN_PEPTIDE_LEN,
            step=PEPTIDE_LEN_STEP,
            help=MIN_PEPTIDE_LEN_HELP,
            key="min_peptide_len",
        )

    with c2:
        max_len = stp.number_input(
            label="Max length",
            min_value=MIN_PEPTIDE_LEN,
            max_value=MAX_PEPTIDE_LEN,
            value=DEFAULT_MAX_PEPTIDE_LEN,
            step=PEPTIDE_LEN_STEP,
            help=MAX_PEPTIDE_LEN_HELP,
            key="max_peptide_len",
        )


    c1, c2 = st.columns(2, vertical_alignment="bottom")

    with c1:
        min_mass = stp.number_input(
            label="Min neutral mass",
            min_value=MIN_PEPTIDE_MASS,
            max_value=MAX_PEPTIDE_MASS,
            value=DEFAULT_MIN_PEPTIDE_MASS,
            step=PEPTIDE_MASS_STEP,
            help=MIN_PEPTIDE_MASS_HELP,
            key="min_peptide_mass",
        )

    with c2:
        max_mass = stp.number_input(
            label="Max neutral mass",
            min_value=MIN_PEPTIDE_MASS,
            max_value=MAX_PEPTIDE_MASS,
            value=DEFAULT_MAX_PEPTIDE_MASS,
            step=PEPTIDE_MASS_STEP,
            help=MAX_PEPTIDE_MASS_HELP,
            key="max_peptide_mass",
        )

    semi_enzymatic = stp.checkbox(
        label="Semi enzymatic?",
        help=SEMI_ENZYMATIC_HELP,
        value=DEFAULT_SEMI_ENZYMATIC,
        key="semi_enzymatic",
    )

    with st.expander("Peptide Properties", expanded=False):

        st.subheader("Charge", divider="grey")
        infer_charge = stp.checkbox(
            label="Infer Charge",
            value=DEFAULT_INFER_CHARGE,
            help=INFER_CHARGE_HELP,
            key="infer_charge",
        )
        min_charge = max_charge = min_mz = max_mz = plus_minus_charge = None
        if infer_charge:
            c1, c2 = st.columns(2, vertical_alignment="bottom")
            with c1:
                min_charge = stp.number_input(
                    label="Min charge",
                    min_value=MIN_CHARGE,
                    max_value=MAX_CHARGE,
                    value=DEFAULT_MIN_CHARGE,
                    step=CHARGE_STEP,
                    help=MIN_CHARGE_HELP,
                    key="min_charge",
                )

                min_mz = stp.number_input(
                    label="Min m/z",
                    min_value=MIN_MZ,
                    max_value=MAX_MZ,
                    value=DEFAULT_MIN_MZ,
                    step=MZ_STEP,
                    help=MIN_MZ_HELP,
                    key="min_mz",
                )

            with c2:
                max_charge = stp.number_input(
                    label="Max charge",
                    min_value=MIN_CHARGE,
                    max_value=MAX_CHARGE,
                    value=DEFAULT_MAX_CHARGE,
                    step=CHARGE_STEP,
                    help=MAX_CHARGE_HELP,
                    key="max_charge",
                )

                max_mz = stp.number_input(
                    label="Max m/z",
                    min_value=MIN_MZ,
                    max_value=MAX_MZ,
                    value=DEFAULT_MAX_MZ,
                    step=MZ_STEP,
                    help=MAX_MZ_HELP,
                    key="max_mz",
                )
            plus_minus_charge = stp.number_input(
                label="Charge +/-",
                min_value=-1,
                max_value=1,
                value=DEFAULT_PLUS_MINUS_CHARGE,
                step=1,
                help=PLUS_MINUS_CHARGE_HELP,
                key="plus_minus_charge",
            )

        st.subheader("Retention Time", divider="grey")
        infer_retention_time = stp.checkbox(
            label="Infer RT",
            value=DEFAULT_INFER_RT,
            help=INFER_RT_HELP,
            key="infer_retention_time",
        )

        retention_time = filter_invalid_rt = None
        if infer_retention_time:
            c1, c2 = st.columns(2, vertical_alignment="bottom")
            with c1:
                retention_time = stp.number_input(
                    label="RT (min)",
                    min_value=0,
                    value=DEFAULT_RT,
                    help=RT_HELP,
                    key="retention_time",
                )
            with c2:
                filter_invalid_rt = stp.checkbox(
                    label="Filter invalid RT",
                    value=DEFAULT_FILTER_INVALID_RT,
                    help=FILTER_INVALID_RT_HELP,
                    key="filter_invalid_rt",
                )

        st.subheader("Ion Mobility", divider="grey")
        infer_ion_mobility = stp.checkbox(
            label="Infer Ion Mobility",
            value=DEFAULT_INFER_ION_MOBILITY,
            help=INFER_ION_MOBILITY_HELP,
            key="infer_ion_mobility",
        )

        st.subheader("Proteotypic-ness", divider="grey")
        infer_proteotypic = stp.checkbox(
            label="Infer Proteotypic-ness",
            value=DEFAULT_INFER_PROTEOTYPIC,
            help=INFER_PROTEOTYPIC_HELP,
            key="infer_proteotypic",
        )
        remove_non_proteotypic = score_threshold = None
        if infer_proteotypic:

            c1, c2 = st.columns(2, vertical_alignment="bottom")
            with c1:
                remove_non_proteotypic = stp.checkbox(
                    label="Filter Non Proteotypic Peptides",
                    value=DEFAULT_REMOVE_NON_PROTEOTYPIC,
                    help=REMOVE_NON_PROTEOTYPIC_HELP,
                )
            with c2:
                score_threshold = stp.number_input(
                    label="Proteotypic Score Threshold",
                    min_value=0.0,
                    max_value=1.0,
                    value=DEFAULT_PROTEOTYPIC_THRESHOLD,
                    step=0.01,
                    help=PROTEOTYPIC_THRESHOLD_HELP,
                    key="proteotypic_score_threshold",
                )

    with st.expander("Static Modifications", expanded=False):
        df = stp.data_editor(
            pd.DataFrame(DEFAULT_STATIC_MODS_DATA),
            use_container_width=True,
            hide_index=True,
            key="static_mods",
            num_rows="dynamic",
            column_config={
                "Amino Acids": stp.column_config.SelectboxColumn(
                    options=list(pt.AMINO_ACIDS),
                    help=AMINO_ACIDS_HELP,
                    default="C",
                ),
                "Mass (Da)": stp.column_config.NumberColumn(
                    format="%.5f",
                    help=MOD_MASS_HELP,
                    default=57.02146,
                ),
            },
        )
        static_modifications = {
            row["Amino Acids"]: row["Mass (Da)"] for _, row in df.iterrows()
        }

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
    )