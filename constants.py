import copy
import os

from peptacular.constants import PROTEASES


def get_env_int(var_name, default):
    return int(os.getenv(var_name, default))


def get_env_float(var_name, default):
    return float(os.getenv(var_name, default))


def get_env_str(var_name, default):
    return os.getenv(var_name, default)


MAX_PROTEIN_LEN = get_env_int('MAX_PROTEIN_LEN', 2_500)
MAX_PROTEIN_INPUT_LENGTH = get_env_int('MAX_PROTEIN_INPUT_LENGTH', 10_000)


MIN_MISSED_CLEAVAGES = 0
MAX_MISSED_CLEAVAGES = get_env_int('MAX_MISSED_CLEAVAGES', 7)
DEFAULT_MISSED_CLEAVAGES = 1
MISSED_CLEAVAGES_STEP = 1

MIN_PEPTIDE_LEN = get_env_int('MIN_PEPTIDE_LEN', 1)
MAX_PEPTIDE_LEN = get_env_int('MAX_PEPTIDE_LEN', 50)
DEFAULT_MIN_PEPTIDE_LEN = 7
DEFAULT_MAX_PEPTIDE_LEN = 30
PEPTIDE_LEN_STEP = 1

MIN_PEPTIDE_MASS = get_env_float('MIN_PEPTIDE_MASS', 0.0)
MAX_PEPTIDE_MASS = get_env_float('MAX_PEPTIDE_MASS', 10_000.0)
DEFAULT_MIN_PEPTIDE_MASS = 500.0
DEFAULT_MAX_PEPTIDE_MASS = 6_000.0
PEPTIDE_MASS_STEP = 100.0

DEFAULT_PROTEASES = ["trypsin/P"]

MASS_TYPE_OPTIONS = ['monoisotopic', 'average']
DEFAULT_MASS_TYPE = 'monoisotopic'
DEFAULT_MASS_TYPE_INDEX = MASS_TYPE_OPTIONS.index(DEFAULT_MASS_TYPE)

MIN_STATIC_MODS = get_env_int('MIN_STATIC_MODS', 0)
MAX_STATIC_MODS = get_env_int('MAX_STATIC_MODS', 10)
DEFAULT_STATIC_MODS = 1
STATIC_MODS_STEP = 1

MIN_VAR_MODS = get_env_int('MIN_VAR_MODS', 0)
MAX_VAR_MODS = get_env_int('MAX_VAR_MODS', 10)
DEFAULT_VAR_MODS = 0
VAR_MOD_STEP = 1

MIN_CHARGE = get_env_int('MIN_CHARGE', 1)
MAX_CHARGE = get_env_int('MAX_CHARGE', 10)
DEFAULT_MIN_CHARGE = 2
DEFAULT_MAX_CHARGE = 5
CHARGE_STEP = 1

MIN_MAX_VAR_MODS = get_env_int('MIN_MAX_VAR_MODS', 0)
MAX_MAX_VAR_MODS = get_env_int('MAX_MAX_VAR_MODS', 3)
DEFAULT_MAX_VAR_MODS = 0
MAX_VAR_MOD_STEP = 1

MIN_MZ = get_env_float('MIN_MZ', 0.0)
MAX_MZ = get_env_float('MAX_MZ', 10_000.0)
DEFAULT_MIN_MZ = 200.0
DEFAULT_MAX_MZ = 1_800.0
MZ_STEP = 100.0

CMAP = 'hot'

DEFAULT_PROTEIN_SEQUENCE = "MASFRLFLLCLAGLVFVSEAGSVGAGEPKCPLMVKVLDAVRGSPAANVGVKVFKKAADETWEPFASGKTSESGELHGLTTEDKFVEGLY" \
                           "KVELDTKSYWKSLGISPFHEFAEVVFTANDSGPRHYTIAALLSPYSYSTTALVSSPKA"

LINK = get_env_str('LINK', 'https://peptidefragmenter.streamlit.app/')

VALID_PROTEASES = copy.deepcopy(PROTEASES)
VALID_PROTEASES.pop('non-specific', None)
VALID_PROTEASES = {k.replace(' ', '_'): v for k, v in VALID_PROTEASES.items()}

BASE_URL = get_env_str('BASE_URL', 'http://localhost:8501')

