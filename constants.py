import copy

from peptacular.constants import PROTEASES

MAX_PROTEIN_LEN = 2_500
MAX_MISSED_CLEAVAGES = 10
MAX_PEPTIDE_LEN = 50
MAX_PEPTIDE_MASS = 10_000

EXAMPLE_PROTEIN = "HRNGEMYACEQEHDKEPHMKIMPHGSGGFFPLVQFGRHFGQLKNLKRPAVHVDTEVLYWCNTRCEFLMWAFDCRIDPRDWGMDHMHCRESRCYASFRG" \
                  "TRGFDNLFYPAKHLEMHGTMISIMQWFQANGDKTLHSTYKFMSPCSGEKRMYQSWKWGEKPRCYSTQHVYCAVDKRMSVWSKCFSQGKALGTKESLNN" \
                  "VDDHHDLKQCVMISSWSTPYCKIPNCAAEQWMETMTMPWDWPPMFIKIVIASDRCVVLHPQLGLHAHGMTRWATTVRKGKIGFYDPGPNMCYWQQQWL" \
                  "FTVAGS"

LINK = 'https://peptidefragmenter.streamlit.app/'

VALID_PROTEASES = copy.deepcopy(PROTEASES)
VALID_PROTEASES.pop('non-specific', None)
