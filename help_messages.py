# Help messages for input options

# Protein input help messages
PROTEIN_ID_HELP = 'A protein accession number / identifier'
PROTEIN_SEQUENCE_HELP = 'An amino acid sequence to digest'

# Protease help messages
PROTEASES_HELP = "The proteases to use for digestion"
CUSTOM_REGEX_HELP = "A custom regular expression to use for digestion. Will be used along with selected proteases"

# Digestion parameters help messages
MISSED_CLEAVAGES_HELP = "Number of missed cleavages to allow during digestion"
MASS_TYPE_HELP = "Mass type to use for calculations"

# Peptide length and mass help messages
MIN_PEPTIDE_LEN_HELP = "Minimum peptide length (inclusive)"
MAX_PEPTIDE_LEN_HELP = "Maximum peptide length (inclusive)"
MIN_PEPTIDE_MASS_HELP = "Minimum peptide neutral mass (inclusive)"
MAX_PEPTIDE_MASS_HELP = "Maximum peptide neutral mass (inclusive)"
SEMI_ENZYMATIC_HELP = "Allow semi enzymatic peptides?"

# Charge help messages
INFER_CHARGE_HELP = "Infer charge of peptides based on +/- 1 of (#Lysines and #Arginines + 1)"
MIN_CHARGE_HELP = "Minimum peptide charge (inclusive)"
MAX_CHARGE_HELP = "Maximum peptide charge (inclusive)"
MIN_MZ_HELP = "Minimum peptide m/z (inclusive)"
MAX_MZ_HELP = "Maximum peptide m/z (inclusive)"
PLUS_MINUS_CHARGE_HELP = "Charge +/- for m/z calculations"

# Retention time help messages
INFER_RT_HELP = "Use retention time model to predict retention times for peptides"
RT_HELP = "Retention time of the peptide in minutes"
FILTER_INVALID_RT_HELP = "Filter out peptides with invalid retention times"

# Ion mobility help message
INFER_ION_MOBILITY_HELP = "Use ion mobility model to predict ion mobility for peptides"

# Proteotypic help messages
INFER_PROTEOTYPIC_HELP = "Use proteotypic model to predict proteotypic peptides"
REMOVE_NON_PROTEOTYPIC_HELP = "Remove peptides that are not proteotypic"
PROTEOTYPIC_THRESHOLD_HELP = "Proteotypic score threshold for filtering peptides"

# Static modifications help messages
AMINO_ACIDS_HELP = "Select amino acids for which to apply the static modification"
MOD_MASS_HELP = "The mass of the modification (in daltons)"
