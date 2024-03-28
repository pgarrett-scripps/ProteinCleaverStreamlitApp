PROTEASE_WIKI = """  
## Protein Cleavage:

Protein cleavage is the biochemical process of breaking proteins down into smaller peptide fragments. This is achieved through the hydrolysis of peptide bonds, which are the chemical links between amino acids. Typically, this process is facilitated by enzymes known as proteases, which expedite the cleavage of these bonds.

## Protease:

A protease is a type of enzyme specifically designed to carry out proteolysis, initiating the breakdown of proteins by cleaving their peptide bonds. Each protease is characterized by its 'cleavage specificity,' a unique ability to identify and cut at particular amino acids or specific sequences of amino acids. This specificity allows for targeted cleavage, influencing the size and sequence of the resulting peptide fragments.

For example Trypsin, the most widely used protease, cleaves after (at the C-Terminus of) Arginine (R) and Lysine (K).

![enter image description here](https://www.bocsci.com/upload/image/trypsin-cleavage-of-peptides.jpg)
  
## Missed Cleavages   

A missed cleavage refers to an instance where a protease does not cut a protein at a spot where it typically should. This oversight can happen for several reasons, including steric hindrance, modifications on the amino acids meant to be cut, or conditions that aren't ideal for the reaction. Steric hindrance, for example, can arise from the protein's complex folding patterns (its secondary, tertiary, or quaternary structures), which might block the protease's access to the intended cleavage site. Furthermore, post-translational modifications—chemical changes to the protein after it's made—near the cleavage site can alter how effectively the protease works. Other important factors affecting the occurrence of missed cleavages include the reaction's pH level, temperature, and the concentration of ions present.

![enter image description here](https://cores.imp.ac.at/fileadmin/additional_pages/core_facilities/protein_chemistry/pix/faq/missed_cleavage1.jpg)

## Semi-Enzymatic Peptides   

Semi-enzymatic peptides result from the enzymatic digestion of proteins where the peptides are cleaved by the protease at only one end. In an ideal digestion process, each peptide generated would exhibit cleavage at both the C-terminus and N-terminus, aside from the first and last peptides in the sequence, which naturally have one end unaltered by the enzyme. Semi-enzymatic peptides deviate from this ideal, showcasing only one end cleaved in alignment with the expected enzymatic action.
 
## Static Modifications

A static modification (also referred to as a fixed modification) is a post-translational modification (PTM) which is  
expected to occur at every instance of the specified residue(s). For example, one of the most common static modifications is the carbamidomethylation of cysteine residues. 

## Variable Modifications   

Contrastingly, variable modifications are those that may or may not occur at every possible site, and these are handled differently in data analysis. Common examples of variable modifications include methionine oxidation and protein N-terminal acetylation.   

## Average vs. Monoisotopic Mass   


  

  
In mass spectrometry, analyte mass can be calculated using two principal methods: monoisotopic mass and average mass. The selection between monoisotopic and average mass depends on the mass spectrometer's resolution. Lower-resolution spectrometers struggle to differentiate between isotopic peaks, leading to the reporting of a peak that represents the average mass of all isotopes of an element present in the molecule. Conversely, high-resolution spectrometers can distinguish individual isotopes, allowing for the precise determination of monoisotopic mass.

### Average Mass   

The average mass, also known as the molecular weight or the molecular mass, is the weighted average of the masses of all isotopes of an element, taking into account their natural abundance. For instance, carbon (C) has two naturally occurring isotopes: C-12 and C-13. The average mass of carbon takes into account the masses and abundances of these two isotopes.  

Average mass is calculated as follows:  
  
```  
Isotope = 12C  
Relative Atomic Mass = 12.0000000(00)  
Isotopic Composition = 0.9893(8)  
```  
  
```  
Isotope = 13C  
Relative Atomic Mass = 13.00335483507(23)  
Isotopic Composition = 0.0107(8)  
```  
  
```  
Average Mass = (12.0 * 0.99) + (13.0 * 0.01) = 12.01 Da  
```  
 
### Monoisotopic Mass   

Monoisotopic mass is the mass of a molecule, ion, or compound calculated using the mass of the most abundant isotope of each element. For carbon, this would be C-12, which has a mass of exactly 12 Da. So, for a peptide or protein, the monoisotopic mass would be calculated by summing the monoisotopic masses of all the individual amino acids.  
"""

HELP ="""
### Introduction to Protein Cleaver
Protein Cleaver is a versatile tool for protein analysis and digestion. It offers a range of features to help users:

1. **Identify Cleavage Sites:** Protein Cleaver assists in locating cleavage sites within protein sequences.

2. **Generate and Filter Peptides:** Users can generate peptides and apply filters to refine the selection.

3. **View Protein Sequence Coverage:** Gain insights into the coverage of a protein sequence by the generated peptides.

4. **Search for Motifs:** Easily search for specific motifs within the resulting peptides.

Protein Cleaver is particularly valuable in understanding why certain peptides might go unnoticed in mass 
spectrometry experiments. Potential reasons include overly strict search parameters or the physical and chemical 
properties of the peptide. To investigate the latter, Protein Cleaver can create a list of peptides that are 
likely to be generated and identified from a given protein sequence, referred to as proteotypic peptides.

Moreover, if users are interested in specific motifs like phosphorylation sites or particular amino acid residues, 
Protein Cleaver can assist in generating a list of peptides that contain these motifs. It also provides a 
visualization of motif coverage across the protein sequence, aiding in comprehensive analysis.

#### Inputting Your Protein Sequence:
- **FASTA Format Compatibility**: You can input your protein sequence in the FASTA format. Don't worry if your 
sequence extends over multiple lines;
- **Accession Number or Protein Identifier**: Alternatively, you can use a protein's accession number (e.g., P04406) 
or its unique protein identifier (e.g., ALBU_HUMAN) to retrieve and analyze the sequence.

#### Specifying Modifications:
- **Static Modifications**: You have the option to apply static modifications globally to all amino acids in the 
sequence.
- **Direct Sequence Modifications**: Modifications can also be included directly in the sequence. Use parentheses for 
standard modifications and square brackets for terminal modifications.
  
  Example: `[-13]MAS(1.2345)FRLFLLCLAGLVFVS[57.0]` indicates specific modifications at both the start and end of the 
  peptide sequence, as well as within it.
  
"""

COLUMN_DESCRIPTIONS = """
### Column Descriptions

**Sequence**: The peptide sequence.

**Start**: The position in the original protein sequence where this peptide begins (0 based indexing).

**End**: The position in the original protein sequence where this peptide ends (0 based indexing).

**Len**: The length of the unmodified peptide sequence, or simply the number of amino acids in the peptide.

**MC**: The number of missed cleavages for this peptide.

**Semi**: Indicates whether the peptide is semi-enzymatic or not.

**NeutralMass**: This is the neutral mass of the peptide, can be either the monoisotopic or average mass depending
on the selected option.

**Charge**: The estimated charge of the peptide, determined by the number of Lysine and Arginine residues + 1.

**Mz**: The m/z value of the peptide, calculated using the neutral mass and charge.

**RT**: The predicted normalized retention time of the peptide, calculated using a XGBoost model.

**IM**: The predicted ion mobility of the peptide, calculated using a XGBoost model.

**Proteotypic**: Indicates whether the peptide is proteotypic or not, determined using a XGBoost model.

**Motifs**: Indicates the number of motif matches within the peptide sequence.
"""

CONTACT = """
If you encounter any issues or have suggestions for improvement, please contact pgarrett@scripps.edu.
This is a work in progress and your feedback is greatly appreciated!
"""

PROTEOTYPIC_MODEL_HELP = """
### Proteotypic Model

#### Overview
This report details the development and evaluation of a machine learning model for predicting whether a peptide 
sequence has been observed (seen) or not. 

#### Data Preparation
- **Seen Peptides**: Extracted from a large number of timstof pro runs, modified sequences were stripped, and length 
constraints (6 to 50 amino acids) were applied.
- **Unseen Peptides**: Derived from a FASTA file of Hela cell proteins, these peptides were generated using a 
trypsin digestion simulation with specific parameters (missed cleavages: 2, semi-cleavage: False, length range: 6-50).
- **Balancing Data**: The dataset was balanced for both the 'Seen' status and sequence length, ensuring an equal number
 of seen and unseen peptides for each length category.

#### Feature Engineering
- **Amino Acid Binning**: Each sequence was represented by the count of each amino acid present, creating 
a feature vector for model input.

#### Model Training and Evaluation
- **Model**: An XGBClassifier was chosen for its high performance and efficiency.
- **Data Split**: The data was split into 80% training and 20% testing.

#### Results
- **Overall Accuracy**: 80.04%
- **Precision and Recall**: 
  - For 'False' (Unseen): Precision - 83%, Recall - 76%
  - For 'True' (Seen): Precision - 78%, Recall - 85%
- **Confusion Matrix**: 
  - True Negatives (Unseen, Correctly Classified): 28,531
  - False Positives (Unseen, Incorrectly Classified): 9,227
  - False Negatives (Seen, Incorrectly Classified): 5,847
  - True Positives (Seen, Correctly Classified): 31,921
- **Sequence Length Accuracy**:
  - 0-10: 74.10%
  - 10-20: 81.39%
  - 20-30: 86.48%
  - 30-40: 89.24%
  - 40-50: 97.92%

"""

RT_MODEL_HELP = """

### RT Model

#### Overview
This report details the development and evaluation of a machine learning model for predicting the retention time of a 
peptide sequence.

#### Data Preparation
- **Peptide Sequences**: Extracted from a large number of timstof pro runs, modified sequences were stripped, and length 
constraints (6 to 50 amino acids) were applied.
- **Alignment**: Retention times were aligned using a hierarchical alignment method
- **Retention Times ('RetTime')**: Normalized by dividing each value by the maximum value in the dataset.

#### Feature Engineering
- **Amino Acid Binning**: Each sequence was represented by the count of each amino acid present, creating 
a feature vector for model input.

#### Model Training and Evaluation
- **Model**: An XGBClassifier was chosen for its high performance and efficiency.
- **Data Split**: The data was split into 80% training and 20% testing.

#### Results
- **Mean Squared Error (MSE)**: 0.00322
- **Pearson Correlation**: 0.970
"""

IM_MODEL_HELP = """
### IM Model

#### Overview
This report details the development and evaluation of a machine learning model for predicting the ion mobility of a 
peptide sequence.

#### Data Preparation
- **Peptide Sequences**: Extracted from a large number of timstof pro runs, modified sequences were stripped, and length 
constraints (6 to 50 amino acids) were applied.
- **Alignment**: Ion mobility values were aligned using a hierarchical alignment method

#### Feature Engineering
- **Amino Acid Binning**: Each sequence was represented by the count of each amino acid present, creating a feature 
vector for model input.
- **Charge State**: The charge state of each peptide was also included as a feature.

#### Model Training and Evaluation
- **Model**: An XGBClassifier was chosen for its high performance and efficiency.
- **Data Split**: The data was split into 80% training and 20% testing.

#### Results
- **Mean Squared Error (MSE)**: 0.001235
- **Pearson Correlation**: 0.9741
"""

MODEL_CODE = """
import pickle
import pandas as pd
from collections import Counter
import numpy as np
import xgboost as xgb

def bin_aa_counts(seq, charge=None):
    amino_acids = 'ACDEFGHIKLMNPQRSTVWY'
    aa_counts = Counter(seq)

    enc = [aa_counts.get(aa, 0) for aa in amino_acids]
    if charge is not None:
        enc.append(charge)
    return enc
    
    
sequences = ['PEPTIDE', 'PEPTIDES']

# RT Model
rt_model = pickle.load(open("rt_model.pkl", "rb"))
preds = rt_model.predict(np.array([bin_aa_counts(seq) for seq in sequences]))

# IM Model
im_model = pickle.load(open("im_model.pkl", "rb"))
preds = im_model.predict(np.array([bin_aa_counts(seq, 2) for seq in sequences]))

# Proteotypic Model
proteotypic_model = pickle.load(open("proteotypic_model.pkl", "rb"))
preds = proteotypic_model.predict(np.array([bin_aa_counts(seq) for seq in sequences]))
"""