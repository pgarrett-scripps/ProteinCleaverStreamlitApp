PROTEASE_WIKI = """  
## Protein Cleavage:  

Protein cleavage is the process by which proteins are broken down into smaller peptides. This occurs through hydrolysis 
of the amide bonds linking amino acids, which is typically catalyzed by an protease.  

## Protease:  

A protease is an enzyme that catalyzes proteolysis and begins the protein cleavage process by breaking down the 
peptide bonds. Different proteases have different 'cleavage specificity', which means they recognize and cleave at 
specific sequences or types of amino acids.  

### How does a Protease work?  

#### 1 - Recognition/Binding:   
In the first step in the protein cleavage process the protease recognizes and binds to a specific amino acid or amino 
acid sequence. This is typically facilitated by the protease's active site, which has a shape and chemical environment 
conducive to binding the target protein.


For example, trypsin binds after Lysine (K) or Arginine (R). Let's use the first Arginine (R) in the following protein 
as an example:
```
Protein:       H-W-P-R-A-T-G-A-K-Y-G-G-L
                     ^
                     |
Protease (Trypsin): -R-
```

#### 2 - Cleavage:   
Once the protease has bound to the target sequence, it catalyzes a hydrolysis reaction, breaking the peptide bond 
between two specific amino acids. This results in the protein being cleaved into smaller peptides or individual amino 
acids.  

```  
Peptide 1: H-W-P-R  
Peptide 2: A-T-G-A-K-Y-G-G-L  
```  

#### 3 - Release:   
After the cleavage has occurred, the smaller peptide sequences are released from the active site, and the protease 
can then go on to catalyze another reaction.  


Since there is still another active site in peptide 2, the protease can continue.  
```
Protease (Trypsin): Ready for next cleavage.

Peptide 2:   A-T-G-A-K-Y-G-G-L
                     ^
                     |
Protease (Trypsin): -K-
```

## Missed Cleavages  

A missed cleavage occurs when a protease fails to cleave a protein at a location where it typically would. Each type 
of protease has specific target amino acid sequences where they are expected to cleave. However, in some cases, 
the protease might not cleave at these specific sites. This event is termed as a missed cleavage.  

### Causes  

Missed cleavages can occur for various reasons such as steric hindrance, modifications on the target residues, or 
suboptimal reaction conditions. Steric hindrance might be caused by the protein's secondary, tertiary, or quaternary 
structures which can prevent the protease from accessing the cleavage site. Post-translational modifications like 
methylation, acetylation, or phosphorylation on the target residues can also influence cleavage efficiency. Additional 
factors such as pH, temperature, and ionic strength of the reaction can also contribute to missed cleavages.

### Impact  

Missed cleavages are of particular importance in proteomics research because they increase the complexity of peptide 
mixtures, thus making protein identification and quantification more challenging. However, information about missed 
cleavages can provide valuable insights into the protein's structure and modifications.

### Example  

Consider the protein sequence,  **H-W-K-R-A-T-K-G-A-L-Y-G-G-L. Digestion** with trypsin, would be expected to yield 
four peptides:

A perfect cleavage scenario would yield four peptides:  

```  
H-W-K  
R  
A-T-K  
G-A-L-Y-G-G-L  
```  

But if the first cleavage site is missed, we would end up with a different set of peptides, with one larger peptide 
instead of two smaller ones. 

```  
H-W-K-R  
A-T-K  
G-A-L-Y-G-G-L  
```  

## Semi-Enzymatic Peptides  

Semi-enzymatic peptides are peptides produced during enzymatic digestion of proteins that have only one cleavage site 
consistent with the specificity of the protease used.  

### Explanation  

During protein digestion, proteases cleave the protein at specific residues or sequences. In a perfect digestion 
scenario, every peptide resulting from this process would have a cleavage site at both the C-terminus and N-terminus 
(with the exception of the first and last peptides in the sequence, which naturally only have one enzymatic terminus).  

However, incomplete digestion or other processes can result in peptides that only have a cleavage-consistent residue 
at one terminus. These are referred to as semi-enzymatic peptides.  

### Importance  

The identification of semi-enzymatic peptides is important for comprehensive proteome coverage. It can provide 
insights into protease efficiency, post-translational modifications, and protein structure.   
### Example  

Consider digestion of the following protein sequence by trypsin: **T-V-K-A-T-R-G-L-I-M**. 

A fully enzymatic digestion would produce these peptides:   
```  
T-V-K  
A-T-R  
G-L-I-M  
```  

However, if a cleavage at the second K is missed and the protease cleaves after M instead, a semi-enzymatic peptide 
is produced: **T-V-K-A-T-R-G-L-I-M**.  

In this semi-enzymatic peptide, the trypsin cleavage site (K) is present only at the N-terminus. The C-terminus ends 
with a methionine, which is not a typical trypsin cleavage site.  

## Static Modifications  

A static modification (also referred to as a fixed modification) is a post-translational modification (PTM) where a 
characteristic molecular group is permanently attached to an amino acid residue in a protein. This modification is 
expected to occur at every instance of the specified residue(s) in every protein analyzed.

For example, one of the most common static modifications is the carbamidomethylation of cysteine residues. In this 
modification, iodoacetamide or iodoacetic acid reacts with the sulfhydryl group on the cysteine residue to form a 
carbamidomethyl group. The mass of the cysteine residue is thus increased by 57.021464 Daltons (Da).  

It is known as a "static" modification because it is assumed to occur 100% of the time at all specified residues. This 
means that during mass spectrometry data analysis, the mass of the specified residues is adjusted to include the mass 
of the modification, and the search algorithm does not consider the possibility of the residue being unmodified.  

Contrastingly, variable modifications are those that may or may not occur at every possible site, and these are handled 
differently in data analysis. Common examples of variable modifications include methionine oxidation and protein 
N-terminal acetylation.  

## Average vs. Monoisotopic Mass  

In the context of mass spectrometry and proteomics, two important concepts related to the mass of atoms, ions, 
molecules, or compounds are average mass and monoisotopic mass.  

### Average Mass  

The average mass, also known as the molecular weight or the molecular mass, is the weighted average of the masses of 
all isotopes of an element, taking into account their natural abundance. For instance, carbon (C) has two naturally 
occurring isotopes: C-12 and C-13. The average mass of carbon takes into account the masses and abundances of these 
two isotopes.  

Similarly, when calculating the average mass of a peptide or a protein, the average masses of all the individual amino 
acids (which again, take into account the different isotopes of all the atoms in the amino acid) are summed.  

### Monoisotopic Mass  

Monoisotopic mass is the mass of a molecule, ion, or compound calculated using the mass of the most abundant isotope of 
each element. For carbon, this would be C-12, which has a mass of exactly 12 Da. So, for a peptide or protein, the 
monoisotopic mass would be calculated by summing the monoisotopic masses of all the individual amino acids.

### Considerations in Proteomics  

The choice between using average mass and monoisotopic mass in a proteomics study depends on the type of mass 
spectrometer and the subsequent data analysis used. Lower resolution mass spectrometers often report average masses, 
whereas higher resolution instruments can discern individual isotopes and thus report monoisotopic masses. Data analysis 
algorithms must use the same type of mass that the instrument reports for accurate identification of peptides and 
proteins.  

Typically, monoisotopic mass is preferred in proteomics as it provides a more exact mass measurement, which can improve 
the accuracy of peptide and protein identifications. However, in some cases, such as with larger proteins or lower 
resolution instruments, it may be more appropriate to use average mass.
"""

HELP = """
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
