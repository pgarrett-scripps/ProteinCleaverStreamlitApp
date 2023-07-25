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


PROTEASE_WIKI = """
## Protein Cleavage:

Protein cleavage is the process by which proteins are broken down into smaller polypeptides or single amino acids. This occurs through a hydrolysis reaction, where a water molecule breaks a peptide bond, the bond linking two amino acids together in a protein.

## Protease:

A protease is an enzyme that conducts proteolysis, i.e., it begins the protein cleavage process by breaking down the peptide bonds. Proteases are crucial for many biological functions including digestion, immune response, cell cycle progression, and protein recycling, among others.

Different proteases have different 'cleavage specificity', which means they recognize and cleave at specific sequences or types of amino acids.

### How does a Protease work?

#### 1 - Recognition/Binding: 

The first step involves the protease recognizing and binding to the specific sequence or type of amino acid where it will cleave. This is typically facilitated by the protease's active site, which has a shape and chemical environment conducive to binding the target protein. Fo


For trypsin, it's after Lysine (K) or Arginine (R). Let's use the first Arginine (R) as an example:
```
Protein:       H-W-P-R-A-T-G-A-K-Y-G-G-L
                     ^
                     |
Protease (Trypsin): -R-
```

#### 2 - Cleavage: 

Once the protease has bound to the target sequence, it catalyzes a hydrolysis reaction, breaking the peptide bond between two specific amino acids. This results in the protein being cleaved into smaller peptides or individual amino acids.

```
Peptide 1: H-W-P-R
Peptide 2: A-T-G-A-K-Y-G-G-L
```

#### 3 - Release: 

After the cleavage has occurred, the smaller peptide sequences are released from the active site, and the protease can then go on to catalyze another reaction.


Since there are still another active site in peptide 2, the protease can continue.
```
Protease (Trypsin): Ready for next cleavage.

Peptide 2:   A-T-G-A-K-Y-G-G-L
                     ^
                     |
Protease (Trypsin): -K-
```

## Missed Cleavages

A missed cleavage occurs when a protease fails to cleave a protein at a location where it typically would. Each type of protease has specific target amino acid sequences where they are expected to cleave. However, in some cases, the protease might not cleave at these specific sites. This event is termed as a missed cleavage.

### Causes

Missed cleavages can occur due to various reasons such as steric hindrance, modifications on the target residues, or suboptimal reaction conditions. Steric hindrance might be caused by the protein's secondary, tertiary, or quaternary structures which can prevent the protease from accessing the cleavage site. Post-translational modifications like methylation, acetylation, or phosphorylation on the target residues can also influence cleavage efficiency. Furthermore, factors like pH, temperature, and ionic strength of the reaction can also contribute to missed cleavages.

### Impact

Missed cleavages are of particular importance in proteomics research. They can increase the complexity of peptide mixtures, thus making protein identification and quantification more challenging. On the other hand, information about missed cleavages can provide valuable insights into the protein's structure and modifications.

### Example

Consider a protein sequence, **H-W-K-R-A-T-K-G-A-L-Y-G-G-L**, and trypsin as the protease which typically cleaves after lysine (K) and arginine (R).

A perfect cleavage scenario would yield four peptides:

```
H-W-K
R
A-T-K
G-A-L-Y-G-G-L
```

But, if a cleavage after the first **K** is missed, we would end up with different peptides:

```
H-W-K-R
A-T-K
G-A-L-Y-G-G-L
```


In the missed cleavage scenario, the first two peptides from the fully cleaved protein are merged into one larger peptide. Understanding and considering these missed cleavage events can be crucial for accurate protein identification and quantification.

## Semi-Enzymatic Peptides

In proteomics, semi-enzymatic peptides are peptides produced during enzymatic digestion of proteins that have only one cleavage site consistent with the specificity of the protease used.

### Explanation

During protein digestion, proteases cleave the protein at specific residues or sequences. For example, trypsin, a commonly used protease in proteomics, typically cleaves proteins at the carboxyl side of lysine and arginine residues. In a perfect digestion scenario, every peptide resulting from this process would have a lysine or arginine residue at both the C-terminus and N-terminus (with the exception of the first and last peptides in the sequence, which naturally only have one enzymatic terminus).

However, incomplete digestion or other processes can result in peptides that only have a cleavage-consistent residue at one terminus. These are referred to as semi-enzymatic peptides. For example, if trypsin were to miss a cleavage site and cleave at a different location, it may result in a peptide with a lysine or arginine only at the N-terminus or C-terminus.

### Importance

The identification of semi-enzymatic peptides is important for comprehensive proteome coverage. It can provide insights into protease efficiency, post-translational modifications, and protein structure. 

### Example

Consider the following protein sequence digested by trypsin: **T-V-K-A-T-R-G-L-I-M**.

A fully enzymatic digestion would produce these peptides: 

```
T-V-K
A-T-R
G-L-I-M
```

However, if a cleavage at the second K is missed and the protease cleaves after M instead, a semi-enzymatic peptide is produced: **T-V-K-A-T-R-G-L-I-M**.

In this semi-enzymatic peptide, the trypsin cleavage site (K) is present only at the N-terminus. The C-terminus ends with a methionine, which is not a typical trypsin cleavage site.

## Static Modifications

In the field of proteomics, a static modification, also referred to as a fixed modification, is a post-translational modification (PTM) where a certain molecular group is permanently attached to an amino acid residue in a protein. This modification is expected to occur at every instance of the specified residue(s) in every protein analyzed.

For example, one of the most common static modifications is the carbamidomethylation of cysteine residues. In this modification, iodoacetamide or iodoacetic acid reacts with the sulfhydryl group on the cysteine residue to form a carbamidomethyl group. The mass of the cysteine residue is thus increased by 57.021464 Daltons (Da).

It is known as a "static" modification because it is assumed to occur 100% of the time at all specified residues. This means that during mass spectrometry data analysis, the mass of the specified residues is adjusted to include the mass of the modification, and the search algorithm does not consider the possibility of the residue being unmodified.

Contrastingly, variable modifications are those that may or may not occur at every possible site, and these are handled differently in data analysis. Common examples of variable modifications include methionine oxidation and protein N-terminal acetylation.

## Average vs. Monoisotopic Mass

In the context of mass spectrometry and proteomics, two important concepts related to the mass of atoms, ions, molecules, or compounds are average mass and monoisotopic mass.

### Average Mass

The average mass, also known as the molecular weight or the molecular mass, is the weighted average of the masses of all isotopes of an element, taking into account their natural abundance. For instance, carbon (C) has two naturally occurring isotopes: C-12 and C-13. The average mass of carbon takes into account the masses and abundances of these two isotopes.

Similarly, when calculating the average mass of a peptide or a protein, the average masses of all the individual amino acids (which again, take into account the different isotopes of all the atoms in the amino acid) are added up.

### Monoisotopic Mass

On the other hand, the monoisotopic mass is the mass of a molecule, ion, or compound calculated using the mass of the most abundant isotope of each element. For carbon, this would be C-12, which has a mass of exactly 12 Da. So, for a peptide or protein, the monoisotopic mass would be calculated by summing up the masses of the individual amino acids, using the monoisotopic mass for each one.

### Considerations in Proteomics

The choice between using average mass and monoisotopic mass in a proteomics study depends on the type of mass spectrometer and the subsequent data analysis used. Lower resolution mass spectrometers often report average masses, whereas higher resolution instruments can discern individual isotopes and thus report monoisotopic masses. Data analysis algorithms must use the same type of mass that the instrument reports for accurate identification of peptides and proteins.

Typically, monoisotopic mass is preferred in proteomics as it provides a more exact mass measurement, which can improve the accuracy of peptide and protein identifications. However, in some cases, such as with larger proteins or lower resolution instruments, it may be more appropriate to use average mass.


"""


HELP ="""

## Column Labels

**PeptideSequence**: The sequence of amino acids that make up a specific peptide generated from the digestion process.

**Start**: The position in the original protein sequence where this peptide begins.

**End**: The position in the original protein sequence where this peptide ends.

**AACount**: This stands for Amino Acid Count, it represents the number of amino acids in the peptide.

**MC**: Missed Cleavages. This is the number of cleavage sites in the peptide that were not cut during the digestion process. In other words, the number of times a protease could have cleaved the peptide further, but did not.

**IsSemi**: This is a boolean (true/false) column indicating whether the peptide is semi-enzymatic. A semi-enzymatic peptide is one where one end of the peptide is a result of enzymatic cleavage, while the other end is not.

**Mass**: This is the mass of the peptide. The mass calculation can be done using either the monoisotopic mass or the average mass of the peptide, depending on what the user selects.

## Statistics:

**Cleavage Sites**: The total number of cleavage sites identified in the protein sequence.

**Total Peptides**: The total number of peptides generated from the digestion process.

**Unique Peptides**: The number of distinct peptide sequences generated from the digestion process.

**Semi Peptides**: The number of semi-enzymatic peptides generated.

**Enzymatic Peptides**: The number of fully enzymatic peptides generated.

**Protein Coverage**: The percentage of the protein sequence that is covered by the generated peptides.
"""
