---

# Protein Cleaver 🔪

Protein Cleaver is a streamlit application designed to compute protease-specific cleavage sites and peptides based on user-provided protein sequences. It also offers functionality to visualize protein coverage, allowing researchers to gain insights into the proteolysis of proteins.

Check it out on Streamlits Community Cloud: https://proteincleaver.streamlit.app/

## Features

- **Input Protein Sequences**: Users can input protein sequences in FASTA format.
- **Visualize Cleavage Sites**: The app highlights cleavage sites directly on the protein sequence.
- **Generate Peptide Data**: Peptide data derived from the cleavage is displayed in tabular format.
- **Visualize Protein Coverage**: The app provides a visual representation of which parts of the protein sequence are covered by the generated peptides.
- **Protein Coverage Plots**: The app offers plots to visualize protein coverage against different criteria such as missed cleavages, peptide length, and peptide mass.

## Usage

1. **Inputting Protein Sequence**:
    - Navigate to the sidebar on the left.
    - Enter your protein sequence in the provided textarea under the label "Protein sequence".
    - The sequence should be in FASTA format.

2. **Configuring the Digestion**:
    - Select the proteases for cleavage.
    - Adjust the maximum number of missed cleavages.
    - Set the minimum and maximum peptide length.
    - Define the peptide mass range.
    - Choose between monoisotopic or average mass types.
    - Add any static modifications if necessary.

3. **View Results**:
    - Navigate through the tabs to view:
        - **Digest**: View the results of the digestion, including cleavage sites, peptides, and associated metrics.
        - **Cleavage**: Visualize the cleavage sites on the protein sequence.
        - **Coverage**: See the protein coverage and associated plots.

4. **Additional Information**:
    - The **Wiki** tab provides additional information on the proteases used.
    - The **Help** tab offers guidance and details on the app's functionality.
