![VBRN Logo](/logo_landscape.png)

# [UVM Proteomics](https://vbrn.org/proteomics/)

### Authors
[Emily Curd](https://scholar.google.com/citations?user=uGHWHbgAAAAJ&hl=en&oi=ao) and James (Jax) Lubkowitz

### Affiliations
[Vermont Biomedical Research Network](https://vbrn.org), [Vermont Integrated Genomics Resource](https://www.med.uvm.edu/vigr/home), [University of Vermont](https://www.uvm.edu), [St. Lawrence University](https://www.stlawu.edu)

### Funding
Research reported in this repository was supported by an Institutional Development Award (IDeA) from the National Institute of General Medical Sciences of the National Institutes of Health under grant number P20GM103449. Its contents are solely the responsibility of the authors and do not necessarily represent the official views of NIGMS or NIH.

This is a R package for the [UVM - VBRN Proteomics Core](https://vbrn.org/proteomics/).  It includes functions for three main types of analyses: Post Translational Modification, Spectral Counting, and Protein Expression.

***Currently Supports or In Development***

Post Translational Modification:
  - Motif Analysis
    
Spectral Counting:
  - Differential protein expression (DEqMS)
    
Protein Expression:
  - Pathway Analysis
  - functional classification

### Installation
Package can be installed from github via
```R
install.packages("devtools")
devtools::install_github("UVM_Proteomics",username="jaxlub")
```
Another download option is to download a zip file from github, unzip it and then run
```R
install.packages("devtools")
source <- devtools:::source_pkg("path/to/package")
install(source)
```

# Post Translational Modification

## Motif Extraction

This function extracts predetermined protein motifs from a master protein in a user specified proteome (e.g. [uniProt mouse proteome](https://www.uniprot.org/proteomes?facets=proteome_type:1&query=(organism_id:10090)). The function identifies the position of the modified amino acid(s) in a master protein, and uses that information as a midpoint to calculate the position range of amino acids in a motif. The blast+ suite is then used to search a proteome and extract the motifs. To do this, a user defined proteome is converted to a blast searchable database (using makeblastdb). The master protein accession is used to search the database (using blastdbcmd), and if present the motif determined by the range calculated above is extracted and appended to the xlsx table.

If multiple master proteins or multiple modified amino acids are identified for a given row of input data, this function will separate that data to search a single master protein / modification at a time resulting in a greater number of rows in the data set.

### Function Contents
```
run_motif_extraction.R
build_blastdb.R
formatter.R
data_pull_and_fill.R
run_blastdbcmd.R
```

### Parameters
This program has 5 parameters.
1. The path to the data file
- An .xlsx peptidel isoforms output file from Proteome Discover (is this correct?) that must includes the following columns information:
  | Modifications_in_Master_Proteins |
  | ---- |
  | P49312 1xPhospho [S308(99.1)] |
  | A0A668KLC6 2xPhospho [T1811(100); S1815(99.2)] |
  | etc... |
__SXXX__ indicates the modified amino acid (in this case the phosphorylated amino acid)
2. The user determined path for the output file
3. The path to ncbi-blast-2.14.0+/bin
- download the [blast+ executables](https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/)
4. The path to Fasta reference file for desired organism [e.g. mouse proteome download](https://rest.uniprot.org/uniprotkb/stream?download=true&format=fasta&query=%28mouse%29%20AND%20%28model_organism%3A10090%29))\
5. The number of amino acids on either side of the modified amino acid to extract from the Master protein (default is num=7)

### The output  
- an .csv file with additional columns
   - accession
   - start
   - end
   - range
   - extracted_sequence

### Example usage
```
#path to the protein sequence for the reference database
ref_protein <- "/Users/eguswa/Documents/Open_UVM_Projects/Proteomics_core/Mouse_ref/protein.faa"

#path to bin in blast suit
ncbi_bin <- "/Applications/ncbi-blast-2.14.0+/bin"

#input data sheet
input_path <- "/Users/eguswa/Documents/Open_UVM_Projects/Proteomics_core/Peptide_extraction/May222023_Eclipse_2022_138_unnormalized_peptideIsoforms_Bioinformatics.xlsx"

#output path
out_dir_path <- "/Users/eguswa/Documents/Open_UVM_Projects/Proteomics_core/Peptide_extraction/May222023_Eclipse_2022_138_with_extracted_sequence.csv"


run_motif_extraction(input_path, output_path, ncbi_bin, ref_protein)
```