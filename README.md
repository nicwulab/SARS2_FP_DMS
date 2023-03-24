## Studying the mutational fitness effects of SARS-CoV-2 spike fusion peptide

### Dependencies
* [Python](https://www.python.org/) (version 3.9)
* [Biopython](https://github.com/biopython/biopython)
* [Seaborn](https://seaborn.pydata.org/)
* [Matplotlib](https://matplotlib.org/)
* [R](https://www.r-project.org/) (version 4.1)
* [PEAR](https://github.com/tseemann/PEAR)

### Input files
* Raw read files in fastq format from NIH SRA database [BioProject PRJNA910585](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA910585)
* [./Fasta/SARS2-FP_with_flank.fa](./Fasta/SARS2-FP_with_flank.fa): Fusion peptide sequences with 21 nt upstream (5' flank)
* [./Fasta/FP_ref.fa](./Fasta/FP_ref.fa): Reference (i.e. wild type) amino acid sequences (primer regions not included)
* [./data/Tree_aa_fitness.csv](./data/Tree_aa_fitness.csv): S mutational fitness estimated from phylogeny by [Bloom & Neher](https://www.biorxiv.org/content/10.1101/2023.01.30.526314v1). The orignal data are from [here](https://github.com/jbloomlab/SARS2-mut-fitness/blob/main/results/aa_fitness/aa_fitness.csv).
* [./data/BA1_DMS_muteffects_observed.csv](./data/BA1_DMS_muteffects_observed.csv): S mutational fitness measured by a pseudovirus system by [Dadonaite et al.](https://www.cell.com/cell/fulltext/S0092-8674(23)00103-4). The orignal data are from [here](https://github.com/dms-vep/SARS-CoV-2_Omicron_BA.1_spike_DMS_mAbs/blob/main/results/muteffects_functional/muteffects_observed.csv)

### Primer design for DMS library construction
1. Generating foward (NNK + internal barcode) and reverse primers (constant)   
``python3 script/lib_primer_design.py``
    - Input file:
      - [./Fasta/SARS2-FP_with_flank.fa](./Fasta/SARS2-FP_with_flank.fa)
    - Output files:
      - [./primer/SARS2-FP_lib_Fprimer_bc.fa](./primer/SARS2-FP_lib_Fprimer_bc.fa)
      - [./primer/SARS2-FP_lib_Rprimer_bc.fa](./primer/SARS2-FP_lib_Rprimer_bc.fa)

2. Generating barcode file   
``python3 script/check_barcode.py``
    - Input files:
      - [./primer/SARS2-FP_lib_Fprimer_bc.fa](./primer/SARS2-FP_lib_Fprimer_bc.fa)
      - [./Fasta/SARS2-FP_with_flank.fa](./Fasta/SARS2-FP_with_flank.fa)
    - Output file:
      - [./data/barcodes.tsv](./data/barcodes.tsv)

### Calculating fitness from DMS data
1. Merge overlapping paired-end reads using [PEAR](https://github.com/tseemann/PEAR)   
``pear -f [FASTQ FILE FOR FORWARD READ] -r [FASTQ FILE FOR FORWARD READ] -o [OUTPUT FASTQ FILE]``   
    - Output files should be placed in the fastq_merged/ folder and named as described in [./doc/filename_merged_fastq.tsv](./doc/filename_merged_fastq.tsv)

2. Counting variants based on nucleotide sequences   
``python3 script/FP_fastq2count.py``   
    - Input files:
      - Merged read files in fastq_merged/ folder
    - Output files:
      - result/FP_DMS_count_nuc.tsv

3. Convert nucleotide sequences to amino acid mutations   
``python3 script/FP_count_nuc2aa.py``   
    - Input files:
      - [./data/barcodes.tsv](./data/barcodes.tsv)
      - [./Fasta/FP_ref.fa](./Fasta/FP_ref.fa)
      - result/FP_DMS_count_nuc.tsv
    - Output files:
      - [./result/FP_DMS_count_aa.tsv](./result/FP_DMS_count_aa.tsv)

4. Compute fitness   
``python3 script/FP_count2fit.py``   
    - Input files:
      - [./result/FP_DMS_count_aa.tsv](./result/FP_DMS_count_aa.tsv)
    - Output file:
      - [./result/FP_DMS_fit.tsv](./result/FP_DMS_fit.tsv)
      - [./result/FP_DMS_fit_by_resi.tsv](./result/FP_DMS_fit_by_resi.tsv)

5. Plot correlation between replicates and compare silent/missense/nonsense   
``Rscript script/plot_QC.R``   
    - Input file:
      - [./result/FP_DMS_fit.tsv](./result/NTD_DMS_fit.tsv)
    - Output files:
      - graph/QC_*.png

6. Plot correlation between fitness measurements in this study and those in previous studies
``python3 script/plot_cor_measures.py``
    - Input file: 
      - [./result/FP_DMS_fit.tsv](./result/FP_DMS_fit.tsv)
      - [./data/BA1_DMS_muteffects_observed.csv]
      - [./data/Tree_aa_fitness.csv]
    - Output files:
      - [./graph/cor_BA1_vs_Calu3.png](./graph/cor_BA1_vs_Calu3.png)
      - [./graph/cor_BA1_vs_E6.png](./graph/cor_BA1_vs_E6.png)
      - [./graph/cor_tree_vs_Calu3.png](./graph/cor_tree_vs_Calu3.png)
      - [./graph/cor_tree_vs_E6.png](./graph/cor_tree_vs_E6.png)

7. Plot heatmap for the fitnss of individual mutations   
``Rscript script/plot_heatmap_fit.R``   
    - Input file:
      - [./result/FP_DMS_fit.tsv](./result/FP_DMS_fit.tsv)
    - Ouput file:
      - graph/FP_fit_heatmap_*.png

8. Plot heatmap for the antibody escape of individual mutations   
``Rscript script/plot_heatmap_escape.R``   
    - Input file:
      - [./result/FP_DMS_fit.tsv](./result/FP_DMS_fit.tsv)
    - Ouput file:
      - graph/FP_escape_*.png
