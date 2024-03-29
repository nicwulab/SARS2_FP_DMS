## Studying the mutational fitness effects of SARS-CoV-2 spike fusion peptide

### Dependencies
* [Python](https://www.python.org/) (version 3.9)
* [Biopython](https://github.com/biopython/biopython)
* [Seaborn](https://seaborn.pydata.org/)
* [Matplotlib](https://matplotlib.org/)
* [R](https://www.r-project.org/) (version 4.1)
* [PEAR](https://github.com/tseemann/PEAR)
* [PyMOL](https://pymol.org/)

### Input files
* Raw read files in fastq format from NIH SRA database [BioProject PRJNA910585](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA910585)
* [./Fasta/SARS2-FP_with_flank.fa](./Fasta/SARS2-FP_with_flank.fa): Fusion peptide sequences with 21 nt upstream (5' flank)
* [./Fasta/FP_ref.fa](./Fasta/FP_ref.fa): Reference (i.e. wild type) amino acid sequences (primer regions not included)
* [./data/Tree_aa_fitness.csv](./data/Tree_aa_fitness.csv): S mutational fitness estimated from phylogeny by [Bloom & Neher](https://www.biorxiv.org/content/10.1101/2023.01.30.526314v1). The orignal data are from [here](https://github.com/jbloomlab/SARS2-mut-fitness/blob/main/results/aa_fitness/aa_fitness.csv).
* [./data/BA1_DMS_muteffects_observed.csv](./data/BA1_DMS_muteffects_observed.csv): S mutational fitness measured by a pseudovirus system by [Dadonaite et al](https://www.cell.com/cell/fulltext/S0092-8674(23)00103-4). The orignal data are from [here](https://github.com/dms-vep/SARS-CoV-2_Omicron_BA.1_spike_DMS_mAbs/blob/main/results/muteffects_functional/muteffects_observed.csv).

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

4. Convert nucleotide sequences to codon variants   
``python3 script/FP_count_nuc2codon.py``   
    - Input files:
      - [./data/barcodes.tsv](./data/barcodes.tsv)
      - [./Fasta/FP_ref.fa](./Fasta/FP_ref.fa)
      - result/FP_DMS_count_nuc.tsv
    - Output files:
      - [./result/FP_DMS_count_codon.tsv](./result/FP_DMS_count_codon.tsv)

5. Compute fitness   
``python3 script/FP_count2fit.py``   
    - Input files:
      - [./result/FP_DMS_count_aa.tsv](./result/FP_DMS_count_aa.tsv)
    - Output file:
      - [./result/FP_DMS_fit.tsv](./result/FP_DMS_fit.tsv)
      - [./result/FP_DMS_fit_by_resi.tsv](./result/FP_DMS_fit_by_resi.tsv)

6. Convert B factor in PDB file into mean fitness value   
``python3 script/convert_Bfactor_to_fit.py``   
    - Input file:
      - [./PDB/7my8.pdb](./PDB/7my8.pdb)
      - [./result/FP_DMS_fit_by_resi.tsv](./result/FP_DMS_fit_by_resi.tsv)
    - Output file: 
      - [./PDB/7my8_fit.pdb](./PDB/7my8_fit.pdb)

### Plotting
1. Plot correlation between replicates and compare silent/missense/nonsense   
``Rscript script/plot_QC.R``   
    - Input file:
      - [./result/FP_DMS_fit.tsv](./result/NTD_DMS_fit.tsv)
    - Output files:
      - graph/QC_*.png

2. Plot correlation between fitness measurements in this study and those in previous studies
``python3 script/plot_cor_measures.py``
    - Input file: 
      - [./result/FP_DMS_fit.tsv](./result/FP_DMS_fit.tsv)
      - [./data/BA1_DMS_muteffects_observed.csv](./data/BA1_DMS_muteffects_observed.csv)
      - [./data/Tree_aa_fitness.csv](./data/Tree_aa_fitness.csv)
    - Output files:
      - graph/cor_*.png

3. Plot heatmap for the fitnss of individual mutations   
``Rscript script/plot_heatmap_fit.R``   
    - Input file:
      - [./result/FP_DMS_fit.tsv](./result/FP_DMS_fit.tsv)
    - Ouput file:
      - graph/FP_fit_heatmap_*.png

4. Plot mean fitness (i.e. mutational tolerance) of individual residue positions   
``Rscript script/plot_mean_fit.R``   
    - Input file:
      - [./result/FP_DMS_fit_by_resi.tsv](./result/FP_DMS_fit_by_resi.tsv)
    - Output file:
      - [./graph/mean_fit.png](./graph/mean_fit.png)

5. Plot meanfitness on structure   
``pymol script/plot_Bfactor_as_fit.pml``
    - Input file:
      - [./PDB/7my8_fit.pdb](./PDB/7my8_fit.pdb)
    - Output file:
      - [./graph/structure_FP.png](./graph/structure_FP.png)

6. Plot heatmap for the antibody escape of individual mutations   
``Rscript script/plot_heatmap_escape.R``   
    - Input file:
      - [./result/FP_DMS_fit.tsv](./result/FP_DMS_fit.tsv)
    - Ouput file:
      - graph/FP_escape_*.png

7. Plot heatmap for the codon variants   
``Rscript script/plot_heatmap_codon_freq.R``   
    - Input file:
      - [./result/FP_DMS_codon_freq.tsv](./result/FP_DMS_codon_freq.tsv)
    - Ouput file:
      - graph/FP_codon_freq_*.png

### Analyze mutation rate
1. Identity the number of amino acid mutations on each merged read    
``python3 script/analyze_mut_rate.py``   
    - Input file:
      - Merged read files in fastq_merged/ folder
      - [./Fasta/FP_ref.fa](./Fasta/FP_ref.fa)
    - Output file:
      - [./result/Lib1_mut_count.tsv](./result/Lib1_mut_count.tsv)
      - [./result/Lib2_mut_count.tsv](./result/Lib2_mut_count.tsv)

2. Plot mutation rate   
``Rscript script/plot_lib_mut_count.R``   
    - Input file:
      - [./result/Lib1_mut_count.tsv](./result/Lib1_mut_count.tsv)
      - [./result/Lib2_mut_count.tsv](./result/Lib2_mut_count.tsv)
    - Output file:
      - [./graph/Mut_count_input_lib.png](./graph/Mut_count_input_lib.png)
