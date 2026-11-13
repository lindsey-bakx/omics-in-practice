# BRB-seq analysis: K562 drug treatment time course

Omics in Practice, Module 2. Group [13]: [Lindsey Bakx]

## Data
- BRB-seq, human K562 cells, 8 timepoints (0h, 6h, 1d, 2d, 3d, 4d, 7d, 8d), 3 replicates each
- Reference genome: GRCh38 (STAR index provided by the course)

## Scripts
| Script | What it does |
|---|---|
| `run_starsolo.sh` | Aligns R2 (cDNA) to the genome and counts UMIs per gene per sample with STARsolo 2.7.9a |
| `convert_matrix.R` | Converts the STARsolo sparse matrix to a gene x sample count table |
| `brbseq_deseq2_analysis.R` | PCA, heatmaps, DESeq2 (each timepoint vs 0h) and volcano plots |

## Main results
- 71.8% of reads mapped uniquely; 95.6% had a valid sample barcode
- Samples cluster by timepoint in the PCA; replicates are highly consistent
- The drug induces an erythroid gene programme (haemoglobin genes, ALAS2, GYPB)
- Gene of interest TBP is not differentially expressed at any timepoint
