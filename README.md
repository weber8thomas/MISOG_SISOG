# Do 5’ regions of human protein-coding genes contain the blueprints for alternative splicing?

In this repository, you will find all the developped scripts and tools used to reperformed the multiple transcript isoform genes (MISOG) VS single transcript isoform genes (SISOG) analysis.

## Environment & requirements

All code was executed under Python 3.7 version.

Following packages are required to run `prepare_refseq.py` and `misog_sisog_analysis.ipynb` files : 
- pandas=1.1.4
- scipy=1.6.1
- numpy=1.18.1
- tqdm=4.48.2
- pyyaml=0.1.7
- matplotlib=3.2.1
- seaborn=0.11.2
- statannot=0.2.3
- jupyterlab=3.0.9


Conda environment can be create with the following : 


## External RefSeq GFF file 

RefSeq used GFF file (GCF_000001405.39) can be retrieved [here](https://ftp.ncbi.nlm.nih.gov/refseq/H_sapiens/annotation/annotation_releases/109.20210514/GCF_000001405.39_GRCh38.p13/GCF_000001405.39_GRCh38.p13_genomic.gff.gz).

## How to rerun analysis ? 

- 0. Create conda environment with : `conda create -c conda-forge -n misogsisog python=3.7 pandas=1.1.4 scipy=1.6.1 numpy=1.18.1 tqdm=4.48.2 yaml=0.1.7 matplotlib=3.2.1 seaborn=0.11.2 statannot=0.2.3 jupyterlab=3.0.9` 
- 1. Activate conda environment with the following : `conda activate misogsisog`
- 2. Run the `python prepare_refseq.py` python file to prepare temporary files for each biological concept (chromosome, gene, mRNA, exon, CDS)
- 3. Run the `misog_sisog_analysis.ipynb` notebook 

## About the notebook

The notebook is built to retrieve each XLSX table, each CSV/Apache parquet file and to produce the visualizations present in the manuscript

## Figures & Supplementary Figures

![Figure1](./Figures/Figure1.png)
**Figure 1 – MISOG and SISOG archetypal gene architecture.**

Major values discussed in the manuscript are summarized graphically (all values are available in the supplementary Tables). The median lengths of each gene component (exon, TER, intron, UTR) are described according to the values observed for the MISOG (red; median number of 3 mRNAs as transcript isoforms) and SISOG (blue). The major MISOG and SISOG discriminant elements (first 5’ introns, 5’ exons and 5’ TER) are indicated by colored labels. Non-discriminant elements are indicated by grey labels.
Gene lengths was calculated from TSS (Transcription Start Site) to TTS (Transcription Termination Site). TER stands for Translated Exonic Regions, the coding part for amino acids of exons.  

![Figure2](./Figures/Figure2.png)
**Figure 2 – Comparative analysis of MISOG and SISOG 5’ and 3’ non-coding regions.**

(A-B) barplots indicate the percentage of MISOG/SISOG transcripts exhibiting from one up to seven 5’ or 3’ UTR exons. (C) Comparison of the MISOG versus SISOG median distances separating TSS from START codon (left) and STOP codon from TTS (right) (all values are available in Table S9).
(A-B) represent the percentage of transcript according to the number of exons in 5’ (A) or 3’ UTR(B).

![Figure S1](./Figures/FigS1.png)
**Figure S1 – Length distribution of gene elements (Intron, Exon, TER and UTR exons) according to ordinal positions.**

All plots represent the length distribution of gene elements according to ordinal position (as illustrated in Fig. S4F). Plots A-D correspond to the first five 5’ elements in the genes (1 up to 5), plots E-H correspond to the last five 3’ elements (-5 to -1). The p-values obtained through Mann-Whitney U tests are displayed above the violin plots. Black p-values correspond to comparisons between MISOG and SISOG while colored p-values represent comparisons between first elements and following ones (MISOG: red, SISOG: blue). Values under each plot represent the numbers of observed elements at the ordinal position for MISOG (red) and SISOG (blue). 

![Figure S2](./Figures/FigS2.png)
**Figure S2 – Cumulative length and number of 5’ and 3’ UTR exons by gene.**

Cumulative total lengths of 5’ (A) and 3’ (B) UTR regions per gene. Number of 5’ and 3’ UTR exons per gene are displayed in (C, D). The p-values obtained through Mann-Whitney U tests are displayed above the boxplots (A, B) and boxenplots (C, D). 

![Figure S3](./Figures/FigS3.png)
**Figure S3 – Percentage of alternative exons and TERs according to their ordinal positions.**

All plots show calculated frequencies of alternative exons (A, B) and alternative TERs (C, D) as illustrated in Fig. S4G. Plots (A, C) correspond to the first five 5’ elements in the genes (1 up to 5), while plots (B,D) correspond to the last five 3’ elements (-5 to -1). The p-values obtained through Mann-Whitney U tests are displayed above the violin plots. 

![Figure S4](./Figures/FigS4.png)
**Figure S4 – Material and Methods schemas**

## How to cite ?

Coming ...
