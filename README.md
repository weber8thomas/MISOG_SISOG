# MISOG VS SISOG Analysis

In this repository, you will find all the developped scripts and tools used to reperformed the MISOG VS SISOG analysis.

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

## How to cite ?

Coming ...