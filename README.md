# MISOG / SISOG analysis

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

`conda create -c conda-forge -n misogsisog python=3.7 pandas=1.1.4 scipy=1.6.1 numpy=1.18.1 tqdm=4.48.2 yaml=0.1.7 matplotlib=3.2.1 seaborn=0.11.2 statannot=0.2.3 jupyterlab=3.0.9` 

## External ENSEMBL files were retrieved with BIOMART 

TODO

## Steps 

Human data processing
1. src/GTEx_correction_by_expression.py
2. src/misog_sisog_analysis_ensembl_gencode.ipynb
3. src/Stats_and_plots.ipynb


Mouse data processing

3. src/download_encode_files.py
4. src/concat_files_ENCODE.py
5. src/process_encode_srs.py


Orthology analysis

6. src/orthologs.ipynb
7. src/Stats_ortho.ipynb
8. src/Visualization_ortho.ipynb

Enrichment analysis

9. src/clustering_enrichment_analysis.ipynb

