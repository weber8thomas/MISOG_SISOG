# General imports
import os
import sys
import pandas as pd
import scipy
pd.options.mode.chained_assignment = None  # default='warn'
import numpy as np
from tqdm import tqdm
sys.path.append('../')
# Other imports
tqdm.pandas()
import yaml
import json

## YAML FILES CONFIG
yaml = yaml.load(open("config/config_files.yaml"), Loader=yaml.FullLoader)
base_dir = yaml['base_directory']
sys.exit('EXIT : Need to specify the base_directory in config file : "conf_files.yaml"') if base_dir == 'TO_CHANGE' else None

genes_raw = pd.read_csv(base_dir + yaml['External']['biomart_ensembl'], compression='gzip', sep='\t')
genes = genes_raw[['Gene stable ID', 'Gene name', 'Chromosome/scaffold name', 'Gene start (bp)', 'Gene end (bp)']].drop_duplicates().sort_values(by='Gene name')
genes['Gene_length'] = genes['Gene end (bp)'] - genes['Gene start (bp)']

# genes['Gene'] = genes['Attributes'].apply(lambda r: [e.replace('ID=gene-', '') for e in r.split(';') if 'ID=' in e][0])
# genes.loc[genes['Gene name'].duplicated(keep=False) == True].head(30)


genes_raw['Chromosome/scaffold name'] = genes_raw['Chromosome/scaffold name'].astype(str)
mrna = genes_raw.loc[~genes_raw['Chromosome/scaffold name'].str.contains('CHR|\.')].copy()
mrna = mrna[['Gene stable ID', 'Transcript stable ID', 'Protein stable ID', 'Transcript start (bp)', 'Transcript end (bp)', 'Transcription start site (TSS)', 'Transcript length (including UTRs and CDS)', 'Transcript support level (TSL)', 'GENCODE basic annotation', 'APPRIS annotation']]

# ! TODO: add script to handle GTEx

gtex_checked = pd.read_csv('/gstock/EXOTIC/data/EXPRESSION/GTEx_V7_transcript_checking.csv.gz', compression='gzip', sep='\t')
gtex_checked['transcript_id'] = gtex_checked['transcript_id'].apply(lambda r : r.split('.')[0])

mrna = mrna.loc[mrna['Transcript stable ID'].isin(gtex_checked.loc[gtex_checked['check'] == True]['transcript_id'].values.tolist())]


miso_siso = mrna[['Gene stable ID', 'Transcript stable ID']].groupby('Gene stable ID')['Transcript stable ID'].nunique().reset_index()
miso_siso.loc[miso_siso['Transcript stable ID'] > 1, 'Miso_siso'] = 'Miso'
miso_siso.loc[miso_siso['Transcript stable ID'] == 1, 'Miso_siso'] = 'Siso'
miso_siso = miso_siso.rename({'Gene stable ID': 'GeneID', 'Transcript stable ID' : 'transcript_count'}, axis=1)
miso_siso.groupby('Miso_siso')['transcript_count'].sum()

merge_raw_gtex = pd.merge(miso_siso, miso_siso_raw, on=['GeneID'])
merge_raw_gtex.loc[merge_raw_gtex['Miso_siso'] == merge_raw_gtex['Raw_Miso_siso'], 'Same'] = True
merge_raw_gtex.loc[merge_raw_gtex['Miso_siso'] != merge_raw_gtex['Raw_Miso_siso'], 'Same'] = False
genes_to_keep = merge_raw_gtex.loc[merge_raw_gtex['Same'] == True]['GeneID'].values.tolist()
genes_to_not_keep = merge_raw_gtex.loc[merge_raw_gtex['Same'] == False]['GeneID'].values.tolist()
