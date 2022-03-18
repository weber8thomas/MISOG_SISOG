from tqdm import tqdm
import pandas as pd
import sys, os
import collections


"""
Small script to concat ENCODE files into a single dataframe to process it easily
5 cols = SRS sequencing
12 cols = LRS sequencing
"""

encode_dl_directory = "/gstock/biolo_datasets/ENCODE/DL/"

dict_df = collections.defaultdict(list)

for file in tqdm(os.listdir(encode_dl_directory)):
    # print(file)
    cols = open(encode_dl_directory + file, "r").readline().strip().split("\t")
    cols_length = len(cols)
    if cols_length == 5:
        dict_df[cols_length].append(pd.read_csv(encode_dl_directory + file, sep="\t"))
    elif cols_length == 12:
        dict_df[cols_length].append(pd.read_csv(encode_dl_directory + file, sep="\t").drop([cols[-1]], axis=1))

if dict_df[5]:
    pd.concat(dict_df[5]).to_csv(
        "/gstock/biolo_datasets/ENCODE/ENCODE_SRS_concat.tsv.gz", compression="gzip", sep="\t", index=False
    )
if dict_df[12]:
    pd.concat(dict_df[12]).to_csv(
        "/gstock/biolo_datasets/ENCODE/ENCODE_LRS_concat.tsv.gz", compression="gzip", sep="\t", index=False
    )
