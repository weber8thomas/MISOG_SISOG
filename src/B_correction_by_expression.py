import math
import os
import pandas as pd
import gzip

pd.options.mode.chained_assignment = None  # default='warn'
import multiprocessing
import parmap
import numpy as np
import collections
from tqdm import tqdm


tqdm.pandas()
from pprint import pprint
import requests
import re
import json
import subprocess
from pandarallel import pandarallel

pandarallel.initialize(nb_workers=30, progress_bar=True)


## YAML FILES CONFIG
import yaml

yaml = yaml.load(open("/home/weber/PycharmProjects/MISOG_SISOG/src/config/config_files.yaml"), Loader=yaml.FullLoader)


class CorrectExpression:
    def __init__(self):

        # if os.path.isfile(path) is False:
        # * 1 Correction by expression
        gtex_corrected = self.launch(yaml["EXPRESSION"]["Final"]["transcript_check"])

        # # * 2 Load RefSeq & remove genes with transcript identical on CDS
        # refseq = pd.read_parquet(yaml["1_GENOMICS"]["Final"]["refseq_cds_with_variable"])
        # utrs = pd.read_parquet(yaml["1_GENOMICS"]["Final"]["refseq_miso_utrs"])
        # refseq = refseq.loc[~refseq["Gene"].isin(utrs.loc[utrs["Nb_combi"] == 1])]

        # # * 3 Remove unexpressed transcripts & recompute frequencies & usage
        # refseq_correct_mrnas = self.correct_refseq_by_expression(
        #     yaml["EXPRESSION"]["Final"]["refseq_corrected_cds_with_variable_test"], refseq, gtex_corrected
        # )
        # tmp = refseq_correct_mrnas.explode("mRNA_exons")
        # print(tmp.Gene.nunique())
        # print(tmp.mRNA_exons.nunique())

        # self.final_corrected_refseq = self.process_corrected_file(path, refseq_correct_mrnas)
        # tmp = self.final_corrected_refseq.explode("mRNA_exons")
        # print(tmp.Gene.nunique())
        # print(tmp.mRNA_exons.nunique())

        # else:
        #     self.final_corrected_refseq = pd.read_parquet(path)
        #     tmp = self.final_corrected_refseq.explode("mRNA_exons")
        #     print(tmp.Gene.nunique())
        #     print(tmp.mRNA_exons.nunique())
        #     print(self.final_corrected_refseq)

    def launch(self, path):

        print("### Launch pipeline correction method / File = {}".format(path))

        if os.path.isfile(path) is False:
            print("# Files don't exist ☒")
            # MAIN FCT
            merge_check = self.merge_check_files(
                path_gtex_reads=yaml["EXPRESSION"]["base_directory"] + yaml["EXPRESSION"]["External"]["transcript_reads"],
                path_gtex_tpm=yaml["EXPRESSION"]["base_directory"] + yaml["EXPRESSION"]["External"]["transcript_tpm"],
            )
            merge_check["transcript_id"] = merge_check["transcript_id"].apply(lambda r: r.split(".")[0])

            # OUTPUT
            merge_check.to_csv(path, compression="gzip", sep="\t")
        else:
            print("# Files exist ✓, Loading ... ")

            merge_check = pd.read_csv(path, compression="gzip", sep="\t")
        print(merge_check)
        return merge_check

    def check_file(self, filepath):

        print("--- Checking GTEx raw files ---")

        # CHECK IF PROCESSING EITHER TPM OR READS
        check_type = "tpm" if "tpm" in filepath else "reads"
        print(check_type)

        o_file = filepath.replace(".gct.gz", "_checked_complete.parquet")
        o_file_lite = filepath.replace(".gct.gz", "_checked_lite.gct.gz")
        print(o_file)

        if os.path.isfile(o_file) is False:
            print("# Files don't exist ☒")

            final_df_list = list()
            for df in tqdm(pd.read_csv(filepath, compression="gzip", sep="\t", skiprows=2, chunksize=1000)):

                # # LOAD FILE
                # df = pd.read_csv(
                #     filepath,
                #     compression="gzip",
                #     sep="\t",
                #     skiprows=2,
                #     # nrows=1000,
                # )

                df["gene_id"] = df["gene_id"].apply(lambda r: r.split(".")[0])

                # MERGE WITH BIOMART TO HAVE GENE SYMBOLS
                # biomart_grch37 = self.load_biomart(biomart_path)
                # df = pd.merge(biomart_grch37, df, on="gene_id")

                # CUTOFFS
                cutoff_reads = 6
                cutoff_tpm = 0.1

                # CHECK FUNCTION
                def check(r, check_type):
                    i = 0
                    if check_type == "reads":
                        cutoff = cutoff_reads
                    elif check_type == "tpm":
                        cutoff = cutoff_tpm
                    for e in r[3:]:
                        if e >= cutoff:
                            i += 1
                    return i

                cols_nb = len(df.columns[2:])

                # PARALLEL APPLY TO CHECK CUTOFF
                df["check_{}_below_cutoff".format(check_type)] = df.apply(lambda r: check(r, check_type), axis=1)
                df["check_{}".format(check_type)] = df["check_{}_below_cutoff".format(check_type)].apply(
                    lambda r: True if r / cols_nb >= 0.2 else False
                )
                # print(df)

                final_df_list.append(df[list(df.columns)[:2] + list(df.columns)[-2:]])

            # OUTPUT COMPLETE & LITE
            # df.to_parquet(o_file, index=False)
            final_df = pd.concat(final_df_list)
            # df = df[list(df.columns)[:2] + list(df.columns)[-2:]]
            print(final_df)
            # exit()
            final_df.to_csv(o_file_lite, index=False, compression="gzip", sep="\t")
        else:

            print("# Files exist ✓")

            final_df = pd.read_csv(o_file_lite, compression="gzip", sep="\t")
            print(final_df)
        return final_df

    def merge_check_files(self, path_gtex_reads, path_gtex_tpm):
        if os.path.isfile(yaml["EXPRESSION"]["Final"]["transcript_check"]) is False:

            # CHECK GTEx FILES TO RETRIEVE TRANSCRIPTS PASSING CUTOFFS
            check_reads_file = self.check_file(path_gtex_reads)
            check_tpm_file = self.check_file(path_gtex_tpm)

            print(check_reads_file)
            print(check_tpm_file)

            # MERGE FILES
            merge_check = pd.merge(
                check_reads_file,
                check_tpm_file,
                on=[
                    "gene_id",
                    "transcript_id",
                ],
            )

            # COMPUTE COLUMN TO RETRIEVE TRANSCRIPTS PASSING BOTH CONDITIONS
            merge_check["check"] = merge_check.apply(
                lambda r: True if set([r["check_reads"], r["check_tpm"]]) == {True} else False,
                axis=1,
            )

            # OUTPUT
            merge_check.to_csv(
                yaml["EXPRESSION"]["Final"]["transcript_check"],
                compression="gzip",
                sep="\t",
                index=False,
            )
        else:
            merge_check = pd.read_csv(
                yaml["EXPRESSION"]["Final"]["transcript_check"],
                compression="gzip",
                sep="\t",
            )
        print(merge_check)
        return merge_check


if __name__ == "__main__":
    c = CorrectExpression()
    # df = c.final_corrected_refseq
    # print(df.loc[df["Share"] == True])
