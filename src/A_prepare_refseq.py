# IMPORTS

import math
import os
import sys
import pandas as pd

# sys.path.append("/home/weber/PycharmProjects/EXOTIC/src")
sys.path.append("/home/weber/PycharmProjects/EXOTIC/clean/src")

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
import seaborn as sns
import json
from utils import utils

from pandarallel import pandarallel

pandarallel.initialize(nb_workers=20, progress_bar=True)

## YAML FILES CONFIG
yaml = utils.load_config_file(config_file="/home/weber/PycharmProjects/gene_isoforms/src/config/config_files.yaml")
base_directory = yaml["base_directory"]


class ProcessRefSeq:
    def __init__(self, path):
        """[Main function to launch steps]

        Arguments:
            path {[str]} -- [Output file path]

        Returns:
            [pd.DataFrame] -- [Final processed refseq dataframe]
        """

        self.base_directory = yaml["base_directory"]

        if os.path.isfile(path) is True:
            print("# Files don't exist ☒")

            # * 0 Load raw file
            refseq_gff = yaml["1_GENOMICS"]["External"]["raw_refseq"]
            refseq_df = pd.DataFrame()
            # refseq_df = self.load_refseq(refseq_gff)

            utils.mkdir(os.path.dirname(self.base_directory + yaml["1_GENOMICS"]["TMP"]["tmp_refseq_chroms"]))

            # * 1 Build tmp files by category
            refseq_df_chroms = self.refseq_chroms_fct(yaml["1_GENOMICS"]["TMP"]["tmp_refseq_chroms"], refseq_df)
            refseq_df_pc_genes = self.refseq_pc_genes_fct(yaml["1_GENOMICS"]["TMP"]["tmp_refseq_pc_genes"], refseq_df)
            refseq_df_mrnas = self.refseq_mrnas_fct(yaml["1_GENOMICS"]["TMP"]["tmp_refseq_mrnas"], refseq_df)
            refseq_df_exons = self.refseq_exons_fct(yaml["1_GENOMICS"]["TMP"]["tmp_refseq_exons"], refseq_df)
            refseq_df_cds = self.refseq_cds_fct(yaml["1_GENOMICS"]["TMP"]["tmp_refseq_cds"], refseq_df)

            # # * 2 Select exons(with UTRs) & coding exons
            # self.select_exons(
            #     yaml["1_GENOMICS"]["TMP"]["tmp_refseq_exons_filtered"],
            #     refseq_df_chroms,
            #     refseq_df_mrnas,
            #     refseq_df_pc_genes,
            #     refseq_df_exons,
            #     "Exon",
            # )

            self.select_exons(
                yaml["1_GENOMICS"]["TMP"]["tmp_refseq_cds_filtered"],
                refseq_df_chroms,
                refseq_df_mrnas,
                refseq_df_pc_genes,
                refseq_df_cds,
                "CDS",
            )
            exit()

            # # * 3 Concat tmp files
            concat_exons_cds = self.concat_exons_cds_fct(
                yaml["1_GENOMICS"]["TMP"]["tmp_refseq_concat_exons_cds_filtered"],
                yaml["1_GENOMICS"]["TMP"]["tmp_refseq_exons_filtered"],
                yaml["1_GENOMICS"]["TMP"]["tmp_refseq_cds_filtered"],
            )

            # # * 4 Groupby and produce final output file
            self.groupby_mrnas = self.groupby_mrnas_fct(yaml["1_GENOMICS"]["Final"]["refseq_processed"], concat_exons_cds)

            # # * 5 Groupby and produce final output file
            self.refseq_cds_with_variable = self.get_shared_exons(
                yaml["1_GENOMICS"]["Final"]["refseq_cds_with_variable"], self.groupby_mrnas
            )

            # # * 6 Groupby and produce final output file
            self.refseq_cds_without_variable = self.remove_variable_fct(
                yaml["1_GENOMICS"]["Final"]["refseq_cds_without_variable"], self.refseq_cds_with_variable
            )
        else:
            print("# Files exist ✓, Loading ... ")

            # print(pd.read_parquet(yaml["1_GENOMICS"]["TMP"]["tmp_refseq_pc_genes"]))
            # exit()

            # # * 3 Concat tmp files
            # concat_exons_cds = self.concat_exons_cds_fct(
            #     yaml["1_GENOMICS"]["TMP"]["tmp_refseq_concat_exons_cds_filtered"],
            #     yaml["1_GENOMICS"]["TMP"]["tmp_refseq_exons_filtered"],
            #     yaml["1_GENOMICS"]["TMP"]["tmp_refseq_cds_filtered"],
            # )

            # # * 4 Groupby and produce final output file
            # self.groupby_mrnas = self.groupby_mrnas_fct(yaml["1_GENOMICS"]["Final"]["refseq_processed"], concat_exons_cds)
            # print(self.groupby_mrnas)
            # exit()

            self.refseq_cds_without_variable = pd.read_parquet(
                self.base_directory + yaml["1_GENOMICS"]["Final"]["refseq_cds_without_variable"]
            )
            print(self.refseq_cds_without_variable)
            self.get_introns_ranges(yaml["1_GENOMICS"]["TMP"]["tmp_refseq_introns"], self.refseq_cds_without_variable)

    @staticmethod
    def load_refseq(path):
        """[Load RefSeq GFF file]

        Arguments:
            path {[str]} -- [Path to the GFF RefSeq file]

        Returns:
            [pd.DataFrame] -- [RefSeq GFF turned into pandas dataframe]
        """
        print("### Load RefSeq / File = {}".format(path))

        refseq_df = pd.read_csv(
            path,
            compression="gzip",
            sep="\t",
            skiprows=9,
            # nrows=10000,
            names=["NC", "RefSeq_validation", "Region_type", "Start", "End", "Score", "Strand", "Phase", "Attributes"],
        )
        refseq_df = refseq_df.dropna(subset=["Start", "End"])
        refseq_df["Start"] = refseq_df["Start"].astype(int)
        refseq_df["End"] = refseq_df["End"].astype(int)
        return refseq_df

    # @staticmethod
    def refseq_chroms_fct(self, path, refseq_df):
        """[Extract chromosomes from RefSeq GFF]

        Arguments:
            path {[str]} -- [Tmp output file]
            refseq_df {[pd.DataFrame]} -- [RefSeq GFF turned into pandas dataframe]

        Returns:
            [pd.DataFrame] -- [RefSeq chromosomes into pandas dataframe]
        """
        print("### Build temp file (chroms part) / File = {}".format(path))

        if os.path.isfile(self.base_directory + path) is False:
            print("# Files don't exist ☒")

            refseq_df_chroms = refseq_df.loc[refseq_df["Region_type"] == "region"]
            index_list = list(refseq_df_chroms.index)

            chroms = [(i, index_list[j + 1] - 1) for j, i in enumerate(index_list) if j < (len(index_list) - 1)]
            refseq_df_chroms = refseq_df_chroms.loc[
                (refseq_df_chroms["NC"].str.contains("NC")) & (refseq_df_chroms["RefSeq_validation"] == "RefSeq")
            ]
            refseq_df_chroms.to_parquet(self.base_directory + path)
        else:
            print("# Files exist ✓, Loading ... ")
            refseq_df_chroms = pd.read_parquet(self.base_directory + path)
        print(refseq_df_chroms)
        print(refseq_df_chroms.shape)
        return refseq_df_chroms

    # @staticmethod
    def refseq_pc_genes_fct(self, path, refseq_df):
        """[Extract protein coding genes from RefSeq GFF]

        Arguments:
            path {[str]} -- [Tmp output file]
            refseq_df {[pd.DataFrame]} -- [RefSeq GFF turned into pandas dataframe]

        Returns:
            [pd.DataFrame] -- [RefSeq protein coding genes into pandas dataframe]
        """
        print("### Build temp file (protein coding genes part) / File = {}".format(path))

        if os.path.isfile(self.base_directory + path) is False:
            print("# Files don't exist ☒")

            refseq_df_pc_genes = refseq_df.loc[
                (refseq_df["Attributes"].str.contains("gene_biotype=protein_coding"))
                & (refseq_df["NC"].str.contains("NC_"))
                # & (refseq_df["RefSeq_validation"].str.contains("BestRefSeq"))
            ]
            refseq_df_pc_genes.to_parquet(self.base_directory + path)
        else:
            print("# Files exist ✓, Loading ... ")
            refseq_df_pc_genes = pd.read_parquet(self.base_directory + path)
        print(refseq_df_pc_genes)
        print(refseq_df_pc_genes.shape)
        return refseq_df_pc_genes

    # @staticmethod
    def refseq_mrnas_fct(self, path, refseq_df):
        """[Extract mRNAs from RefSeq GFF]

        Arguments:
            path {[str]} -- [Tmp output file]
            refseq_df {[pd.DataFrame]} -- [RefSeq GFF turned into pandas dataframe]

        Returns:
            [pd.DataFrame] -- [RefSeq mRNAs into pandas dataframe]
        """
        print("### Build temp file (mRNAs part) / File = {}".format(path))

        if os.path.isfile(self.base_directory + path) is False:
            print("# Files don't exist ☒")

            refseq_df_mrna = refseq_df.loc[
                (refseq_df["Attributes"].str.contains("NM_")) & (refseq_df["Region_type"] == "mRNA") & (refseq_df["NC"].str.contains("NC_"))
            ]
            refseq_df_mrna.to_parquet(self.base_directory + path)
        else:
            print("# Files exist ✓, Loading ... ")
            refseq_df_mrna = pd.read_parquet(self.base_directory + path)
        print(refseq_df_mrna)
        print(refseq_df_mrna.shape)
        return refseq_df_mrna

    # @staticmethod
    def refseq_exons_fct(self, path, refseq_df):
        """[Extract exons including UTRs from RefSeq GFF]

        Arguments:
            path {[str]} -- [Tmp output file]
            refseq_df {[pd.DataFrame]} -- [RefSeq GFF turned into pandas dataframe]

        Returns:
            [pd.DataFrame] -- [RefSeq exons including UTRs into pandas dataframe]
        """
        print("### Build temp file (Exons (with UTRs) part) / File = {}".format(path))

        if os.path.isfile(self.base_directory + path) is False:
            print("# Files don't exist ☒")

            refseq_df_exons = refseq_df.loc[
                (refseq_df["Attributes"].str.contains("exon-NM"))
                & (refseq_df["Region_type"] == "exon")
                & (refseq_df["NC"].str.contains("NC_"))
            ]
            refseq_df_exons.to_parquet(self.base_directory + path)
        else:
            print("# Files exist ✓, Loading ... ")
            refseq_df_exons = pd.read_parquet(self.base_directory + path)
        print(refseq_df_exons)
        print(refseq_df_exons.shape)
        return refseq_df_exons

    # @staticmethod
    def refseq_cds_fct(self, path, refseq_df):
        """[Extract coding exons (CDS) from RefSeq GFF]

        Arguments:
            path {[str]} -- [Tmp output file]
            refseq_df {[pd.DataFrame]} -- [RefSeq GFF turned into pandas dataframe]

        Returns:
            [pd.DataFrame] -- [RefSeq coding exons (CDS) into pandas dataframe]
        """
        print("### Build temp file (coding exons part) / File = {}".format(path))

        if os.path.isfile(self.base_directory + path) is False:
            print("# Files don't exist ☒")

            refseq_df_cds = refseq_df.loc[
                (refseq_df["Attributes"].str.contains("NP_")) & (refseq_df["Region_type"] == "CDS") & (refseq_df["NC"].str.contains("NC_"))
            ]
            refseq_df_cds.to_parquet(self.base_directory + path)
        else:
            print("# Files exist ✓, Loading ... ")
            refseq_df_cds = pd.read_parquet(self.base_directory + path)
        print(refseq_df_cds)
        print(refseq_df_cds.shape)
        return refseq_df_cds

if __name__ == "__main__":

    c = ProcessRefSeq(base_directory + yaml["1_GENOMICS"]["Final"]["refseq_cds_without_variable"])
    # print(c.groupby_mrnas)
