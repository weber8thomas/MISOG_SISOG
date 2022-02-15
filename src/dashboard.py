from datetime import timedelta
from pathlib import Path
from time import sleep

import numpy as np
import pandas as pd
import plotly_express as px
import seaborn as sns
import matplotlib.pyplot as plt
import streamlit as st

st.set_page_config(layout="wide")


@st.cache
def load_data():
    exons = pd.read_parquet("/gstock/GeneIso/V2/Exons.parquet")
    cds = pd.read_parquet("/gstock/GeneIso/V2/CDS.parquet")
    five_UTR = pd.read_parquet("/gstock/GeneIso/V2/5_UTR.parquet")
    three_UTR = pd.read_parquet("/gstock/GeneIso/V2/3_UTR.parquet")
    introns = pd.read_parquet("/gstock/GeneIso/V2/Introns.parquet")
    return exons, cds, five_UTR, three_UTR, introns


"""
# Test
"""

exons, cds, five_UTR, three_UTR, introns = load_data()

exons = exons.loc[exons["Exon_nb"] <= 5].reset_index(drop=True)
exons_to_display = exons.head(50)
exons_to_display

cds = cds.loc[cds["CDS_nb"] <= 5].reset_index(drop=True)

introns = introns.loc[introns["Introns_nb"] <= 5].reset_index(drop=True)
introns["Length"] = introns["Introns"].apply(lambda r: int(r.split("-")[1]) - int(r.split("-")[0]))


miso = exons.loc[exons["Miso_siso"] == "Miso"].GeneID.unique().tolist()[:5]
siso = exons.loc[exons["Miso_siso"] == "Siso"].GeneID.unique().tolist()[:5]
option = st.multiselect("Gene ?", exons.GeneID.unique().tolist(), miso + siso)


col1, col2 = st.columns(2)

with col1:

    """
    # Data exploration
    ## Barplots
    """

    fig1 = plt.figure()
    sns.boxplot(
        data=exons.loc[exons["GeneID"].isin(option)].sort_values(by="Miso_siso"),
        x="Exon_nb",
        y="Exon_length",
        hue="Miso_siso",
        showfliers=False,
    )

    fig1


with col2:

    """
    # Data2 exploration
    ## Barplots2
    """

    fig2 = plt.figure()
    sns.boxplot(
        data=introns.loc[introns["GeneID"].isin(option)].sort_values(by="Miso_siso"),
        x="Introns_nb",
        y="Length",
        hue="Miso_siso",
        showfliers=False,
    )

    fig2

st.columns(1)


"""
# Data3 exploration
## Barplots2
"""

stats = (
    exons.drop_duplicates(subset=["GeneID", "Chromosome/scaffold name", "Exon region start (bp)", "Exon region end (bp)"])
    .groupby("Miso_siso")["Exon_length"]
    .describe()
)
stats
