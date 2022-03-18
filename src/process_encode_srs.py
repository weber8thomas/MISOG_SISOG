import pandas as pd
from tqdm import tqdm

tqdm.pandas()

# Load concat file SRS
df = pd.read_csv("/gstock/biolo_datasets/ENCODE/ENCODE_SRS_concat.tsv.gz", compression="gzip", sep="\t")

# Process files
df[["ENST", "ENSG", "VEGAT", "VEGAG", "transcript_id", "GeneID", "transcript_length", "transcript_biotype", ""]] = df[
    "target_id"
].str.split("|", expand=True)
df = df.loc[df["ENSG"].isna() != True]

df["ENST"] = df["ENST"].apply(lambda r: r.split(".")[0])
df["ENSG"] = df["ENSG"].apply(lambda r: r.split(".")[0])


def filter_tpm(tpm):
    return tpm.loc[tpm > 0.1].shape[0] / tpm.shape[0]


tpm_transcripts_ratio = df.groupby("ENST")["tpm"].progress_apply(filter_tpm)
print(tpm_transcripts_ratio)

tpm_transcripts_ratio = tpm_transcripts_ratio.rename("TPM_ratio").reset_index()
tpm_transcripts_ratio.to_csv(
    "/gstock/biolo_datasets/ENCODE/ENCODE_SRS_TPM_summary.tsv.gz", compression="gzip", sep="\t", index=False
)

# df = pd.merge(df, tpm_transcripts_ratio.rename('TPM_ratio').reset_index(), on='ENST')
# df