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

## External files

### GTEx expression data

GTEx V8 human expression data files were recovered from : https://www.gtexportal.org/home/datasets


Following files were used : 
- [Transcript read counts](https://storage.googleapis.com/gtex_analysis_v8/rna_seq_data/GTEx_Analysis_2017-06-05_v8_RSEMv1.3.0_transcript_expected_count.gct.gz)
- [Transcript TPMs](https://storage.googleapis.com/gtex_analysis_v8/rna_seq_data/GTEx_Analysis_2017-06-05_v8_RSEMv1.3.0_transcript_tpm.gct.gz)

To get tissue information, use the following mapping [TSV file](https://storage.googleapis.com/gtex_analysis_v8/annotations/GTEx_Analysis_v8_Annotations_SampleAttributesDS.txt)

Location on CSTB server : `/gstock/biolo_datasets/variation/benchmark/Databases/GTEX/V8/RNA_Seq_data/`


### ENSEMBL files were retrieved with BIOMART 

#### Human level

Ensembl gene/mRNA/exons definitions were retrieved through [BIOMART](https://www.ensembl.org/biomart/martview) in 2 times.

1. First file for all human genes/mRNA with the following criteria : 
- **Database :** Ensembl Genes 105
- **Dataset :** Human genes (GRCh38.p13)
- **Filters:**
    - Gene type: protein_coding
    - Transcript type: protein_coding
    - GENCODE basic annotation: Only
- **Attributes: (Features category)**
  - GENE properties: 
    - Gene stable ID
    - Transcript stable ID
    - Protein stable ID
    - Gene description
    - Chromosome/scaffold name
    - Gene start (bp)
    - Gene end (bp)
    - Strand
    - Transcript start (bp)
    - Transcript end (bp)
    - Transcription start site (TSS)
    - Transcript length (including UTRs and CDS)
    - Transcript support level (TSL)
    - GENCODE basic annotation
    - APPRIS annotation
    - Ensembl Canonical
    - RefSeq match transcript (MANE Select)
    - RefSeq match transcript (MANE Plus Clinical)
    - Gene name
    - Source of gene name
    - Transcript name
    - Source of transcript name
    - Gene % GC content
    - Gene type
    - Transcript type
    - Source (gene)
    - Source (transcript)
  - EXTERNAL properties: 
    - CCDS ID *(External reference)*
    - UniProtKB isoform ID *(External reference)*

2. Second file for human exons 
- **Database :** Ensembl Genes 105
- **Dataset :** Human genes (GRCh38.p13)
- **Filters:**
    - Gene type: protein_coding
    - Transcript type: protein_coding
    - GENCODE basic annotation: Only
- **Attributes (Structure category)** 
  - GENE properties:
    - Gene stable ID
    - Transcript stable ID
    - Exon region start (bp)
    - Exon region end (bp)
    - 5' UTR start
    - 5' UTR end
    - 3' UTR start
    - 3' UTR end
    - CDS Length
    - Strand
    - Transcription start site (TSS)
    - Transcript length (including UTRs and CDS)
    - Transcript start (bp)
    - Transcript end (bp)
    - Gene start (bp)
    - Gene end (bp)
    - Chromosome/scaffold name
  - EXON properties : 
    - CDS start
    - CDS end
    - Constitutive exon
    - Exon rank in transcript
    - Start phase
    - End phase
    - cDNA coding start
    - cDNA coding end
    - Genomic coding start
    - Genomic coding end
    - Exon stable ID


File 1 (gene/mRNA level) local path : `/gstock/GeneIso/data/External/ENSEMBL/mart_export.txt.gz`

File 2 (exons level) local path : `/gstock/GeneIso/data/External/ENSEMBL/mart_export_exons.txt.gz`

#### Orthologs

Orthologs were retrieved in BIOMART with the following criteria : 

- **Database :** Ensembl Genes 105
- **Dataset :** Human genes (GRCh38.p13)
- **Filters:**
    - Gene type: protein_coding
    - Transcript type: protein_coding
    - GENCODE basic annotation: Only
- **Attributes (Homologues category)**
  - GENE properties:
      - Gene stable ID
      - Gene name
      - Source of gene name
  - ORTHOLOGUES (K-O) : 
    - Mouse orthologues:
      - ***ALL*** properties 

- File location : `/gstock/GeneIso/data/External/ENSEMBL/mmusculus_orthologs.txt.gz`

#### Mouse level

To retrieve genes/mRNA/exons from *Mus Musculus*, same procedure was done as the one for [Human level](####Human-level) but with the dataset : **Mouse genes (GRCm39)**.

### ENCODE mouse expression data

Mouse expression data were retrieved through ENCODE project.

Following ENCODE parameters were used to retrieve files :

- *searchTerm=transcript+quantification*
- *type=File*
- *file_format=tsv*
- *assembly=mm10*
- *status=released*
- *audit.NOT_COMPLIANT.category!=insufficient+read+depth*
- *audit.WARNING.category!=low+replicate+concordance*
- *audit.WARNING.category!=low+read+depth*
- *audit.WARNING.category!=missing+analysis_step_run*

Manual observation of files structure according to nb of cols gave the following analysis:

- 123 - 12 cols ??? LRS
- 291 - 17 cols ??? Gene TPM level ??? not usable
- 730 - 19 cols ??? Gene TPM level ??? not usable
- 535 - 5 cols ??? SRS with TPM at transcript level
    - Prot coding
    - Same cutoff than GTEx

Only SRS with TPM level were finally used, with same filtering criteria than GTEx (see GTEx part).



Download list : `/gstock/biolo_datasets/ENCODE/file_report_2022_3_2_10h_59m.tsv`

Associate URL query in the file header : 
```
https://www.encodeproject.org/report/?type=File&searchTerm=transcript+quantification&output_type=transcript+quantifications&file_format=tsv&assembly=mm10&status=released&award.project=ENCODE&lab.title=ENCODE+Processing+Pipeline&audit.NOT_COMPLIANT.category%21=insufficient+read+depth&audit.WARNING.category%21=low+read+depth&audit.WARNING.category%21=low+replicate+concordance&field=%40id&field=dataset&field=href&limit=all
```

Expression files downloaded are located here : `/gstock/biolo_datasets/ENCODE/DL`



## Steps 

Human data processing
1. src/GTEx_correction_by_expression.py

Take as input the config file  (`src/config/config_files.yaml`) as the following : `python GTEx_correction_by_expression.py config/config_files.yaml`
Need to specify transcript read counts & TPM file as well as output file in config file  like the following : 
```
EXPRESSION:
  base_directory: /gstock/biolo_datasets/variation/benchmark/Databases/GTEX/V8/RNA_Seq_data/
  External:
    transcript_reads: GTEx_Analysis_2017-06-05_v8_RSEMv1.3.0_transcript_expected_count.gct.gz
    transcript_tpm: GTEx_Analysis_2017-06-05_v8_RSEMv1.3.0_transcript_tpm.gct.gz
  Final:
    transcript_check: /gstock/EXOTIC/data/EXPRESSION/GTEx_V8_transcript_checking.gct.gz

```

2. src/misog_sisog_analysis_ensembl_gencode.ipynb

Jupyter NB use config file with gene/mRNA & exon level files to compute intermediate files located in `/gstock/GeneIso/V2` :
- genes = "/gstock/GeneIso/V2/Genes.parquet"
- mrna = "/gstock/GeneIso/V2/mRNA.parquet"
- exons = "/gstock/GeneIso/V2/Exons.parquet"
- cds = "/gstock/GeneIso/V2/CDS.parquet"
- five_UTR = "/gstock/GeneIso/V2/5_UTR.parquet"
- three_UTR = "/gstock/GeneIso/V2/3_UTR.parquet"
- introns = "/gstock/GeneIso/V2/Introns.parquet"


3. src/Stats_and_plots.ipynb

Intermediate files reused for make plots & statistics


Mouse data processing

4. src/download_encode_files.py

Used file as argument : `/gstock/biolo_datasets/ENCODE/file_report_2022_3_2_10h_59m.tsv`

5. src/concat_files_ENCODE.py

Concat SRS file into a single dataframe.

6. src/process_encode_srs.py

Identification of transcripts that respect expression filters.

Orthology analysis

7. src/orthologs.ipynb

Identification of closest orthologs in MM

8. src/Stats_ortho.ipynb
9. src/Visualization_ortho.ipynb

Enrichment analysis

10. src/clustering_enrichment_analysis.ipynb

