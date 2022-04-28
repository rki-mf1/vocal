# VOCAL: Variant Of Concern ALert and prioritization 

The goal of VOCAL is to detect sc-2 emerging variants from collected bases of genomes, before their annotation by phylogenetic analysis.
It does so by parsing sc2 genomes and detecting amino acids mutations in the spike proteins that can be associated with a phenotypic change. The phenotypic changes are annotated according to the knowledge accumulated on previous variants. Owing to the limited size of the genome, convergent evolution is expected to take place. 

# Quick start VOCAL

Note: * warning !! VOCAL supports mac OS and Linux system only

## Installation

clone this repository:
```
git clone https://github.com/rki-mf1/vocal.git
```

You can easly install all dependencies with conda:
```
cd vocal
conda create -n vocal -f env/environment.yml
conda activate vocal
```

## Running Vocal
... in two steps.

### Step1: Annotate mutations in the Spike protein

```bash
python vocal/vocal.py -i test-data/sample-test.fasta -o results/variant_table.tsv
```
This creates by default a `variant_table.tsv` file with all mutations. 

⚠️**Note**: when `VOCAL` is run without option, it realigns each query sequence to the reference Wuhan sequence NC_045512 using the pairwise alignment function in the biopython library.
 
🐌 SLOW ??:  The alignment option in vocal uses a biopython pairwise aligner and can be relatively slow. It is thus recommended to first generate an alignment file of all the sequences before running vocal annotation of the mutations.
The alignment file (in PSL format) can be created using the tool `pblat` that can be downloaded [here](https://icebert.github.io/pblat/) or simply installed through our provided conda environment.

👀 Thus, we can use precomputed whole-genome alignments of the fasta file as a PSL file ( `--PSL` option) to improve alignment speed.

**To generate a PSL file with alignments**

Example command to run VOCAL with a PSL file;
```bash
pblat ref.fasta output.psl input.fasta -threads=32
```
```bash
python vocal/vocal.py -i sequences.fasta --psl output.psl -o variant_table.tsv
```

Then (continue from # 1): Annotate mutation phenotypes
```bash
python vocal/Mutations2Function.py -i variant_table.tsv -a data/table_cov2_mutations_annotation.tsv -o variants_with_phenotypes.tsv 
```
By default, this will create the consolidated table ("`variants_with_phenotypes.tsv`") of mutations with phenotype annotation.

### Step2:  Detect/Alert emerging variants

```bash
Rscript --vanilla "vocal/Script_VOCAL_unified.R" \
-v data/ \
-f variants_with_phenotypes.tsv \
-o results/ 
```

in case we want to include meta data
```bash
Rscript --vanilla "vocal/Script_VOCAL_unified.R" \
-v data/ \
-f variants_with_phenotypes.tsv \
-f meta.tsv \
-o results/ 
```

meta data must have



| Syntax      | Description |
| ----------- | ----------- |
| -f        | File corresponding to the variants (output from Mutations2Function.py *variants_with_phenotypes.tsv)     |
| -a        | File containing metadata on the samples        |
| -v        | Directory path where VOCAL database is stored (files concerned: ECDC_assigned_variants.csv and escape_data_bloom_lab.csv and filiation_information)        |
| -o        | Output directory        |
| -file_annot_latest        | File with latest lineage annotation from Desh system [optional argument]      |
| --lineage_column        | Column name of reporting LINEAGE information in metadata file  |
| --date_column        | Column name of reporting sampling date information in metadata file   |
| --geoloc_column        | Column name of geolocalisation information in metadata file      |
| --id_column  | Column name of the sample ID in metadata file   |

--------

# How to interprete result.

Vocal output an alert level in 5 different colours which can be classified into 3 ratings.

| Alert color      | Description | Impact | 
| ----------- | ----------- | ----------- |
| Pink | Variant is known as VOC/VOI and containing MOC or new mutations.   | HIGH |
| Red | Not VOC/VOI but contain high MOC or ROI, and a new matuation (likely to cause a problem/ new dangerous).  | HIGH |
| Orange | Variant contains moderately muations, or also possibly consider them either VUM or De-escalated variant.   | MODERATE |
| Lila | Mostly harmless variant (near-zero mutation size for MOC or ROI). | LOW |
| Grey | No evidence of impact (either no MOC or no ROI).     | LOW |

# VOCAL Ecosystem

## Report

```bash
python  vocal/Reporter.py  \
        -s vocal-alerts-samples-all.csv \
        -c vocal-alerts-clusters-summaries-all.csv \
        -o vocal-report.html \
        --from_date $VOCAL_LOOKBACK \
        --to_date $VOCAL_SELECTED_DATE
```

## Selector
Selector tool offer two functionalites to aid processing in VOCAL.
* A. 'select-desh' command

    ⚠️**Note**: This tool works only for the autopilot system.!!⚠️

    This tool is used to select a fasta file with a given date range, and it will output a result in a merged fasta file.  
    Two important arguments are required. 

    * '`-s`', A start date.
    * '`-e`', Number of day we want to look back.

    For example;  we will select 2021-11-08 and use a 5-day lookback (`-s 2021-11-08 -e 5`). This will give results; 2021-11-08 2021-11-07 2021-11-06 2021-11-05 2021-11-04 2021-11-03 2021-11-02

    Example command;
    ```bash
    python Selector.py select-desh -i ~/wissdaten/MF1_SC2_Backup/sc2/incoming/desh -o output/ -s 2021-11-08 -e 35
    ```
* B. 'convert-covSonar' command

    This tool convert covSonar output (from match command) into Vocal format which can be used in detection.
    
    Requreid argument;
    * `-i`, Input file from covSonar output (*match command)
    * `-o`, Output file (e.g., variants_with_phenotype_sc2-global.tsv)')
    * `-a`, Vocal DB, please use "sc2-vocal/data/table_cov2_mutations_annotation.tsv",

    Optional argument;
    *  `--cpus`, Number of cpus to use (default: 1)
  
    Example command;
    ```bash
    python vocal/Selector.py  convert-covSonar -i covSonar.match.tsv -a data/table_cov2_mutations_annotation.tsv -o table_cov2_mutations_annotation.tsv --cpus 80
    ```

## Update the VOCAL's Knowledge base 

This tool is used to update VOCAL annotation database. Currently, we only support 
* [NetaZuckerman](https://github.com/NetaZuckerman/covid19/blob/master/mutationsTable.xlsx) source
* []

We can provide `--online` tag, the program  will try to download the latest files from sources.

Example command;

To create a new annotation file.
```bash
python  vocal/update.vocalDB.py --online  -o table_cov2_mutations_annotation.tsv
```

If we want to update our table_cov2_mutations_annotation.tsv, just provide the old table_cov2_mutations_annotation.tsv in -i tag.
```bash
python vocal/update.vocalDB.py --online -i data/table_cov2_mutations_annotation.tsv  \
--out_file data/table_cov2_mutations_annotation.tsv
```

If we want to update our table_cov2_mutations_annotation.tsv without download from source (manually download)
```bash
python vocal/update.vocalDB.py -i data/table_cov2_mutations_annotation.tsv  \
-n ../vocal-test-annotation_DB/NetaZuckerman_mutationsTable.xlsx  \
-t ../vocal-test-annotation_DB/Tabelle_VOC-PCR-Finder.xlsx  \
--out_file data/table_cov2_mutations_annotation.tsv
```

## Update the Parent-Child Lineage.

We use parent-child information from [Pango-designation](https://github.com/cov-lineages/pango-designation/) Github. The script requires `lineags.csv` and `alias_key.json` files from Pango-designation Github to reconstruct the information.

Example command;
```bash
python update.lineages.py -l lineags.csv -a alias_key.json -o lineages.all.tsv
```
or automatically download  
```bash
python update.lineages.py --online -o lineages.all.tsv
```
