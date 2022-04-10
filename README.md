# VOCAL: Variant Of Concern ALert and prioritization 

The goal of Vocal is to detect sc-2 emerging variants from collected bases of genomes, before their annotation by phylogenetic analysis.
It does so by parsing sc2 genomes and detecting amino acids mutations in the spike proteins that can be associated with a phenotypic change. The phenotypic changes are annotated according to the knowledge accumulated on previous variants. Owing to the limited size of the genome, convergent evolution is expected to take place. 

# Quick start VOCAL

Note*: The executable package (`conda install VOCAL` or `pip install VOCAL`) will be avaialable soon. 

## Installation

clone this repository to your computer:
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
# 1)
python vocal.py -i sequences.fasta 
```
This creates by default a `variant_table.tsv` file with all mutations. 

⚠️**Note**: when `vocal` is run without option, it realigns each query sequence to the reference Wuhan sequence NC_045512 using the pairwise alignment function in the biopython library.

Alternatively you can use precomputed whole genome alignments of the fasta file as a PSL file with the `--PSL` option. (Improve alignment speed)

**To generate a PSL file with alignments**

The alignment option in vocal uses biopython pairwise aligner and can be relatively slow. It is thus recommended to first generate an alignment file of all the sequences before running vocal annotation of the mutations.
The alignment file (in PSL format) can be created using the tool `pblat` that can be downloaded [here](https://icebert.github.io/pblat/) or simply installed through conda:

Example command to run vocal with a PSL file;
```bash
pblat ref.fasta output.psl input.fasta -threads=32

python vocal.py -i sequences.fasta --psl output.psl -o vocal.tsv

```

```bash
# 2) Then (continue from # 1): Annotate mutation phenotypes
python Mutations2Function.py -i variant_table.tsv -a data/table_cov2_mutations_annotation.csv -o variants_with_phenotypes.tsv 
```
By default, this will create the consolidated table `variants_with_phenotypes.tsv` of mutations with phenotype annotation.

### Step2:  Detect/Alert emerging variants

```bash
# 3)
Rscript --vanilla "vocal/Script_VOCAL_unified.R" \
-v data/ \
-f variants_with_phenotypes.tsv \
-o /results/ \
-a meta-information.tsv \
--lineage_column lineage \
--date_column date \
--id_column accession

```

| Syntax      | Description |
| ----------- | ----------- |
| -f        | File corresponding to the variants (output from Mutations2Function.py *variants_with_phenotypes.tsv)     |
| -a        | File containing metadata on the samples        |
| -v        | Directory path where Vocal database is stored (files concerned: ECDC_assigned_variants.csv and escape_data_bloom_lab.csv and filiation_information)        |
| -o        | Output directory        |
| -file_annot_latest        | File with latest lineage annotation from Desh system [optional argument]      |
| --lineage_column        | Column name of reporting LINEAGE information in metadata file  |
| --date_column        | Column name of reporting sampling date information in metadata file   |
| --geoloc_column        | Column name of geolocalisation information in metadata file      |
| --id_column  | Column name of the sample ID in metadata file   |

--------
## How to interprete result.

Vocal output an alert level in 5 different colours which can be classified into 3 ratings.

| Alert color      | Description | Impact | 
| ----------- | ----------- | ----------- |
| Pink | Variant is known as VOC/VOI and containing MOC or new mutations.   | HIGH |
| Red | Not VOC/VOI but contain high MOC or ROI, and a new matuation (likely to cause a problem/ new dangerous).  | HIGH |
| Orange | Variant contains moderately muations, or also possibly consider them either VUM or De-escalated variant.   | MODERATE |
| Lila | Mostly harmless variant (near-zero mutation size for MOC or ROI). | LOW |
| Grey | No evidence of impact (either no MOC or no ROI).     | LOW |
