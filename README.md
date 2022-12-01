<div id="top"></div>

<div align="center">
<h1 align="center"> VOCAL: Variant Of Concern ALert and prioritization </h1>
</div>
The goal of VOCAL is to detect sc-2 emerging variants from collected bases of genomes, before their annotation by phylogenetic analysis.
It does so by parsing sc2 genomes and detecting amino acids mutations in the spike proteins that can be associated with a phenotypic change. The phenotypic changes are annotated according to the knowledge accumulated on previous variants. Owing to the limited size of the genome, convergent evolution is expected to take place. 

# Documentation.

<a href="https://rki-mf1.github.io/vocal-doc/"><strong>Explore the docs »</strong></a>

# Getting Started

⚠️**Note**: 🔌 Right now, VOCAL tested on Linux and Mac system only 💻 

## Installation

clone this repository:
```
git clone https://github.com/rki-mf1/vocal.git
```

You can easly install all dependencies with conda:
```
cd vocal
conda env create -n vocal -f environment.yml
conda activate vocal
```

## Running Vocal
... in three steps.

### Step1: Annotate mutations in the Spike protein

```bash
python vocal/vocal.py -i test-data/sample-test.fasta -o results/variant_table.tsv
```
This creates by default a `variant_table.tsv` file with all mutations. 

⚠️**Note**: when `VOCAL` is run without option, it realigns each query sequence to the reference Wuhan sequence NC_045512 using the pairwise alignment function in the biopython library.
 
🐌 SLOW ??:  The alignment option in vocal uses a biopython pairwise aligner and can be relatively slow. It is thus recommended to first generate an alignment file of all the sequences before running vocal annotation of the mutations.
The alignment file (in PSL format) can be created using the tool `pblat` that can be downloaded [here](https://icebert.github.io/pblat/) or simply installed through our provided conda environment.

👀 Thus, if we want to use precomputed whole-genome alignments of the fasta file as a PSL file ( `--PSL` option) to improve alignment speed please see the below section, otherwise please continue to **step2**.

**To generate a PSL file with alignments**

Example command to generate PSL format.
```bash
pblat test-data/ref.fna test-data/sample-test.fasta -threads=4 results/output.psl
```

To run VOCAL with a PSL file;
```bash
python vocal/vocal.py -i test-data/sample-test.fasta --PSL results/output.psl -o results/variant_table.tsv
```

### Step2: Annotate mutation phenotypes

```bash
python vocal/Mutations2Function.py -i results/variant_table.tsv -a data/table_cov2_mutations_annotation.tsv -o results/variants_with_phenotypes.tsv 
```
By default, this step will create the consolidated table ("`variants_with_phenotypes.tsv`") of mutations with phenotype annotation. 

### Step3: Detect/Alert emerging variants

```bash
Rscript --vanilla "vocal/Script_VOCAL_unified.R" \
-f results/variants_with_phenotypes.tsv \
-o results/ 
```

in case we want to include metadata file, use (-a)
```bash
Rscript --vanilla "vocal/Script_VOCAL_unified.R" \
-f results/variants_with_phenotypes.tsv \
-a test-data/meta.tsv \
-o results/ 
```
⚠️**Note**: meta data must have these information
* ID column (match with sample ID in FASTA file)
* LINEAGE column (e.g., B.1.1.7, BA.1)
* SAMPLING DATE column (the date that a sample was collected) (format YYYY-mm-dd)

Finally, we can easily generate report into HTML format at the end of the analysis.

```bash
python  vocal/Reporter.py  \
        -s results/vocal-alerts-samples-all.csv \
        -c results/vocal-alerts-clusters-summaries-all.csv \
        -o results/vocal-report.html 
```

Please visit <a href="https://rki-mf1.github.io/vocal-doc/"><strong>Explore the docs »</strong></a>

# How to interprete result.

Vocal output an alert level in 5 different colours which can be classified into 3 ratings.

| Alert color      | Description | Impact | 
| ----------- | ----------- | ----------- |
| Pink | Variant is known as VOC/VOI and containing MOC or new mutations.   | HIGH |
| Red | Not VOC/VOI but contain high MOC or ROI, and a new matuation (likely to cause a problem/ new dangerous).  | HIGH |
| Orange | Variant contains moderately muations, or also possibly consider them either VUM or De-escalated variant.   | MODERATE |
| Lila | Mostly harmless variant (near-zero mutation size for MOC or ROI). | LOW |
| Grey | No evidence of impact (either no MOC or no ROI).     | LOW |

# Contact

Did you find a bug?🐛 Suggestion/Feedback/Feature request?👨‍💻 please visit [GitHub Issues](https://github.com/rki-mf1/vocal/issues)

For business inquiries or professional support requests 🍺 please contact 
Dr. Hölzer, Martin(<HoelzerM@rki.de>) or Dr. Richard, Hugues (<RichardH@rki.de>)

# Acknowledgments

* Original Idea: SC2 Evolution Working group at the Robert Koch Institute in Berlin

* Funding: Supported by the European Centers for Disease Control [grant number ECDC GRANT/2021/008 ECD.12222].



