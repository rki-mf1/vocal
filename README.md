<div id="top"></div>

<div align="center">
<h1 align="center"> VOCAL: Variant Of Concern ALert and prioritization </h1>
</div>
The goal of VOCAL is to detect SARS-CoV-2 emerging variants from collected bases of genomes, before their annotation by phylogenetic analysis.
It does so by parsing SARS-CoV-2 genomes and detecting amino acids mutations in the spike proteins that can be associated with a phenotypic change. The phenotypic changes are annotated according to the knowledge accumulated on previous variants. Owing to the limited size of the genome, convergent evolution is expected to take place. 

# Documentation.

<a href="https://rki-mf1.github.io/vocal-doc/"><strong>Explore the docs »</strong></a>

# Getting Started

⚠️**Note**: 🔌 Right now, VOCAL tested on Linux system only 💻 

## Quick Installation

To run the pipeline, you need to have `conda` and `Nextflow` installed and set up.
All other dependencies will be installed over `conda` in the pipeline.

To install `conda`, use the following bash commands if you are working on **Linux**:
```bash
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
bash Miniconda3-latest-Linux-x86_64.sh
```

To install `conda`, use the following bash commands if you are working on **Mac**:
```bash
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-MacOSX-arm64.sh
bash Miniconda3-latest-MacOSX-arm64.sh
```

Then, `Nextflow` an be installed over `conda`:
```bash
conda create -n nextflow -c bioconda nextflow
conda activate nextflow
```

The VOCAL repository can be cloned from Git:
```bash
git clone https://github.com/rki-mf1/vocal.git
```

### Call help

```bash
nextflow run main.nf --help
```

## Running VOCAL

```bash
nextflow run main.nf  \
     -profile conda,local \
     --fasta 'test/sample-test.fasta'
```

⚠️**Note**: When `VOCAL` is run without option `--psl`, it realigns each query sequence to the reference sequence Wuhan NC_045512 using the pairwise alignment function in the biopython library.

### With pblat
 
🐌 SLOW?: The alignment option in VOCAL uses a biopython pairwise aligner and can be relatively slow. It is thus recommended to first generate an alignment file of all the sequences before running VOCAL annotation of the mutations. The alignment file (in PSL format) can be created using the tool `pblat` by adding the option `--psl`.

```bash
nextflow run main.nf  \
     -profile conda,local \
     --fasta 'test/sample-test.fasta' \
     --psl
```

### With metadata file

```bash
nextflow run main.nf  \
     -profile conda,local \
     --fasta 'test/sample-test.fasta' \
     --metadata 'test/meta.tsv' \
     --psl
```

⚠️**Note**: Metadata must have these information
* ID column (match with sample ID in FASTA file)
* LINEAGE column (e.g., B.1.1.7, BA.1)
* SAMPLING DATE column (the date that a sample was collected) (format YYYY-MM-DD)

## Parameter list

```
fasta                    REQUIRED! The path to the fasta file with the sequences for VOCAL.
                         [ default: '' ]
metadata                 The path to a metadate file for the sequences.
                         [ default: '' ]
psl                      Run process with pblat alignment.
                         [ default: false ]
strict                   Run process with strict alert levels (without orange).
                         [ default: 'n' ]
```

# How to interprete result.

VOCAL output an alert level in four different colours which can be classified into 3 ratings.

| Alert color      | Description | Impact | 
| ----------- | ----------- | ----------- |
| Pink | Variant is known as VOC/VOI and containing MOC or new mutations.   | HIGH |
| Red | Not VOC/VOI but contain high MOC or ROI, and a new matuation (likely to cause a problem/ new dangerous).  | HIGH |
| Orange | Variant contains moderately muations, or also possibly consider them either VUM or De-escalated variant.   | MODERATE |
| Grey | Near-zero mutation size for MOC or ROI or either no MOC or no ROI.     | LOW |

# Contact

Did you find a bug?🐛 Suggestion/Feedback/Feature request?👨‍💻 please visit [GitHub Issues](https://github.com/rki-mf1/vocal/issues)

For business inquiries or professional support requests 🍺 please contact 
Dr. Martin Hölzer(<HoelzerM@rki.de>) or Dr. Hugues Richard(<RichardH@rki.de>)

# Acknowledgments

* Original Idea: SC2 Evolution Working group - M. von Kleist, S. Clavignac-Spencer

* Funding: Supported by the European Centers for Disease Control [grant number ECDC GRANT/2021/008 ECD.12222].



