# Information about external datasets used for phenotype annotation

## Hackmd dataset
#### *`table_cov2_mutations_annotation.tsv`*

The [2021-09 version](data/2021-09/table_cov2_mutations_annotation.tsv) is an in house ressource developed by people from MF1, P3 and P5 group at the Robert Koch Institute (Berlin, Germany).

The [2022-11 version](data/2022-11/table_cov2_mutations_annotation.tsv) was built with an in house Mutation Of Concern table and the results of a the [sc2-mutation-frequency-calculator](https://github.com/rki-mf1/sc2-mutation-frequency-calculator) for the Lineage Defining Mutations on a dataset with sequences from October 2021 to November 2022.


## Antibody escape mutations by deep mutational screening
#### *`escape_data_bloom_lab.csv`*

This list was curated by the Bloom lab and accessible at [this adress](https://jbloomlab.github.io/SARS2_RBD_Ab_escape_maps/) and in CSV format at [this adress](https://raw.githubusercontent.com/jbloomlab/SARS2_RBD_Ab_escape_maps/main/processed_data/escape_data.csv).


## List of Variants of Concerns and Variants of Interest 
#### *`ECDC_assigned_variants.csv`*

We use the list from the ECDC and accessible [here](https://www.ecdc.europa.eu/en/covid-19/variants-concern) or 
[here in table format](https://github.com/erikalmecdc/ecdc_virology/blob/main/assigned_variants.csv).

The last commit to the table from ECDC was on September 8, 2021. Afterwards, we updated the table manually with the changelog that can be found [here](https://www.ecdc.europa.eu/sites/default/files/documents/Variants%20page%20changelog_8.pdf). 
We added the children lineages for all of the VOC/VOI/VUM to the table if the parent lineage was listed as VOC/VOI/VUM to be in line with the definitions at that time. This will be no longer necessary for datasets covering timepoints after [March 15th, 2023 as the definitions were changed by WHO](https://www.who.int/news/item/16-03-2023-statement-on-the-update-of-who-s-working-definitions-and-tracking-system-for-sars-cov-2-variants-of-concern-and-variants-of-interest).


## List of Lineages
#### *`lineage.all.tsv`*

To get a updated list of lineages, run [`update.lineages.py`](bin/update.lineages.py) on the needed version of the [alias_key.json from the pango-designation repository](https://github.com/cov-lineages/pango-designation/blob/master/pango_designation/alias_key.json).
