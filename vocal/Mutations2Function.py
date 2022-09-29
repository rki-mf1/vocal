"""

Annotate the list of mutations found in a dataframe

@hrichard

"""

import argparse

from data_loader_tmp import DF_SPIKE_ROI
import pandas as pd


def aggregateMutationTable(
    df_annot,
    group_cols=["gene", "amino acid", "type"],
    ID_col="ID",
    comment_col="comment",
    sep=";",
):
    """
    Aggregates the information from a table of mutations according to group_cols and merges the
    informations from infos_cols
    """
    df_aggregated = df_annot.groupby(group_cols).agg(
        ID_list=pd.NamedAgg(
            column=ID_col, aggfunc=lambda x: sep.join([str(a) for a in x.tolist()])
        ),
        comment_list=pd.NamedAgg(
            column=comment_col, aggfunc=lambda x: sep.join([str(a) for a in x.tolist()])
        ),
    )
    return df_aggregated.reset_index()


def merge_variants_annotation(
    df_mutations,
    df_annot,
    mut_merge=["target_gene", "aa_pattern"],
    annot_merge=["gene", "amino acid"],
):
    """
    Merge a df of mutations with a df of lineages and VOC.
    Reports as well if the mutation is in a region of concern
    by concatenating the table DF_SPIKE_ROI in data_loader_tmp
    !!! There has to be a type column in each of the DataFrame
    """
    df_merge_position = pd.merge(
        df_mutations,
        df_annot,
        left_on=mut_merge,
        right_on=annot_merge,
        how="left",
        validate="many_to_many",
    )
    df_merge_position = df_merge_position.fillna(value={"type": "NotAnnotated"}).rename(
        columns={"ID_list": "infos"}
    )
    ##TODO and check
    df_merge_region = pd.merge(
        df_mutations,
        DF_SPIKE_ROI,
        left_on=["target_gene", "aa_pos_ref_start"],
        right_on=["gene", "aa_position"],
        how="inner",
    ).rename(columns={"functional domain": "infos"})

    return pd.concat([df_merge_position, df_merge_region], ignore_index=True)


def annotateVariantTable(df_variants, annot_file):
    """
    DataFrame * annot_file -> DataFrame
    Annotate the variants in the DataFrame using the table of annotations given in annot_file
    the mutations in annot_file are aggregated
    """
    if annot_file.endswith(".csv"):
        df_annot = pd.read_csv(annot_file)
    elif annot_file.endswith(".tsv"):
        df_annot = pd.read_csv(annot_file, sep="\t")
    else:
        raise TypeError("Not recognized file type")
    df_agg_annot = aggregateMutationTable(df_annot)
    df_variant_with_annot = merge_variants_annotation(df_variants, df_agg_annot)
    return df_variant_with_annot


def main():
    parser = argparse.ArgumentParser(
        description="Mutations2Functions: Intersect variant table with a set of sequences",
        epilog="Usage: python Mutations2Functions.py -i mutations_table.tsv -p mutation_phenotypes.tsv -o [variants_with_phenotypes.tsv]",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument(
        "-i", "--input", required=True, help="table of mutations (produced by Vocal)"
    )
    parser.add_argument(
        "-a",
        "--annotation",
        default="table_cov2_mutations_annotation.csv",
        help="Table with information about lineage defining mutation and Variants Of Concern",
    )
    parser.add_argument(
        "-o",
        "--output",
        default="variants_with_phenotypes.tsv",
        help="output table with the set of variants detected over all sequences, with annotation",
    )
    parser.add_argument(
        "-L",
        "--largetable",
        action="store_true",
        help="write all columns for the result table (there may be redundancies)",
    )

    args = parser.parse_args()

    tab_file = args.input
    annot_file = args.annotation

    out_table = args.output
    # Reading in the files

    df_variants = pd.read_csv(tab_file, sep="\t")
    df_variant_with_annot = annotateVariantTable(df_variants, annot_file).sort_values(
        ["ID", "aa_pos_ref_start", "variant_type", "type"]
    )

    columns_save = [
        "ID",
        "target_gene",
        "aa_pattern",
        "nt_pattern",
        "aa_pos_ref_start",
        "variant_type",
        "variant_size",
        "type",
        "infos",
    ]

    if args.largetable:
        df_variant_with_annot.to_csv(out_table, sep="\t", index=False)
    else:
        df_variant_with_annot[columns_save].to_csv(out_table, sep="\t", index=False)
    ##End


if __name__ == "__main__":
    main()
