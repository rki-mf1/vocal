import argparse as ap
import datetime
from multiprocessing import Pool
import os
import re
import sys
import timeit as ti

from Mutations2Function import aggregateMutationTable
from Mutations2Function import merge_variants_annotation
import numpy as np
import pandas as pd
from utils.utility import ROOT_DIR
from utils.utility import S_proteinseq

print(ROOT_DIR)
deletion_df = pd.read_csv(
    os.path.join(ROOT_DIR, "data/fixformat/new_fix_deletion.csv"), header=0
)
deletion_df = deletion_df.drop(deletion_df[deletion_df.vocal_format == "?"].index)
deletion_df["consonar_aa_patern"] = deletion_df["consonar_aa_patern"].str.strip()


def is_file_empty(file_path):
    """Check if file is empty by confirming if its size is 0 bytes"""
    # Check if file exist and it is empty
    if not os.path.exists(file_path):
        print("file doesn't exist")
        return True
    if os.stat(file_path).st_size == 0:
        print("file is empty")
        return True
    return False


#
def fasta_writer(df, out_file):
    out_dir = os.path.dirname(out_file)
    if not os.path.exists(out_dir):
        # Create a new directory because it does not exist
        os.makedirs(out_dir)

    # Exporting to fa
    # adding '>' for accessions
    df["Accession"] = ">" + df["Accession"]
    df.to_csv(out_file, sep="\n", index=None, header=None)
    print("Write fasta file:", out_file)


#
def fasta_reader(file):
    print("Read:", file)
    if is_file_empty(file):
        return pd.DataFrame()
    fasta_df = pd.read_csv(file, sep=">", lineterminator=">", header=None)
    fasta_df[["Accession", "Sequence"]] = fasta_df[0].str.split("\n", 1, expand=True)
    fasta_df.drop(0, axis=1, inplace=True)
    fasta_df["Sequence"] = fasta_df["Sequence"].replace("\n", "", regex=True)
    return fasta_df


# use pandas instead BioPython
def scan_data(in_dir, start_date, end_date):
    date_list = pd.date_range(end_date, start_date).format(
        formatter=lambda x: x.strftime("%Y-%m-%d")
    )
    # print(date_list)
    df = pd.concat(
        (fasta_reader(os.path.join(in_dir, i, "day.fasta")) for i in date_list),
        ignore_index=True,
    )
    return df


def main(args):
    start_date = args.start_date.strftime("%Y-%m-%d")
    end_date = (args.start_date - datetime.timedelta(days=args.look_back)).strftime(
        "%Y-%m-%d"
    )
    input_dir = args.in_dir
    out_file = args.out_file
    print("Start:", start_date, "End:", end_date)

    df = scan_data(input_dir, start_date, end_date)
    if len(df) > 0:
        fasta_writer(df, out_file)
    else:
        # create empty fasta file
        with open(out_file, "w") as document:
            pass


# Convert deletion to Vocal deletion format
# Find variant_size, variant_type
def fix_format(_df):
    temp = re.compile("([*a-zA-Z]+)([0-9]+)([*a-zA-Z]+)")  # "R214REPE" "R256E"
    temp_del = re.compile("([*a-zA-Z]+)([0-9]+)-([0-9]+)([*a-zA-Z]+)")  # HV69-70del
    # for  row in tqdm(_df.itertuples(index=True), total=_df.shape[0],disable=None):
    for row in _df.itertuples(index=True):
        idx = row.Index
        try:
            if "del:" in row.aa_profile:  # Deletion
                _split_list = row.aa_profile.split(":")
                delete_pos = int(_split_list[1])
                delete_size = int(_split_list[2])
                new_format = deletion_fixAAformat(row.aa_profile, row.lineage)
                if not new_format:
                    new_format = (
                        S_proteinseq[delete_pos - 1 : (delete_pos + delete_size - 1)]
                        + str(delete_pos)
                        + "-"
                        + str(delete_pos + delete_size - 1)
                        + "del"
                    )

                    _df.loc[idx, "new_aa_profile"] = new_format
                    _df.loc[idx, "variant_type"] = "D"
                    _df.loc[idx, "variant_size"] = delete_size
                    _df.loc[idx, "aa_pos_ref_start"] = delete_pos
                else:
                    _pattern = temp_del.match(new_format).groups()
                    if len(_pattern) == 4:
                        _df.loc[idx, "new_aa_profile"] = new_format
                        _df.loc[idx, "variant_type"] = "D"
                        _df.loc[idx, "variant_size"] = len(_pattern[0])
                        _df.loc[idx, "aa_pos_ref_start"] = _pattern[1]

                    else:
                        print("Found unexpected Deletion format, stop!!!")
                        raise
            else:
                _pattern = temp.match(row.aa_profile).groups()
                if len(_pattern) == 3:
                    _ref = _pattern[0]
                    _alt = _pattern[2]
                    _pos = _pattern[1]
                    if len(_ref) < len(_alt):  # insertion
                        new_format = _alt + str(_pos) + "-" + str(_pos) + "ins"
                        _df.loc[idx, "new_aa_profile"] = new_format
                        _df.loc[idx, "variant_type"] = "I"
                        _df.loc[idx, "variant_size"] = len(set(_alt) - set(_ref))
                        _df.loc[idx, "aa_pos_ref_start"] = int(_pos)
                    else:  # SNP
                        _df.loc[idx, "new_aa_profile"] = row.aa_profile
                        _df.loc[idx, "variant_type"] = "M"
                        _df.loc[
                            idx, "variant_size"
                        ] = 1  #  works for now, in future,we might need to chage how we calculate it.
                        _df.loc[idx, "aa_pos_ref_start"] = int(_pos)
                else:  # Unsupported case
                    _df.loc[idx, "new_aa_profile"] = row.aa_profile
                    _df.loc[idx, "variant_type"] = "U"
                    _df.loc[idx, "variant_size"] = 0
                    _df.loc[idx, "aa_pos_ref_start"] = 0
                    print("Warning: Found unsupported format")
                    print(row)
        except Exception:
            print(row)
            raise
    _df["variant_size"] = _df["variant_size"].astype(np.int8)
    _df["aa_pos_ref_start"] = _df["aa_pos_ref_start"].astype(np.int32)
    return _df


def deletion_fixAAformat(aa_pattern, lineage):
    # pattern = nt_pattern.split(' ')
    # skip ? #
    result = deletion_df[
        (deletion_df["consonar_aa_patern"].str.contains(aa_pattern))
        & (deletion_df["ID"] == lineage.strip())
    ]

    if (len(result)) == 0:
        # print('False')
        return False
    else:
        # print('Detect lineage to be fixed:',lineage)
        # print('New format:',result['vocal_format'].values[0])
        return result["vocal_format"].values[0]


def parallelize_dataframe(df, func, num_cores):
    _tmp_lis = np.array_split(df, num_cores)

    # with Pool(processes=num_cores) as pool:
    #    res = pool.starmap(func, zip_items)
    pool = Pool(num_cores)
    full_df_list = pool.map(func, _tmp_lis)

    # finish all tasks
    pool.close()
    pool.join()
    # print("Tmp result: ", full_paht_list)
    _df = pd.concat(full_df_list)
    return _df


def main_covSonar(args):
    annot_file = args.annotation
    covSonar_output = args.input
    output_file = args.output
    num_cores = args.cpus
    df_annot = pd.read_csv(annot_file, sep="\t")

    if annot_file.endswith(".csv"):
        df_annot = pd.read_csv(annot_file)
    elif annot_file.endswith(".tsv"):
        df_annot = pd.read_csv(annot_file, sep="\t")
    else:
        raise TypeError("Not recognized file type")

    df_agg_annot = aggregateMutationTable(df_annot)

    # check input
    df_variants = pd.read_csv(covSonar_output, sep="\t")
    want_to_join = df_variants[["accession", "lineage"]]
    # Start to clear format

    print(
        "These samples will not be processed due to no AA profile present:",
        df_variants[df_variants.aa_profile.isnull()]["accession"].to_list(),
    )
    df_variants = df_variants[
        df_variants.aa_profile.notnull()
    ]  # remove null AA profile
    b = pd.DataFrame(
        df_variants.aa_profile.str.split(" ").tolist(), index=df_variants.accession
    ).stack()
    b = b.reset_index()[[0, "accession"]]  # aa_profile variable is currently labeled 0
    b.columns = ["aa_profile", "accession"]  # renaming aa_profile
    b["aa_profile"] = b["aa_profile"].str.strip()  # clear all space
    b["target_gene"] = b["aa_profile"].str.split(":").str[0]
    b = b[b["target_gene"] == "S"]  # only S gene
    b.reset_index(drop=True, inplace=True)
    b["aa_profile"] = b["aa_profile"].str.split(":", 1).str[1]
    b = b.merge(want_to_join, on="accession")

    _new_df = parallelize_dataframe(b, fix_format, num_cores)
    # _new_df=fix_format(b)
    _new_df.rename(
        columns={
            "aa_profile": "covSonar_aa_profile",
            "new_aa_profile": "aa_pattern",
            "accession": "ID",
        },
        inplace=True,
    )
    _new_df["nt_pattern"] = ""
    # print(len(_new_df))
    # Final
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
    # s["ID", "target_gene", "aa_pattern", "aa_pos_ref_start","variant_type", "variant_size", "type", "infos",'covSonar_aa_profile']
    df_variant_with_annot = merge_variants_annotation(
        _new_df, df_agg_annot
    ).sort_values(["ID", "aa_pos_ref_start", "variant_type", "type"])
    df_variant_with_annot = df_variant_with_annot[columns_save]
    # df_variant_with_annot.drop_duplicates( inplace=True)
    df_variant_with_annot[columns_save].to_csv(output_file, sep="\t", index=False)
    # print(len(df_variant_with_annot))


if __name__ == "__main__":
    yesterday = (datetime.date.today() - datetime.timedelta(days=1)).strftime(
        "%Y-%m-%d"
    )

    parser = ap.ArgumentParser(prog="Selector.py", description="")
    subparsers = parser.add_subparsers(
        help="Select mutation profile from different sources and transform it to Vocal input"
    )
    subparsers.dest = "tool"
    subparsers.required = True

    # print("This tool works for autopilot system only. It uses for selecting a seqeunce with a given timeframe \n")

    parser_select_covsonar = subparsers.add_parser(
        "convert-covSonar",
        help="convert covSonar output (from match command) into VOCAL format which can be used in detection.",
    )
    parser_select_covsonar.add_argument(
        "-i",
        "--input",
        required=True,
        help="Input file from covSonar output (*match command)",
    )
    parser_select_covsonar.add_argument(
        "-o",
        "--output",
        required=True,
        default="variants_with_phenotype_sc2-global.tsv",
        help="Output file  (e.g., variants_with_phenotype_sc2-global.tsv)",
    )
    parser_select_covsonar.add_argument(
        "--cpus",
        metavar="int",
        help="number of cpus to use (default: 1)",
        type=int,
        default=1,
    )
    parser_select_covsonar.add_argument(
        "-a",
        "--annotation",
        required=True,
        default="data/table_cov2_mutations_annotation.tsv",
        help="Table with information about lineage defining mutation and Variants Of Concern (data/table_cov2_mutations_annotation.tsv)",
    )

    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit()
    args = parser.parse_args()
    t1 = ti.default_timer()
    if args.tool == "convert-covSonar":
        main_covSonar(args)
    else:
        sys.exit("please input a correct command")
    t2 = ti.default_timer()
    print("Processing time: {} seconds".format(round(t2 - t1, 2)))
