#!/usr/bin/python
# Maintainer: KongkitimanonK

import argparse as ap
import ast
import os
import re
import shutil
import sys
import timeit as ti

import numpy as np
import pandas as pd
import requests
from utils.utility import get_current_date
from utils.utility import ROOT_DIR
from utils.utility import update_version, S_proteinseq

# URL
dict_annotation_db = {
    "NetaZuckerman": "https://github.com/NetaZuckerman/covid19/blob/master/mutationsTable.xlsx?raw=true",
    "SC2-Variants": "https://raw.githubusercontent.com/3dgiordano/SARS-CoV-2-Variants/main/data/variants.csv",
}

deletion_df = pd.read_csv(
    os.path.join(ROOT_DIR, "data/fixformat/new_fix_deletion.csv"), header=0
)
deletion_df = deletion_df.drop(deletion_df[deletion_df.vocal_format == "?"].index)


def deletion_fixAAformat(nt_pattern, aa_pattern, lineage):
    # skip ? #
    result = deletion_df.loc[
        (deletion_df["nucleotide"] == nt_pattern)
        & (deletion_df["amino acid"] == aa_pattern)
        & (deletion_df["ID"] == lineage)
    ]
    if (len(result)) == 0:
        return False
    else:
        # print('Detect lineage',lineage)
        # print('old format', aa_pattern, nt_pattern)
        # print('new format',result['vocal_format'].values[0])
        return result["vocal_format"].values[0]


def create_vocal_nuc_format(_frame):
    _array_type = _frame["mutation type"].unique()
    pattern = ""
    if len(_array_type) == 1:
        _type = _array_type[0]
        start_pos = _frame["position"].iloc[0]  # first row as start position
        _size = len(_frame)
        if "deletion" in _type:  # deletion
            pattern = "del:" + str(start_pos) + ":" + str(_size)
        elif _type == "insertion":  # insertion
            pattern = "ins:" + str(start_pos) + ":" + str(_size)
        else:
            print(_frame)
            print("something wrong in create_vocal_format")
            raise
    else:
        print(_frame)
        print("something wrong in create_vocal_format")
        raise
    return pattern


def NetaZuckerman(input_NetaZuckerman):
    listOfMutation = [
        "SNP",
        "deletion",
        "SNP_silent",
        "extragenic",
        "deletion_frameshift",
        "SNP_stop",
        "insertion",
    ]
    nut_pattern = re.compile("([*a-zA-Z]+)([0-9]+)")
    # FilterNetaZuckerman_dfs only S Gene
    NetaZuckerman_dfs = pd.read_excel(
        input_NetaZuckerman, sheet_name=0, converters={"Position": int}
    )
    _df = NetaZuckerman_dfs[NetaZuckerman_dfs["protein"] == "S"]
    _df = _df.reset_index(drop=True)
    # Fix variant name
    _df["variantname"] = _df["variantname"].str.strip()
    _df["variantname"] = _df["variantname"].str.replace("-", " ")
    _df["variantname"] = _df["variantname"].str.split(" ").str[0]

    _df["nuc sub"] = ""
    _df["variant_vocal"] = ""
    _df["comment"] = ""
    _df["type"] = "LineageDefiningMutation"
    for i, row in _df.iterrows():  # fix format
        if (
            row["mutation type"] == "deletion"
            or row["mutation type"] == "deletion_frameshift"
        ):
            _pattern = nut_pattern.match(row["variant"]).groups()
            if len(_pattern) == 2:
                new_reps = create_vocal_nuc_format(
                    _df[
                        (_df["variantname"] == row["variantname"])
                        & (_df["variant"] == row["variant"])
                    ]
                )
                vocal_format = deletion_fixAAformat(
                    new_reps, row["variant"], row["variantname"]
                )
                if vocal_format is not False:
                    _df.at[i, "variant_vocal"] = vocal_format
                else:
                    new_pattern = (
                        _pattern[0] + str(_pattern[1]) + "-" + str(_pattern[1]) + "del"
                    )
                    _df.at[i, "variant_vocal"] = new_pattern
                _df.at[i, "nuc sub"] = new_reps
            else:
                print(_pattern)
                print("something wrong with this deletion, please reach us asap")
        elif "SNP" in row["mutation type"]:  # if SNP
            if row["nuc sub"] == "" or pd.isnull(row["nuc sub"]):
                _df.at[i, "variant_vocal"] = row["variant"]
                _df.at[i, "nuc sub"] = row["nucsub"]
        elif "insertion" in row["mutation type"]:  # if insertion
            _pattern = nut_pattern.match(row["variant"]).groups()
            if len(_pattern) == 2:
                new_pattern = (
                    _pattern[0] + str(_pattern[1]) + "-" + str(_pattern[1]) + "ins"
                )
                _df.at[i, "variant_vocal"] = new_pattern
                _df.at[i, "nuc sub"] = create_vocal_nuc_format(
                    _df[
                        (_df["variantname"] == row["variantname"])
                        & (_df["variant"] == row["variant"])
                    ]
                )
            else:
                print(_pattern)
                print("something wrong with this insertion, please reach us asap")
        elif row["mutation type"] not in listOfMutation:
            print(row["mutation type"])
            print("We are not support this mutation, please reach us asap")
    # drop unnecessary column and change column name
    _df = _df.drop(
        [
            "position",
            "reference",
            "mutation",
            "variant",
            "mutation type",
            "annotation",
            "varname",
            "nucsub",
        ],
        axis=1,
    )
    _df.rename(
        columns={
            "variantname": "ID",
            "protein": "gene",
            "nuc sub": "nucleotide",
            "variant_vocal": "amino acid",
        },
        inplace=True,
    )
    # deduplicate we only need one record for both deletion and insertion
    _df = _df.drop_duplicates(subset=["ID", "amino acid"])

    return _df


def match_star_lineage(df, path_lineage):
    lineage_df = pd.read_csv(path_lineage, sep="\t", header=0)
    df.rename(columns={"Pango lineage": "Pango_lineage"}, inplace=True)
    select_list = df[df["Pango_lineage"].str.contains(r"\*")]

    for row in select_list.itertuples():
        _status = row.Status
        to_search = row.Pango_lineage
        lst_dict = []
        _selected_lins = lineage_df[lineage_df["lineage"].str.contains(to_search)][
            "lineage"
        ].to_list()
        for x in _selected_lins:
            lst_dict.append({"Pango_lineage": x, "Status": _status, "Label": x})
        df = df.append(lst_dict)

    df = df[~df["Pango_lineage"].str.contains(r"\*")]
    df.rename(columns={"Pango_lineage": "Pango lineage"}, inplace=True)
    return df


def write_to_file(df, file_path, format="\t"):
    # if not os.path.exists(out_dir):
    # Create a new directory because it does not exist
    #    os.makedirs(out_dir)
    df.to_csv(file_path, sep=format, index=False)


def write_download_to_file(content, out_dir, filename, extenion):
    # extenion = ".xlsx"
    # if filename == "":
    #   extenion

    file_path = os.path.join(out_dir, filename + extenion)

    if not os.path.exists(out_dir):
        # Create a new directory because it does not exist
        os.makedirs(out_dir)

    with open(file_path, "wb") as f:
        f.write(content)


def start_download(tmp_dir, what_to_download):
    print("Start Download", what_to_download)
    if what_to_download == "NetaZuckerman":
        extension = ".xlsx"
    elif what_to_download == "SC2-Variants":
        extension = ".csv"
    else:
        extension = ""

    for key, url in dict_annotation_db.items():
        if what_to_download == key:
            r = requests.get(url)
            if r.status_code == 200:
                filename = key
                write_download_to_file(r.content, tmp_dir, filename, extension)
            else:
                print("Can not download, status code is: " + str(r.status_code))
                raise requests.exceptions.HTTPError
    return os.path.join(tmp_dir, filename + extension)


def main_lineage_status(args, key):
    _online = args.online
    _our_annotation_path = args.input_annotation_file
    tmp_dir = args.tmp_directory
    output = args.out_file

    # if online
    if _online:
        _SC2_variant_path = start_download(tmp_dir, key)
    else:
        _SC2_variant_path = args.lineage_status_table
    sc2_anno_df = pd.read_csv(_SC2_variant_path, header=0)
    # Read our annotation file
    if os.path.isfile(_our_annotation_path):
        print("Read our annotation file")
        our_anno_df = pd.read_csv(_our_annotation_path, header=0)
        # print(our_anno_df)
    else:
        print("Not Found:", _our_annotation_path, ", we create new one")
        our_anno_df = pd.DataFrame(
            columns=[
                "Label",
                "Status",
                "EcdcLabel",
                "Date added",
                "Comment",
                "definition_flexible",
                "mutation",
                "label",
                "Pango lineage",
            ]
        )

    # convert the other status into VOCAL format
    # WHO/CDC/ECDC/PHE Variants of Concern (VOC) -VOC
    # WHO/CDC/ECDC Variants of Interest (VOI) - VOI
    # WHO Alerts for Further Monitoring (AFM) - Monitoring
    # CDC Variants Being Monitored (VBM) - Monitoring
    # ECDC Variants Under Monitoring (VUM) - Monitoring
    # PHE Variants Under Investigation (VUI) - Monitoring
    sc2_anno_df["type"].replace([np.nan, "FMV"], "De-escalated", inplace=True)
    sc2_anno_df["type"].replace(
        ["AFM", "VBM", "VUM", "VUI"], "Monitoring", inplace=True
    )
    sc2_anno_df["type"].replace(["VOC-SUM"], "VOC", inplace=True)
    # change column names

    sc2_anno_df.rename(
        columns={"pango": "Pango lineage", "type": "Status"}, inplace=True
    )
    sc2_anno_df.drop(columns=["label", "interest"], inplace=True)
    sc2_anno_df["Label"] = sc2_anno_df["Pango lineage"]
    # Forloop update
    df = pd.concat([our_anno_df, sc2_anno_df])
    df.drop_duplicates(["Label"], keep="last", inplace=True)

    # Fix lineage
    if args.lineage_list:
        df = match_star_lineage(df, args.lineage_list)
    else:
        pass

    print("Write Result")
    write_to_file(df, output, ",")
    version_file = os.path.join(ROOT_DIR, "data/.version")
    update_version("ECDC assigned variants", get_current_date(), version_file)


def main_mutation(args, key):
    _online = args.online
    _our_annotation_path = args.input_annotation_file
    tmp_dir = args.tmp_directory
    output = args.out_file

    if not os.path.exists(tmp_dir):
        os.makedirs(tmp_dir)

    # if online
    if _online:
        _NetaZuckerman_path = start_download(tmp_dir, key)
    else:
        _NetaZuckerman_path = args.netazuckerman

    # Read our annotation file
    if os.path.isfile(_our_annotation_path):
        print("Read our annotation file")
        our_anno_df = pd.read_csv(_our_annotation_path, sep="\t", header=0)
        print(our_anno_df)
        our_anno_df.drop(
            ["section", "REF", "date of detection", "significance"],
            axis=1,
            inplace=True,
            errors="ignore",
        )
    else:
        print("Not Found:", _our_annotation_path, ", create new one")
        our_anno_df = pd.DataFrame(
            columns=["gene", "amino acid", "nucleotide", "type", "ID", "comment"]
        )

    # Updating the SC2 mutation from NetaZuckerman
    if os.path.isfile(_NetaZuckerman_path):
        print("Update the SC2 mutation from NetaZuckerman ")
        final_NetaZuckerman_dfs_full_table = NetaZuckerman(_NetaZuckerman_path)
        print(final_NetaZuckerman_dfs_full_table)
        # update_version()
        version_file = os.path.join(ROOT_DIR, "data/.version")
        update_version(
            "SARS-CoV-2 mutation information", get_current_date(), version_file
        )
    else:
        print("Not Found:", _NetaZuckerman_path)
        final_NetaZuckerman_dfs_full_table = pd.DataFrame(
            columns=["gene", "amino acid", "nucleotide", "type", "ID", "comment"]
        )

    # Updating MOCs
    # currently and unfortunetly, we don't support for update MOCs

    # Write result
    final_our_anno_df = pd.concat(
        [our_anno_df, final_NetaZuckerman_dfs_full_table],
        axis=0,
        ignore_index=True,
        sort=False,
    )
    final_our_anno_df = final_our_anno_df.drop_duplicates(
        subset=["amino acid", "type", "ID"], keep="first"
    )  # remove duplicated rows
    final_our_anno_df.nucleotide = final_our_anno_df.nucleotide.fillna("?")
    write_to_file(final_our_anno_df, output)

def main_mutation_PS(args):
    _our_annotation_path = args.input_annotation_file
    _our_mutation_table_ps = args.mutation_table_ps
    output = args.out_file
    df = pd.read_csv(_our_annotation_path, sep="\t", header=0)
    df2 = pd.read_csv(_our_mutation_table_ps, sep="\t", header=0)
    _tmp_df = df2[(df2.FEL < 0.05) | (df2.MEME < 0.05)].reset_index(drop=True)
    _tmp_df["Position"] = _tmp_df.Position.apply(int)

    # map REF
    for index, row in _tmp_df.iterrows():
        _tmp_df.loc[index, "ref"] = S_proteinseq[row.Position - 1]

    _tmp_df.AA_Counts = _tmp_df.AA_Counts.astype("str").apply(
        lambda x: ast.literal_eval(x)
    )
    _final_df = pd.DataFrame(
        columns=["gene", "amino acid", "nucleotide", "type", "ID", "comment"]
    )
    lst_dict = []
    for index, row in _tmp_df.iterrows():
        for alt, value in row.AA_Counts.items():
            if row.ref == alt:
                continue
            patter = row.ref + str(row.Position) + alt

            lst_dict.append(
                {
                    "gene": row.Gene,
                    "ID": "",
                    "comment": "",
                    "nucleotide": "?",
                    "type": "PositiveSelection",
                    "amino acid": patter,
                }
            )
            # print("The key and value are ({}) = ({})".format(key, value))
    _final_df = _final_df.append(lst_dict)
    _1_to_go_df = pd.concat([df, _final_df], ignore_index=True)
    _1_to_go_df = _1_to_go_df[_1_to_go_df.type != "MutationOfConcern"]
    _1_to_go_df.to_csv(output, sep="\t", index=False)

if __name__ == "__main__":
    parser = ap.ArgumentParser(
        description=" This tool is used to update DBs of the VOCAL",
        epilog="please visit https://github.com/rki-mf1/sc2-vocal/wiki for the full detail",
    )
    subparsers = parser.add_subparsers(help="sc2-mutation, sc2-lineage-status")
    subparsers.dest = "tool"
    subparsers.required = True

    general_parser = ap.ArgumentParser(add_help=False)
    general_parser.add_argument("-o", "--out_file", required=True, help="Output name")
    general_parser.add_argument(
        "-tmp",
        "--tmp_directory",
        default=".tmp.vocal",
        help="A directory used to hold temporary files, and the directory will automatically cleaned up upon exiting [default: ../sc2-vocal.tmp/ ]",
    )
    general_parser.add_argument(
        "-i",
        "--input_annotation_file",
        default="",
        required=False,
        help="If original annotation file is provided (e.g., table_cov2_mutations_annotation.tsv), it will update the given annotation file otherwise it will create a new one",
    )
    general_parser.add_argument(
        "-a",
        "--online",
        help="Download from annotation sources (always get a latest data from source)",
        action="store_true",
    )

    parser_sc2_mutation = subparsers.add_parser(
        "sc2-mutation", parents=[general_parser], help="update defined mutation."
    )
    parser_sc2_mutation.add_argument(
        "-n",
        "--mutation-table",
        required=False,
        help="please visit https://github.com/NetaZuckerman/covid19/blob/master/mutationsTable.xlsx (NetaZuckerman_mutationsTable.xlsx)",
        default="mutationsTable.xlsx",
    )
    parser_sc2_PS_mutation = subparsers.add_parser(
        "sc2-mutation-ps",
        parents=[general_parser],
        help="update positive selection site.",
    )
    parser_sc2_PS_mutation.add_argument(
        "-ps",
        "--mutation-table-ps",
        required=True,
        help="please see data/PS-desh-jan-march.csv",
    )
    parser_sc2_lineage_status = subparsers.add_parser(
        "sc2-lineage-status", parents=[general_parser], help="update lineage status."
    )
    parser_sc2_lineage_status.add_argument(
        "-l",
        "--lineage-status-table",
        required=False,
        help="please visit https://github.com/3dgiordano/SARS-CoV-2-Variants/blob/main/data/variants.csv (variants.csv)",
        default="variants.csv",
    )
    parser_sc2_lineage_status.add_argument(
        "-L",
        "--lineage_list",
        required=False,
        help="please find lineage.all.tsv in data directory",
    )
    if len(sys.argv) == 1:
        parser.print_help()
        print("Example of usage:")
        print(
            "python vocal/update.vocalDB.py sc2-mutation -i data/table_cov2_mutations_annotation.tsv --online --out_file data/table_cov2_mutations_annotation.new.tsv \n"
        )
        print(
            "python vocal/update.vocalDB.py sc2-mutation -i data/table_cov2_mutations_annotation.tsv -n data/NetaZuckerman_mutationsTable.xlsx  --out_file data/table_cov2_mutations_annotation.new.tsv \n"
        )
        print(
            "python vocal/update.vocalDB.py sc2-mutation-ps -i data/table_cov2_mutations_annotation.tsv -ps data/PS-desh-jan-march.csv --out_file data/table_cov2_mutations_annotation.PS.tsv \n"
        )
        print(
            "python vocal/update.vocalDB.py sc2-lineage-status --online -o data/ECDC_assigned_variants.new.csv -i data/ECDC_assigned_variants.csv -L data/lineage.all.tsv"
        )
        print("python vocal/update.vocalDB.py -h")
        sys.exit()
    args = parser.parse_args()

    tmp_dir = args.tmp_directory
    if not os.path.exists(tmp_dir):
        os.makedirs(tmp_dir)

    # print(args)
    t1 = ti.default_timer()
    if args.tool == "sc2-mutation":
        main_mutation(args, "NetaZuckerman")
    elif args.tool == "sc2-lineage-status":
        main_lineage_status(args, "SC2-Variants")
    elif args.tool == "sc2-mutation-ps":
        main_mutation_PS(args)

    t2 = ti.default_timer()
    # clean tmp directory
    try:
        shutil.rmtree(tmp_dir)
    except OSError as e:
        print("Error: %s - %s." % (e.filename, e.strerror))

    print("Processing time: {} seconds".format(round(t2 - t1, 2)))
