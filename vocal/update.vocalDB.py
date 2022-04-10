#!/usr/bin/python
# Maintainer: KongkitimanonK

import argparse as ap
import timeit as ti
import sys
import os
import datetime
import pandas as pd
import requests
import shutil
import re
from utils.utility import ROOT_DIR

# URL
dict_annotation_db = {
    "NetaZuckerman": "https://github.com/NetaZuckerman/covid19/blob/master/mutationsTable.xlsx?raw=true",
    "RKI-VOC-PCR-Finder": "https://www.rki.de/DE/Content/InfAZ/N/Neuartiges_Coronavirus/DESH/Tabelle_VOC-PCR-Finder.xlsx?__blob=publicationFile",
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
        input_NetaZuckerman, sheet_name="bodek  (1)", converters={"Position": int}
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
                if vocal_format != False:
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


def write_to_file(df, file_path):
    # if not os.path.exists(out_dir):
    # Create a new directory because it does not exist
    #    os.makedirs(out_dir)
    df.to_csv(file_path, sep="\t", index=False)


def write_download_to_file(content, out_dir, filename):
    extenion = ".xlsx"
    # if filename == "":
    #   extenion

    file_path = os.path.join(out_dir, filename + extenion)

    if not os.path.exists(out_dir):
        # Create a new directory because it does not exist
        os.makedirs(out_dir)

    with open(file_path, "wb") as f:
        f.write(content)


def start_download(tmp_dir):
    for key, url in dict_annotation_db.items():
        print(key)
        r = requests.get(url)
        if r.status_code == 200:
            filename = key
            write_download_to_file(r.content, tmp_dir, filename)
        else:
            print("Can not download, status code is: " + str(r.status_code))
            raise requests.exceptions.HTTPError


def main(args):
    _online = args.online
    _our_annotation_path = args.input_annotation_file
    tmp_dir = args.tmp_directory
    output = args.out_file

    if not os.path.exists(tmp_dir):
        os.makedirs(tmp_dir)

    # if online
    if _online:
        start_download(tmp_dir)
        _NetaZuckerman_path = os.path.join(tmp_dir, "NetaZuckerman")
        _RKI_VOC_PCR_Finder_path = os.path.join(tmp_dir, "RKI-VOC-PCR-Finder")
        pass
    else:
        _NetaZuckerman_path = args.netazuckerman
        _RKI_VOC_PCR_Finder_path = args.tabelle_VOC_PCR_Finder

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
    else:
        print("Not Found:", _NetaZuckerman_path)
        final_NetaZuckerman_dfs_full_table = pd.DataFrame(
            columns=["gene", "amino acid", "nucleotide", "type", "ID", "comment"]
        )

    # will add lineages UPDATER.py
    # if os.path.isfile(_RKI_VOC_PCR_Finder_path):
    #    new_RKI_VOC_PCR_Finder_dfs = RKI_VOC_FINDER(_RKI_VOC_PCR_Finder_path)
    # else:
    #    print('Not Found:',_RKI_VOC_PCR_Finder_path)
    #    new_RKI_VOC_PCR_Finder_dfs = pd.DataFrame(columns=['gene','amino acid','nucleotide','type','ID',"comment"])

    # Updating MOCs
    ## currently and unfortunetly, we don't support for update MOCs

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

    # clean tmp directory
    try:
        shutil.rmtree(tmp_dir)
    except OSError as e:
        print("Error: %s - %s." % (e.filename, e.strerror))


if __name__ == "__main__":
    parser = ap.ArgumentParser(
        description=" This tool used for update lineage-defining mutations of the Vocal's knowledge base only",
    )
    parser.add_argument(
        "-tmp",
        "--tmp_directory",
        default="../sc2-vocal.tmp/",
        help="A directory used to hold temporary files, and the directory will automatically cleaned up upon exiting [default: ../sc2-vocal.tmp/ ]",
    )
    parser.add_argument(
        "-i",
        "--input_annotation_file",
        required=False,
        help="If original annotation file is provided (table_cov2_mutations_annotation.tsv, tsv format), it will update the given annotation file, otherwise it will create a new one",
        default="../data/table_cov2_mutations_annotation.tsv",
    )
    parser.add_argument(
        "-a",
        "--online",
        help="Download from annotation sources instead using a given path in command line (always get a latest data from source)",
        action="store_true",
    )
    parser.add_argument(
        "-n",
        "--netazuckerman",
        required=False,
        help="please visit https://github.com/NetaZuckerman/covid19/blob/master/mutationsTable.xlsx (NetaZuckerman_mutationsTable.xlsx)",
        default="../data/table_cov2_mutations_annotation.csv",
    )
    parser.add_argument(
        "-t",
        "--tabelle_VOC_PCR_Finder",
        required=False,
        help=" please visit https://www.rki.de/DE/Content/InfAZ/N/Neuartiges_Coronavirus/DESH/Tabelle_VOC-PCR-Finder.xlsx?__blob=publicationFile (Tabelle_VOC-PCR-Finder.xlsx)",
        default="../data/table_cov2_mutations_annotation.csv",
    )

    parser.add_argument(
        "-o",
        "--out_file",
        required=True,
        help="Ouput in tsv format",
        default="../data/table_cov2_mutations_annotation.tsv",
    )

    if len(sys.argv) == 1:
        parser.print_help()
        print('Example of usage:')
        print('python update.vocalDB.py -i ../data/table_cov2_mutations_annotation.tsv --online --out_file ../data/table_cov2_mutations_annotation.tsv \n')
        print('python update.vocalDB.py -i data/table_cov2_mutations_annotation.tsv -n ../vocal-test-annotation_DB/NetaZuckerman_mutationsTable.xlsx  --out_file data/table_cov2_mutations_annotation.1.tsv \n')
        print('python update.vocalDB.py -h')
        sys.exit()
    args = parser.parse_args()
    # print(args)
    t1 = ti.default_timer()
    main(args)
    t2 = ti.default_timer()
    print("Processing time: {} seconds".format(round(t2 - t1), 2))


def NetaZuckerman_old(input_NetaZuckerman):
    NetaZuckerman_dfs = pd.read_excel(
        input_NetaZuckerman, sheet_name=None, converters={"Position": int}
    )
    all_sheets = []
    # construct one big dataframe from multiple sheets
    for name, sheet in NetaZuckerman_dfs.items():
        sheet["ID"] = name
        all_sheets.append(sheet)
    NetaZuckerman_dfs_full_table = pd.concat(all_sheets)
    NetaZuckerman_dfs_full_table.reset_index(inplace=True, drop=True)
    NetaZuckerman_dfs_full_table = NetaZuckerman_dfs_full_table.rename(
        columns={"protein": "gene"}
    )
    NetaZuckerman_dfs_full_table = NetaZuckerman_dfs_full_table.drop_duplicates(
        subset=["gene", "Mutation type", "varname", "ID"], keep="first"
    ).reset_index(drop=True)
    # NetaZuckerman_dfs_full_table[NetaZuckerman_dfs_full_table['Mutation type'].notna()]
    # print(NetaZuckerman_dfs_full_table)
    listOfMutation = [
        "SNP",
        "deletion",
        "SNP_silent",
        "extragenic",
        "deletion_frameshift",
        "SNP_stop",
    ]
    for i, row in NetaZuckerman_dfs_full_table.iterrows():  # fix format

        if (
            row["Mutation type"] == "deletion"
            or row["Mutation type"] == "deletion_frameshift"
        ):
            linage = row["ID"].strip().replace("-", " ").split(" ")[0]
            if row["nuc sub"] == "" or pd.isnull(
                row["nuc sub"]
            ):  # temp solution for now
                if (
                    row["ID"] == "AY.121- Mediator (Delta D)"
                    or row["ID"] == "AY.121 (Delta D based)"
                ):
                    if row["gene"] == "S":
                        new_reps = "del:22029:6"
                    elif row["gene"] == "ORF8":
                        new_reps = "del:28248:6"

                    ## fix AA deletion format....
                    vocal_format = deletion_fixAAformat(
                        new_reps, row["variant"], linage
                    )
                    if vocal_format != False:
                        NetaZuckerman_dfs_full_table.at[i, "variant"] = vocal_format

                    NetaZuckerman_dfs_full_table.at[i, "nuc sub"] = new_reps

                else:
                    print(row)
                    print("cannot fix this missing nuc sub at > ", row["ID"])
            else:
                # print(row['ID'], row["Position"], row["nuc sub"])
                new_reps = (
                    "del:"
                    + str(row["Position"])
                    + ":"
                    + str(sum(c.isalpha() for c in row["nuc sub"]))
                )
                # print(new_reps, " ",row["ID"] )
                ## fix AA deletion format....
                vocal_format = deletion_fixAAformat(new_reps, row["variant"], linage)
                if vocal_format != False:
                    NetaZuckerman_dfs_full_table.at[i, "variant"] = vocal_format

                NetaZuckerman_dfs_full_table.at[i, "nuc sub"] = new_reps

        elif "SNP" in row["Mutation type"]:  # if SNP
            if row["nuc sub"] == "" or pd.isnull(row["nuc sub"]):
                new_reps = row["Reference"] + str(row["Position"]) + row["Mutation"]
                # print(new_reps, " ",row["ID"] )
                NetaZuckerman_dfs_full_table.at[i, "nuc sub"] = new_reps
        elif row["Mutation type"] not in listOfMutation:
            print(row["Mutation type"])

    NetaZuckerman_dfs_full_table["type"] = "LineageDefiningMutation"
    NetaZuckerman_dfs_full_table.rename(
        columns={"variant": "amino acid", "nuc sub": "nucleotide"}, inplace=True
    )
    NetaZuckerman_dfs_full_table = NetaZuckerman_dfs_full_table[
        ["gene", "amino acid", "nucleotide", "ID", "type"]
    ]
    # remove any other except 'S' gene
    NetaZuckerman_dfs_full_table = NetaZuckerman_dfs_full_table[
        NetaZuckerman_dfs_full_table["gene"] == "S"
    ]
    NetaZuckerman_dfs_full_table["comment"] = ""

    # clean duplicate again, keep the first one
    NetaZuckerman_dfs_full_table = NetaZuckerman_dfs_full_table.drop_duplicates(
        subset=["gene", "amino acid", "nucleotide", "ID"], keep="first"
    ).reset_index(drop=True)

    # Clean ID tag
    # we have seperate PDI and non-PDI lineage into two dataframes then, we process each one of them individually
    # For other lineage except PDI*
    _df_1 = NetaZuckerman_dfs_full_table.loc[
        ~NetaZuckerman_dfs_full_table["ID"].str.contains("PDI*", regex=True)
    ]
    _df_1 = _df_1.reset_index(drop=True)
    _df_1["ID"] = _df_1["ID"].str.strip()  # remove white space
    _df_1["ID"] = _df_1["ID"].str.replace("-", " ")
    _df_1["ID"] = _df_1["ID"].str.split(" ").str[0]
    # For PDI <---- we have a little problem, we dont know much about them
    _df_PDI = NetaZuckerman_dfs_full_table.loc[
        NetaZuckerman_dfs_full_table["ID"].str.contains("PDI*", regex=True)
    ]
    _df_PDI = _df_PDI.reset_index(drop=True)
    final_NetaZuckerman_dfs_full_table = _df_1.append(
        _df_PDI, sort=False, ignore_index=True
    )
    return final_NetaZuckerman_dfs_full_table


def RKI_VOC_FINDER(input_RKI_VOC_PCR_Finder):
    RKI_VOC_PCR_Finder_dfs = pd.read_excel(
        input_RKI_VOC_PCR_Finder, sheet_name="VOC PCR Finder", index_col=0
    )
    RKI_VOC_PCR_Finder_dfs = RKI_VOC_PCR_Finder_dfs.dropna(axis=1, how="all")
    # drop last column
    RKI_VOC_PCR_Finder_dfs.drop(
        columns=RKI_VOC_PCR_Finder_dfs.columns[-1], axis=1, inplace=True
    )

    # replace '-' with Zero
    RKI_VOC_PCR_Finder_dfs = RKI_VOC_PCR_Finder_dfs.replace("-", 0)

    RKI_VOC_PCR_Finder_dfs.loc[
        :,
        (
            RKI_VOC_PCR_Finder_dfs.iloc[
                1:,
            ]
            != 0
        ).any(axis=0),
    ]
    new_RKI_VOC_PCR_Finder_dfs = pd.DataFrame(
        columns=["gene", "amino acid", "nucleotide", "type", "ID", "comment"]
    )

    for column in RKI_VOC_PCR_Finder_dfs:
        series_row = RKI_VOC_PCR_Finder_dfs[column]
        # main_lineage = "None"
        for index, value in series_row.items():
            # print(f"Index : {index}, Value : {value}")
            if index == "Main lineage":
                continue
            elif index == "Number of sequences detected":
                continue
            else:  # for all amino acid mutation
                if value == 0:
                    continue  # skip for zero count
                else:
                    new_RKI_VOC_PCR_Finder_dfs = new_RKI_VOC_PCR_Finder_dfs.append(
                        {
                            "amino acid": index,
                            "gene": "?",
                            "nucleotide": "?",
                            "ID": column,
                            "type": "LineageDefiningMutation",
                        },
                        ignore_index=True,
                    )
    new_RKI_VOC_PCR_Finder_dfs["ID"] = (
        new_RKI_VOC_PCR_Finder_dfs["ID"].str.split(" ").str[1]
    )
    new_RKI_VOC_PCR_Finder_dfs["ID"] = new_RKI_VOC_PCR_Finder_dfs[
        "ID"
    ].str.strip()  # remove white space
    return new_RKI_VOC_PCR_Finder_dfs
