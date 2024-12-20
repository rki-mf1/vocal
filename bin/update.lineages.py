#!/usr/bin/python
# Maintainer: KongkitimanonK
# The method originally came from https://github.com/cov-lineages/pango-designation.
# We just adapt and change some parts to be used in covSonar and VOCAL.
import os
import sys
import pandas as pd
import argparse
import datetime
from tempfile import  mkdtemp
import json
import requests
import shutil
from utils.utility import update_version, get_current_date, ROOT_DIR


class Aliasor:
    def __init__(self, alias_file):

        aliases = pd.read_json(alias_file)

        self.alias_dict = {}
        for column in aliases.columns:
            if column.startswith("X"):
                self.alias_dict[column] = column
            else:
                self.alias_dict[column] = aliases[column][0]

        self.alias_dict["A"] = "A"
        self.alias_dict["B"] = "B"

        self.realias_dict = {v: k for k, v in self.alias_dict.items()}

    def compress(self, name):
        name_split = name.split(".")
        # print(name_split)
        if len(name_split) < 5:
            return name
        letter = self.realias_dict[".".join(name_split[0:4])]
        if len(name_split) == 5:
            # print('len5:'+letter + '.' + name_split[4])
            return letter + "." + name_split[4]
        else:
            # print('len6:'+letter + '.' + ".".join(name_split[4:]))
            return letter + "." + ".".join(name_split[4:])

    def uncompress(self, name):
        name_split = name.split(".")
        # print(name_split)
        letter = name_split[0]
        unaliased = self.alias_dict[letter]
        if len(name_split) == 1:
            return name
        if len(name_split) == 2:
            # print('len2:'+unaliased + '.' + name_split[1])
            return unaliased + "." + name_split[1]
        else:
            # print('len3:'+unaliased + '.' + ".".join(name_split[1:]))
            return unaliased + "." + ".".join(name_split[1:])


def lts(lineage):
    items = []
    for item in lineage.split("."):
        item_string = str(item)
        items.append((5 - len(item)) * "0" + item_string)
    return "".join(items)


def download_source(tmp_dir):
    lineages_url = "https://raw.githubusercontent.com/cov-lineages/pango-designation/master/lineages.csv"
    alias_key_url = "https://raw.githubusercontent.com/cov-lineages/pango-designation/master/pango_designation/alias_key.json"
    lineag = os.path.join(tmp_dir, "lineags.csv")
    alias_key = os.path.join(tmp_dir, "alias_key.json")
    print("Download lineags")
    url_content = requests.get(lineages_url).content
    csv_file = open(lineag, "wb")
    csv_file.write(url_content)
    csv_file.close()
    print("Download alias_key")
    items = requests.get(alias_key_url)
    data = items.json()
    with open(alias_key, "w") as f:
        json.dump(data, f)
    return alias_key, lineag


def process_lineage(alias_key_path, lineages_path, output):
    print("Calculate all lineages")
    
    # handle duplicate values
    with open(alias_key_path) as f:
		# load json objects to dictionaries
        data_dict = json.load(f)

    for k, v in data_dict.items():
        if type(v) is list:
            data_dict[k] = list(set(v))
	# rewrite the json
    with open(alias_key_path ,'w') as nf:
        json.dump(data_dict, nf)

    aliasor = Aliasor(alias_key_path)
    df_lineages = pd.read_csv(lineages_path)
    lineages = df_lineages.lineage.unique()
    #%%
    uncompressed_lineages = []
    sorted_lineages = []
    print("Calculate parent-child relationship")
    for ch in map(aliasor.uncompress, lineages):
        uncompressed_lineages.append(ch)
    uncompressed_lineages.sort(key=lts)
    for ch in map(aliasor.compress, uncompressed_lineages):
        sorted_lineages.append(ch)
    #%%
    print("To output")
    df = pd.DataFrame()
    for _id in lineages:
        sub_lineage_char = aliasor.realias_dict.get(_id)
        sub_lineage_list = []

        for name_ in sorted_lineages:
            letter = name_.split(".")[0]
            if sub_lineage_char == letter:
                sub_lineage_list.append(name_)

        sub_lineage_list = list(filter((_id).__ne__, sub_lineage_list))
        if len(sub_lineage_list):
            df = df.append(
                {"lineage": _id, "sublineage": ",".join(sub_lineage_list)},
                ignore_index=True,
            )
        else:
            df = df.append({"lineage": _id, "sublineage": "none"}, ignore_index=True)
    df.to_csv(output, sep="\t", index=False)

    # update_version()
    version_file = os.path.join(ROOT_DIR, "data/.version")
    update_version(
        "Lineage and sublineage relationship", get_current_date(), version_file
    )


def main(args):
    if args.online:
        tmp_dirname = mkdtemp(prefix=".tmp_")
        print("tmp. directory:", tmp_dirname)
        alias_key, lineage = download_source(tmp_dirname)
        process_lineage(alias_key, lineage, args.output)

        if os.path.isdir(tmp_dirname):
            shutil.rmtree(tmp_dirname)
    else:
        process_lineage(args.alias_key, args.lineages, args.output)
    print("Done")


if __name__ == "__main__":
    today = datetime.date.today().strftime("%Y-%m-%d")
    parser = argparse.ArgumentParser(
        description="This script is used for creating a lineage/sublineage relationship file.(https://github.com/cov-lineages/pango-designation)",
        epilog=""" Usage example: \n python Lineages.UPDATER.py --online """,
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument(
        "-l",
        "--lineages",
        default="lineages.csv",
        help="Download from https://github.com/cov-lineages/pango-designation/blob/master/lineages.csv?raw=true",
    )
    parser.add_argument(
        "-a",
        "--alias_key",
        default="alias_key.json",
        help="Download from https://github.com/cov-lineages/pango-designation/blob/master/pango_designation/alias_key.json",
    )
    parser.add_argument(
        "--online",
        help="Download from annotation sources instead using a given path in command line (If this is enable, it always get latest data from source regardless of a given -l and -a)",
        action="store_true",
    )
    parser.add_argument(
        "-o",
        "--output",
        default="lineage.all.{}.tsv".format(today),
        help="output all lineages/sublineages as a tab delimited file",
    )
    args = parser.parse_args()
    if len(sys.argv) == 1:
        print('Usage example: python update.lineages.py -l lineags.csv -a alias_key.json -o lineages.all.tsv')
        print('python update.lineages.py --online')
        print('python update.lineages.py -h')
    else:
        main(args)
