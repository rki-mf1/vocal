#!/usr/bin/python
# Maintainer: KongkitimanonK
# The method originally came from https://github.com/cov-lineages/pango-designation.
# We just adapt and change some parts to be used in covSonar and VOCAL.
import argparse
import datetime
import json
import os
import shutil
import sys
from tempfile import mkdtemp

import pandas as pd
import requests
from utils.utility import get_current_date
from utils.utility import ROOT_DIR
from utils.utility import update_version

try:  # noqa: C901
    from pango_aliasor.aliasor import Aliasor
except ModuleNotFoundError:  # pragma: no cover
    print(
        "Dependency `pango_aliasor` missing, please install using `pip install pango_aliasor`"
    )
    print("Fall back to original Aliasor...")

    class Aliasor:
        def __init__(self, alias_file=None):
            import json

            if alias_file is None:
                import importlib.resources

                with importlib.resources.open_text(
                    "pango_designation", "alias_key.json"
                ) as file:
                    file = json.load(file)

            else:
                with open(alias_file) as file:
                    file = json.load(file)

            self.alias_dict = {}
            for column in file.keys():
                if type(file[column]) is list or file[column] == "":
                    self.alias_dict[column] = column
                else:
                    self.alias_dict[column] = file[column]

            self.realias_dict = {v: k for k, v in self.alias_dict.items()}

        def compress(self, name):
            name_split = name.split(".")
            levels = len(name_split) - 1
            num_indirections = (levels - 1) // 3
            if num_indirections <= 0:
                return name
            alias = ".".join(name_split[0 : (3 * num_indirections + 1)])
            ending = ".".join(name_split[(3 * num_indirections + 1) :])
            return self.realias_dict[alias] + "." + ending

        def uncompress(self, name):
            name_split = name.split(".")
            letter = name_split[0]
            try:
                unaliased = self.alias_dict[letter]
            except KeyError:
                return name
            if len(name_split) == 1:
                return name
            if len(name_split) == 2:
                return unaliased + "." + name_split[1]
            else:
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
    aliasor = Aliasor(alias_key_path)
    df_lineages = pd.read_csv(lineages_path)
    lineages = df_lineages.lineage.unique()

    uncompressed_lineages = []
    sorted_lineages = []

    # Calculating parent-child relationship
    cleanedlineages = [x for x in lineages if str(x) != "nan"]
    uncompressed_lineages = list(map(aliasor.uncompress, cleanedlineages))
    uncompressed_lineages.sort(key=lts)
    sorted_lineages = list(map(aliasor.compress, uncompressed_lineages))

    _final_list = []
    for _id in sorted_lineages:
        alias_lineage_char = aliasor.uncompress(_id)
        sub_lineage_list = []
        row_dict = {}
        # print(_id, '=',alias_lineage_char)

        for name_ in uncompressed_lineages:  # fetch all lineage again
            root = ""
            for index, letter in enumerate(name_.split(".")):
                if index != 0:
                    letter = root + "." + letter
                root = letter
                if letter == alias_lineage_char:
                    sub_lineage_list.append(aliasor.compress(name_))
        # remove root lineage
        sub_lineage_list.remove(_id)
        if len(sub_lineage_list) > 0:
            row_dict["lineage"] = _id
            row_dict["sublineage"] = ",".join(sub_lineage_list)
        else:
            row_dict["lineage"] = _id
            row_dict["sublineage"] = "none"
        _final_list.append(row_dict)

    df = pd.DataFrame.from_dict(_final_list, orient="columns")
    df = df.sort_values(by=["lineage"])
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
        print(
            "Usage example: python update.lineages.py -l lineags.csv -a alias_key.json -o lineages.all.tsv"
        )
        print("python update.lineages.py --online")
        print("python update.lineages.py -h")
    else:
        main(args)
