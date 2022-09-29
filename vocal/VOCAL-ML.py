import argparse
import os
import pickle
import sys

import pandas as pd


def get_feature(feature_PS_mode=False):
    if feature_PS_mode:
        list_of_x_cols = [
            "s_moc_M",
            "s_pm_M",
            "s_roi_M",
            "sum_of_avg.normalized_site_total_escape_M",
        ]
    else:
        list_of_x_cols = [
            "s_moc_M",
            "s_moc_D",
            "s_roi_M",
            "s_pm_M",
            "sum_of_avg.normalized_site_total_escape_M",
        ]
    return list_of_x_cols


DT_PS = "vocal/models/DT.best_estimator.PS.pkl"
with open(DT_PS, "rb") as f:
    model = pickle.load(f)


def main(arge):
    print("VOCAL - Auto mode: Prepare and Predict data")
    input = arge.input
    output = arge.output
    _df = pd.read_csv(input)
    list_of_x_cols = get_feature(True)
    if arge.cluster_file:  # change column name
        _df.rename(
            columns={
                "s_moc_M.avg": "s_moc_M",
                "s_moc_D.avg": "s_moc_D",
                "s_roi_M.avg": "s_roi_M",
                "s_pm_M.avg": "s_pm_M",
            },
            inplace=True,
        )
    if not _df.empty:
        X_test = _df[list_of_x_cols]
        y_pred_test = model.predict(X_test)
        _df.insert(loc=2, column="ML_concern", value=y_pred_test)
        _df["ML_concern"] = _df["ML_concern"].map({1: "yes", 0: "no"})
    else:
        _df.insert(loc=2, column="ML_concern", value="")
    if arge.cluster_file:  # change them back
        _df.rename(
            columns={
                "s_moc_M": "s_moc_M.avg",
                "s_moc_D": "s_moc_D.avg",
                "s_roi_M": "s_roi_M.avg",
                "s_pm_M": "s_pm_M.avg",
            },
            inplace=True,
        )
    _df.to_csv(output, index=False)
    print("Done !!")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="VOCAL automode (V.Alpha)",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument(
        "-i",
        "--input",
        required=True,
        help="We use VOCAL to annotate positive selection feature. The input must be a sample file (e.g., vocal-alerts-samples-all.PS.csv)",
    )
    parser.add_argument(
        "-o",
        "--output",
        required=True,
        help="Output of this prediction",
    )
    parser.add_argument(
        "--cluster-file",
        help="In case we want to predict a cluster file (vocal_alerts_clusters_summaries_all.csv)",
        action="store_true",
    )
    args = parser.parse_args()
    if len(sys.argv) == 1:
        print("Usage example: ")
        print(
            "python vocal/VOCAL-ML.py -i vocal-alerts-samples-all.PS.csv  -o vocal-alerts-samples-all.pred.PS.csv"
        )
    else:
        main(args)
