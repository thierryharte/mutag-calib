#!/usr/bin/env python3

import os
import subprocess
import argparse
import pandas as pd

from allowed_categories import ALLOWED_CATEGORIES


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("base_dir", help="Base directory containing datacards")
    parser.add_argument("--csv-all-results", help="Output summary CSV file", default="ALL_FIT_RESULTS.csv")
    parser.add_argument("--split-c", action="store_true",
                        help="datacards have SF_c split into SF_c_pass (frozen) and SF_c_fail")
    parser.add_argument("--c-freeze", action="store_true",
                        help="if true, feeze all SFc factors")
    parser.add_argument("--only-nominal", action="store_true", help="Only fit the central value. Can be relevant if internalised tau21 is used.")
    args = parser.parse_args()

    BASE_DIR = args.base_dir
    SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
    RUN_FIT = os.path.join(SCRIPT_DIR, "fit_diagnostics.py")
    EXTRACT = os.path.join(SCRIPT_DIR, "extract_fit_results.py")
    output_file = args.csv_all_results
    TAU21_NOMINAL = "0p30"

    extra_args = ["--split-c"] if args.split_c else []
    if args.c_freeze:
        extra_args.append("--c-freeze")

    summary_rows = []
    failed = []

    for year in sorted(os.listdir(BASE_DIR)):
        year_path = os.path.join(BASE_DIR, year)
        if not os.path.isdir(year_path):
            continue

        for category in sorted(os.listdir(year_path)):
            if category not in ALLOWED_CATEGORIES:
                continue

            category_path = os.path.join(year_path, category)

            for cut in sorted(os.listdir(category_path)):
                if args.only_nominal and TAU21_NOMINAL not in cut:
                    continue
                cut_path = os.path.join(category_path, cut)
                if not os.path.isdir(cut_path):
                    continue

                print(f"\n[INFO] {year} | {category} | {cut}")

                # remove stale results so a failed fit is not picked up from a previous run
                csv_file = os.path.join(cut_path, "fitResults.csv")
                if os.path.isfile(csv_file):
                    os.remove(csv_file)

                # 1) Run FitDiagnostics, 2) Extract results
                # a failing category is recorded and skipped instead of stopping the loop
                try:
                    for script in (RUN_FIT, EXTRACT):
                        subprocess.run(
                            ["python3", script] + extra_args,
                            cwd=cut_path,
                            check=True
                        )
                except subprocess.CalledProcessError as err:
                    print(f"[ERROR] {os.path.basename(err.cmd[1])} failed, skipping {year} | {category} | {cut}")
                    failed.append(cut_path)
                    continue

                # 3) Summary rows from single results files
                if os.path.isfile(csv_file):
                    df = pd.read_csv(csv_file)
                    summary_rows.append(df)

    # save global summary
    if summary_rows:
        summary_df = pd.concat(summary_rows, ignore_index=True)
        if os.path.isfile(output_file):
            summary_df.to_csv(output_file, index=False, mode='a', header=False)
        else:
            summary_df.to_csv(output_file, index=False, mode='a')

        print(f"\n[OK] Global summary saved in {output_file}")
    else:
        print("\n[WARN] None result collected")

    if failed:
        print(f"\n[WARN] {len(failed)} fit(s) failed, see fitDiagnostics.*.log in:")
        for path in failed:
            print(f"  {path}")


if __name__ == "__main__":
    main()
