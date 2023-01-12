#!/usr/bin/python3

"""combine_cov.py: Combine multiple coverage files of samples from CoverM.

"""

import glob
import os
import subprocess

import pandas as pd

__author__ = "Vijini Mallawaarachchi"
__copyright__ = "Copyright 2022, Phables Project"
__license__ = "MIT"
__type__ = "Support Script"
__maintainer__ = "Vijini Mallawaarachchi"
__email__ = "viji.mallawaarachchi@gmail.com"


def main():

    # Get arguments
    # -----------------------

    covpath = snakemake.params.covpath
    output_path = snakemake.params.output

    # Validate inputs
    # ---------------------------------------------------

    # Handle for missing trailing forwardslash in output folder path
    if output_path[-1:] != "/":
        output_path = output_path + "/"

    # Create output folder if it does not exist
    if not os.path.isdir(output_path):
        subprocess.run("mkdir -p " + output_path, shell=True)

    # Get coverage values from samples
    # ---------------------------------------------------

    # Get coverage files
    cov_files = glob.glob(f"{covpath}/*.tsv")

    final_df = pd.DataFrame()

    for file in cov_files:
        df = pd.read_csv(file, sep="\t", header=0)

        if final_df.empty:
            final_df = df
        else:
            final_df = pd.concat(
                [final_df, df[list(df.columns)[1]]], axis=1, join="inner"
            )

    print(f"Dataframe shape: {final_df.shape}")

    # Save dataframe to file
    final_df.to_csv(output_path + "coverage.tsv", sep="\t", index=False)
    print(f"The combined coverage values can be found at {output_path}coverage.tsv")

    # Exit program
    # --------------

    print("Thank you for using combine_cov!")


if __name__ == "__main__":
    main()
