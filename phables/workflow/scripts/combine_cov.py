#!/usr/bin/python3

"""combine_cov.py: Combine multiple coverage files of samples.

"""

import glob
import logging
import os
import subprocess

import pandas as pd

__author__ = "Vijini Mallawaarachchi"
__copyright__ = "Copyright 2023, Phables Project"
__license__ = "MIT"
__type__ = "Support Script"
__maintainer__ = "Vijini Mallawaarachchi"
__email__ = "viji.mallawaarachchi@gmail.com"


def main():
    # Get arguments
    # -----------------------

    covpath = snakemake.params.covpath
    output_path = snakemake.params.output
    log = snakemake.params.log

    # Setup logger
    # ----------------------------------------------------------------------

    logger = logging.getLogger("combine_cov")
    logger.setLevel(logging.DEBUG)
    logging.captureWarnings(True)
    formatter = logging.Formatter("%(asctime)s - %(levelname)s - %(message)s")
    consoleHeader = logging.StreamHandler()
    consoleHeader.setFormatter(formatter)
    consoleHeader.setLevel(logging.INFO)
    logger.addHandler(consoleHeader)

    # Setup output path for log file
    if log is None:
        fileHandler = logging.FileHandler(f"{output_path}/combine_cov.log")
    else:
        fileHandler = logging.FileHandler(f"{log}")

    fileHandler.setLevel(logging.DEBUG)
    fileHandler.setFormatter(formatter)
    logger.addHandler(fileHandler)

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
        logger.info(f"Reading file {file}")
        df = pd.read_csv(file, sep="\t", header=0)

        if final_df.empty:
            final_df = df
        else:
            final_df = pd.concat(
                [final_df, df[list(df.columns)[1]]], axis=1, join="inner"
            )

    logger.info(f"Dataframe shape: {final_df.shape}")

    # Save dataframe to file
    final_df.to_csv(output_path + "coverage.tsv", sep="\t", index=False)
    logger.info(
        f"The combined coverage values can be found at {output_path}coverage.tsv"
    )

    # Exit program
    # --------------

    logger.info("Thank you for using combine_cov!")


if __name__ == "__main__":
    main()
