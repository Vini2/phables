#!/usr/bin/python3

"""format_koverage_results.py: Format koverage results.

"""

import logging
import os
import subprocess
from collections import defaultdict

import pandas as pd
from Bio import SeqIO

__author__ = "Vijini Mallawaarachchi"
__copyright__ = "Copyright 2023, Phables Project"
__license__ = "MIT"
__type__ = "Support Script"
__maintainer__ = "Vijini Mallawaarachchi"
__email__ = "viji.mallawaarachchi@gmail.com"


def main():
    # Get arguments
    # -----------------------

    samples_file = snakemake.params.samples_file
    koverage_tsv = snakemake.params.koverage_tsv
    seq_file = snakemake.params.seq_file
    info_file = snakemake.params.info_file
    output_path = snakemake.params.output_path
    log = snakemake.params.log

    # Setup logger
    # ----------------------------------------------------------------------

    logger = logging.getLogger("format_coverage")
    logger.setLevel(logging.DEBUG)
    logging.captureWarnings(True)
    formatter = logging.Formatter("%(asctime)s - %(levelname)s - %(message)s")
    consoleHeader = logging.StreamHandler()
    consoleHeader.setFormatter(formatter)
    consoleHeader.setLevel(logging.INFO)
    logger.addHandler(consoleHeader)

    # Setup output path for log file
    if log is None:
        fileHandler = logging.FileHandler(f"{log}")
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

    # Get sample-wise genome coverage stats
    # ----------------------------------------------------------------------

    # Log inputs
    logger.info(f"Samples file: {samples_file}")
    logger.info(f"Koverage results: {koverage_tsv}")
    logger.info(f"Output path: {output_path}")

    # Get sample names
    mysamples = [s.split("\t")[0] for s in open(samples_file, "r")]
    logger.debug(mysamples)

    # Initialise dataframe
    df_read_counts = pd.DataFrame(columns=["contig_phables"] + mysamples)
    df_rpkm = pd.DataFrame(columns=["contig_phables"] + mysamples)
    df_mean_cov = pd.DataFrame(columns=["contig_phables"] + mysamples)

    # Get coverage stats of genomes in each sample
    read_counts = defaultdict(lambda: defaultdict(list))
    rpkm = defaultdict(lambda: defaultdict(list))
    mean_cov = defaultdict(lambda: defaultdict(list))

    with open(koverage_tsv, "r") as mf:
        for line in mf.readlines()[1:]:
            strings = line.strip().split("\t")
            read_counts[strings[1]][strings[0]] = int(float(strings[2]))
            rpkm[strings[1]][strings[0]] = float(strings[4])
            mean_cov[strings[1]][strings[0]] = float(strings[7])

    # Add records to dataframe
    counter = 0
    for genome in read_counts:
        read_counts_row = read_counts[genome]
        read_counts_row["contig_phables"] = genome
        read_counts_row = dict(read_counts_row)
        read_counts_row_df = pd.DataFrame(read_counts_row, index=[counter])
        df_read_counts = pd.concat([df_read_counts, read_counts_row_df])

        rpkm_row = rpkm[genome]
        rpkm_row["contig_phables"] = genome
        rpkm_row = dict(rpkm_row)
        rpkm_row_df = pd.DataFrame(rpkm_row, index=[counter])
        df_rpkm = pd.concat([df_rpkm, rpkm_row_df])

        mean_cov_row = mean_cov[genome]
        mean_cov_row["contig_phables"] = genome
        mean_cov_row = dict(mean_cov_row)
        mean_cov_row_df = pd.DataFrame(mean_cov_row, index=[counter])
        df_mean_cov = pd.concat([df_mean_cov, mean_cov_row_df])

        counter += 1

    # Save dataframe to file
    df_read_counts.to_csv(
        f"{output_path}sample_genome_read_counts.tsv", sep="\t", index=False
    )
    df_rpkm.to_csv(f"{output_path}sample_genome_rpkm.tsv", sep="\t", index=False)
    df_mean_cov.to_csv(
        f"{output_path}sample_genome_mean_coverage.tsv", sep="\t", index=False
    )

    logger.info(
        f"Raw read counts mapped to resolved genomes can be found in {output_path}sample_genome_read_counts.tsv"
    )
    logger.info(
        f"RPKM values of resolved genomes can be found in {output_path}sample_genome_rpkm.tsv"
    )
    logger.info(
        f"Estimated mean read depth of resolved genomes can be found in {output_path}sample_genome_mean_coverage.tsv"
    )

    # Make sequence information file
    with open(info_file, "w") as myfile:
        myfile.write(f"contig_phables_name\tlength\tcontig_or_phables\n")
        for index, record in enumerate(SeqIO.parse(seq_file, "fasta")):
            if "phage_comp" in record.id:
                myfile.write(f"{record.id}\t{len(record.seq)}\tphables\n")
            else:
                myfile.write(f"{record.id}\t{len(record.seq)}\tcontig\n")

    logger.info(f"Sequence information file can be found in {info_file}")

    # Exit program
    # --------------

    logger.info("Thank you for using format_koverage_results!")


if __name__ == "__main__":
    main()
