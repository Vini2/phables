#!/usr/bin/python3

"""gfa2fasta.py: Obtain the sequences corresponding to edges in the Flye and Miniasm assembly graphs in FASTA format.

The assembly graph file of Flye (assembly_graph.gfa) should be provided as inputs.

"""

import logging
import os
import re
import subprocess
import sys

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

__author__ = "Vijini Mallawaarachchi"
__copyright__ = "Copyright 2023, Phables Project"
__license__ = "MIT"
__type__ = "Support Script"
__maintainer__ = "Vijini Mallawaarachchi"
__email__ = "viji.mallawaarachchi@gmail.com"


def main():
    # Get arguments
    # -----------------------

    assembly_graph_file = snakemake.params.graph
    output_path = snakemake.params.output
    log = snakemake.params.log
    prefix = ""

    # Setup logger
    # ----------------------------------------------------------------------

    logger = logging.getLogger("gfa2fasta")
    logger.setLevel(logging.DEBUG)
    logging.captureWarnings(True)
    formatter = logging.Formatter("%(asctime)s - %(levelname)s - %(message)s")
    consoleHeader = logging.StreamHandler()
    consoleHeader.setFormatter(formatter)
    consoleHeader.setLevel(logging.INFO)
    logger.addHandler(consoleHeader)

    # Setup output path for log file
    if log is None:
        fileHandler = logging.FileHandler(f"{output_path}/gfa2fasta.log")
    else:
        fileHandler = logging.FileHandler(f"{log}")

    fileHandler.setLevel(logging.DEBUG)
    fileHandler.setFormatter(formatter)
    logger.addHandler(fileHandler)

    # Check assembly graph file
    if not os.path.isfile(assembly_graph_file):
        logger.error(
            "Failed to open the assembly graph file. Please make sure to provife the .gfa file."
        )
        logger.info("Exiting gfa2fasta.py...\nBye...!\n")
        sys.exit(1)

    # Check if output folder exists
    # ---------------------------------------------------

    # Handle for missing trailing forwardslash in output folder path
    if output_path[-1:] != "/":
        output_path = f"{output_path}/"

    # Create output folder if it does not exist
    if not os.path.isdir(output_path):
        subprocess.run("mkdir -p " + output_path, shell=True)

    # Get the sequences corresponding to edges of the graph.
    # ---------------------------------------------------

    logger.info("Obtaining edge sequences")

    sequenceset = []

    with open(assembly_graph_file) as file:
        line = file.readline()

        while line != "":
            if "S" in line:
                strings = line.split("\t")

                record = SeqRecord(
                    Seq(re.sub("[^GATC]", "", str(strings[2]).upper())),
                    id=str(strings[1]),
                    name=str(strings[1]),
                    description="",
                )

                sequenceset.append(record)

            line = file.readline()

    logger.info("Writing edge sequences to FASTA file")

    with open(f"{output_path}{prefix}edges.fasta", "w") as output_handle:
        SeqIO.write(sequenceset, output_handle, "fasta")

    logger.info(
        f"The FASTA file with unitig sequences can be found at {output_handle.name}"
    )

    # Exit program
    # --------------

    logger.info("Thank you for using gfa2fasta!")


if __name__ == "__main__":
    main()
