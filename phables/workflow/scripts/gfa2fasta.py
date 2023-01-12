#!/usr/bin/python3

"""gfa2fasta.py: Obtain the sequences corresponding to edges in the Flye and Miniasm assembly graphs in FASTA format.

The assembly graph file of Flye (assembly_graph.gfa) should be provided as inputs.

"""

import os
import re
import subprocess
import sys

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

__author__ = "Vijini Mallawaarachchi"
__copyright__ = "Copyright 2022, Phables Project"
__license__ = "MIT"
__type__ = "Support Script"
__maintainer__ = "Vijini Mallawaarachchi"
__email__ = "viji.mallawaarachchi@gmail.com"


def main():

    # Get arguments
    # -----------------------

    assembler = snakemake.params.assembler
    assembler_name = ""
    assembly_graph_file = snakemake.params.graph
    output_path = snakemake.params.output
    prefix = ""

    # Check assembly graph file
    if not os.path.isfile(assembly_graph_file):
        print("\nFailed to open the assembly graph file.")
        print("Exiting gfa2fasta.py...\nBye...!\n")
        sys.exit(1)

    # Check assembler type
    if assembler.lower() == "flye":
        assembler_name = "Flye"
    elif assembler.lower() == "miniasm":
        assembler_name = "Miniasm"

    # Check if output folder exists
    # ---------------------------------------------------

    # Handle for missing trailing forwardslash in output folder path
    if output_path[-1:] != "/":
        output_path = output_path + "/"

    # Create output folder if it does not exist
    if not os.path.isdir(output_path):
        subprocess.run("mkdir -p " + output_path, shell=True)

    # Get the sequences corresponding to edges of the graph.
    # ---------------------------------------------------

    print("\nObtaining edge sequences")

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

    print("\nWriting edge sequences to FASTA file")

    if assembler.lower() == "flye":
        final_file = "edges.fasta"
    elif assembler.lower() == "miniasm":
        final_file = "unitigs.fasta"

    with open(output_path + prefix + final_file, "w") as output_handle:
        SeqIO.write(sequenceset, output_handle, "fasta")

    print(
        "\nThe FASTA file with",
        assembler_name,
        "sequences can be found at",
        output_handle.name,
    )

    # Exit program
    # --------------

    print("\nThank you for using gfa2fasta!\n")


if __name__ == "__main__":
    main()
