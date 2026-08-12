import logging
import time

import piqtree
from cogent3 import load_aligned_seqs

__author__ = "Vijini Mallawaarachchi"
__copyright__ = "Copyright 2022, Phables Project"
__license__ = "MIT"
__version__ = "2.0.0"
__maintainer__ = "Vijini Mallawaarachchi"
__email__ = "viji.mallawaarachchi@gmail.com"
__status__ = "Stable Release"


def main():

    # Get arguments
    # ----------------------------------------------------------------------

    aligned = snakemake.params.aligned
    output = snakemake.params.output
    nthreads = int(snakemake.params.nthreads)
    seed = int(snakemake.params.seed)
    log = snakemake.params.log

    # Setup logger
    # ----------------------------------------------------------------------

    logger = logging.getLogger(f"build tree")
    logger.setLevel(logging.DEBUG)
    logging.captureWarnings(True)
    formatter = logging.Formatter("%(asctime)s - %(levelname)s - %(message)s")
    consoleHeader = logging.StreamHandler()
    consoleHeader.setFormatter(formatter)
    consoleHeader.setLevel(logging.INFO)
    logger.addHandler(consoleHeader)

    # Setup output path for log file
    fileHandler = logging.FileHandler(f"{log}")
    fileHandler.setLevel(logging.DEBUG)
    fileHandler.setFormatter(formatter)
    logger.addHandler(fileHandler)

    logger.info(f"Starting phylogenetic tree reconstruction with IQ-TREE from piqtree")
    logger.info(f"Input aligned fasta file: {aligned}")
    logger.info(f"Output tree file: {output}")
    logger.info(f"Number of threads: {nthreads}")
    logger.info(f"Random seed: {seed}")
    logger.info(f"Using piqtree version: {piqtree.__version__}")

    # Load Sequences
    logger.info(f"Loading aligned sequences from {aligned}")
    aln = load_aligned_seqs(aligned, moltype="dna")

    # Find the best model using ModelFinder
    logger.info(f"Finding the best model using ModelFinder")
    results = piqtree.model_finder(aln, rand_seed=seed, num_threads=nthreads)

    # Reconstruct a phylogenetic tree with IQ-TREE
    logger.info(f"Building phylogenetic tree using model: {results.best_bic}")
    tree = piqtree.build_tree(
        aln, results.best_bic, rand_seed=seed, num_threads=nthreads
    )

    # Write tree to newick file
    tree.write(output)
    logger.info(f"Tree file in newick format written to {output}")


if __name__ == "__main__":
    main()
