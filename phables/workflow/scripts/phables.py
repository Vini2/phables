#!/usr/bin/env python3

import logging
import time

from phables_utils import (
    component_utils,
    edge_graph_utils,
    gene_utils,
    long_utils,
    short_utils,
)
from phables_utils.coverage_utils import (
    get_junction_pe_coverage,
    get_sub_path_coverage,
    get_unitig_coverage,
)
from phables_utils.output_utils import (
    init_files,
    write_component_info,
    write_component_phrog_info,
    write_res_genome_info,
    write_unitigs,
)

__author__ = "Vijini Mallawaarachchi"
__copyright__ = "Copyright 2022, Phables Project"
__license__ = "MIT"
__version__ = "1.3.2"
__maintainer__ = "Vijini Mallawaarachchi"
__email__ = "viji.mallawaarachchi@gmail.com"
__status__ = "Stable Release"


# Phables main code
# ----------------------------------------------------------------------


def main():
    # Get arguments
    # ----------------------------------------------------------------------
    graph = snakemake.params.graph
    coverage = snakemake.params.coverage
    bampath = snakemake.params.bampath
    hmmout = snakemake.params.hmmout
    phrogs = snakemake.params.phrogs
    minlength = int(snakemake.params.minlength)
    mincov = int(snakemake.params.mincov)
    compcount = int(snakemake.params.compcount)
    maxpaths = int(snakemake.params.maxpaths)
    mgfrac = float(snakemake.params.mgfrac)
    evalue = float(snakemake.params.evalue)
    seqidentity = float(snakemake.params.seqidentity)
    covtol = float(snakemake.params.covtol)
    alpha = float(snakemake.params.alpha)
    longreads = bool(snakemake.params.longreads)
    prefix = snakemake.params.prefix
    output = snakemake.params.output
    nthreads = int(snakemake.params.nthreads)
    log = snakemake.params.log

    # Setup logger
    # ----------------------------------------------------------------------

    logger = logging.getLogger(f"phables {__version__}")
    logger.setLevel(logging.DEBUG)
    logging.captureWarnings(True)
    formatter = logging.Formatter("%(asctime)s - %(levelname)s - %(message)s")
    consoleHeader = logging.StreamHandler()
    consoleHeader.setFormatter(formatter)
    consoleHeader.setLevel(logging.INFO)
    logger.addHandler(consoleHeader)

    # Setup output path for log file
    if log is None:
        fileHandler = logging.FileHandler(f"{output}/phables.log")
    else:
        fileHandler = logging.FileHandler(f"{log}")
    fileHandler.setLevel(logging.DEBUG)
    fileHandler.setFormatter(formatter)
    logger.addHandler(fileHandler)

    logger.info(
        "Welcome to Phables: from fragmented assemblies to high-quality bacteriophage genomes."
    )

    logger.info(f"Input arguments: ")
    logger.info(f"Assembly graph file: {graph}")
    logger.info(f"Unitig coverage file: {coverage}")
    logger.info(f"BAM files path: {bampath}")
    logger.info(f"Unitig .hmmout file: {hmmout}")
    logger.info(f"Unitig phrog annotations file: {phrogs}")
    logger.info(f"Minimum length of unitigs to consider: {minlength}")
    logger.info(f"Minimum coverage of paths to output: {mincov}")
    logger.info(f"Minimum unitig count to consider a component: {compcount}")
    logger.info(f"Maximum number of paths to resolve for a component: {maxpaths}")
    logger.info(f"Length threshold to consider single copy marker genes: {mgfrac}")
    logger.info(f"Maximum e-value for phrog annotations: {evalue}")
    logger.info(f"Minimum sequence identity for phrog annotations: {seqidentity}")
    logger.info(f"Coverage tolerance for extending subpaths: {covtol}")
    logger.info(f"Coverage multipler for flow interval modelling: {alpha}")
    logger.info(f"Input long reads: {longreads}")
    logger.info(f"Prefix for genome identifiers: {prefix}")
    logger.info(f"Number of threads to use: {nthreads}")
    logger.info(f"Output folder: {output}")

    if prefix is None or prefix == "":
        prefix = ""
    else:
        prefix = f"{prefix}_"

    start_time = time.time()

    # Init files
    # ----------------------------------------------------------------------
    init_files(output)

    # Get assembly graph
    # ----------------------------------------------------------------------
    (
        assembly_graph,
        oriented_links,
        link_overlap,
        unitig_names,
        unitig_names_rev,
        graph_unitigs,
        self_looped_nodes,
        edges_lengths,
    ) = edge_graph_utils.build_assembly_graph(graph)

    logger.info(
        f"Total number of vertices in the assembly graph: {len(assembly_graph.vs)}"
    )
    logger.info(
        f"Total number of links in the assembly graph: {len(assembly_graph.es)}"
    )

    # Get single unitigs
    # ----------------------------------------------------------------------
    circular = edge_graph_utils.get_circular(self_looped_nodes, graph_unitigs)

    # Get unitigs with bacterial single copy marker genes
    # ----------------------------------------------------------------------
    smg_unitigs = gene_utils.get_smg_unitigs(hmmout, mgfrac)

    # Get unitigs with PHROGs
    # ----------------------------------------------------------------------
    unitig_phrogs, phrog_dict = gene_utils.get_phrog_unitigs(
        phrogs, evalue, seqidentity
    )

    # Get components with viral components
    # ----------------------------------------------------------------------
    pruned_vs, comp_phrogs, likely_complete = component_utils.get_components(
        assembly_graph,
        unitig_names,
        smg_unitigs,
        unitig_phrogs,
        circular,
        edges_lengths,
        minlength,
        phrog_dict,
    )
    logger.info(f"Total number of components found: {len(pruned_vs)}")

    # Get unitig coverages
    # ----------------------------------------------------------------------

    unitig_coverages = get_unitig_coverage(coverage)

    # Resolve genomes
    # ----------------------------------------------------------------------

    # If long reads are provided
    if longreads:
        logger.info(f"Long reads provided")

        # Get sub path coverages
        sub_path_cov = edge_graph_utils.get_all_sub_paths(assembly_graph, unitig_names)
        sub_path_cov = get_sub_path_coverage(sub_path_cov, bampath, output)

        # Resolve genomes
        (
            resolved_edges,
            all_resolved_paths,
            all_components,
            cycle_components,
            linear_components,
            resolved_components,
            resolved_linear,
            single_unitigs,
            resolved_cyclic,
            case1_found,
            case1_resolved,
            case2_found,
            case2_resolved,
            case3_found,
            case3_resolved,
            phage_like_edges,
            all_phage_like_edges,
            unresolved_phage_like_edges,
        ) = long_utils.resolve_long(
            assembly_graph,
            pruned_vs,
            unitig_names,
            unitig_names_rev,
            self_looped_nodes,
            graph_unitigs,
            minlength,
            link_overlap,
            unitig_coverages,
            compcount,
            oriented_links,
            sub_path_cov,
            likely_complete,
            alpha,
            mincov,
            covtol,
            maxpaths,
            prefix,
            output,
            nthreads,
        )

    # Else default to short reads
    else:
        logger.info(f"Short reads provided")

        # Get junction pe coverages
        junction_pe_coverage = get_junction_pe_coverage(bampath, output)

        # Resolve genomes
        (
            resolved_edges,
            all_resolved_paths,
            all_components,
            cycle_components,
            linear_components,
            resolved_components,
            resolved_linear,
            single_unitigs,
            resolved_cyclic,
            case1_found,
            case1_resolved,
            case2_found,
            case2_resolved,
            case3_found,
            case3_resolved,
            phage_like_edges,
            all_phage_like_edges,
            unresolved_phage_like_edges,
        ) = short_utils.resolve_short(
            assembly_graph,
            pruned_vs,
            unitig_names,
            unitig_names_rev,
            self_looped_nodes,
            graph_unitigs,
            minlength,
            link_overlap,
            unitig_coverages,
            compcount,
            oriented_links,
            junction_pe_coverage,
            likely_complete,
            alpha,
            mincov,
            covtol,
            maxpaths,
            prefix,
            output,
            nthreads,
        )

    # Log final summary information
    # ----------------------------------------------------------------------
    logger.info(f"Total number of cyclic components found: {len(cycle_components)}")
    logger.info(f"Total number of cyclic components resolved: {len(resolved_cyclic)}")
    logger.info(f"Single unitigs identified: {len(single_unitigs)}")
    logger.info(f"Total number of linear components found: {len(linear_components)}")
    logger.info(f"Total number of linear components resolved: {len(resolved_linear)}")
    logger.info(
        f"Total number of cyclic components found including single unitigs: {len(cycle_components) + len(single_unitigs)}"
    )
    logger.info(
        f"Total number of components resolved: {len(single_unitigs)+len(resolved_cyclic)+len(resolved_linear)}"
    )
    logger.info(f"Case 1 (resolved/found): {len(case1_resolved)}/{len(case1_found)}")
    logger.info(f"Case 2 (resolved/found): {len(case2_resolved)}/{len(case2_found)}")
    logger.info(f"Case 3 (resolved/found): {len(case3_resolved)}/{len(case3_found)}")
    logger.info(f"Total number of genomes resolved: {len(all_resolved_paths)}")

    if len(all_resolved_paths) == 0:
        logger.info(f"No genomes were resolved.")
    else:
        logger.info(f"Resolved genomes can be found in {output}/resolved_paths.fasta")

    # Write edges to file
    # ----------------------------------------------------------------------

    write_unitigs(
        phage_like_edges, unitig_names, graph_unitigs, "phage_like_edges", output
    )
    write_unitigs(
        all_phage_like_edges,
        unitig_names,
        graph_unitigs,
        "all_phage_like_edges",
        output,
    )
    write_unitigs(resolved_edges, unitig_names, graph_unitigs, "resolved_edges", output)
    write_unitigs(
        unresolved_phage_like_edges,
        unitig_names,
        graph_unitigs,
        "unresolved_phage_like_edges",
        output,
    )

    # Record path information
    # ----------------------------------------------------------------------

    filename = write_res_genome_info(all_resolved_paths, output)
    if len(all_resolved_paths) > 0:
        logger.info(f"Resolved genome information can be found in {output}/{filename}")

    # Record component information
    # ----------------------------------------------------------------------

    filename = write_component_info(all_components, output)
    if len(all_components) > 0:
        logger.info(
            f"Resolved component information can be found in {output}/{filename}"
        )

    filename = write_component_phrog_info(resolved_components, comp_phrogs, output)
    if len(resolved_components) > 0:
        logger.info(
            f"PHROGs found in resolved components can be found in {output}/{filename}"
        )

    # Get elapsed time
    # ----------------------------------------------------------------------

    # Determine elapsed time
    elapsed_time = time.time() - start_time

    # Print elapsed time for the process
    logger.info(f"Elapsed time: {elapsed_time} seconds")

    # Exit program
    # ----------------------------------------------------------------------

    logger.info("Thank you for using Phables!")


if __name__ == "__main__":
    main()
