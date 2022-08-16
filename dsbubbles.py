#!/usr/bin/env python3

import logging
import subprocess
import sys
import time

import click
import networkx as nx
from igraph import *
from tqdm import tqdm

from dsbubbles_utils import (component_utils, edge_graph_utils, edge_utils,
                             gene_utils)
from dsbubbles_utils.genome_utils import GenomeComponent, GenomePath

__author__ = "Vijini Mallawaarachchi"
__copyright__ = "Copyright 2022, ds-bubbles Project"
__license__ = "MIT"
__version__ = "0.1"
__maintainer__ = "Vijini Mallawaarachchi"
__email__ = "vijini.mallawaarachchi@flinders.edu.au"
__status__ = "Development"

MAX_VAL = sys.maxsize
FASTA_LINE_LEN = 60

# Sample command
# -------------------------------------------------------------------
# python dsbubbles.py  -g /path/to/assembly_graph.gfa
#                      -c /path/to/assembly.fasta
#                      -p /path/to/assembly_info.txt
#                      -o /path/to/output_folder
# -------------------------------------------------------------------


# Setup arguments
# ----------------------------------------------------------------------


@click.command()
@click.option(
    "--graph",
    "-g",
    required=True,
    help="path to the assembly graph file",
    type=click.Path(exists=True),
)
@click.option(
    "--contigs",
    "-c",
    required=True,
    help="path to the contigs file",
    type=click.Path(exists=True),
)
@click.option(
    "--paths",
    "-p",
    required=True,
    help="path to the contig paths file",
    type=click.Path(exists=True),
)
@click.option(
    "--hmmout",
    "-hm",
    required=True,
    help="path to the contig .hmmout file",
    type=click.Path(exists=True),
)
@click.option(
    "--phrogs",
    "-ph",
    required=True,
    help="path to the contig phrog annotations file",
    type=click.Path(exists=True),
)
@click.option(
    "--minlength",
    "-ml",
    default=2000,
    required=False,
    help="minimum length of circular contigs to consider",
    type=int,
)
@click.option(
    "--biglength",
    "-bl",
    default=10000,
    required=False,
    help="minimum length of a big circular contig",
    type=int,
)
@click.option(
    "--pathdiff",
    "-pd",
    default=2000,
    required=False,
    help="length difference threshold to filter paths of a component",
    type=int,
)
@click.option(
    "--mgfrac",
    "-mgf",
    default=0.5,
    required=False,
    help="length threshold to consider single copy marker genes",
    type=float,
)
@click.option(
    "--alignscore",
    "-as",
    default=90,
    required=False,
    help="minimum alignment score (%) for phrog annotations",
    type=float,
)
@click.option(
    "--seqidentity",
    "-si",
    default=0.5,
    required=False,
    help="minimum sequence identity for phrog annotations",
    type=float,
)
@click.option(
    "--degree",
    "-d",
    default=10,
    required=False,
    help="minimum in/out degree of nodes in a component to be complex",
    type=int,
)
@click.option(
    "--output",
    "-o",
    required=True,
    help="path to the output folder",
    type=click.Path(exists=True),
)
def main(
    graph,
    contigs,
    paths,
    hmmout,
    phrogs,
    minlength,
    biglength,
    pathdiff,
    mgfrac,
    alignscore,
    seqidentity,
    degree,
    output,
):

    """ds-bubbles: Resolve bacteriophage genomes from viral bubbles in metagenomic data."""

    # Setup logger
    # ----------------------------------------------------------------------

    logger = logging.getLogger("dsbubbles 0.1")
    logger.setLevel(logging.DEBUG)
    logging.captureWarnings(True)
    formatter = logging.Formatter("%(asctime)s - %(levelname)s - %(message)s")
    consoleHeader = logging.StreamHandler()
    consoleHeader.setFormatter(formatter)
    consoleHeader.setLevel(logging.INFO)
    logger.addHandler(consoleHeader)

    # Setup output path for log file
    fileHandler = logging.FileHandler(f"{output}/dsbubbles.log")
    fileHandler.setLevel(logging.DEBUG)
    fileHandler.setFormatter(formatter)
    logger.addHandler(fileHandler)

    logger.info(
        "Welcome to ds-bubbles: Resolve bacteriophage genomes from viral bubbles in metagenomic data."
    )

    logger.info(f"Input arguments: ")
    logger.info(f"Assembly graph file: {graph}")
    logger.info(f"Contigs file: {contigs}")
    logger.info(f"Contig paths file: {paths}")
    logger.info(f"Contig .hmmout file: {hmmout}")
    logger.info(f"Contig phrog annotations file: {phrogs}")
    logger.info(f"Minimum length of contigs to consider: {minlength}")
    logger.info(f"Minimum length of a big circular contig: {biglength}")
    logger.info(f"Length difference threshold to filter paths of a component: {pathdiff}")
    logger.info(f"Length threshold to consider single copy marker genes: {mgfrac}")
    logger.info(f"Minimum alignment score (%) for phrog annotations: {alignscore}")
    logger.info(f"Minimum sequence identity for phrog annotations: {seqidentity}")
    logger.info(
        f"Minimum in/out degree of nodes in a component to be complex: {degree}"
    )
    logger.info(f"Output folder: {output}")

    start_time = time.time()

    # Get assembly graph
    # ----------------------------------------------------------------------
    (
        assembly_graph,
        edge_list,
        contig_names,
        contig_names_rev,
        graph_contigs,
        edge_depths,
        self_looped_nodes,
        edges_lengths,
    ) = edge_graph_utils.build_assembly_graph(graph)

    logger.info(
        f"Total number of vertices in the assembly graph: {len(assembly_graph.vs)}"
    )
    logger.info(
        f"Total number of links in the assembly graph: {len(assembly_graph.es)}"
    )

    # Get circular contigs
    # ----------------------------------------------------------------------
    circular = edge_utils.get_circular(paths)

    # Get contigs with bacterial single copy marker genes
    # ----------------------------------------------------------------------
    smg_contigs, contig_smgs = gene_utils.get_smg_contigs(hmmout, mgfrac)

    # Get contigs with PHROGs
    # ----------------------------------------------------------------------
    contig_phrogs = gene_utils.get_phrog_contigs(phrogs, alignscore, seqidentity)

    # Get components with viral bubbles
    # ----------------------------------------------------------------------
    pruned_vs = component_utils.get_components(
        assembly_graph,
        contig_names,
        smg_contigs,
        contig_phrogs,
        circular,
        edges_lengths,
        minlength,
    )
    logger.info(f"Total number of components found: {len(pruned_vs)}")

    # Resolve genomes
    # ----------------------------------------------------------------------

    resolved_edges = set()

    all_resolved_paths = []

    all_components = []

    for my_count in tqdm(pruned_vs, desc="Resolving components"):

        my_genomic_paths = []

        logger.debug(f"my_count: {my_count}")

        pruned_graph = assembly_graph.subgraph(pruned_vs[my_count])

        if len(pruned_vs[my_count]) > 1:

            has_long_circular = False

            # Check if the bubble is complex
            is_complex_graph = False

            for node in pruned_vs[my_count]:

                in_degree_node = assembly_graph.degree(node, mode="in")
                out_degree_node = assembly_graph.degree(node, mode="out")

                if out_degree_node > degree or in_degree_node > degree:
                    is_complex_graph = True

                contig_name = contig_names[node]

                if (
                    contig_name in self_looped_nodes
                    and len(graph_contigs[contig_names[node]]) > biglength
                ):
                    has_long_circular = True
                    cycle_number = 1

                    path_string = graph_contigs[contig_name]

                    genome_path = GenomePath(
                        f"phage_comp_{my_count}_cycle_{cycle_number}",
                        [pruned_vs[my_count][0]],
                        path_string,
                        coverage,
                        len(graph_contigs[contig_name]),
                    )
                    my_genomic_paths.append(genome_path)

            if not has_long_circular:

                # Create Directed Graph
                G = nx.DiGraph()

                # Add a list of nodes:
                named_vertex_list = pruned_graph.vs()["label"]
                G.add_nodes_from(named_vertex_list)

                my_edges = []

                for edge in pruned_graph.es:
                    source_vertex_id = edge.source
                    target_vertex_id = edge.target
                    my_edges.append((source_vertex_id, target_vertex_id))

                # Add a list of edges:
                G.add_edges_from(my_edges)

                path_lengths = []
                path_coverages = []

                mypath_strings = []

                if is_complex_graph:
                    logger.debug("Complex_graph")
                    my_cycles = nx.cycle_basis(G.to_undirected())
                else:
                    logger.debug("Simple graph")
                    my_cycles = list(nx.simple_cycles(G))

                cycle_number = 1

                # Return a list of cycles described as a list of nodes
                for cycle in my_cycles:

                    if len(cycle) > 1:

                        # cycle_family = set()

                        # has_cycles.append(my_count)

                        total_length = 0

                        coverage = MAX_VAL

                        node_order = []

                        path_string = ""

                        for node in cycle:

                            contig_name = contig_names[pruned_vs[my_count][node]]

                            # if contig_name in edge_family:
                            #     cycle_family.add(edge_family[contig_name])
                            # else:
                            #     cycle_family.add("Unclassified")

                            # if pruned_vs[my_count][node] in [2228, 2231, 2234]:
                            #     print(my_count)
                            # print(node, ":", pruned_vs[my_count][node], contig_names[pruned_vs[my_count][node]], ":", edge_depths[contig_names[pruned_vs[my_count][node]]], end=" --> ")

                            resolved_edges.add(contig_names[pruned_vs[my_count][node]])

                            node_order.append(contig_names[pruned_vs[my_count][node]])
                            total_length += len(graph_contigs[contig_name])
                            path_string += graph_contigs[contig_name]

                            if coverage > edge_depths[contig_name]:
                                coverage = edge_depths[contig_name]

                        # print("\ntotal_length:", total_length)
                        path_lengths.append(total_length)
                        # print("path coverage:", coverage)
                        path_coverages.append(coverage)

                        mypath_strings.append(path_string)

                        # Save path

                        if len(path_string) > 0 and coverage != MAX_VAL:

                            genome_path = GenomePath(
                                f"phage_comp_{my_count}_cycle_{cycle_number}",
                                node_order,
                                path_string,
                                coverage,
                                total_length,
                            )
                            my_genomic_paths.append(genome_path)

                            # viral_comp_length["phage_comp_"+str(my_count)+"_cycle_"+str(cycle_number)] = len(path_string)

                            # if len(cycle_family) == 1:
                            #     viral_comp_family["phage_comp_"+str(my_count)+"_cycle_"+str(cycle_number)] = list(cycle_family)[0]

                        cycle_number += 1


        else:

            contig_name = contig_names[pruned_vs[my_count][0]]

            resolved_edges.add(contig_name)

            # if contig_name in edge_family:
            #     cycle_family.add(edge_family[contig_name])
            # else:
            #     cycle_family.add("Unclassified")

            path_string = graph_contigs[contig_name]

            cycle_number = 1

            genome_path = GenomePath(
                f"phage_comp_{my_count}_cycle_{cycle_number}",
                [pruned_vs[my_count][0]],
                path_string,
                coverage,
                len(graph_contigs[contig_name]),
            )
            my_genomic_paths.append(genome_path)

            # viral_comp_length["viral_comp_"+str(my_count)+"_path_"+str(cycle_number)] = len(path_string)

            # has_cycles.append(my_count)

            # contig_name = contig_names[pruned_vs[my_count][0]]

            # if contig_name in edge_family:
            #     viral_comp_family["viral_comp_"+str(my_count)+"_path_"+str(cycle_number)] = edge_family[contig_name]


        my_genomic_paths.sort(key=lambda x: x.length, reverse=True)

        final_genomic_paths = []

        if len(my_genomic_paths) > 1:

            # Get component stats
            graph_degree = assembly_graph.degree(pruned_vs[my_count], mode="all")
            in_degree = assembly_graph.degree(pruned_vs[my_count], mode="in")
            out_degree = assembly_graph.degree(pruned_vs[my_count], mode="out")

            path_lengths = []
            path_coverages = []

            len_dif_threshold = pathdiff

            prev_length = my_genomic_paths[0].length

            for genomic_path in my_genomic_paths:
                current_len_dif = abs(prev_length - genomic_path.length)
                # cycle_nodes = set(genomic_path.node_order)

                if current_len_dif < len_dif_threshold:

                    path_is_subset = False

                    for final_path in final_genomic_paths:
                        if (set(genomic_path.node_order).issubset(set(final_path.node_order))):
                            path_is_subset = True
                            break

                    if not path_is_subset:
                        prev_length = genomic_path.length

                        logger.debug(f"{genomic_path.id}\t{genomic_path.length}")
                        path_lengths.append(genomic_path.length)
                        path_coverages.append(genomic_path.coverage)
                        final_genomic_paths.append(genomic_path)
                        all_resolved_paths.append(genomic_path)
                else:
                    break

            genome_comp = GenomeComponent(
                f"phage_comp_{my_count}",
                len(pruned_vs[my_count]),
                len(final_genomic_paths),
                max(graph_degree),
                max(in_degree),
                max(out_degree),
                sum(graph_degree) / len(graph_degree),
                sum(in_degree) / len(in_degree),
                sum(out_degree) / len(out_degree),
                pruned_graph.density(loops=False),
                max(path_lengths),
                min(path_lengths),
                max(path_lengths)/min(path_lengths),
                path_lengths[path_coverages.index(max(path_coverages))],
                path_lengths[path_coverages.index(min(path_coverages))],
                path_lengths[path_coverages.index(max(path_coverages))]/path_lengths[path_coverages.index(min(path_coverages))],
                max(path_coverages),
                min(path_coverages),
                max(path_coverages)/min(path_coverages)
            )
            all_components.append(genome_comp)

        else:
            for genomic_path in my_genomic_paths:
                final_genomic_paths.append(genomic_path)
                all_resolved_paths.append(genomic_path)
                logger.debug(f"{genomic_path.id}\t{genomic_path.length}")

        # Write to file

        with open(f"{output}/resolved_paths.fasta", "a+") as myfile:

            for genomic_path in final_genomic_paths:

                myfile.write(f">{genomic_path.id}\n")

                chunks = [
                    genomic_path.path[i : i + FASTA_LINE_LEN]
                    for i in range(0, genomic_path.length, FASTA_LINE_LEN)
                ]

                for chunk in chunks:
                    myfile.write(f"{chunk}\n")

        output_genomes_path = f"{output}/resolved_phages"
        subprocess.run("mkdir -p " + output_genomes_path, shell=True)

        for genomic_path in final_genomic_paths:

            with open(
                f"{output}/resolved_phages/{genomic_path.id}.fasta", "w+"
            ) as myfile:

                myfile.write(">" + genomic_path.id + "\n")

                chunks = [
                    genomic_path.path[i : i + FASTA_LINE_LEN]
                    for i in range(0, genomic_path.length, FASTA_LINE_LEN)
                ]

                for chunk in chunks:
                    myfile.write(chunk + "\n")

    logger.info(f"Total number of genomes resolved: {len(all_resolved_paths)}")
    logger.info(f"Resolved genomes can be found in {output}/resolved_paths.fasta")

    # Record path information
    # ----------------------------------------------------------------------

    with open(f"{output}/resolved_genome_info.txt", "w") as myfile:
        myfile.write(f"Path\tCoverage\tLength\tNode order\n")
        for genomic_path in all_resolved_paths:
            myfile.write(
                f"{genomic_path.id}\t{genomic_path.coverage}\t{genomic_path.length}\t{genomic_path.node_order}\n"
            )

    logger.info(
        f"Resolved genome information can be found in {output}/resolved_genome_info.txt"
    )

    # Record component information
    # ----------------------------------------------------------------------

    with open(f"{output}/resolved_component_info.txt", "w") as myfile:
        myfile.write(f"Component\tNumber of nodes\tNumber of paths\tMaximum degree\t")
        myfile.write(f"Maximum in degree\tMaximum out degree\tAverage degree\t")
        myfile.write(f"Average in degree\tAverage out degree\tDensity\t")
        myfile.write(f"Maximum path length\tMinimum path length\tLength ratio (long/short)\t")
        myfile.write(f"Maximum coverage path length\tMinimum coverage path length\tLength ratio (highest cov/lowest cov)\t")
        myfile.write(f"Maximum coverage\tMinimum coverage\tCoverage ratio (highest/lowest)\n")

        for component in all_components:
            myfile.write(f"{component.id}\t{component.n_nodes}\t{component.n_paths}\t{component.max_degree}\t")
            myfile.write(f"{component.max_in_degree}\t{component.max_out_degree}\t{component.avg_degree}\t")
            myfile.write(f"{component.avg_in_degree}\t{component.avg_out_degree}\t{component.density}\t")
            myfile.write(f"{component.max_path_length}\t{component.min_path_length}\t{component.min_max_len_ratio}\t")
            myfile.write(f"{component.max_cov_path_length}\t{component.min_cov_path_length}\t{component.min_max_cov_len_ratio}\t")
            myfile.write(f"{component.max_cov}\t{component.min_cov}\t{component.min_max_cov_ratio}\n")

    logger.info(
        f"Resolved component information can be found in {output}/resolved_component_info.txt"
    )

    # Get elapsed time
    # ----------------------------------------------------------------------

    # Determine elapsed time
    elapsed_time = time.time() - start_time

    # Print elapsed time for the process
    logger.info("Elapsed time: " + str(elapsed_time) + " seconds")

    # Exit program
    # ----------------------------------------------------------------------

    logger.info("Thank you for using ds-bubbles!")


if __name__ == "__main__":
    main()
