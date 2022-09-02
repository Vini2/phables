#!/usr/bin/env python3

import logging
import subprocess
import sys
import time

import click
import networkx as nx
import numpy as np
from igraph import *
from tqdm import tqdm

from phables_utils import (component_utils, edge_graph_utils, edge_utils,
                             gene_utils)
from phables_utils.genome_utils import GenomeComponent, GenomePath
from phables_utils.CycFlowDec import CycFlowDec

__author__ = "Vijini Mallawaarachchi"
__copyright__ = "Copyright 2022, PhaBles Project"
__license__ = "MIT"
__version__ = "0.1"
__maintainer__ = "Vijini Mallawaarachchi"
__email__ = "vijini.mallawaarachchi@flinders.edu.au"
__status__ = "Development"

MAX_VAL = sys.maxsize
FASTA_LINE_LEN = 60
REPEAT_MIN_LENGTH = 1000

# Sample command
# -------------------------------------------------------------------
# python phables.py  -g /path/to/assembly_graph.gfa
#                    -c /path/to/assembly.fasta
#                    -p /path/to/assembly_info.txt
#                    -o /path/to/output_folder
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
    "--coverage",
    "-cov",
    required=True,
    help="path to the coverage file",
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
    "--mincov",
    "-mcov",
    default=1,
    required=False,
    help="minimum coverage of paths to output",
    type=int,
)
@click.option(
    "--compcount",
    "-cc",
    default=500,
    required=False,
    help="maximum contig count to consider a component",
    type=int,
)
@click.option(
    "--mgfrac",
    "-mgf",
    default=0.2,
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
    default=0.3,
    required=False,
    help="minimum sequence identity for phrog annotations",
    type=float,
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
    coverage,
    minlength,
    mincov,
    compcount,
    mgfrac,
    alignscore,
    seqidentity,
    output,
):

    """PhaBles: Resolve bacteriophage genomes from phage bubbles in metagenomic data."""

    # Setup logger
    # ----------------------------------------------------------------------

    logger = logging.getLogger("phables 0.1")
    logger.setLevel(logging.DEBUG)
    logging.captureWarnings(True)
    formatter = logging.Formatter("%(asctime)s - %(levelname)s - %(message)s")
    consoleHeader = logging.StreamHandler()
    consoleHeader.setFormatter(formatter)
    consoleHeader.setLevel(logging.INFO)
    logger.addHandler(consoleHeader)

    # Setup output path for log file
    fileHandler = logging.FileHandler(f"{output}/phables.log")
    fileHandler.setLevel(logging.DEBUG)
    fileHandler.setFormatter(formatter)
    logger.addHandler(fileHandler)

    logger.info(
        "Welcome to PhaBles: Resolve bacteriophage genomes from phage bubbles in metagenomic data."
    )

    logger.info(f"Input arguments: ")
    logger.info(f"Assembly graph file: {graph}")
    logger.info(f"Contigs file: {contigs}")
    logger.info(f"Contig paths file: {paths}")
    logger.info(f"Contig .hmmout file: {hmmout}")
    logger.info(f"Contig phrog annotations file: {phrogs}")
    logger.info(f"Contig coverage file: {coverage}")
    logger.info(f"Minimum length of contigs to consider: {minlength}")
    logger.info(f"Minimum coverage of paths to output: {mincov}")
    logger.info(f"minimum contig count to consider a component: {compcount}")
    logger.info(f"Length threshold to consider single copy marker genes: {mgfrac}")
    logger.info(f"Minimum alignment score (%) for phrog annotations: {alignscore}")
    logger.info(f"Minimum sequence identity for phrog annotations: {seqidentity}")
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

    # Get contig coverages
    # ----------------------------------------------------------------------

    contig_coverages = {}

    with open(coverage, "r") as myfile:
        for line in myfile.readlines():

            if not line.startswith("Contig"):
                strings = line.strip().split()

                contig_name = strings[0].replace("contig", "edge")

                coverage_sum = sum([float(x) for x in strings[1:]])

                contig_coverages[contig_name] = coverage_sum

    # Resolve genomes
    # ----------------------------------------------------------------------

    resolved_edges = set()

    all_resolved_paths = []

    all_components = []

    cycle_components = set()
    resolved_components = set()
    circular_contigs = set()

    for my_count in tqdm(pruned_vs, desc="Resolving components"):

        my_genomic_paths = []

        candidate_nodes = pruned_vs[my_count]

        pruned_graph = assembly_graph.subgraph(candidate_nodes)

        comp_self_looped_nodes = [x for x in pruned_vs[my_count] if x in self_looped_nodes]

        # if 1352 not in candidate_nodes:
        #     continue

        logger.debug(f"my_count: {my_count}")
        

        logger.debug(f"number of contigs: {len(candidate_nodes)}")
        logger.debug(f"{candidate_nodes}")

        
        if len(candidate_nodes) == 2:

            all_self_looped = True

            for node in candidate_nodes:
                if node not in self_looped_nodes:
                    all_self_looped = False
            
            if all_self_looped:

                cycle_components.add(my_count)

                contig1 = ""
                contig2 = ""

                for edge in pruned_graph.es:
                    source_vertex_id = edge.source
                    target_vertex_id = edge.target

                    if source_vertex_id != target_vertex_id:
                        contig1 = candidate_nodes[source_vertex_id]
                        contig2 = candidate_nodes[target_vertex_id]

                if contig1 != "" and contig2 != "":

                    contig1_name = contig_names[contig1]
                    contig2_name = contig_names[contig2]

                    resolved_edges.add(contig1_name)
                    resolved_edges.add(contig2_name)

                    path_string = graph_contigs[contig1_name] + graph_contigs[contig2_name]

                    cycle_number = 1

                    genome_path = GenomePath(
                        f"phage_comp_{my_count}_cycle_{cycle_number}",
                        [contig1_name, contig2_name],
                        [contig1, contig2],
                        path_string,
                        min(contig_coverages[contig1_name], contig_coverages[contig2_name]),
                        len(path_string),
                    )
                    my_genomic_paths.append(genome_path)
                    resolved_components.add(my_count)


        elif len(candidate_nodes) > 2 and len(candidate_nodes) <= compcount:

            is_complex = False
            cycles_found = False

            all_degrees = []
            all_clean_degrees = []

            for node in candidate_nodes:

                # Get degree without self-looping nodes
                in_neighbours = assembly_graph.neighbors(node, mode="in")
                out_neighbours = assembly_graph.neighbors(node, mode="out")
                in_degree_node = len([x for x in in_neighbours])
                out_degree_node = len([x for x in out_neighbours])
                clean_in_degree_node = len([x for x in in_neighbours if x not in self_looped_nodes])
                clean_out_degree_node = len([x for x in out_neighbours if x not in self_looped_nodes])

                all_clean_degrees.append(clean_in_degree_node)
                all_clean_degrees.append(clean_out_degree_node)

                all_degrees.append(in_degree_node)
                all_degrees.append(out_degree_node)

            clean_max_degree = max(all_clean_degrees)
            component_max_degree = max(all_degrees)

            logger.debug(f"clean_max_degree: {clean_max_degree}")
            logger.debug(f"component_max_degree: {component_max_degree}")

            # if component_max_degree >= degree:
            #     is_complex = True
            #     logger.debug(f"Complex component")

            # Create Directed Graph
            G=nx.DiGraph()

            my_edges = []

            for edge in pruned_graph.es:
                source_vertex_id = edge.source
                target_vertex_id = edge.target
                my_edges.append((source_vertex_id, target_vertex_id))

            # Add a list of edges:
            G.add_edges_from(my_edges)

            try:
                nx.find_cycle(G, orientation="original")
                cycles_found = True
                cycle_components.add(my_count)

            except nx.exception.NetworkXNoCycle as e:
                logger.debug(f"{e}")

            if cycles_found:

                # Create Directed Graph
                G=nx.DiGraph()

                for i in range(len(candidate_nodes)):
                    for j in range(len(candidate_nodes)):
                        
                        if (i,j) in my_edges:

                            consider_edge = False

                            if candidate_nodes[i] not in self_looped_nodes and candidate_nodes[j] not in self_looped_nodes:
                                consider_edge = True

                            elif len(comp_self_looped_nodes) == 1:
                                consider_edge = True

                            elif candidate_nodes[i] in self_looped_nodes and candidate_nodes[j] not in self_looped_nodes:
                                if len(graph_contigs[contig_names[candidate_nodes[i]]]) > REPEAT_MIN_LENGTH:
                                    consider_edge = True

                            elif candidate_nodes[i] not in self_looped_nodes and candidate_nodes[j] in self_looped_nodes:
                                if len(graph_contigs[contig_names[candidate_nodes[j]]]) > REPEAT_MIN_LENGTH:
                                    consider_edge = True

                            if consider_edge:
                                # min_cov = min(edge_depths[contig_names[candidate_nodes[i]]], edge_depths[contig_names[candidate_nodes[j]]])

                                cov_1 = MAX_VAL
                                cov_2 = MAX_VAL

                                if contig_names[candidate_nodes[i]] in contig_coverages:
                                    cov_1 = contig_coverages[contig_names[candidate_nodes[i]]]
                                if contig_names[candidate_nodes[j]] in contig_coverages:
                                    cov_2 = contig_coverages[contig_names[candidate_nodes[j]]]

                                min_cov = min([cov_1, cov_2])
                                G.add_edge(i, j, weight=min_cov)
                                # G.add_edge(j, i, weight=min_cov)
                        elif (j,i) in my_edges:
                        #     G.add_edge(j, i, weight=min_cov)
                            continue
                        else:
                            G.add_edge(i, j, weight=0)

                # for edge in G.edges.data("weight", default=1):
                #     print(edge)

                cycle_number = 1

                cycles_iterated = False

                try:
                    F = np.zeros([G.number_of_nodes(),G.number_of_nodes()])

                    for edge in G.edges.data("weight", default=1):
                    #     print(edge)
                        F[edge[0],edge[1]] = edge[2]

                    n_points = 10
                    init_steps = 20
                    step_add = 10
                    run_post = 2
                    state = 0
                    tol = 1E-5

                    myCycFlowDec = CycFlowDec(F,state,tol)
                    myCycFlowDec.run(init_steps-run_post,run_post)
                    steps = [init_steps]
                    MREs = [myCycFlowDec.calc_MRE(0)]

                    j = 1
                    while j < n_points:
                        myCycFlowDec.run(step_add-run_post,run_post)
                        steps.append(steps[-1] + step_add)
                        MREs.append(myCycFlowDec.calc_MRE(0))
                        j += 1

                    my_cc = 0

                    logger.debug(f"Number of cycles found: {len(myCycFlowDec.cycles.keys())}")

                    for cycle in myCycFlowDec.cycles.keys():

                        cycle_len = sum([len(graph_contigs[contig_names[candidate_nodes[node]]]) for node in cycle])

                        if myCycFlowDec.cycles[cycle] >= mincov and cycle_len >= minlength:

                            my_cc += 1

                            path_string = ""
                            total_length = 0
                            coverage_val = myCycFlowDec.cycles[cycle]

                            for node in cycle:
                                contig_name = contig_names[candidate_nodes[node]]
                                path_string += graph_contigs[contig_name]
                                total_length += len(graph_contigs[contig_name])

                            genome_path = GenomePath(
                                f"phage_comp_{my_count}_cycle_{cycle_number}",
                                [contig_names[candidate_nodes[x]] for x in cycle],
                                [candidate_nodes[x] for x in cycle],
                                path_string,
                                coverage_val,
                                total_length,
                            )
                            my_genomic_paths.append(genome_path)

                            cycle_number += 1

                    cycles_iterated = True
                
                except ZeroDivisionError as e:
                    logger.debug(f"Encountered {e} Could not find paths")
                    # is_complex = True

                if cycles_iterated and len(myCycFlowDec.cycles.keys()) == 0:
                    logger.debug(f"No paths were found. Switching to detect cycles")
                    is_complex = True

            # if is_complex:

            #     # Create Directed Graph
            #     G = nx.DiGraph()

            #     my_edges = []

            #     for edge in pruned_graph.es:
            #         source_vertex_id = edge.source
            #         target_vertex_id = edge.target
            #         my_edges.append((source_vertex_id, target_vertex_id))

            #     # Add a list of edges:
            #     G.add_edges_from(my_edges)

            #     logger.debug("Complex_graph")
            #     my_cycles = nx.cycle_basis(G.to_undirected())

            #     path_lengths = []
            #     path_coverages = []

            #     # Return a list of cycles described as a list of nodes
            #     for cycle in my_cycles:

            #         if len(cycle) > 1:

            #             total_length = 0

            #             coverage = MAX_VAL

            #             node_order = []
            #             node_id_order = []

            #             path_string = ""

            #             for node in cycle:

            #                 contig_name = contig_names[candidate_nodes[node]]

            #                 node_order.append(contig_names[candidate_nodes[node]])
            #                 node_id_order.append(candidate_nodes[node])
            #                 total_length += len(graph_contigs[contig_name])
            #                 path_string += graph_contigs[contig_name]

            #                 if contig_name in contig_coverages:
            #                     if coverage > contig_coverages[contig_name]:
            #                         coverage = contig_coverages[contig_name]

            #             path_lengths.append(total_length)
            #             path_coverages.append(coverage)

            #             # Save path

            #             if len(path_string) > 0 and coverage != MAX_VAL and coverage > 0:

            #                 genome_path = GenomePath(
            #                     f"phage_comp_{my_count}_cycle_{cycle_number}",
            #                     node_order,
            #                     node_id_order,
            #                     path_string,
            #                     coverage,
            #                     total_length,
            #                 )
            #                 my_genomic_paths.append(genome_path)

            #             cycle_number += 1


        elif len(candidate_nodes) == 1:

            contig_name = contig_names[candidate_nodes[0]]

            resolved_edges.add(contig_name)

            path_string = graph_contigs[contig_name]

            cycle_number = 1

            genome_path = GenomePath(
                f"phage_comp_{my_count}_cycle_{cycle_number}",
                [contig_names[candidate_nodes[0]]],
                [candidate_nodes[0]],
                path_string,
                contig_coverages[contig_name],
                len(graph_contigs[contig_name]),
            )
            my_genomic_paths.append(genome_path)
            resolved_components.add(my_count)

            circular_contigs.add(my_count)

        # single_in_out_nodes_visited = set()
        visited_count = {}

        my_genomic_paths.sort(key=lambda x: x.length, reverse=True)

        final_genomic_paths = []

        if len(my_genomic_paths) > 1:

            # Get component stats
            graph_degree = assembly_graph.degree(candidate_nodes, mode="all")
            in_degree = assembly_graph.degree(candidate_nodes, mode="in")
            out_degree = assembly_graph.degree(candidate_nodes, mode="out")

            path_lengths = []
            path_coverages = []

            prev_length = my_genomic_paths[0].length

            for genomic_path in my_genomic_paths:
                current_len_dif = abs(prev_length - genomic_path.length)

                # if current_len_dif < len_dif_threshold:

                path_is_subset = False

                for final_path in final_genomic_paths:
                    if set(genomic_path.node_order).issubset(
                        set(final_path.node_order)
                    ):
                        path_is_subset = True
                        break

                for path_node in genomic_path.node_id_order:
                    if path_node in visited_count:
                        if (
                            visited_count[path_node]
                            >= clean_max_degree
                        ):
                            path_is_subset = True

                if not path_is_subset:
                    prev_length = genomic_path.length

                    logger.debug(f"{genomic_path.id}\t{genomic_path.length}")
                    path_lengths.append(genomic_path.length)
                    path_coverages.append(genomic_path.coverage)
                    final_genomic_paths.append(genomic_path)
                    all_resolved_paths.append(genomic_path)

                    for path_node in genomic_path.node_id_order:
                        resolved_edges.add(contig_names[path_node])
                        # if path_node in single_in_out_nodes:
                        #     single_in_out_nodes_visited.add(path_node)

                        if path_node in visited_count:
                            visited_count[path_node] += 1
                        else:
                            visited_count[path_node] = 1

                else:
                    break

            if len(final_genomic_paths) > 1:

                genome_comp = GenomeComponent(
                    f"phage_comp_{my_count}",
                    len(candidate_nodes),
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
                    max(path_lengths) / min(path_lengths),
                    path_lengths[path_coverages.index(max(path_coverages))],
                    path_lengths[path_coverages.index(min(path_coverages))],
                    path_lengths[path_coverages.index(max(path_coverages))]
                    / path_lengths[path_coverages.index(min(path_coverages))],
                    max(path_coverages),
                    min(path_coverages),
                    max(path_coverages) / min(path_coverages),
                )
                all_components.append(genome_comp)

            if len(final_genomic_paths) > 0:
                resolved_components.add(my_count)

        else:
            for genomic_path in my_genomic_paths:
                final_genomic_paths.append(genomic_path)
                all_resolved_paths.append(genomic_path)
                logger.debug(f"{genomic_path.id}\t{genomic_path.length}")
                resolved_components.add(my_count)

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

                myfile.write(f">{genomic_path.id}\n")

                chunks = [
                    genomic_path.path[i : i + FASTA_LINE_LEN]
                    for i in range(0, genomic_path.length, FASTA_LINE_LEN)
                ]

                for chunk in chunks:
                    myfile.write(f"{chunk}\n")

    resolved_cyclic = len(cycle_components) - (len(resolved_components) - len(circular_contigs))

    logger.info(f"Total number of cyclic components found: {len(cycle_components)}")
    logger.info(f"Total number of cyclic components resolved: {resolved_cyclic}")
    logger.info(f"Circular contigs identified: {len(circular_contigs)}")
    logger.info(f"Total number of cyclic components found: {len(cycle_components) + len(circular_contigs)}")
    logger.info(f"Total number of components resolved: {len(resolved_components)}")
    logger.info(f"Total number of genomes resolved: {len(all_resolved_paths)}")
    logger.info(f"Resolved genomes can be found in {output}/resolved_paths.fasta")

    with open(f"{output}/resolved_edges.fasta", "w+") as myfile:
        for edge in resolved_edges:

            edge_seq = graph_contigs[edge]

            myfile.write(">" + edge + "\n")

            chunks = [
                edge_seq[i : i + FASTA_LINE_LEN]
                for i in range(0, len(edge_seq), FASTA_LINE_LEN)
            ]

            for chunk in chunks:
                myfile.write(chunk + "\n")

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
        myfile.write(f"Component\t")
        myfile.write(f"Number of nodes\t")
        myfile.write(f"Number of paths\t")
        myfile.write(f"Maximum degree\t")
        myfile.write(f"Maximum in degree\t")
        myfile.write(f"Maximum out degree\t")
        myfile.write(f"Average degree\t")
        myfile.write(f"Average in degree\t")
        myfile.write(f"Average out degree\t")
        myfile.write(f"Density\t")
        myfile.write(f"Maximum path length\t")
        myfile.write(f"Minimum path length\t")
        myfile.write(f"Length ratio (long/short)\t")
        myfile.write(f"Maximum coverage path length\t")
        myfile.write(f"Minimum coverage path length\t")
        myfile.write(f"Length ratio (highest cov/lowest cov)\t")
        myfile.write(f"Maximum coverage\t")
        myfile.write(f"Minimum coverage\t")
        myfile.write(f"Coverage ratio (highest/lowest)\n")

        for component in all_components:
            myfile.write(f"{component.id}\t")
            myfile.write(f"{component.n_nodes}\t")
            myfile.write(f"{component.n_paths}\t")
            myfile.write(f"{component.max_degree}\t")
            myfile.write(f"{component.max_in_degree}\t")
            myfile.write(f"{component.max_out_degree}\t")
            myfile.write(f"{component.avg_degree}\t")
            myfile.write(f"{component.avg_in_degree}\t")
            myfile.write(f"{component.avg_out_degree}\t")
            myfile.write(f"{component.density}\t")
            myfile.write(f"{component.max_path_length}\t")
            myfile.write(f"{component.min_path_length}\t")
            myfile.write(f"{component.min_max_len_ratio}\t")
            myfile.write(f"{component.max_cov_path_length}\t")
            myfile.write(f"{component.min_cov_path_length}\t")
            myfile.write(f"{component.min_max_cov_len_ratio}\t")
            myfile.write(f"{component.max_cov}\t")
            myfile.write(f"{component.min_cov}\t")
            myfile.write(f"{component.min_max_cov_ratio}\n")

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

    logger.info("Thank you for using phables!")


if __name__ == "__main__":
    main()
