#!/usr/bin/env python3

import math
import multiprocessing
import pickle
from concurrent.futures import ProcessPoolExecutor
import logging
import sys
import time

import networkx as nx
from phables_utils import edge_graph_utils, flow_utils
from phables_utils.genome_utils import GenomeComponent, GenomePath
from phables_utils.output_utils import write_path, write_path_fasta
from tqdm import tqdm

MAX_VAL = sys.maxsize
LEN_THRESHOLD = 0.95

# Create logger
logger = logging.getLogger("phables 2.0.0")

# Upper bound on the number of tasks handed to the pool in
# resolve_short_parallel. In practice this means ONE COMPONENT PER TASK for any
# realistic assembly, which is deliberate -- see the chunking comment there.
# The cap only exists so a pathological graph with hundreds of thousands of
# components cannot turn into hundreds of thousands of pool submissions.
MAX_MFD_TASKS = 50000


def resolve_short(
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
):
    """
    Resolve genomic paths using short reads
    """

    resolved_edges = set()

    all_resolved_paths = []

    all_components = []

    cycle_components = set()
    linear_components = set()
    resolved_components = set()
    resolved_linear = set()
    single_unitigs = set()
    resolved_cyclic = set()

    case1_found = set()
    case1_resolved = set()
    case2_found = set()
    case2_resolved = set()
    case3_found = set()
    case3_resolved = set()

    phage_like_edges = set()
    all_phage_like_edges = set()

    unresolved_phage_like_edges = set()

    for my_count in tqdm(pruned_vs, desc="Resolving components"):
        component_time_start = time.time()

        my_genomic_paths = []

        original_candidate_nodes = pruned_vs[my_count]

        candidate_nodes = pruned_vs[my_count]

        pruned_graph = assembly_graph.subgraph(candidate_nodes)

        has_cycles = False

        logger.debug(f"my_count: {my_count}")

        logger.debug(f"number of unitigs: {len(candidate_nodes)}")
        logger.debug(f"{candidate_nodes}")

        all_phage_like_edges = all_phage_like_edges.union(set(candidate_nodes))

        comp_all_edges = set(set(candidate_nodes))
        comp_resolved_edges = set()

        in_degree = []
        out_degree = []

        case_name = ""

        # Case 2 components
        if len(candidate_nodes) == 2:
            all_self_looped = False
            one_circular = False

            if (
                unitig_names[candidate_nodes[0]] in self_looped_nodes
                and unitig_names[candidate_nodes[1]] in self_looped_nodes
            ):
                all_self_looped = True
            else:
                if unitig_names[candidate_nodes[0]] in self_looped_nodes:
                    one_circular = True
                    all_self_looped = False
                if unitig_names[candidate_nodes[1]] in self_looped_nodes:
                    one_circular = True
                    all_self_looped = False

            unitig1 = ""
            unitig2 = ""

            for edge in pruned_graph.es:
                source_vertex_id = edge.source
                target_vertex_id = edge.target

                if source_vertex_id != target_vertex_id:
                    unitig1 = candidate_nodes[source_vertex_id]
                    unitig2 = candidate_nodes[target_vertex_id]

            unitig1_name = unitig_names[unitig1]
            unitig2_name = unitig_names[unitig2]

            unitig1_len = len(str(graph_unitigs[unitig1_name]))
            unitig2_len = len(str(graph_unitigs[unitig2_name]))

            if unitig1 != "" and unitig2 != "":
                if all_self_looped:
                    # Case 2 - both are circular
                    case_name = "case2_circular"

                    case2_found.add(my_count)

                    cycle_components.add(my_count)

                    phage_like_edges = phage_like_edges.union(set(candidate_nodes))
                    comp_resolved_edges = comp_resolved_edges.union(
                        set(candidate_nodes)
                    )

                    unitig_to_consider = -1
                    unitig_name = ""

                    repeat_unitig = -1
                    repeat_unitig_name = ""

                    if unitig1_len > unitig2_len and unitig1_len > minlength:
                        unitig_to_consider = unitig1
                        unitig_name = unitig1_name
                        repeat_unitig = unitig2
                        repeat_unitig_name = unitig2_name
                    elif unitig2_len > unitig1_len and unitig2_len > minlength:
                        unitig_to_consider = unitig2
                        unitig_name = unitig2_name
                        repeat_unitig = unitig1
                        repeat_unitig_name = unitig1_name

                    if unitig_to_consider != -1:
                        logger.debug(
                            f"Case 2 component: {unitig1_name} is {unitig1_len} bp long and {unitig2_name} is {unitig2_len} bp long."
                        )
                        cycle_number = 1
                        resolved_edges.add(unitig_to_consider)
                        resolved_edges.add(repeat_unitig)
                        path_string = (
                            str(graph_unitigs[repeat_unitig_name])
                            + str(
                                graph_unitigs[unitig_name][
                                    link_overlap[(repeat_unitig, unitig_to_consider)] :
                                ]
                            )
                            + str(
                                graph_unitigs[repeat_unitig_name][
                                    link_overlap[(unitig_to_consider, repeat_unitig)] :
                                ]
                            )
                        )
                        logger.debug(
                            f"Terminal repeat detected is {repeat_unitig_name}"
                        )

                        genome_path = GenomePath(
                            id=f"{prefix}phage_comp_{my_count}_cycle_{cycle_number}",
                            bubble_case=case_name,
                            node_order=[
                                f"{repeat_unitig_name}+",
                                f"{unitig_name}+",
                                f"{repeat_unitig_name}-",
                            ],
                            node_order_human=f"{repeat_unitig_name}:fwd,{unitig_name}:fwd,{repeat_unitig_name}:rev",
                            node_id_order=[
                                repeat_unitig,
                                unitig_to_consider,
                                repeat_unitig,
                            ],
                            path=path_string,
                            coverage=int(unitig_coverages[unitig_name]),
                            length=len(path_string),
                            gc=(path_string.count("G") + path_string.count("C"))
                            / len(path_string)
                            * 100,
                        )
                        my_genomic_paths.append(genome_path)
                        resolved_components.add(my_count)
                        resolved_cyclic.add(my_count)
                        case2_resolved.add(my_count)

                # Case 2 - only one is circular
                elif one_circular:
                    case_name = "case2_linear"

                    case2_found.add(my_count)

                    cycle_components.add(my_count)

                    phage_like_edges = phage_like_edges.union(set(candidate_nodes))
                    comp_resolved_edges = comp_resolved_edges.union(
                        set(candidate_nodes)
                    )

                    unitig_to_consider = -1
                    unitig_name = ""

                    repeat_unitig = -1
                    repeat_unitig_name = ""

                    if (
                        unitig1_len > unitig2_len
                        and unitig1_len > minlength
                        and unitig2_name in self_looped_nodes
                    ):
                        unitig_to_consider = unitig1
                        unitig_name = unitig1_name
                        repeat_unitig = unitig2
                        repeat_unitig_name = unitig2_name
                    elif (
                        unitig2_len > unitig1_len
                        and unitig2_len > minlength
                        and unitig1_name in self_looped_nodes
                    ):
                        unitig_to_consider = unitig2
                        unitig_name = unitig2_name
                        repeat_unitig = unitig1
                        repeat_unitig_name = unitig1_name

                    if unitig_to_consider != -1:
                        logger.debug(
                            f"Case 2 component: {unitig1_name} is {unitig1_len} bp long and {unitig2_name} is {unitig2_len} bp long."
                        )
                        cycle_number = 1
                        resolved_edges.add(unitig_to_consider)
                        resolved_edges.add(repeat_unitig)

                        # Get repeat count
                        repeat_count = max(
                            int(
                                unitig_coverages[repeat_unitig_name]
                                / unitig_coverages[unitig_name]
                            ),
                            1,
                        )
                        logger.debug(f"Repeat count: {repeat_count}")

                        path_string = (
                            str(
                                graph_unitigs[unitig_name][
                                    link_overlap[(repeat_unitig, unitig_to_consider)] :
                                ]
                            )
                            + str(
                                graph_unitigs[repeat_unitig_name][
                                    link_overlap[(unitig_to_consider, repeat_unitig)] :
                                ]
                            )
                            * repeat_count
                        )
                        logger.debug(
                            f"Terminal repeat detected is {repeat_unitig_name}"
                        )

                        # Format path node order
                        path_with_repeats = [f"{unitig_name}+"] + [
                            f"{repeat_unitig_name}+" for x in range(repeat_count)
                        ]

                        repeat_order = f"{repeat_unitig_name}:fwd," * repeat_count
                        path_with_repeats_human = (
                            f"{unitig_name}:fwd,{repeat_order[:-1]}"
                        )
                        node_id_order_with_repeats = [unitig_to_consider] + [
                            repeat_unitig for x in range(repeat_count)
                        ]

                        genome_path = GenomePath(
                            id=f"{prefix}phage_comp_{my_count}_cycle_{cycle_number}",
                            bubble_case=case_name,
                            node_order=path_with_repeats,
                            node_order_human=path_with_repeats_human,
                            node_id_order=node_id_order_with_repeats,
                            path=path_string,
                            coverage=int(unitig_coverages[unitig_name]),
                            length=len(path_string),
                            gc=(path_string.count("G") + path_string.count("C"))
                            / len(path_string)
                            * 100,
                        )
                        my_genomic_paths.append(genome_path)
                        resolved_components.add(my_count)
                        resolved_cyclic.add(my_count)
                        case2_resolved.add(my_count)

        # Case 3 components
        elif len(candidate_nodes) > 2 and len(candidate_nodes) <= compcount:
            case_name = "case3_circular"

            # Create initial directed graph with coverage values
            # ----------------------------------------------------------------------
            G_edge = nx.DiGraph()

            my_counter = 0

            node_indices = {}
            node_indices_rev = {}

            cycle_edges = {}

            clean_node_count = 0

            max_comp_cov = -1

            for vertex in pruned_graph.vs["id"]:
                unitig_name = unitig_names[vertex]

                # Find the maximum coverage within the component
                if (
                    unitig_name in unitig_coverages
                    and unitig_coverages[unitig_name] > max_comp_cov
                ):
                    max_comp_cov = unitig_coverages[unitig_name]

                if unitig_name not in self_looped_nodes:
                    clean_node_count += 1

                for node in oriented_links[unitig_name]:
                    consider_edge = False

                    if not (
                        unitig_name in self_looped_nodes and node in self_looped_nodes
                    ):
                        consider_edge = True

                    if consider_edge:
                        cov_1 = MAX_VAL
                        cov_2 = MAX_VAL

                        if unitig_name in unitig_coverages:
                            cov_1 = unitig_coverages[unitig_name]
                        if node in unitig_coverages:
                            cov_2 = unitig_coverages[node]

                        if min([cov_1, cov_2]) != 0:
                            min_cov = min([cov_1, cov_2])
                        else:
                            min_cov = max([cov_1, cov_2])

                        for edge in oriented_links[unitig_name][node]:
                            cycle_edges[(unitig_name + edge[0], node + edge[1])] = int(
                                min_cov
                            )

            logger.debug(f"clean_node_count: {clean_node_count}")

            for cedge in cycle_edges:
                G_edge.add_edge(cedge[0], cedge[1], weight=cycle_edges[cedge])

            two_comp = sorted(nx.weakly_connected_components(G_edge), key=len)
            logger.debug(f"No. of weakly connected components: {len(two_comp)}")

            if len(two_comp) >= 2:
                G_edge.remove_nodes_from(list(two_comp[0]))

            try:
                cycles_found = nx.find_cycle(G_edge, orientation="original")
                if len(cycles_found) > 0:
                    has_cycles = True
            except nx.exception.NetworkXNoCycle:
                logger.debug(f"No cycles found in component {my_count}")

            if has_cycles:
                logger.debug(
                    f"Potentially cycles can be detected in component {my_count}."
                )

                # Remove dead-ends (nodes with no incoming or no outgoing edges)
                # ----------------------------------------------------------------------
                dead_ends_to_remove = edge_graph_utils.remove_dead_ends(G_edge)

                if len(dead_ends_to_remove) > 0:
                    for node in dead_ends_to_remove:
                        node_id = unitig_names_rev[node[:-1]]
                        if node_id in candidate_nodes:
                            candidate_nodes.remove(node_id)

                    G_edge.remove_nodes_from(dead_ends_to_remove)

                    logger.debug(f"Dead-ends found and removed: {dead_ends_to_remove}")

                # Identify source/sink vertex
                # ----------------------------------------------------------------------

                source_sink_candidates = flow_utils.get_source_sink_circular(
                    G_edge, graph_unitigs, minlength, self_looped_nodes
                )

                source_sink = 0

                logger.debug(f"Original candidate_nodes: {candidate_nodes}")
                logger.debug(
                    f"Identified candidate source_sinks from BFS: {source_sink_candidates}"
                )

                if len(source_sink_candidates) > 0:
                    # Identify the longest source/sink vertex
                    max_length = -1
                    max_length_st_vertex = -1

                    for vertex in source_sink_candidates:
                        if len(graph_unitigs[vertex[:-1]]) > max_length:
                            max_length = len(graph_unitigs[vertex[:-1]])
                            max_length_st_vertex = vertex

                    source_sink = unitig_names_rev[max_length_st_vertex[:-1]]
                    logger.debug(
                        f"Identified source_sink from BFS: {source_sink}, {max_length_st_vertex}"
                    )

                    candidate_nodes.remove(source_sink)
                    candidate_nodes.insert(0, source_sink)
                    logger.debug(f"Ordered candidate_nodes: {candidate_nodes}")

                else:
                    logger.debug(f"No source/sink node detected")
                    continue

                # Create refined directed graph for flow network
                # ----------------------------------------------------------------------
                G = nx.DiGraph()

                for u, v, cov in G_edge.edges(data=True):
                    if u not in node_indices:
                        node_indices[u] = my_counter
                        node_indices_rev[my_counter] = u
                        my_counter += 1
                    if v not in node_indices:
                        node_indices[v] = my_counter
                        node_indices_rev[my_counter] = v
                        my_counter += 1

                    logger.debug(f"Edge: {u}, {v}, {cov['weight']}")

                    G.add_edge(node_indices[u], node_indices[v], weight=cov["weight"])

                # Get connections and degree information
                # ----------------------------------------------------------------------
                in_degree = []
                out_degree = []

                for node in list(G.nodes):
                    if node_indices_rev[node][:-1] not in self_looped_nodes:
                        clean_indegree = len(
                            [
                                x
                                for x in G.predecessors(node)
                                if node_indices_rev[x][:-1] not in self_looped_nodes
                            ]
                        )
                        in_degree.append(clean_indegree)

                        clean_outdegree = len(
                            [
                                x
                                for x in G.successors(node)
                                if node_indices_rev[x][:-1] not in self_looped_nodes
                            ]
                        )
                        out_degree.append(clean_outdegree)

                degrees = in_degree + out_degree

                if len(degrees) == 0:
                    logger.debug(
                        f"Skipping component as no clean connections were found"
                    )
                    continue

                # Create flow network
                # ----------------------------------------------------------------------
                network_edges = list()

                edge_list_indices = {}

                subpaths = {}
                subpath_count = 0

                visited_edges = []

                logger.debug(f"G_edge.nodes: {list(G_edge.nodes)}")
                logger.debug(f"G_edge.edges: {G_edge.edges(data=True)}")

                for u, v, cov in G_edge.edges(data=True):
                    u_name = unitig_names_rev[u[:-1]]
                    v_name = unitig_names_rev[v[:-1]]

                    u_index = candidate_nodes.index(u_name)
                    v_index = candidate_nodes.index(v_name)

                    edge_list_indices[u_index] = u
                    edge_list_indices[v_index] = v

                    juction_cov = junction_pe_coverage[(u[:-1], v[:-1])]

                    if v_index == 0:
                        final_vertex = len(candidate_nodes)
                    else:
                        final_vertex = v_index

                    if (u_index, final_vertex) not in visited_edges and (
                        final_vertex,
                        u_index,
                    ) not in visited_edges:
                        # Get coverage interval
                        cov_lower_bound = cov["weight"]
                        cov_upper_bound = int(max_comp_cov * alpha)

                        logger.debug(
                            f"({v}, {u}), {juction_cov}, {cov_lower_bound}, {cov_upper_bound}"
                        )

                        if juction_cov == 0:
                            network_edges.append(
                                (u_index, final_vertex, 0, cov_upper_bound)
                            )
                        else:
                            network_edges.append(
                                (
                                    u_index,
                                    final_vertex,
                                    cov_lower_bound,
                                    cov_upper_bound,
                                )
                            )

                        visited_edges.append((u_index, final_vertex))

                        # Add subpaths
                        if juction_cov >= mincov:
                            logger.debug(f"Adding subpath {[u_index, final_vertex]}")
                            subpaths[subpath_count] = [u_index, final_vertex]
                            subpath_count += 1

                            # Extend subpaths using coverages of successors and predecessors
                            # -----------------------------------

                            # Extend subpath using coverages of predecessors
                            for u_pred in G_edge.predecessors(u):
                                u_pred_name = unitig_names_rev[u_pred[:-1]]
                                u_pred_index = candidate_nodes.index(u_pred_name)
                                u_pred_cov = unitig_coverages[u_pred[:-1]]
                                u_cov = unitig_coverages[u[:-1]]

                                if (
                                    final_vertex != 0
                                    and u_index != 0
                                    and u_pred_index != final_vertex
                                ):
                                    if (
                                        abs(min(u_pred_cov, u_cov) - cov["weight"])
                                        < covtol
                                    ):
                                        subpaths[subpath_count] = [
                                            u_pred_index,
                                            u_index,
                                            final_vertex,
                                        ]
                                        logger.debug(
                                            f"Extending subpath based on predecessor coverage {[u_pred_index, u_index, final_vertex]}"
                                        )
                                        subpath_count += 1

                            # Extend subpath using coverages of successors
                            for v_succ in G_edge.successors(v):
                                v_succ_name = unitig_names_rev[v_succ[:-1]]
                                v_succ_index = candidate_nodes.index(v_succ_name)
                                v_succ_cov = unitig_coverages[v_succ[:-1]]
                                v_cov = unitig_coverages[v[:-1]]

                                if (
                                    v_succ_index != 0
                                    and u_index != 0
                                    and final_vertex != 0
                                    and final_vertex != len(candidate_nodes)
                                    and v_succ_index != u_index
                                ):
                                    if (
                                        abs(min(v_succ_cov, v_cov) - cov["weight"])
                                        < covtol
                                    ):
                                        subpaths[subpath_count] = [
                                            u_index,
                                            final_vertex,
                                            v_succ_index,
                                        ]
                                        logger.debug(
                                            f"Extending subpath based on successor coverage {[u_index, final_vertex, v_succ_index]}"
                                        )
                                        subpath_count += 1

                        else:
                            # Extend subpaths of l=3 based on paired-end reads
                            # aligned to successors and predecessors
                            # -----------------------------------

                            # Extend subpath using coverages of predecessors
                            for u_pred in G_edge.predecessors(u):
                                if junction_pe_coverage[(u_pred[:-1], v[:-1])] > 0:
                                    u_pred_name = unitig_names_rev[u_pred[:-1]]
                                    u_pred_index = candidate_nodes.index(u_pred_name)
                                    if (
                                        final_vertex != 0
                                        and u_index != 0
                                        and u_pred_index != final_vertex
                                    ):
                                        subpaths[subpath_count] = [
                                            u_pred_index,
                                            u_index,
                                            final_vertex,
                                        ]
                                        logger.debug(
                                            f"Extending subpath {[u_pred_index, u_index, final_vertex]}"
                                        )
                                        subpath_count += 1

                            # Extend subpath using coverages of successors
                            for v_succ in G_edge.successors(v):
                                if junction_pe_coverage[(u[:-1], v_succ[:-1])] > 0:
                                    v_succ_name = unitig_names_rev[v_succ[:-1]]
                                    v_succ_index = candidate_nodes.index(v_succ_name)
                                    if (
                                        v_succ_index != 0
                                        and u_index != 0
                                        and final_vertex != 0
                                        and final_vertex != len(candidate_nodes)
                                        and v_succ_index != u_index
                                    ):
                                        subpaths[subpath_count] = [
                                            u_index,
                                            final_vertex,
                                            v_succ_index,
                                        ]
                                        logger.debug(
                                            f"Extending subpath {[u_index, final_vertex, v_succ_index]}"
                                        )
                                        subpath_count += 1

                logger.debug(f"edge_list_indices: {edge_list_indices}")
                logger.debug(f"subpaths: {subpaths}")

                # Create flow network and run MFD-ILP
                # ----------------------------------------------------------------------
                G_mfd = {
                    "Nodes": len(list(G_edge.nodes)),
                    "list of edges": network_edges,
                    "subpaths": subpaths,
                }
                logger.debug(f"G_mfd: {G_mfd}")
                solution_paths = flow_utils.solve_mfd(G_mfd, maxpaths, output, nthreads)
                logger.debug(f"Number of paths found: {len(solution_paths)}")

                cycle_components.add(my_count)
                case3_found.add(my_count)

                # Iterate through solution paths
                # ----------------------------------------------------------------------
                if len(solution_paths) != 0:
                    phage_like_edges = phage_like_edges.union(
                        set(original_candidate_nodes)
                    )

                    cycle_number = 1

                    for solution_path in solution_paths:
                        coverage_val = solution_paths[solution_path]["weight"]

                        # Filter path by coverage
                        if coverage_val >= mincov:
                            logger.debug(
                                f"Path {cycle_number} coverage: {coverage_val}"
                            )

                            # Create graph for path
                            G_path = nx.DiGraph()

                            # Fill graph with data
                            G_path.add_edges_from(solution_paths[solution_path]["path"])
                            logger.debug(
                                f"solution path: {solution_paths[solution_path]['path']}"
                            )

                            if 0 in list(G_path.nodes):
                                # Get all simple paths from node 0 to last node
                                try:
                                    candidate_paths = list(
                                        nx.all_simple_paths(
                                            G_path, 0, len(candidate_nodes)
                                        )
                                    )

                                    if len(candidate_paths) > 0:
                                        logger.debug(
                                            f"candidate_paths: {candidate_paths[0]}"
                                        )

                                        # Get mapped unitigs in order from the flow network
                                        path_order = []
                                        for path_edge in candidate_paths[0]:
                                            if path_edge != len(candidate_nodes):
                                                path_order.append(
                                                    edge_list_indices[path_edge]
                                                )

                                        logger.debug(f"path_order: {path_order}")

                                        # Get the order of unitigs in path
                                        path_string = ""
                                        total_length = 0

                                        previous_edge = 0

                                        for nodeid in range(len(path_order)):
                                            node = path_order[nodeid]
                                            unitig_name = node[:-1]

                                            if node.endswith("+"):
                                                unitig_seq = str(
                                                    graph_unitigs[unitig_name]
                                                )
                                            else:
                                                unitig_seq = str(
                                                    graph_unitigs[
                                                        unitig_name
                                                    ].reverse_complement()
                                                )

                                            # If first node in path
                                            if previous_edge == 0:
                                                path_string += unitig_seq
                                                total_length += len(unitig_seq)

                                            else:
                                                trimmed_seq = unitig_seq[
                                                    link_overlap[
                                                        (previous_edge, node)
                                                    ] :
                                                ]
                                                path_string += trimmed_seq
                                                total_length += len(trimmed_seq)

                                            previous_edge = node

                                        # Format genomic path
                                        path_node_order_human = ""

                                        for c in path_order:
                                            if c.endswith("+"):
                                                path_node_order_human += (
                                                    f"{c[:-1]}:fwd,"
                                                )
                                            else:
                                                path_node_order_human += (
                                                    f"{c[:-1]}:rev,"
                                                )

                                        path_node_order_human = path_node_order_human[
                                            :-1
                                        ]

                                        # Create GenomePath object with path details
                                        genome_path = GenomePath(
                                            id=f"{prefix}phage_comp_{my_count}_cycle_{cycle_number}",
                                            bubble_case=case_name,
                                            node_order=[x for x in path_order],
                                            node_order_human=path_node_order_human,
                                            node_id_order=[
                                                unitig_names_rev[x[:-1]]
                                                for x in path_order
                                            ],
                                            path=path_string,
                                            coverage=int(coverage_val),
                                            length=total_length,
                                            gc=(
                                                path_string.count("G")
                                                + path_string.count("C")
                                            )
                                            / len(path_string)
                                            * 100,
                                        )
                                        my_genomic_paths.append(genome_path)
                                        logger.debug(f"total_length: {total_length}")

                                        cycle_number += 1

                                except nx.exception.NodeNotFound:
                                    logger.debug(
                                        f"Could not resolve a continuous path."
                                    )

                    logger.debug(f"Number of paths selected: {cycle_number-1}")

                    if cycle_number > 1:
                        resolved_components.add(my_count)
                        resolved_cyclic.add(my_count)
                        case3_resolved.add(my_count)

                else:
                    logger.debug(f"No paths detected")
                    continue

            else:
                logger.debug(f"No cycles detected. Found a complex linear component.")

                case_name = "case3_linear"

                linear_components.add(my_count)

                # Identify source/sink vertex
                # ----------------------------------------------------------------------

                source_candidates, sink_candidates = flow_utils.get_source_sink_linear(
                    G_edge, self_looped_nodes
                )

                logger.debug(f"Original candidate_nodes: {candidate_nodes}")
                logger.debug(f"Identified candidate sources: {source_candidates}")
                logger.debug(f"Identified candidate sinks: {sink_candidates}")

                if len(source_candidates) > 0 and len(sink_candidates) > 0:
                    source_node_indices = [
                        unitig_names_rev[x[:-1]] for x in source_candidates
                    ]
                    sink_node_indices = [
                        unitig_names_rev[x[:-1]] for x in sink_candidates
                    ]

                    # Create refined directed graph for flow network
                    # ----------------------------------------------------------------------
                    G = nx.DiGraph()

                    for u, v, cov in G_edge.edges(data=True):
                        if u not in node_indices:
                            node_indices[u] = my_counter
                            node_indices_rev[my_counter] = u
                            my_counter += 1
                        if v not in node_indices:
                            node_indices[v] = my_counter
                            node_indices_rev[my_counter] = v
                            my_counter += 1

                        logger.debug(f"Edge: {u}, {v}, {cov['weight']}")

                        G.add_edge(
                            node_indices[u], node_indices[v], weight=cov["weight"]
                        )

                    # Get connections and degree information
                    # ----------------------------------------------------------------------
                    in_degree = []
                    out_degree = []

                    for node in list(G.nodes):
                        if node_indices_rev[node][:-1] not in self_looped_nodes:
                            clean_indegree = len(
                                [
                                    x
                                    for x in G.predecessors(node)
                                    if node_indices_rev[x][:-1] not in self_looped_nodes
                                ]
                            )
                            in_degree.append(clean_indegree)

                            clean_outdegree = len(
                                [
                                    x
                                    for x in G.successors(node)
                                    if node_indices_rev[x][:-1] not in self_looped_nodes
                                ]
                            )
                            out_degree.append(clean_outdegree)

                    degrees = in_degree + out_degree

                    if len(degrees) == 0:
                        logger.debug(
                            f"Skipping component as no clean connections were found"
                        )
                        continue

                    # Create flow network
                    # ----------------------------------------------------------------------
                    network_edges = list()

                    edge_list_indices = {}

                    subpaths = {}
                    subpath_count = 0

                    visited_edges = []

                    logger.debug(f"G_edge.nodes: {list(G_edge.nodes)}")
                    logger.debug(f"G_edge.edges: {G_edge.edges(data=True)}")

                    for u, v, cov in G_edge.edges(data=True):
                        u_name = unitig_names_rev[u[:-1]]
                        v_name = unitig_names_rev[v[:-1]]

                        u_index = candidate_nodes.index(u_name) + 1
                        v_index = candidate_nodes.index(v_name) + 1

                        edge_list_indices[u_index] = u
                        edge_list_indices[v_index] = v

                        juction_cov = junction_pe_coverage[(u[:-1], v[:-1])]

                        if (u_index, v_index) not in visited_edges and (
                            v_index,
                            u_index,
                        ) not in visited_edges:
                            # Get coverage interval
                            cov_lower_bound = cov["weight"]
                            cov_upper_bound = int(max_comp_cov * alpha)

                            logger.debug(
                                f"({v}, {u}), ({u_index}, {v_index}) {juction_cov}, {cov_lower_bound}, {cov_upper_bound}"
                            )

                            if juction_cov == 0:
                                network_edges.append(
                                    (u_index, v_index, 0, cov_upper_bound)
                                )
                            else:
                                network_edges.append(
                                    (
                                        u_index,
                                        v_index,
                                        cov_lower_bound,
                                        cov_upper_bound,
                                    )
                                )

                            visited_edges.append((u_index, v_index))

                            # Add subpaths
                            if juction_cov >= mincov:
                                logger.debug(f"Adding subpath {[u_index, v_index]}")
                                subpaths[subpath_count] = [u_index, v_index]
                                subpath_count += 1

                                # Extend subpaths using coverages of successors and predecessors
                                # -----------------------------------

                                # Extend subpath using coverages of predecessors
                                for u_pred in G_edge.predecessors(u):
                                    u_pred_name = unitig_names_rev[u_pred[:-1]]
                                    u_pred_index = (
                                        candidate_nodes.index(u_pred_name) + 1
                                    )
                                    u_pred_cov = unitig_coverages[u_pred[:-1]]
                                    u_cov = unitig_coverages[u[:-1]]

                                    if (
                                        (v_index - 1) not in source_node_indices
                                        and (u_index - 1) not in source_node_indices
                                        and u_pred_index != v_index
                                    ):
                                        if (
                                            abs(min(u_pred_cov, u_cov) - cov["weight"])
                                            < covtol
                                        ):
                                            subpaths[subpath_count] = [
                                                u_pred_index,
                                                u_index,
                                                v_index,
                                            ]
                                            logger.debug(
                                                f"Extending subpath based on predecessor coverage {[u_pred_index, u_index, v_index]}"
                                            )
                                            subpath_count += 1

                                # Extend subpath using coverages of successors
                                for v_succ in G_edge.successors(v):
                                    v_succ_name = unitig_names_rev[v_succ[:-1]]
                                    v_succ_index = (
                                        candidate_nodes.index(v_succ_name) + 1
                                    )
                                    v_succ_cov = unitig_coverages[v_succ[:-1]]
                                    v_cov = unitig_coverages[v[:-1]]

                                    if (
                                        (v_succ_index - 1) not in source_node_indices
                                        and (u_index - 1) not in source_node_indices
                                        and (v_index - 1) not in source_node_indices
                                        and (v_index - 1) not in sink_node_indices
                                        and v_succ_index != u_index
                                    ):
                                        if (
                                            abs(min(v_succ_cov, v_cov) - cov["weight"])
                                            < covtol
                                        ):
                                            subpaths[subpath_count] = [
                                                u_index,
                                                v_index,
                                                v_succ_index,
                                            ]
                                            logger.debug(
                                                f"Extending subpath based on successor coverage {[u_index, v_index, v_succ_index]}"
                                            )
                                            subpath_count += 1

                            else:
                                # Extend subpaths of l=3 based on paired-end reads
                                # aligned to successors and predecessors
                                # -----------------------------------

                                # Extend subpath using coverages of predecessors
                                for u_pred in G_edge.predecessors(u):
                                    if junction_pe_coverage[(u_pred[:-1], v[:-1])] > 0:
                                        u_pred_name = unitig_names_rev[u_pred[:-1]]
                                        u_pred_index = (
                                            candidate_nodes.index(u_pred_name) + 1
                                        )
                                        if (
                                            (v_index - 1) not in source_node_indices
                                            and (u_index - 1) not in source_node_indices
                                            and u_pred_index != v_index
                                        ):
                                            subpaths[subpath_count] = [
                                                u_pred_index,
                                                u_index,
                                                v_index,
                                            ]
                                            logger.debug(
                                                f"Extending subpath {[u_pred_index, u_index, v_index]}"
                                            )
                                            subpath_count += 1

                                # Extend subpath using coverages of successors
                                for v_succ in G_edge.successors(v):
                                    if junction_pe_coverage[(u[:-1], v_succ[:-1])] > 0:
                                        v_succ_name = unitig_names_rev[v_succ[:-1]]
                                        v_succ_index = (
                                            candidate_nodes.index(v_succ_name) + 1
                                        )
                                        if (
                                            (v_succ_index - 1)
                                            not in source_node_indices
                                            and (u_index - 1) not in source_node_indices
                                            and (v_index - 1) not in source_node_indices
                                            and (v_index - 1) not in sink_node_indices
                                            and v_succ_index != u_index
                                        ):
                                            subpaths[subpath_count] = [
                                                u_index,
                                                v_index,
                                                v_succ_index,
                                            ]
                                            logger.debug(
                                                f"Extending subpath {[u_index, v_index, v_succ_index]}"
                                            )
                                            subpath_count += 1

                    # Add common start to source links
                    for source_v in source_candidates:
                        source_node_index = (
                            candidate_nodes.index(unitig_names_rev[source_v[:-1]]) + 1
                        )
                        source_node_cov = unitig_coverages[source_v[:-1]]
                        cov_upper_bound = int(max_comp_cov * alpha)

                        network_edges.append(
                            (
                                0,
                                source_node_index,
                                source_node_cov,
                                cov_upper_bound,
                            )
                        )

                        subpaths[subpath_count] = [0, source_node_index]
                        subpath_count += 1

                    # Add common sink to end links
                    for sink_v in sink_candidates:
                        sink_node_index = (
                            candidate_nodes.index(unitig_names_rev[sink_v[:-1]]) + 1
                        )
                        sink_node_cov = unitig_coverages[sink_v[:-1]]
                        cov_upper_bound = int(max_comp_cov * alpha)

                        network_edges.append(
                            (
                                sink_node_index,
                                len(candidate_nodes) + 1,
                                sink_node_cov,
                                cov_upper_bound,
                            )
                        )

                        subpaths[subpath_count] = [
                            sink_node_index,
                            len(candidate_nodes) + 1,
                        ]
                        subpath_count += 1

                    logger.debug(f"edge_list_indices: {edge_list_indices}")
                    logger.debug(f"subpaths: {subpaths}")

                    # Create flow network and run MFD-ILP
                    # ----------------------------------------------------------------------
                    G_mfd = {
                        "Nodes": len(list(G_edge.nodes)),
                        "list of edges": network_edges,
                        "subpaths": subpaths,
                    }
                    logger.debug(f"G_mfd: {G_mfd}")
                    solution_paths = flow_utils.solve_mfd(
                        G_mfd, maxpaths, output, nthreads
                    )
                    logger.debug(f"Number of paths found: {len(solution_paths)}")

                    case3_found.add(my_count)

                    # Iterate through solution paths
                    # ----------------------------------------------------------------------
                    if len(solution_paths) != 0:
                        phage_like_edges = phage_like_edges.union(
                            set(original_candidate_nodes)
                        )

                        cycle_number = 1

                        for solution_path in solution_paths:
                            coverage_val = solution_paths[solution_path]["weight"]

                            # Filter path by coverage
                            if coverage_val >= mincov:
                                logger.debug(
                                    f"Path {cycle_number} coverage: {coverage_val}"
                                )

                                # Create graph for path
                                G_path = nx.DiGraph()

                                # Fill graph with data
                                G_path.add_edges_from(
                                    solution_paths[solution_path]["path"]
                                )
                                logger.debug(
                                    f"solution path: {solution_paths[solution_path]['path']}"
                                )

                                if 0 in list(G_path.nodes):
                                    # Get all simple paths from node 0 to last node
                                    try:
                                        candidate_paths = list(
                                            nx.all_simple_paths(
                                                G_path, 0, len(candidate_nodes) + 1
                                            )
                                        )

                                        if len(candidate_paths) > 0:
                                            logger.debug(
                                                f"candidate_paths: {candidate_paths[0]}"
                                            )

                                            # Get mapped unitigs in order from the flow network
                                            path_order = []
                                            for path_edge in candidate_paths[0]:
                                                if not (
                                                    path_edge == 0
                                                    or path_edge
                                                    == len(candidate_nodes) + 1
                                                ):
                                                    path_order.append(
                                                        edge_list_indices[path_edge]
                                                    )

                                            logger.debug(f"path_order: {path_order}")

                                            # Get the order of unitigs in path
                                            path_string = ""
                                            total_length = 0

                                            previous_edge = 0

                                            for nodeid in range(len(path_order)):
                                                node = path_order[nodeid]
                                                unitig_name = node[:-1]

                                                if node.endswith("+"):
                                                    unitig_seq = str(
                                                        graph_unitigs[unitig_name]
                                                    )
                                                else:
                                                    unitig_seq = str(
                                                        graph_unitigs[
                                                            unitig_name
                                                        ].reverse_complement()
                                                    )

                                                # If first node in path
                                                if previous_edge == 0:
                                                    path_string += unitig_seq
                                                    total_length += len(unitig_seq)

                                                else:
                                                    trimmed_seq = unitig_seq[
                                                        link_overlap[
                                                            (previous_edge, node)
                                                        ] :
                                                    ]
                                                    path_string += trimmed_seq
                                                    total_length += len(trimmed_seq)

                                                previous_edge = node

                                            # Format genomic path
                                            path_node_order_human = ""

                                            for c in path_order:
                                                if c.endswith("+"):
                                                    path_node_order_human += (
                                                        f"{c[:-1]}:fwd,"
                                                    )
                                                else:
                                                    path_node_order_human += (
                                                        f"{c[:-1]}:rev,"
                                                    )

                                            path_node_order_human = (
                                                path_node_order_human[:-1]
                                            )

                                            # Create GenomePath object with path details
                                            genome_path = GenomePath(
                                                id=f"{prefix}phage_comp_{my_count}_cycle_{cycle_number}",
                                                bubble_case=case_name,
                                                node_order=[x for x in path_order],
                                                node_order_human=path_node_order_human,
                                                node_id_order=[
                                                    unitig_names_rev[x[:-1]]
                                                    for x in path_order
                                                ],
                                                path=path_string,
                                                coverage=int(coverage_val),
                                                length=total_length,
                                                gc=(
                                                    path_string.count("G")
                                                    + path_string.count("C")
                                                )
                                                / len(path_string)
                                                * 100,
                                            )
                                            my_genomic_paths.append(genome_path)
                                            logger.debug(
                                                f"total_length: {total_length}"
                                            )

                                            cycle_number += 1

                                    except nx.exception.NodeNotFound:
                                        logger.debug(
                                            f"Could not resolve a continuous path."
                                        )

                        logger.debug(f"Number of paths selected: {cycle_number-1}")

                        if cycle_number > 1:
                            resolved_components.add(my_count)
                            resolved_linear.add(my_count)
                            case3_resolved.add(my_count)

                    else:
                        logger.debug(f"No paths detected")
                        continue

        # Case 1 components - single unitigs
        elif len(candidate_nodes) == 1:
            unitig_name = unitig_names[candidate_nodes[0]]

            if unitig_name in self_looped_nodes or likely_complete[my_count]:
                case1_found.add(my_count)

                if unitig_name in self_looped_nodes:
                    case_name = "case1_circular"
                else:
                    case_name = "case1_linear"

                resolved_edges.add(candidate_nodes[0])
                comp_resolved_edges.add(candidate_nodes[0])

                path_string = str(graph_unitigs[unitig_name])

                cycle_number = 1

                # Create GenomePath object with path details
                genome_path = GenomePath(
                    id=f"{prefix}phage_comp_{my_count}_cycle_{cycle_number}",
                    bubble_case=case_name,
                    node_order=[unitig_names[candidate_nodes[0]]],
                    node_order_human=f"{unitig_names[candidate_nodes[0]]}:fwd",
                    node_id_order=[candidate_nodes[0]],
                    path=path_string,
                    coverage=int(unitig_coverages[unitig_name]),
                    length=len(graph_unitigs[unitig_name]),
                    gc=(path_string.count("G") + path_string.count("C"))
                    / len(path_string)
                    * 100,
                )
                my_genomic_paths.append(genome_path)
                resolved_components.add(my_count)
                single_unitigs.add(my_count)
                case1_resolved.add(my_count)

                phage_like_edges = phage_like_edges.union(set(candidate_nodes))

        # Record final paths for the component
        # ----------------------------------------------------------------------

        # Order resolved paths in descending order of length
        my_genomic_paths.sort(key=lambda x: (x.length, x.coverage), reverse=True)

        final_genomic_paths = []
        visited_nodes = set()
        comp_resolved_paths = set()

        frac_unitigs = 1
        n_paths = 0

        if len(my_genomic_paths) > 0:
            # Get the degree of the component
            graph_degree = assembly_graph.degree(original_candidate_nodes)

            path_lengths = []
            path_coverages = []

            largest_length = my_genomic_paths[0].length

            # Filter genomic paths
            for genomic_path in my_genomic_paths:
                passed = False

                if genomic_path.length > largest_length * LEN_THRESHOLD:
                    passed = True

                if case_name == "case3_linear":
                    passed = True

                path_node_order_string = ",".join(genomic_path.node_order)

                if path_node_order_string in comp_resolved_paths:
                    passed = False

                if passed:
                    logger.debug(
                        f"{genomic_path.id}\t{genomic_path.length}\t{genomic_path.coverage}"
                    )
                    logger.debug(f"{genomic_path.node_order}")
                    path_lengths.append(genomic_path.length)
                    path_coverages.append(genomic_path.coverage)
                    final_genomic_paths.append(genomic_path)
                    visited_nodes = visited_nodes.union(set(genomic_path.node_order))
                    comp_resolved_paths.add(path_node_order_string)
                    n_paths += 1

                    for path_node in genomic_path.node_id_order:
                        comp_resolved_edges.add(path_node)

            frac_unitigs = len(visited_nodes) / len(original_candidate_nodes)

            resolved_edges = resolved_edges.union(comp_resolved_edges)

            logger.debug(f"frac_unitigs: {frac_unitigs}")

            # Filter components
            if (
                len(final_genomic_paths) > 1
                and len(in_degree) > 0
                and len(out_degree) > 0
            ):
                coverage_frac = (
                    max(path_coverages) / min(path_coverages)
                    if min(path_coverages) > 0
                    else 1
                )

                # Create GenomeComponent object with component details
                genome_comp = GenomeComponent(
                    f"{prefix}phage_comp_{my_count}",
                    len(original_candidate_nodes),
                    len(final_genomic_paths),
                    max(graph_degree),
                    min(graph_degree),
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
                    coverage_frac,
                    frac_unitigs,
                )
                all_components.append(genome_comp)

            # Linear paths get the same length floor single-unitig components already
            # get in component_utils.get_components (edges_lengths[unitig] > minlength,
            # there called cicular_len) -- multi-unitig linear paths had no length gate
            # at all, so any component with a single phage-hallmark hit produced a
            # "resolved genome" regardless of how short the MFD-resolved path actually
            # came out. Circular paths are exempt: a closed cycle is itself strong
            # completeness evidence independent of length, matching why the
            # single-unitig gate only ever applied to that case in the first place.
            final_genomic_paths = [
                p
                for p in final_genomic_paths
                if not p.bubble_case.endswith("_linear") or p.length > minlength
            ]

            if len(final_genomic_paths) > 0:
                resolved_components.add(my_count)
                all_resolved_paths += final_genomic_paths
                component_elapsed_time = time.time() - component_time_start
                logger.debug(
                    f"Elapsed time to resolve component {my_count} with {len(original_candidate_nodes)} nodes: {component_elapsed_time} seconds"
                )

        else:
            # single unitigs
            for genomic_path in my_genomic_paths:
                final_genomic_paths.append(genomic_path)
                all_resolved_paths.append(genomic_path)
                logger.debug(f"{genomic_path.id}\t{genomic_path.length}")
                resolved_components.add(my_count)

        # Get unresolved edges
        unresolved_edges = comp_all_edges.difference(comp_resolved_edges)
        unresolved_phage_like_edges = unresolved_phage_like_edges.union(
            unresolved_edges
        )
        logger.debug(f"Unresolved edges in comp {my_count}: {unresolved_edges}")

        # Write genome path to file
        # ----------------------------------------------------------------------
        write_path(final_genomic_paths, output)
        write_path_fasta(final_genomic_paths, f"{output}/resolved_phages")

    return (
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
    )


# Set in the PARENT before the pool is created; forked workers inherit it
# through the process image, so none of it is ever serialised.
#
# This was previously passed as ProcessPoolExecutor(initargs=(kwargs,)), which
# pickles the whole thing once per worker. That is fine for a small graph and
# fatal for a real one: on a 485k-vertex assembly the payload is the igraph
# object plus every unitig's sequence in graph_unitigs -- gigabytes, serialised
# eight times through a pipe. It killed the pool seconds after startup, and
# multiprocessing additionally cannot send a single object larger than ~2GB at
# all. Inheriting by fork moves that cost to zero, and copy-on-write means the
# eight workers share the pages rather than each holding a full copy.
_WORKER_KWARGS = None


def _resolve_short_chunk(chunk):
    """Worker entry point: one slice of components, nothing else.

    Module-level (not a closure) so it is picklable by ProcessPoolExecutor --
    though only the small `chunk` is ever pickled; the bulk inputs arrive via
    _WORKER_KWARGS, inherited from the parent.
    """
    return resolve_short(pruned_vs=chunk, **_WORKER_KWARGS)


def resolve_short_parallel(
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
    workers=1,
):
    """
    resolve_short, with components processed in parallel across processes.

    Deliberately implemented WITHOUT touching resolve_short's ~1400-line body:
    that function is already parameterised by the set of components to process
    (`pruned_vs`) and already returns every accumulator it builds, so splitting
    the components into contiguous chunks, running one chunk per worker, and
    merging the returned tuples is equivalent to running it once over all of
    them. The alternative -- extracting the loop body into a per-component
    function -- would mean restructuring 1400 lines of nested branching for no
    additional benefit.

    This is only sound because no component's logic depends on another's
    results. Verified against the body: every touch of a shared accumulator is
    a pure `X.add(...)` / `X = X.union(...)` / `X.append(...)`, and there is not
    one conditional or membership test against them anywhere in the loop --
    per-component decisions use the loop-local `comp_*` sets instead. If that
    ever stops being true, chunking silently changes results, so it is worth
    re-checking before extending this.

    Chunks are CONTIGUOUS and merged in chunk order, so `all_resolved_paths`
    ends up in exactly the same order as the sequential run. That matters
    because downstream naming numbers the genomes by position: a different
    order would rename every genome without changing the biology, which is a
    nasty kind of non-reproducibility.

    Falls back to a plain sequential call for workers <= 1 or a single chunk,
    so the default path is byte-for-byte the code that ran before.
    """
    kwargs = dict(
        assembly_graph=assembly_graph,
        unitig_names=unitig_names,
        unitig_names_rev=unitig_names_rev,
        self_looped_nodes=self_looped_nodes,
        graph_unitigs=graph_unitigs,
        minlength=minlength,
        link_overlap=link_overlap,
        unitig_coverages=unitig_coverages,
        compcount=compcount,
        oriented_links=oriented_links,
        junction_pe_coverage=junction_pe_coverage,
        likely_complete=likely_complete,
        alpha=alpha,
        mincov=mincov,
        covtol=covtol,
        maxpaths=maxpaths,
        prefix=prefix,
        output=output,
        # Each worker is its own process, so the MILP solver inside it should
        # not also try to use every core -- that would oversubscribe by
        # workers x nthreads. Profiling showed solver threads make no
        # measurable difference anyway (model construction dominates), so 1 is
        # the right per-worker value rather than a compromise.
        nthreads=1,
    )

    keys = list(pruned_vs)
    if workers <= 1 or len(keys) <= 1:
        return resolve_short(pruned_vs=pruned_vs, **{**kwargs, "nthreads": nthreads})

    # Everything below has to cross a process boundary, so it all has to pickle.
    # Checked up front rather than discovered when the pool starts: a failure
    # there aborts the whole phables run, and losing a completed assembly to a
    # performance optimisation is a bad trade. Hit for real -- oriented_links
    # was a defaultdict built with a lambda, which cannot be pickled (fixed at
    # source in edge_graph_utils._oriented_links_inner), and it took out a real
    # 120-component run. Falling back to the sequential path keeps that a slow
    # run instead of a failed one.
    try:
        pickle.dumps(kwargs)
    except Exception as e:
        logger.warning(
            f"Cannot run components in parallel -- some input is not picklable "
            f"({type(e).__name__}: {e}). Falling back to sequential; the result "
            f"is unaffected, only the runtime."
        )
        return resolve_short(pruned_vs=pruned_vs, **{**kwargs, "nthreads": nthreads})

    workers = min(workers, len(keys))

    # Workers must INHERIT the bulk inputs rather than be sent them. On a real
    # assembly those inputs are gigabytes (see _WORKER_KWARGS), and pickling
    # them per worker is both ruinously slow and subject to multiprocessing's
    # ~2GB per-object ceiling. Fork gives the children the parent's memory
    # image directly, at no serialisation cost and, thanks to copy-on-write,
    # very little extra memory.
    #
    # If fork is unavailable (Windows; macOS defaults to spawn but can still
    # fork explicitly) there is no cheap way to hand over that much data, so
    # this runs sequentially rather than attempting a copy that would either
    # fail outright or exhaust the node's memory.
    try:
        mp_context = multiprocessing.get_context("fork")
    except ValueError:
        logger.warning(
            "The 'fork' start method is unavailable, so component-level "
            "parallelism would have to copy the whole assembly graph to every "
            "worker. Running sequentially instead; the result is unaffected, "
            "only the runtime."
        )
        return resolve_short(pruned_vs=pruned_vs, **{**kwargs, "nthreads": nthreads})

    # ONE COMPONENT PER TASK. Component cost is heavily skewed and which
    # components are expensive is not known in advance, so any chunk larger than
    # one lets a single task draw several expensive components while the rest of
    # the pool sits idle -- and everything queued behind them in that same chunk
    # waits too.
    #
    # This used to be ceil(len(keys) / (workers * 4)). Measured on a real
    # 424,819-vertex assembly (SRR12983578: 5,307 components reaching the
    # solver, median 3 nodes, max 197 after --compcount drops 15 larger ones),
    # comparing the worst task's share of total work:
    #
    #   cost model    chunks of 9      one per task     speedup of the critical path
    #   size              2.0%             0.7%              2.35x
    #   size^2           12.9%             6.8%              1.89x
    #   size^3           24.0%            14.1%              1.70x
    #
    # One component per task hits the hard floor in every model -- the floor
    # being the single largest component, which no split can subdivide. It is
    # also what an interleaved assignment achieves, without interleaving's
    # problem: `executor.map` preserves input order, and results are merged in
    # that order to stay identical to the sequential path (genome numbering
    # included), so reordering the input would have to be undone at the merge.
    # Keeping tasks in key order and simply making them smaller changes the
    # schedule without touching the result.
    #
    # The per-task cost is a dict of one component crossing the process
    # boundary; the bulk inputs are inherited by fork and never sent (see
    # _WORKER_KWARGS), so more tasks cost essentially nothing extra.
    size = max(1, math.ceil(len(keys) / MAX_MFD_TASKS))
    chunks = [
        {k: pruned_vs[k] for k in keys[i : i + size]} for i in range(0, len(keys), size)
    ]

    logger.info(
        f"Resolving {len(keys)} components across {workers} worker process(es) "
        f"in {len(chunks)} chunk(s)"
    )

    # Published to the module namespace BEFORE the pool exists, so every forked
    # child sees it. Cleared afterwards so the parent does not keep a second
    # reference to the graph alive for the rest of the run.
    global _WORKER_KWARGS
    _WORKER_KWARGS = kwargs
    try:
        with ProcessPoolExecutor(
            max_workers=workers, mp_context=mp_context
        ) as executor:
            # .map preserves input order, which is what keeps the merged
            # all_resolved_paths identical to the sequential run.
            results = list(executor.map(_resolve_short_chunk, chunks))
    except Exception as e:
        # A worker dying (OOM above all -- eight processes touching a large
        # graph can outgrow the node) surfaces here as BrokenProcessPool. Better
        # to spend the extra wall-clock than to lose an assembly that already
        # cost hours, so this retries the whole thing sequentially.
        logger.warning(
            f"Parallel component resolution failed ({type(e).__name__}: {e}). "
            f"Falling back to sequential -- the result is unaffected, only the "
            f"runtime. If this is memory, lower --mfd-workers."
        )
        return resolve_short(pruned_vs=pruned_vs, **{**kwargs, "nthreads": nthreads})
    finally:
        _WORKER_KWARGS = None

    # Field 1 (all_resolved_paths) and field 2 (all_components) are lists and
    # are extended; every other field is a set and is unioned. Merged in chunk
    # order for the reason given above.
    merged = list(results[0])
    for result in results[1:]:
        for i, value in enumerate(result):
            if isinstance(value, list):
                merged[i] = merged[i] + value
            else:
                merged[i] = merged[i] | value
    return tuple(merged)
