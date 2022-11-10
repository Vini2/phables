import networkx as nx

from .FD_Inexact import SolveInstances


def get_source_sink(G_edge, graph_contigs, minlength, self_looped_nodes):
    """
    Identify source/sink vertex
    """

    source_sink_candidates = []

    for node in list(G_edge.nodes):

        contig_name = node[:-1]

        if (
            contig_name not in self_looped_nodes
            and len(graph_contigs[contig_name]) > minlength
        ):

            bfs_successors = dict(enumerate(nx.bfs_successors(G_edge, node)))

            last_layer = list(bfs_successors.keys())[-1]

            node_is_st = True

            for item in bfs_successors[last_layer][1]:

                if item[:-1] not in self_looped_nodes:

                    item_successors = list(G_edge.successors(item))

                    if (
                        len(item_successors) > 0
                        and list(G_edge.successors(item))[0] != node
                    ):
                        node_is_st = False
                        break
                    if len(item_successors) == 0:
                        node_is_st = False

            if len(bfs_successors[last_layer][1]) == 0:
                node_is_st = False

            if node_is_st:
                source_sink_candidates.append(node)

    return source_sink_candidates


def solve_mfd(G, max_paths, output):
    """
    Get paths by solving MFD
    """

    listOfGraphs = {}
    listOfGraphs[0] = G

    outputfile = f"{output}/results_MFD.txt"
    recordfile = f"{output}/results_MFD_details.txt"

    solution_paths = SolveInstances(listOfGraphs, max_paths, outputfile, recordfile)

    return solution_paths
