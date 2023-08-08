import networkx as nx

from .FD_Inexact import SolveInstances


def get_source_sink_circular(G_edge, graph_unitigs, minlength, self_looped_nodes):
    """
    Identify source/sink vertex for circular components
    """

    source_sink_candidates = []

    for node in list(G_edge.nodes):
        unitig_name = node[:-1]

        if (
            unitig_name not in self_looped_nodes
            and len(graph_unitigs[unitig_name]) > minlength
        ):
            # Get BFS layers
            bfs_layers = dict(enumerate(nx.bfs_layers(G_edge, node)))

            # Get last later
            last_layer = list(bfs_layers.keys())[-1]

            node_is_st = True

            # Check if successors of those in last_layer is same as the node
            for item in bfs_layers[last_layer]:
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

            if len(bfs_layers[last_layer]) == 0:
                node_is_st = False

            if node_is_st:
                source_sink_candidates.append(node)

    return source_sink_candidates


def get_source_sink_linear(G_edge, self_looped_nodes):
    """
    Identify source/sink vertex for linear components
    """

    source_candidates = []
    sink_candidates = []

    for node in list(G_edge.nodes):
        unitig_name = node[:-1]

        if unitig_name not in self_looped_nodes:
            indegree = len([x for x in G_edge.predecessors(node)])
            outdegree = len([x for x in G_edge.successors(node)])
            if indegree > 0 and outdegree == 0:
                sink_candidates.append(node)
            elif indegree == 0 and outdegree > 0:
                source_candidates.append(node)

    return source_candidates, sink_candidates


def solve_mfd(G, max_paths, output, nthreads):
    """
    Get paths by solving MFD
    """

    listOfGraphs = {}
    listOfGraphs[0] = G

    outputfile = f"{output}/results_MFD.txt"
    recordfile = f"{output}/results_MFD_details.txt"

    solution_paths = SolveInstances(
        listOfGraphs, max_paths, outputfile, recordfile, nthreads
    )

    return solution_paths
