import networkx as nx

from .FD_Inexact import SolveInstances

def get_source_sink(G_edge, graph_contigs, minlength, self_looped_nodes, contig_names_rev):

    source_sink_candidates = []

    for node in list(G_edge.nodes):
    
        contig_name = node[:-1]

        if contig_name not in self_looped_nodes and len(graph_contigs[contig_name]) > minlength:

            bfs_successors = dict(enumerate(nx.bfs_successors(G_edge, node)))
            print(node, "-->", bfs_successors)
            
            last_layer = list(bfs_successors.keys())[-1]

            node_is_st = True

            for item in bfs_successors[last_layer][1]:

                if item[:-1] not in self_looped_nodes:

                    item_successors = list(G_edge.successors(item))
                    print(item, item_successors)

                    if len(item_successors) > 0 and list(G_edge.successors(item))[0] != node:
                        node_is_st = False
                        break
                    if len(item_successors) == 0:
                        node_is_st = False

            if len(bfs_successors[last_layer][1]) == 0:
                node_is_st = False

            if node_is_st:
                source_sink_candidates.append(contig_names_rev[contig_name])

    return source_sink_candidates


def solve_mfd(G, output):

    listOfGraphs = {}
    listOfGraphs[0] = G

    print(listOfGraphs)

    outputfile = output+"results_gurobi_subpath.txt"
    recordfile = output+"results_gurobi_subpath_details.txt"

    solution_paths = SolveInstances(listOfGraphs, outputfile, recordfile)

    return solution_paths