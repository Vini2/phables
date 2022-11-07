from Bio import SeqIO
from Bio.Seq import Seq
from igraph import *
from collections import defaultdict


class BidirectionalError(Exception):
    """Must set a unique value in a BijectiveMap."""

    def __init__(self, value):
        self.value = value
        msg = 'The value "{}" is already in the mapping.'
        super().__init__(msg.format(value))


class BidirectionalMap(dict):
    """Invertible map."""

    def __init__(self, inverse=None):
        if inverse is None:
            inverse = self.__class__(inverse=self)
        self.inverse = inverse

    def __setitem__(self, key, value):
        if value in self.inverse:
            raise BidirectionalError(value)

        self.inverse._set_item(value, key)
        self._set_item(key, value)

    def __delitem__(self, key):
        self.inverse._del_item(self[key])
        self._del_item(key)

    def _del_item(self, key):
        super().__delitem__(key)

    def _set_item(self, key, value):
        super().__setitem__(key, value)


def get_edge_lengths(edge_file):

    contig_lengths = {}

    for index, record in enumerate(SeqIO.parse(edge_file, "fasta")):
        contig_lengths[record.id] = len(record.seq)

    return contig_lengths


def get_links(assembly_graph_file):

    node_count = 0

    graph_contigs = {}

    edge_depths = {}

    edges_lengths = {}

    oriented_links = defaultdict(lambda: defaultdict(list))

    links = []

    my_map = BidirectionalMap()

    # Get links from .gfa file
    with open(assembly_graph_file) as file:

        for line in file.readlines():

            # Identify lines with link information
            if line.startswith("L"):

                strings = line.split("\t")

                link1 = strings[1]
                link2 = strings[3]

                link1_orientation = strings[2]
                link2_orientation = strings[4]
                # read_count = int(strings[6].split(":")[-1])

                link = []
                link.append(link1)
                link.append(link2)
                links.append(link)

                if link1 != link2:
                    if link1_orientation == "+" and link2_orientation == "+":
                        oriented_links[link1][link2].append(("+","+"))
                        oriented_links[link2][link1].append(("-","-"))
                    elif link1_orientation == "-" and link2_orientation == "-":
                        oriented_links[link1][link2].append(("-","-"))
                        oriented_links[link2][link1].append(("+","+"))
                    elif link1_orientation == "+" and link2_orientation == "-":
                        oriented_links[link1][link2].append(("+","-"))
                        oriented_links[link2][link1].append(("+","-"))
                    elif link1_orientation == "-" and link2_orientation == "+":
                        oriented_links[link1][link2].append(("-","+"))
                        oriented_links[link2][link1].append(("-","+"))
                    
                # if link1_orientation == "+" and link2_orientation == "+":
                #     link = []
                #     link.append(link1)
                #     link.append(link2)
                #     links.append(link)
                # elif link1_orientation == "-" and link2_orientation == "-":
                #     link = []
                #     link.append(link2)
                #     link.append(link1)
                #     links.append(link)
                # else:
                #     link = []
                #     link.append(link1)
                #     link.append(link2)
                #     links.append(link)
                

            elif line.startswith("S"):

                strings = line.strip().split()

                my_map[node_count] = strings[1]

                graph_contigs[strings[1]] = Seq(strings[2])

                depth = int(strings[3].split(":")[-1])
                edge_depths[strings[1]] = depth
                edges_lengths[strings[1]] = len(strings[2])

                node_count += 1

            line = file.readline()

    return node_count, graph_contigs, links, oriented_links, my_map, edge_depths, edges_lengths


def get_graph_edges(links, contig_names_rev):

    self_looped_nodes = []

    edge_list = []
    # weights = []
    # weights_dict = {}

    # Iterate links
    for link in links:
        # Remove self loops
        if link[0] != link[1]:
            # Add edge to list of edges
            edge_list.append((contig_names_rev[link[0]], contig_names_rev[link[1]]))
            # weights.append(link[2])
            # weights_dict[(contig_names_rev[link[0]], contig_names_rev[link[1]])] = link[2]
            # weights_dict[(contig_names_rev[link[1]], contig_names_rev[link[0]])] = link[2]
        else:
            self_looped_nodes.append(link[0])

    return edge_list, self_looped_nodes


def build_assembly_graph(assembly_graph_file):

    (
        node_count,
        graph_contigs,
        links,
        oriented_links,
        contig_names,
        edge_depths,
        edges_lengths,
    ) = get_links(assembly_graph_file)

    # Get reverse mapping of contig identifiers
    contig_names_rev = contig_names.inverse

    # Create graph
    assembly_graph = Graph(directed=False)

    # Add vertices
    assembly_graph.add_vertices(node_count)
    # print("Total number of contigs available: " + str(len(list(assembly_graph.vs))))

    # Name vertices with contig identifiers
    for i in range(node_count):
        assembly_graph.vs[i]["id"] = i
        assembly_graph.vs[i]["name"] = contig_names[i]
        assembly_graph.vs[i]["label"] = contig_names[i] + "\nID:" + str(i)

    edge_list, self_looped_nodes = get_graph_edges(
        links=links, contig_names_rev=contig_names_rev
    )

    # print(len(edge_list), edge_list)

    # Add edges to the graph
    assembly_graph.add_edges(edge_list)
    # assembly_graph.es['weight'] = weights
    # assembly_graph.es['label'] = weights
    # assembly_graph.es['name'] = weights

    # Simplify the graph
    assembly_graph.simplify(multiple=True, loops=False, combine_edges=None)

    # print("Total number of edges in the assembly graph: " + str(len(list(assembly_graph.es))))

    return (
        assembly_graph,
        oriented_links,
        contig_names,
        contig_names_rev,
        graph_contigs,
        self_looped_nodes,
        edges_lengths,
    )


def get_circular(paths):
    
    circular = {}

    with open(paths, "r") as myfile:

        for line in myfile.readlines():
            if not line.startswith("#"):
                strings = line.strip().split()

                if strings[3] == "Y":
                    contig_name = strings[0].replace("contig", "edge")
                    contig_length = int(strings[1])
                    circular[contig_name] = contig_length

    return circular


def remove_dead_ends(G_edge):

    has_dead_ends = True

    while has_dead_ends:

        to_remove = []

        for node in list(G_edge.nodes):
            if not (G_edge.in_degree(node) > 0 and G_edge.out_degree()(node)) > 0:
                to_remove.append(node)
        

        if len(to_remove) > 0:
            G_edge.remove_nodes_from(to_remove)
            print("remove_nodes_from", to_remove, has_dead_ends)
        else:
            has_dead_ends = False

    return G_edge