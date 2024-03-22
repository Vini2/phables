import copy
import logging
from collections import defaultdict

from Bio import SeqIO
from Bio.Seq import Seq
from igraph import Graph

# Create logger
logger = logging.getLogger("phables 1.3.2")


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


def get_unitig_lengths(edge_file):
    """
    Get length of the unitigs
    """

    unitig_lengths = {}

    for index, record in enumerate(SeqIO.parse(edge_file, "fasta")):
        unitig_lengths[record.id] = len(record.seq)

    return unitig_lengths


def get_links(assembly_graph_file):
    """
    Get links from the assembly graph
    """

    node_count = 0
    graph_contigs = {}
    edges_lengths = {}
    oriented_links = defaultdict(lambda: defaultdict(list))
    link_overlap = defaultdict(int)
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
                overlap = int(strings[5].strip()[:-1])

                link = []
                link.append(link1)
                link.append(link2)
                links.append(link)

                if link1 != link2:
                    if link1_orientation == "+" and link2_orientation == "+":
                        oriented_links[link1][link2].append(("+", "+"))
                        link_overlap[(f"{link1}+", f"{link2}+")] = overlap
                        oriented_links[link2][link1].append(("-", "-"))
                        link_overlap[(f"{link2}-", f"{link1}-")] = overlap
                    elif link1_orientation == "-" and link2_orientation == "-":
                        oriented_links[link1][link2].append(("-", "-"))
                        link_overlap[(f"{link1}-", f"{link2}-")] = overlap
                        oriented_links[link2][link1].append(("+", "+"))
                        link_overlap[(f"{link2}+", f"{link1}+")] = overlap
                    elif link1_orientation == "+" and link2_orientation == "-":
                        oriented_links[link1][link2].append(("+", "-"))
                        link_overlap[(f"{link1}+", f"{link2}-")] = overlap
                        oriented_links[link2][link1].append(("+", "-"))
                        link_overlap[(f"{link2}+", f"{link1}-")] = overlap
                    elif link1_orientation == "-" and link2_orientation == "+":
                        oriented_links[link1][link2].append(("-", "+"))
                        link_overlap[(f"{link1}-", f"{link2}+")] = overlap
                        oriented_links[link2][link1].append(("-", "+"))
                        link_overlap[(f"{link2}-", f"{link1}+")] = overlap

            elif line.startswith("S"):
                strings = line.strip().split()
                my_map[node_count] = strings[1]
                graph_contigs[strings[1]] = Seq(strings[2])
                edges_lengths[strings[1]] = len(strings[2])
                node_count += 1

            line = file.readline()

    return (
        node_count,
        graph_contigs,
        links,
        oriented_links,
        link_overlap,
        my_map,
        edges_lengths,
    )


def get_graph_edges(links, contig_names_rev):
    """
    Returns the edges of the assembly graph
    """

    self_looped_nodes = []

    edge_list = []

    # Iterate links
    for link in links:
        # Remove self loops
        if link[0] != link[1]:
            # Add edge to list of edges
            edge_list.append((contig_names_rev[link[0]], contig_names_rev[link[1]]))
        else:
            self_looped_nodes.append(link[0])

    return edge_list, self_looped_nodes


def build_assembly_graph(assembly_graph_file):
    """
    Build the assembly graph
    """

    (
        node_count,
        graph_contigs,
        links,
        oriented_links,
        link_overlap,
        contig_names,
        edges_lengths,
    ) = get_links(assembly_graph_file)

    # Get reverse mapping of contig identifiers
    contig_names_rev = contig_names.inverse

    # Create graph
    assembly_graph = Graph(directed=False)

    # Add vertices
    assembly_graph.add_vertices(node_count)

    # Name vertices with contig identifiers
    for i in range(node_count):
        assembly_graph.vs[i]["id"] = i
        assembly_graph.vs[i]["name"] = contig_names[i]
        assembly_graph.vs[i]["label"] = contig_names[i] + "\nID:" + str(i)

    edge_list, self_looped_nodes = get_graph_edges(
        links=links, contig_names_rev=contig_names_rev
    )

    # Add edges to the graph
    assembly_graph.add_edges(edge_list)

    # Simplify the graph
    assembly_graph.simplify(multiple=True, loops=False, combine_edges=None)

    return (
        assembly_graph,
        oriented_links,
        link_overlap,
        contig_names,
        contig_names_rev,
        graph_contigs,
        self_looped_nodes,
        edges_lengths,
    )


def get_circular(self_looped_nodes, graph_unitigs):
    """
    Get circular unitigs
    """

    circular = {}

    for unitig in self_looped_nodes:
        circular[unitig] = len(str(graph_unitigs[unitig]))

    # with open(paths, "r") as myfile:

    #     for line in myfile.readlines():
    #         if not line.startswith("#"):
    #             strings = line.strip().split()

    #             if strings[3] == "Y":
    #                 contig_name = strings[0].replace("contig", "edge")
    #                 contig_length = int(strings[1])
    #                 circular[contig_name] = contig_length

    return circular


def remove_dead_ends(G_edge):
    """
    Remove dead-ends from the component
    """

    new_G = copy.deepcopy(G_edge)

    has_dead_ends = True

    dead_ends_to_remove = []

    while has_dead_ends:
        to_remove = []

        for node in list(new_G.nodes):
            if not (new_G.in_degree(node) > 0 and new_G.out_degree()(node)) > 0:
                to_remove.append(node)

        if len(to_remove) > 0:
            new_G.remove_nodes_from(to_remove)
            logger.debug(f"Removing dead-ends: {to_remove}")
        else:
            has_dead_ends = False

        dead_ends_to_remove += to_remove

    return set(dead_ends_to_remove)


def get_all_sub_paths(assembly_graph, unitig_names):
    """
    Get all sub paths of length 2 and 3
    """

    sub_paths = defaultdict(int)

    for v in range(assembly_graph.vcount()):
        # Get all paths starting from vertex 'v' of length exactly 2
        paths_from_v = assembly_graph.get_all_simple_paths(v, cutoff=2)

        for path in paths_from_v:
            if len(path) == 3:  # Length 3 means 3 vertices (2 edges)
                node1 = unitig_names[path[0]]
                node2 = unitig_names[path[1]]
                node3 = unitig_names[path[2]]
                sub_paths[tuple([node1, node2, node3])] = 0

            elif len(path) == 2:  # Length 2 means 2 vertices (1 edge)
                node1 = unitig_names[path[0]]
                node2 = unitig_names[path[1]]
                sub_paths[tuple([node1, node2])] = 0

    return sub_paths
