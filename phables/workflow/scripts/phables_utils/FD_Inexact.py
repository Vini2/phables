#!/usr/bin/env python3

# Source: https://github.com/algbio/flowpaths

import itertools
import logging
import math

import flowpaths as fp
import networkx as nx

__author__ = "Vijini Mallawaarachchi"
__copyright__ = "Copyright 2026"
__license__ = "MIT"
__version__ = "0.5.0"
__maintainer__ = "Vijini Mallawaarachchi"
__email__ = "viji.mallawaarachchi@gmail.com"
__status__ = "Development"


# Create logger
logger = logging.getLogger(__name__)


def read_input(graphfile, number_subpath):
    trip_data = open(graphfile, "r").read().split("\n")
    i = 0
    listOfGraphs = {}
    k = 0

    while True:
        if "#" in trip_data[i]:
            i = i + 1
            N = int(trip_data[i])
            edges = list()
            subpaths = {}
            while True:
                i = i + 1
                if "#" in trip_data[i]:
                    break
                if "" == trip_data[i]:
                    break
                if "subpaths" in trip_data[i]:
                    for j in range(0, number_subpath):
                        i = i + 1
                        line = trip_data[i].split(" ")
                        subpaths[j] = line[0 : len(line) - 1]
                    i = i + 4
                if i >= len(trip_data) - 1:
                    break
                line = trip_data[i].split(" ")
                edges.append((line[0], line[1], line[2], line[3]))
            G = {"Nodes": N, "list of edges": edges, "subpaths": subpaths}
            listOfGraphs[k] = G
            k += 1
        if i >= len(trip_data) - 1:
            break

    return listOfGraphs


# FD-Subpath-Inexact-Flowpaths
# --------------------------------------------
class InexactFlowDecomposition(fp.AbstractPathModelDAG):
    def __init__(
        self,
        G,
        lb,
        ub,
        num_paths,
        subpath_constraints=None,
        threads=1,
    ):
        self.G = fp.stDAG(G)
        self.lb = lb
        self.ub = ub
        self._solution = None

        trusted_edges_for_safety = self.G.get_non_zero_flow_edges(flow_attr=self.lb)

        super().__init__(
            self.G,
            num_paths,
            subpath_constraints=subpath_constraints or [],
            optimization_options={
                "trusted_edges_for_safety": trusted_edges_for_safety,
            },
            solver_options={
                "threads": threads,
                "log_to_console": "false",
            },
        )

        self.create_solver_and_paths()
        self._encode_flow_intervals()

    def _encode_flow_intervals(self):
        maximum_allowed_path_weight = max(
            data.get(self.ub, 0) for _, _, data in self.G.edges(data=True)
        )

        self.path_weights_vars = self.solver.add_variables(
            self.path_indexes,
            name_prefix="w",
            lb=0,
            ub=maximum_allowed_path_weight,
            var_type="integer",
        )
        self.pi_vars = self.solver.add_variables(
            self.edge_indexes,
            name_prefix="pi",
            lb=0,
            ub=maximum_allowed_path_weight,
        )

        for u, v, data in self.G.edges(data=True):
            if (u, v) in self.G.source_sink_edges:
                continue

            for i in range(self.k):
                self.solver.add_binary_continuous_product_constraint(
                    binary_var=self.edge_vars[(u, v, i)],
                    continuous_var=self.path_weights_vars[i],
                    product_var=self.pi_vars[(u, v, i)],
                    lb=0,
                    ub=maximum_allowed_path_weight,
                    name=f"product_u={u}_v={v}_i={i}",
                )

            edge_weight = self.solver.quicksum(
                self.pi_vars[(u, v, i)] for i in range(self.k)
            )
            self.solver.add_constraint(
                edge_weight >= data[self.lb],
                name=f"lowerbound_u={u}_v={v}",
            )
            self.solver.add_constraint(
                edge_weight <= data[self.ub],
                name=f"upperbound_u={u}_v={v}",
            )

    def get_solution(self):
        self.check_is_solved()

        if self._solution is not None:
            return self._solution

        weights = self.solver.get_values(self.path_weights_vars)
        self._solution = {
            "paths": self.get_solution_paths(),
            "weights": [weights[i] for i in range(self.k)],
        }

        return self._solution

    def get_objective_value(self):
        return self.solver.get_objective_value()

    def get_lowerbound_k(self):
        weight_function = {
            (u, v): 1
            for u, v, edge_data in self.G.edges(data=True)
            if edge_data.get(self.lb, 0) > 0
        }
        return self.G.compute_max_edge_antichain(weight_function=weight_function)

    def is_valid_solution(self):
        if not self.is_solved():
            return False

        solution = self.get_solution()
        for u, v, edge_data in self.G.base_graph.edges(data=True):
            observed_weight = 0
            for path, weight in zip(solution["paths"], solution["weights"]):
                path_edges = set(itertools.pairwise(path))
                if (u, v) in path_edges:
                    observed_weight += weight

            if (
                observed_weight < edge_data[self.lb]
                or observed_weight > edge_data[self.ub]
            ):
                return False

        return True


def _get_subpath_constraints(subpaths, edges, to_flowpaths_node=str):
    subpath_constraints = []
    for subpath in subpaths.values():
        subpath_edges = list(itertools.pairwise(subpath))
        if subpath_edges and all(edge in edges for edge in subpath_edges):
            subpath_constraints.append(
                [
                    (to_flowpaths_node(edge[0]), to_flowpaths_node(edge[1]))
                    for edge in subpath_edges
                ]
            )

    return subpath_constraints


def _path_to_edges(path, from_flowpaths_node=lambda node: node):
    restored_path = [from_flowpaths_node(node) for node in path]
    return list(itertools.pairwise(restored_path))


def flowMultipleDecomposition(data, K, nthreads):
    graph = data["graph"]
    node_labels = {node: str(node) for node in graph.nodes}
    original_node_labels = {
        flowpaths_node: node for node, flowpaths_node in node_labels.items()
    }
    flowpaths_graph = nx.relabel_nodes(graph, node_labels, copy=True)
    subpath_constraints = _get_subpath_constraints(
        data["subpaths"], data["edges"], node_labels.__getitem__
    )

    try:
        model = InexactFlowDecomposition(
            flowpaths_graph,
            lb="flow_low",
            ub="flow_up",
            num_paths=K,
            subpath_constraints=subpath_constraints,
            threads=nthreads,
        )
        model.solve()

        if model.is_solved():
            solution = model.get_solution()
            data["message"] = "solved"
            data["runtime"] = model.solve_statistics.get(
                f"milp_solve_time_for_num_paths_{K}", 0
            )
            data["weights"] = solution["weights"]
            data["solution"] = [
                _path_to_edges(path, original_node_labels.__getitem__)
                for path in solution["paths"]
            ]
        else:
            data["message"] = "unsolved"
            data["runtime"] = model.solve_statistics.get(
                f"milp_solve_time_for_num_paths_{K}", 0
            )
    except (ValueError, RuntimeError, AttributeError) as e:
        data["message"] = "unsolved"
        logger.debug(f"Flowpaths could not solve the MFD instance with K={K}: {e}")

    return data


def get_lowerbound_k(data):
    """
    Smallest number of paths any valid decomposition of this component could
    possibly use.

    The K search below tries K = 1, 2, 3, ... until one is feasible, rebuilding
    the whole MILP each time (K is structural to the model -- every variable is
    indexed by it -- so it genuinely cannot be reused, and neither flowpaths nor
    HiGHS-via-flowpaths offers a warm start). Every attempt below the true answer
    is therefore a complete model build that can only ever come back infeasible,
    and model construction is ~95% of the cost of an attempt. Starting the search
    at a proven lower bound skips exactly those wasted builds.

    Two bounds, both taken from flowpaths' own MinFlowDecomp.get_lowerbound_k:

    - the graph's WIDTH: the minimum number of paths needed to cover every edge.
      Any decomposition must cover every edge carrying flow, so it cannot use
      fewer paths than this.
    - ceil(log2(number of distinct flow values)): k paths can produce at most
      2^k distinct subset sums, so k must be at least log2 of the number of
      distinct values that have to be represented.

    Both are lower bounds, so the maximum of them is too, and starting there can
    never skip a feasible smaller K -- the search still returns the same first
    feasible K, just without the doomed attempts before it. Verified on synthetic
    components: identical K and identical path sets with and without this.

    Returns 1 (i.e. the original behaviour) if anything about the computation
    fails. A lower bound is an optimisation, not a correctness requirement, so it
    must never be the reason a component stops resolving.
    """
    try:
        graph = data["graph"]
        if graph.number_of_edges() == 0:
            return 1

        # flowpaths' stDAG wants string node ids, same relabelling
        # flowMultipleDecomposition does before building its model.
        node_labels = {node: str(node) for node in graph.nodes}
        flowpaths_graph = nx.relabel_nodes(graph, node_labels, copy=True)
        lowerbound = fp.stDAG(flowpaths_graph).get_width()

        # Only edges that must carry flow constrain the decomposition; an edge
        # whose lower bound is 0 need not be covered at all.
        distinct_flows = {
            int(data["flows_low"][edge])
            for edge in data["edges"]
            if data["flows_low"].get(edge, 0) > 0
        }
        if len(distinct_flows) > 1:
            lowerbound = max(lowerbound, math.ceil(math.log2(len(distinct_flows))))

        return max(1, lowerbound)
    except Exception as e:
        logger.debug(f"Could not compute a lower bound on K ({e}); starting from 1")
        return 1


def FD_Algorithm(data, max_paths, nthreads):
    solutionWeights = 0
    solutionSet = 0

    # See get_lowerbound_k: K < lowerbound is provably infeasible, and each such
    # attempt costs a full MILP build. Capped at max_paths so a bound above the
    # user's --maxpaths doesn't turn into a search over an empty range with
    # different semantics -- that case has no solution within max_paths either
    # way, and this keeps the "attempt at least one K" behaviour identical.
    lowerbound = min(get_lowerbound_k(data), max_paths)
    if lowerbound > 1:
        logger.debug(
            f"Starting the K search at {lowerbound} rather than 1 "
            f"({lowerbound - 1} provably-infeasible attempt(s) skipped)"
        )

    for i in range(lowerbound, max_paths + 1):
        data = flowMultipleDecomposition(data, i, nthreads)
        if data["message"] == "solved":
            solutionSet = data["solution"]
            solutionWeights = data["weights"]
            break

    # Get solution paths and weights
    solution_paths = {}

    if solutionSet != 0:
        for i in range(0, len(solutionSet)):
            solution_paths[i] = {"weight": solutionWeights[i], "path": solutionSet[i]}
            # print("W:",solutionWeights[i], solutionSet[i])

    return data, solution_paths


def SolveInstances(Graphs, max_paths, outfile, recfile, nthreads):
    # Not populated with per-path debug output -- that was already true in the
    # original Gurobi-based implementation (the file handles were opened but
    # never written to). Still created here because workflow/test_phables.smk
    # declares both as Snakemake rule outputs; a declared output missing from
    # disk after the rule runs is a hard Snakemake failure, not a warning.
    open(outfile, "w").close()
    open(recfile, "w").close()

    for s in range(0, 1):
        f_low = {}
        f_up = {}
        Edges = set()
        V = set()
        listOfEdges = Graphs[s]["list of edges"]

        for k in range(0, len(listOfEdges)):
            a, b, c, d = listOfEdges[k]
            Edges.add((a, b))
            V.add(a)
            V.add(b)
            f_low[a, b] = int(float(c))
            f_up[a, b] = int(float(d))

        # creation of graphs
        # creation of graphs
        G = nx.DiGraph()
        G.add_nodes_from(V)
        for edge in Edges:
            G.add_edge(edge[0], edge[1], flow_low=f_low[edge], flow_up=f_up[edge])

        # creation of adjacent matrix
        AD_in = {}
        AD_out = {}

        for v in V:
            setAdj = set()
            for i, j in list(G.out_edges(v)):
                if i != v:
                    setAdj.add(i)
                if j != v:
                    setAdj.add(j)

            AD_out[v] = list(setAdj)

            setAdj = set()
            for i, j in list(G.in_edges(v)):
                if i != v:
                    setAdj.add(i)
                if j != v:
                    setAdj.add(j)

            AD_in[v] = list(setAdj)

        # calculating source, sinks and max flows
        S = [x for x in G.nodes() if G.out_degree(x) >= 1 and G.in_degree(x) == 0]
        D = [x for x in G.nodes() if G.out_degree(x) == 0 and G.in_degree(x) >= 1]
        maxW = max(f_up.values())

        # definition of data

        data = {
            "edges": Edges,
            "flows_low": f_low,
            "flows_up": f_up,
            "vertices": V,
            "graph": G,
            "Kmax": len(Edges),
            "weights": {},
            "sources": S,
            "targets": D,
            "message": {},
            "solution": 0,
            "maxFlow": maxW,
            "adj_in": AD_in,
            "adj_out": AD_out,
            "subpaths": Graphs[s]["subpaths"],
            # Never read by anything -- the K search takes its starting point
            # from get_lowerbound_k(), which computes a real bound per component
            # rather than assuming a constant. Left in place because this dict is
            # passed around wholesale and removing a key is a wider change than
            # it looks; flagged so it isn't mistaken for the live lower bound.
            "minK": 2,
            "runtime": 0,
        }

        data, solution_paths = FD_Algorithm(data, max_paths, nthreads)

    return solution_paths
