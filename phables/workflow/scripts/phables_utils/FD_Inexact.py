#!/usr/bin/env python3

# Source: https://github.com/algbio/flowpaths

import itertools
import json
import logging
import math
import os
import time

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
# The same named logger every other phables_utils module uses. phables.py
# attaches its DEBUG file handler to "phables 2.0.0"; a logger named __name__
# ("phables_utils.FD_Inexact") is not a child of that, so with the previous
# line every message from this module -- including per-K attempt timings --
# went to the unconfigured root logger and never appeared in
# phables_output.log. That is why a 9-hour flow decomposition on Setonix left
# no record of which K it was stuck on.
logger = logging.getLogger("phables 2.0.0")

# Runtime knobs, set once per run by phables.py via configure(). Module globals
# rather than new positional parameters on resolve_short()/resolve_long():
# those signatures are ~20 arguments deep across two 1400-line functions and
# two call sites each, and --mfd-workers children are forked, so a module
# global set in the parent is inherited by every worker for free.
_MFD_TIME_LIMIT = float("inf")   # per-attempt HiGHS time limit in seconds
_MFD_DUMP_SLOW_S = None          # dump the instance when a component exceeds this
_MFD_DUMP_DIR = None


def configure(time_limit=None, dump_slow_s=None, dump_dir=None):
    """Set the per-run MFD knobs (see the globals above). Called by phables.py."""
    global _MFD_TIME_LIMIT, _MFD_DUMP_SLOW_S, _MFD_DUMP_DIR
    if time_limit is not None and float(time_limit) > 0:
        _MFD_TIME_LIMIT = float(time_limit)
    if dump_slow_s is not None:
        _MFD_DUMP_SLOW_S = float(dump_slow_s)
        _MFD_DUMP_DIR = dump_dir


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
        allow_empty_paths=False,
        time_limit=None,
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
                # True only for the feasibility ORACLE in FD_Algorithm: "does any
                # decomposition with <= K paths exist?". Paths may then be empty,
                # which makes feasibility exactly monotone in K by construction.
                # Never True for a model whose solution is returned.
                "allow_empty_paths": allow_empty_paths,
            },
            solver_options={
                "threads": threads,
                "log_to_console": "false",
                **({"time_limit": float(time_limit)} if time_limit is not None and time_limit != float("inf") else {}),
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


def flowMultipleDecomposition(data, K, nthreads, allow_empty_paths=False):
    graph = data["graph"]
    t_build = time.perf_counter()
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
            allow_empty_paths=allow_empty_paths,
            time_limit=_MFD_TIME_LIMIT,
        )
        data["build_time"] = time.perf_counter() - t_build
        model.solve()
        # HiGHS status string: kOptimal / kInfeasible / kTimeLimit / ... The
        # search below needs to tell "proved infeasible" from "ran out of time".
        try:
            data["status"] = model.solver.get_model_status()
        except Exception:
            data["status"] = "unknown"

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
        data["status"] = "error"
        data.setdefault("build_time", time.perf_counter() - t_build)
        logger.debug(f"Flowpaths could not solve the MFD instance with K={K}: {e}")

    return data


def get_lowerbound_k(data):
    """
    A VALID lower bound on the number of paths in any decomposition of this
    component: the maximum antichain of edges that must carry flow.

    Every edge whose lower bound is positive has to lie on at least one path,
    and two edges that no single source-to-sink path can both traverse need two
    different paths. So the largest set of pairwise-incomparable positive edges
    is a lower bound on K. Nothing else is assumed about the flows.

    What this replaces, and why (both measured on 17,189 real components dumped
    from Setonix logs, checked against the original K = 1, 2, 3, ... ladder):

    - stDAG.get_width(), the minimum path cover of ALL edges. Edges with a zero
      lower bound need not be covered at all, and 28% of real edges have one
      (junction coverage 0), so this over-counted in 12% of components and made
      the search START ABOVE the true minimum -- returning a 2-path
      decomposition where 1 path suffices, and so on. 75 of 99 wrong answers.
    - ceil(log2(#distinct flow values)). Valid for exact flows (k paths give at
      most 2^k distinct subset sums) but NOT for intervals: lows of 20, 21, 22
      under an upper bound of 30 are all met by one path of weight 22, yet the
      term claims K >= 2. 36 wrong answers on their own.

    Together those made 7.7% of real components come back with more paths than
    the minimum -- while saving no measurable time (1,084 s vs 1,098 s for the
    plain ladder over the same set). The bound here never skips a feasible K:
    verified identical to the from-1 ladder on all 1,293 ground-truth networks.

    Returns 1 on any failure. A bound is an optimisation; it must never be the
    reason a component does not resolve.
    """
    try:
        graph = data["graph"]
        required = [
            (u, v) for (u, v) in data["edges"] if data["flows_low"].get((u, v), 0) > 0
        ]
        if graph.number_of_edges() == 0 or not required:
            # No edge has to carry flow. Guarded EXPLICITLY: flowpaths treats an
            # empty weight_function as "use the defaults", i.e. weight 1 on every
            # edge, which would silently turn this back into the invalid all-edges
            # width. Found on a real component whose five edges all had low = 0.
            return 1
        node_labels = {node: str(node) for node in graph.nodes}
        st = fp.stDAG(nx.relabel_nodes(graph, node_labels, copy=True))
        weight = {(node_labels[u], node_labels[v]): 1 for (u, v) in required}
        width = st.compute_max_edge_antichain(get_antichain=False, weight_function=weight)
        return max(1, int(width))
    except Exception as e:
        logger.debug(f"Could not compute a lower bound on K ({e}); starting from 1")
        return 1


def FD_Algorithm(data, max_paths, nthreads):
    """
    Find the smallest K in [1, max_paths] for which an inexact flow
    decomposition exists, and return it. The answer is identical to trying
    K = 1, 2, 3, ... in turn; the differences are in how much work is skipped.

    HOW REAL COMPONENTS BEHAVE (1,131 of the largest, re-solved locally):
      - 86% of resolvable components are feasible at the very first K tried;
        the rest need one or two more. Feasible solves are cheap (~0.03 s).
      - 58% of the largest components have NO decomposition within max_paths.
        For those the old ladder tried every K up to max_paths, and proving
        infeasibility gets exponentially slower with K: 0.5 s, 3 s, 30 s, 40 s,
        then hours. 91% of all solver time went into components that ended
        with no result. That is the 9-hour tail.

    SO:
      1. Valid lower bound above max_paths -> provably unresolvable, no solve.
      2. Try K = lower bound. Solved -> done (the common case, no extra cost).
      3. Otherwise ask ONE question before climbing: is K = max_paths feasible?
         Feasibility is monotone in K for this model -- a solution with K paths
         extends to K+1 by appending a zero-weight duplicate of its last path
         (the lexicographic symmetry breaking is non-strict, and the safety
         fixings never touch more than Kmin paths) -- so "no" at max_paths is a
         proof that every K in between is infeasible too. Confirmed on 780/780
         resolved real networks. No -> unresolvable, ladder skipped. Yes ->
         climb as before, and if the climb reaches max_paths the probe's own
         solution IS the answer (same model, same solver, same seed), so it is
         reused rather than solved twice. The smallest feasible K and the model
         at it are unchanged either way, so the returned decomposition is
         identical to the plain ladder's.

         NOT done with allow_empty_paths=True, although that is the textbook
         way to make the question monotone: on a real component that variant
         needed 450 s to prove K=10 infeasible where the standard model needed
         58 s -- empty paths add symmetry and loosen the relaxation.
      4. A per-attempt time limit (configure(time_limit=...)) turns a
         would-be-hours proof into a bounded wait. A timeout is NOT a proof of
         infeasibility, so the search stops there and reports the component
         unresolved rather than continue to a K it cannot call minimal. Off by
         default: with no limit every answer is exact.
    """
    n_edges = len(data["edges"])
    t0 = time.perf_counter()
    attempts = []            # (K, status, build s, solve s, oracle?)
    solution_set, solution_weights = 0, 0

    def attempt(target, K, oracle=False):
        # `oracle` only tags the attempt in the log ("K10*"); the model is the
        # standard one -- see the docstring for why not allow_empty_paths.
        target = flowMultipleDecomposition(target, K, nthreads)
        attempts.append((K, target.get("status", "?"), target.get("build_time", 0.0),
                         target.get("runtime", 0.0), oracle))
        return target

    def limit_text():
        # A timeout can also come from a solver-level limit set outside
        # configure(), in which case _MFD_TIME_LIMIT is still inf.
        return (f"the {_MFD_TIME_LIMIT:g}s time limit" if _MFD_TIME_LIMIT != float("inf")
                else "the solver time limit")

    def summary(outcome):
        detail = " ".join(
            f"K{K}{'*' if oracle else ''}:{str(st).replace('k', '', 1)}/{solve:.1f}s"
            for K, st, build, solve, oracle in attempts
        )
        logger.info(
            f"MFD [{n_edges} edges, lb={lowerbound}]: {outcome} after {len(attempts)} "
            f"attempt(s), {time.perf_counter() - t0:.1f}s total{(' -- ' + detail) if detail else ''}"
        )

    lowerbound = get_lowerbound_k(data)
    if lowerbound > max_paths:
        summary(f"unresolvable: lower bound {lowerbound} > --maxpaths {max_paths}, no solve needed")
        data["mfd_elapsed"] = time.perf_counter() - t0
        return data, {}

    data = attempt(data, lowerbound)
    if data["message"] == "solved":
        solution_set, solution_weights = data["solution"], data["weights"]
        summary(f"solved with K={lowerbound}")
    elif data.get("status") == "kTimeLimit":
        logger.warning(
            f"MFD [{n_edges} edges]: K={lowerbound} hit {limit_text()}; "
            f"cannot prove it infeasible, so the component is left unresolved"
        )
        summary("unresolved (time limit)")
    else:
        resolvable = True
        probe = None
        if lowerbound < max_paths:
            # Standard model on a shallow copy of the data dict: same constraints
            # the ladder would build at K = max_paths, so a feasible answer here
            # can stand in for that final rung.
            probe = attempt(dict(data), max_paths, oracle=True)
            if probe["message"] != "solved":
                resolvable = False
                if probe.get("status") == "kTimeLimit":
                    logger.warning(
                        f"MFD [{n_edges} edges]: feasibility check at K={max_paths} hit {limit_text()}; component left unresolved"
                    )
                    summary("unresolved (time limit at the K=max probe)")
                else:
                    summary(f"unresolvable within --maxpaths {max_paths}; "
                            f"K={lowerbound + 1}..{max_paths - 1} skipped")
        if resolvable:
            for K in range(lowerbound + 1, max_paths):
                data = attempt(data, K)
                if data["message"] == "solved":
                    solution_set, solution_weights = data["solution"], data["weights"]
                    summary(f"solved with K={K}")
                    break
                if data.get("status") == "kTimeLimit":
                    logger.warning(
                        f"MFD [{n_edges} edges]: K={K} hit {limit_text()}; "
                        f"component left unresolved"
                    )
                    summary("unresolved (time limit)")
                    break
            else:
                if probe is not None and probe["message"] == "solved":
                    # Every K below max_paths is now proven infeasible, so the
                    # minimum is max_paths -- and the probe already solved exactly
                    # that model.
                    solution_set, solution_weights = probe["solution"], probe["weights"]
                    summary(f"solved with K={max_paths} (reusing the feasibility probe)")
                else:
                    # lowerbound == max_paths: the single attempt above failed.
                    summary(f"unresolvable within --maxpaths {max_paths}")

    solution_paths = {}
    if solution_set != 0:
        for i in range(0, len(solution_set)):
            solution_paths[i] = {"weight": solution_weights[i], "path": solution_set[i]}

    data["mfd_elapsed"] = time.perf_counter() - t0
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

        # Optional: keep the exact instance of any component that took longer
        # than configure(dump_slow_s=...). This is what makes a pathological
        # component reproducible off the cluster -- the network as handed to the
        # solver, in the same format read_input()/SolveInstances consume.
        if _MFD_DUMP_SLOW_S is not None and data.get("mfd_elapsed", 0) > _MFD_DUMP_SLOW_S and _MFD_DUMP_DIR:
            try:
                os.makedirs(_MFD_DUMP_DIR, exist_ok=True)
                fname = os.path.join(
                    _MFD_DUMP_DIR,
                    f"mfd_{len(Edges)}edges_{data['mfd_elapsed']:.0f}s_pid{os.getpid()}_{int(time.time()*1000)}.json",
                )
                with open(fname, "w") as fh:
                    json.dump({"Nodes": Graphs[s]["Nodes"],
                               "list of edges": [list(e) for e in listOfEdges],
                               "subpaths": {str(k): list(v) for k, v in Graphs[s]["subpaths"].items()},
                               "elapsed_s": data["mfd_elapsed"], "paths": len(solution_paths)}, fh)
                logger.info(f"MFD: slow component ({data['mfd_elapsed']:.0f}s) dumped to {fname}")
            except Exception as e:
                logger.debug(f"could not dump slow MFD instance: {e}")

    return solution_paths
