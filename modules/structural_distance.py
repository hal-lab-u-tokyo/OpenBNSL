# modules/structural_distance.py

import numpy as np
import networkx as nx
from pgmpy.base import PDAG
from itertools import combinations


def _assert_same_nodes(gt: PDAG, pred: PDAG):
    """Assert that the ground truth and predicted PDAGs have the same nodes."""
    gt_nodes = set(gt.nodes())
    pr_nodes = set(pred.nodes())
    if gt_nodes != pr_nodes:
        missing = gt_nodes - pr_nodes
        extra = pr_nodes - gt_nodes
        raise ValueError(
            f"Node sets differ. missing_in_pred={missing}, extra_in_pred={extra}"
        )


def _split_edge_sets(G: PDAG):
    """Split edges of a PDAG into directed and undirected sets."""
    dir_edges = set()
    undir_edges = set()
    for u, v in G.edges():
        if G.has_edge(v, u):
            if u != v:
                undir_edges.add(frozenset((u, v)))
        else:
            if u != v:
                dir_edges.add((u, v))
    return dir_edges, undir_edges


def _skeleton_from(dir_edges, undir_edges):
    """Create a skeleton graph from directed and undirected edge sets."""
    skel = set(undir_edges)
    for u, v in dir_edges:
        skel.add(frozenset((u, v)))
    return skel


def PDAG2CPDAG(pdag: PDAG) -> PDAG:
    """
    Compute the completed partially directed acyclic graph (CPDAG) of a given PDAG.
    Args:
        pdag (PDAG): The input PDAG to convert.
    Returns:
        PDAG: The completed PDAG.
    """
    cpdag = PDAG()
    cpdag.add_nodes_from(pdag.nodes)
    cpdag.add_edges_from(pdag.edges)
    reversedag = pdag.reverse(copy=True)
    cpdag.add_edges_from(reversedag.edges)

    pdag_removed = PDAG()
    pdag_removed.add_nodes_from(pdag.nodes)
    pdag_removed.add_edges_from(pdag.edges)
    for X, Y in combinations(pdag.nodes(), 2):
        if pdag_removed.has_edge(X, Y) and pdag_removed.has_edge(Y, X):
            pdag_removed.remove_edge(Y, X)
            pdag_removed.remove_edge(X, Y)

    for X in pdag_removed.nodes:
        if pdag_removed.in_degree(X) > 1:
            parents = list(pdag_removed.predecessors(X))
            for i in range(len(parents)):
                for j in range(i + 1, len(parents)):
                    Y, Z = parents[i], parents[j]
                    if (not pdag.has_edge(Y, Z)) and (not pdag.has_edge(Z, Y)):
                        if cpdag.has_edge(X, Y) and cpdag.has_edge(Y, X):
                            cpdag.remove_edge(X, Y)
                        if cpdag.has_edge(X, Z) and cpdag.has_edge(Z, X):
                            cpdag.remove_edge(X, Z)
    return cpdag


def structural_errors(ground_truth_pdag: PDAG, predicted_graph: PDAG) -> dict:
    """
    Compute structural errors between two PDAGs.
    Args:
        ground_truth_pdag (PDAG): The ground truth PDAG.
        predicted_graph (PDAG): The predicted PDAG.
    Returns:
        dict: A dictionary containing the structural errors:
        - SHD: Structural Hamming Distance (= EE + ME + DE)
        - ME: Missing Edge          (ground truth has edge, prediction does not)
        - EE: Extra Edge            (prediction has edge, ground truth does not)
        - DE: Directional Error     (= ED + MD + RD)
        - ED: Extra Direction       (undirected in ground truth, directed in prediction)
        - MD: Missing Direction     (directed in ground truth, undirected in prediction)
        - RD: Reversed Direction    (X -> Y in ground truth, Y -> X in prediction)
    """
    gt = PDAG2CPDAG(ground_truth_pdag)
    pr = PDAG2CPDAG(predicted_graph)

    _assert_same_nodes(gt, pr)

    dir_true, undir_true = _split_edge_sets(gt)
    dir_pred, undir_pred = _split_edge_sets(pr)

    skel_true = _skeleton_from(dir_true, undir_true)
    skel_pred = _skeleton_from(dir_pred, undir_pred)

    ME = len(skel_true - skel_pred)
    EE = len(skel_pred - skel_true)

    ED = 0
    for e in undir_true:
        u, v = tuple(e)
        ed_flag = ((u, v) in dir_pred) ^ ((v, u) in dir_pred)  # 片方だけ
        if ed_flag:
            ED += 1

    MD = 0
    for u, v in dir_true:
        if frozenset((u, v)) in undir_pred:
            MD += 1

    RD = 0
    for u, v in dir_true:
        if (v, u) in dir_pred:
            RD += 1

    DE = ED + MD + RD
    SHD = ME + EE + DE

    return {"SHD": SHD, "ME": ME, "EE": EE, "DE": DE, "ED": ED, "MD": MD, "RD": RD}
