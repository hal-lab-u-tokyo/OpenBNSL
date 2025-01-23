import pandas as pd
import numpy as np
import networkx as nx
from pgmpy.base import DAG
from pgmpy.base import PDAG
from itertools import combinations
import argparse
import os
from modules.visualize_graph import display_graph_info as show

def retrieve_adjacency_matrix(graph, order_nodes=None, weight=False):
    """Retrieve the adjacency matrix from the nx.DiGraph or numpy array."""
    if isinstance(graph, np.ndarray):
        return graph
    elif isinstance(graph, nx.DiGraph):
        if order_nodes is None:
            order_nodes = graph.nodes()
        if not weight:
            return np.array(nx.adjacency_matrix(graph, order_nodes, weight=None).todense())
        else:
            return np.array(nx.adjacency_matrix(graph, order_nodes).todense())
    else:
        raise TypeError("Only networkx.DiGraph and np.ndarray (adjacency matrixes) are supported.")


def SHD(target, pred):
    true_labels = retrieve_adjacency_matrix(target)
    predictions = retrieve_adjacency_matrix(pred, target.nodes() 
                                            if isinstance(target, nx.DiGraph) else None)
    diff = np.abs(true_labels - predictions)
    diff = diff + diff.transpose()
    diff[diff > 1] = 1  # Ignoring the double edges.
    return int(np.sum(diff)/2)

def missing_edge(target, pred):
    true_labels = retrieve_adjacency_matrix(target)
    predictions = retrieve_adjacency_matrix(pred, target.nodes() 
    if isinstance(target, nx.DiGraph) else None)
    diff = true_labels - predictions - predictions.transpose()
    diff[diff < 0] = 0  # Ignoring other types of error.
    diff = diff + diff.transpose()
    diff[diff > 1] = 1  # Ignoring the double edges.
    return int(np.sum(diff)/2)

def extra_edge(target, pred):
    true_labels = retrieve_adjacency_matrix(target)
    predictions = retrieve_adjacency_matrix(pred, target.nodes() 
    if isinstance(target, nx.DiGraph) else None)
    diff = predictions - true_labels - true_labels.transpose()
    diff[diff < 0] = 0  # Ignoring other types of error.
    diff = diff + diff.transpose()
    diff[diff > 1] = 1  # Ignoring the double edges.
    return int(np.sum(diff)/2)

def extra_direction(target, pred):
    true_labels = retrieve_adjacency_matrix(target)
    predictions = retrieve_adjacency_matrix(pred, target.nodes() 
    if isinstance(target, nx.DiGraph) else None)
    count = 0
    for x in range(true_labels.shape[0]):
        for y in range(true_labels.shape[1]):
            if true_labels[x, y] == 1 and true_labels[y, x] == 1 and predictions[x, y] == 1 and predictions[y, x] == 0:
                count += 1
    return count

def missing_direction(target, pred):
    true_labels = retrieve_adjacency_matrix(target)
    predictions = retrieve_adjacency_matrix(pred, target.nodes() 
    if isinstance(target, nx.DiGraph) else None)
    count = 0
    for x in range(true_labels.shape[0]):
        for y in range(true_labels.shape[1]):
            if true_labels[x, y] == 1 and true_labels[y, x] == 0 and predictions[x, y] == 1 and predictions[y, x] == 1:
                count += 1
    return count

def reversed_direction(target, pred):
    true_labels = retrieve_adjacency_matrix(target)
    predictions = retrieve_adjacency_matrix(pred, target.nodes() 
    if isinstance(target, nx.DiGraph) else None)
    count = 0
    for x in range(true_labels.shape[0]):
        for y in range(true_labels.shape[1]):
            if true_labels[x, y] == 1 and true_labels[y, x] == 0 and predictions[x, y] == 0 and predictions[y, x] == 1:
                count += 1
    return count

def directional_error(target, pred):
    return SHD(target, pred) - missing_edge(target, pred) - extra_edge(target, pred)

def structural_errors(target, pred):
    """
    Conpute Structural Hamming Distance (SHD), Extra Edge (EE), Missing Edge (ME), Directional Error (DE),ED (Extra Direction), MD (Missing Direction), RD(Reversed Direction)
    SHD = EE + ME + DE
    DE = ED + MD + RD
    true    predicted
    no      -           SHD = 1  EE = 1
    no      ->          SHD = 1  EE = 1
    ->      <-          SHD = 1  DE = 1 (Reversed Direction)
    ->      no          SHD = 1  ME = 1
    ->      -           SHD = 1  DE = 1 (Missing Direction)
    -       no          SHD = 1  ME = 1
    -       ->          SHD = 1  DE = 1 (Extra Direction)

    Parameters
    ----------
    target : nx.DiGraph
        true CPDAG
    pred : nx.DiGraph
        predicted CPDAG
    Returns
    ------
    errors : list
    [SHD, missing_edge, extra_edge, directional_error, extra_direction, missing_direction, reversed_direction]
    """
    errors = []
    errors.append(SHD(target, pred))
    errors.append(missing_edge(target, pred))
    errors.append(extra_edge(target, pred))
    errors.append(directional_error(target, pred))
    errors.append(extra_direction(target, pred))
    errors.append(missing_direction(target, pred))
    errors.append(reversed_direction(target, pred))
    return errors


def PDAG2CPDAG(pdag): 
    """
    Compute the completed partially directed acyclic graph (CPDAG) of a given PDAG.
    
    Parameters
    ----------
    pdag : nx.DiGraph
        input PDAG

    Returns
    -------
    cpdag : nx.DiGraph
        output CPDAG
    """
    cpdag = PDAG()
    cpdag.add_nodes_from(pdag.nodes)
    
    #make skeleton
    cpdag.add_edges_from(pdag.edges)
    reversedag = pdag.reverse(copy = True)
    cpdag.add_edges_from(reversedag.edges)

    # vstructuredag = PDAG()
    # vstructuredag.add_nodes_from(pdag.nodes)

    #pdag_removed = pdag.copy() #remove all undirected edges
    pdag_removed = PDAG()
    pdag_removed.add_nodes_from(pdag.nodes)
    pdag_removed.add_edges_from(pdag.edges)
    for pair in list(combinations(pdag.nodes(), 2)):
            X, Y = pair
            if pdag_removed.has_edge(X, Y) and pdag_removed.has_edge(Y, X): 
                pdag_removed.remove_edge(Y, X)
                pdag_removed.remove_edge(X, Y)
    for X in pdag_removed.nodes:
        if pdag_removed.in_degree(X) > 1:
            for Y in pdag_removed.predecessors(X): #for every parent of V-structure fix it
                for Z in pdag_removed.predecessors(X):
                    if Y != Z and not pdag.has_edge(Y, Z) and not pdag.has_edge(Z, Y):
                        if cpdag.has_edge(Y, X) and cpdag.has_edge(X, Y):
                            cpdag.remove_edge(X, Y)
                        if cpdag.has_edge(Z, X) and cpdag.has_edge(X, Z):
                            cpdag.remove_edge(X, Z)
                #cpdag.remove_edge(X, Y)


    #             vstructuredag.add_edge(Y, X)
    
    # for X in pdag.nodes:
    #     for Y in pdag.nodes:
    #         if cpdag.has_edge(X, Y) and cpdag.has_edge(Y, X):
    #             if vstructuredag.in_degree(X) > 0 and vstructuredag.in_degree(Y) == 0:
    #                 vstructuredag.add_edge(X, Y)
    #                 cpdag.remove_edge(Y, X)
    #             elif vstructuredag.in_degree(Y) > 0 and vstructuredag.in_degree(X) == 0:
    #                 vstructuredag.add_edge(Y, X)
    #                 cpdag.remove_edge(X, Y)

    return cpdag

def DAG2CPDAG(dag): 
    """
    Compute the completed partially directed acyclic graph (CPDAG) of a given DAG.
    
    Parameters
    ----------
    dag : nx.DiGraph
        input DAG

    Returns
    -------
    cpdag : nx.DiGraph
        output CPDAG
    """
    cpdag = PDAG()
    cpdag.add_nodes_from(dag.nodes)
    
    #make skeleton
    cpdag.add_edges_from(dag.edges)
    reversedag = dag.reverse(copy = True)
    cpdag.add_edges_from(reversedag.edges)

    # vstructuredag = PDAG()
    # vstructuredag.add_nodes_from(dag.nodes)

    for X in dag.nodes:
        if dag.in_degree(X) > 1:
            for Y in dag.predecessors(X): #for every parent of V-structure fix it
                for Z in dag.predecessors(X):
                    if Y != Z and not dag.has_edge(Y, Z) and not dag.has_edge(Z, Y):
                        if cpdag.has_edge(Y, X) and cpdag.has_edge(X, Y):
                            cpdag.remove_edge(X, Y)
                        if cpdag.has_edge(Z, X) and cpdag.has_edge(X, Z):
                            cpdag.remove_edge(X, Z)
                # vstructuredag.add_edge(Y, X)
    
    # for X in dag.nodes:
    #     for Y in dag.nodes:
    #         if cpdag.has_edge(X, Y) and cpdag.has_edge(Y, X):
    #             if vstructuredag.in_degree(X) > 0 and vstructuredag.in_degree(Y) == 0:
    #                 vstructuredag.add_edge(X, Y)
    #                 cpdag.remove_edge(Y, X)
    #             elif vstructuredag.in_degree(Y) > 0 and vstructuredag.in_degree(X) == 0:
    #                 vstructuredag.add_edge(Y, X)
    #                 cpdag.remove_edge(X, Y)

    return cpdag