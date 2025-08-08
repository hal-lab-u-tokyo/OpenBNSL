import numpy as np
import networkx as nx
from pgmpy.base import PDAG
from itertools import combinations


def _retrieve_adjacency_matrix(
    graph: PDAG, order_nodes=None, weight=None
) -> np.ndarray:
    """
    Retrieve the adjacency matrix of a PDAG.
    Parameters
    ----------
    graph : PDAG
        The PDAG object.
    order_nodes : list, optional
        The order of nodes in the adjacency matrix. If None, the order of original graph is used.
    weight : str, optional
        The edge attribute to use as weight. If None, the adjacency matrix is unweighted.
    Returns
    -------
    np.ndarray
        The adjacency matrix of the PDAG.
    """

    if not isinstance(graph, PDAG):
        raise TypeError("graph must be a PDAG object")
    if order_nodes is None:
        order_nodes = graph.nodes()
    return np.array(nx.adjacency_matrix(graph, order_nodes, weight=weight).todense())


def structural_hamming_distance(ground_truth_pdag: PDAG, predicted_graph: PDAG) -> int:
    """
    Compute the structural Hamming distance between two PDAGs.
    Parameters
    ----------
    ground_truth_pdag : PDAG
        The ground_truth_pdag PDAG (ground truth).
    predicted_graph : PDAG
        The predicted PDAG.
    Returns
    -------
    int
        The structural Hamming distance between the ground_truth_pdag and predicted PDAGs.
    """

    true_labels = _retrieve_adjacency_matrix(ground_truth_pdag)
    predictions = _retrieve_adjacency_matrix(
        predicted_graph,
        (
            ground_truth_pdag.nodes()
            if isinstance(ground_truth_pdag, nx.DiGraph)
            else None
        ),
    )
    diff = np.abs(true_labels - predictions)
    diff = diff + diff.transpose()
    diff[diff > 1] = 1  # Ignoring the double edges.
    return int(np.sum(diff) / 2)


def missing_edge(ground_truth_pdag: PDAG, predicted_graph: PDAG) -> int:
    """
    Compute the missing edge between two PDAGs.
    Parameters
    ----------
    ground_truth_pdag : PDAG
        The ground_truth_pdag PDAG (ground truth).
    predicted_graph : PDAG
        The predicted PDAG.
    Returns
    -------
    int
        The missing edge error between the ground_truth_pdag and predicted PDAGs.
    """
    true_labels = _retrieve_adjacency_matrix(ground_truth_pdag)
    predictions = _retrieve_adjacency_matrix(
        predicted_graph,
        (
            ground_truth_pdag.nodes()
            if isinstance(ground_truth_pdag, nx.DiGraph)
            else None
        ),
    )
    diff = true_labels - predictions - predictions.transpose()
    diff[diff < 0] = 0  # Ignoring other types of error.
    diff = diff + diff.transpose()
    diff[diff > 1] = 1  # Ignoring the double edges.
    return int(np.sum(diff) / 2)


def extra_edge(ground_truth_pdag: PDAG, predicted_graph: PDAG) -> int:
    """
    Compute the extra edge between two PDAGs.
    Parameters
    ----------
    ground_truth_pdag : PDAG
        The ground_truth_pdag PDAG (ground truth).
    predicted_graph : PDAG
        The predicted PDAG.
    Returns
    -------
    int
        The extra edge error between the ground_truth_pdag and predicted PDAGs.
    """

    true_labels = _retrieve_adjacency_matrix(ground_truth_pdag)
    predictions = _retrieve_adjacency_matrix(
        predicted_graph,
        (
            ground_truth_pdag.nodes()
            if isinstance(ground_truth_pdag, nx.DiGraph)
            else None
        ),
    )
    diff = predictions - true_labels - true_labels.transpose()
    diff[diff < 0] = 0  # Ignoring other types of error.
    diff = diff + diff.transpose()
    diff[diff > 1] = 1  # Ignoring the double edges.
    return int(np.sum(diff) / 2)


def extra_direction(ground_truth_pdag: PDAG, predicted_graph: PDAG) -> int:
    """
    Compute the extra direction between two PDAGs.
    Parameters
    ----------
    ground_truth_pdag : PDAG
        The ground_truth_pdag PDAG (ground truth).
    predicted_graph : PDAG
        The predicted PDAG.
    Returns
    -------
    int
        The extra direction error between the ground_truth_pdag and predicted PDAGs.
    """

    true_labels = _retrieve_adjacency_matrix(ground_truth_pdag)
    predictions = _retrieve_adjacency_matrix(
        predicted_graph,
        (
            ground_truth_pdag.nodes()
            if isinstance(ground_truth_pdag, nx.DiGraph)
            else None
        ),
    )
    count = 0
    for x in range(true_labels.shape[0]):
        for y in range(true_labels.shape[1]):
            if (
                true_labels[x, y] == 1
                and true_labels[y, x] == 1
                and predictions[x, y] == 1
                and predictions[y, x] == 0
            ):
                count += 1
    return count


def missing_direction(ground_truth_pdag: PDAG, predicted_graph: PDAG) -> int:
    """
    Compute the missing direction between two PDAGs.
    Parameters
    ----------
    ground_truth_pdag : PDAG
        The ground_truth_pdag PDAG (ground truth).
    predicted_graph : PDAG
        The predicted PDAG.
    Returns
    -------
    int
        The missing direction error between the ground_truth_pdag and predicted PDAGs.
    """

    true_labels = _retrieve_adjacency_matrix(ground_truth_pdag)
    predictions = _retrieve_adjacency_matrix(
        predicted_graph,
        (
            ground_truth_pdag.nodes()
            if isinstance(ground_truth_pdag, nx.DiGraph)
            else None
        ),
    )
    count = 0
    for x in range(true_labels.shape[0]):
        for y in range(true_labels.shape[1]):
            if (
                true_labels[x, y] == 1
                and true_labels[y, x] == 0
                and predictions[x, y] == 1
                and predictions[y, x] == 1
            ):
                count += 1
    return count


def reversed_direction(ground_truth_pdag: PDAG, predicted_graph: PDAG) -> int:
    """
    Compute the reversed direction error between two PDAGs.
    Parameters
    ----------
    ground_truth_pdag : PDAG
        The ground_truth_pdag PDAG (ground truth).
    predicted_graph : PDAG
        The predicted PDAG.
    Returns
    -------
    int
        The reversed direction error between the ground_truth_pdag and predicted PDAGs.
    """

    true_labels = _retrieve_adjacency_matrix(ground_truth_pdag)
    predictions = _retrieve_adjacency_matrix(
        predicted_graph,
        (
            ground_truth_pdag.nodes()
            if isinstance(ground_truth_pdag, nx.DiGraph)
            else None
        ),
    )
    count = 0
    for x in range(true_labels.shape[0]):
        for y in range(true_labels.shape[1]):
            if (
                true_labels[x, y] == 1
                and true_labels[y, x] == 0
                and predictions[x, y] == 0
                and predictions[y, x] == 1
            ):
                count += 1
    return count


def directional_error(ground_truth_pdag: PDAG, predicted_graph: PDAG) -> int:
    """
    Compute the directional error between two PDAGs.
    Parameters
    ----------
    ground_truth_pdag : PDAG
        The ground_truth_pdag PDAG (ground truth).
    predicted_graph : PDAG
        The predicted PDAG.
    Returns
    -------
    int
        The directional error between the ground_truth_pdag and predicted PDAGs.
    """

    return (
        structural_hamming_distance(ground_truth_pdag, predicted_graph)
        - missing_edge(ground_truth_pdag, predicted_graph)
        - extra_edge(ground_truth_pdag, predicted_graph)
    )


# TODO remove this function
def _structural_errors(ground_truth_pdag, predicted_graph):
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
    ground_truth_pdag : nx.DiGraph
        true CPDAG
    predicted_graph : nx.DiGraph
        predicted CPDAG
    Returns
    ------
    errors : list
    [SHD, missing_edge, extra_edge, directional_error, extra_direction, missing_direction, reversed_direction]
    """
    errors = []
    errors.append(
        structural_hamming_distance(ground_truth_pdag, predicted_graph)
    )  # SHD = EE + ME + DE
    errors.append(missing_edge(ground_truth_pdag, predicted_graph))
    errors.append(extra_edge(ground_truth_pdag, predicted_graph))
    errors.append(
        directional_error(ground_truth_pdag, predicted_graph)
    )  # DE = ED + MD + RD
    errors.append(extra_direction(ground_truth_pdag, predicted_graph))
    errors.append(missing_direction(ground_truth_pdag, predicted_graph))
    errors.append(reversed_direction(ground_truth_pdag, predicted_graph))
    return errors


def structural_errors(ground_truth_pdag: PDAG, predicted_graph: PDAG) -> dict:
    """
    Compute structural errors between two PDAGs.
    Parameters
    ----------
    ground_truth_pdag : PDAG
        The ground_truth_pdag PDAG (ground truth).
    predicted_graph : PDAG
        The predicted PDAG.

    Returns
    -------
    errors : dict
        A dictionary containing the structural errors:
        - SHD: Structural Hamming Distance (= EE + ME + DE)
        - ME: Missing Edge
        - EE: Extra Edge
        - DE: Directional Error (= ED + MD + RD)
        - ED: Extra Direction
        - MD: Missing Direction
        - RD: Reversed Direction
    """
    # Convert PDAG to CPDAG
    ground_truth_pdag = PDAG2CPDAG(ground_truth_pdag)
    predicted_graph = PDAG2CPDAG(predicted_graph)

    errors = {}
    errors["SHD"] = structural_hamming_distance(ground_truth_pdag, predicted_graph)
    errors["ME"] = missing_edge(ground_truth_pdag, predicted_graph)
    errors["EE"] = extra_edge(ground_truth_pdag, predicted_graph)
    errors["DE"] = directional_error(ground_truth_pdag, predicted_graph)
    errors["ED"] = extra_direction(ground_truth_pdag, predicted_graph)
    errors["MD"] = missing_direction(ground_truth_pdag, predicted_graph)
    errors["RD"] = reversed_direction(ground_truth_pdag, predicted_graph)
    return errors


def PDAG2CPDAG(pdag: PDAG) -> PDAG:
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
    cpdag.add_edges_from(pdag.edges)
    reversedag = pdag.reverse(copy=True)
    cpdag.add_edges_from(reversedag.edges)

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
            for Y in pdag_removed.predecessors(
                X
            ):  # for every parent of V-structure fix it
                for Z in pdag_removed.predecessors(X):
                    if Y != Z and not pdag.has_edge(Y, Z) and not pdag.has_edge(Z, Y):
                        if cpdag.has_edge(Y, X) and cpdag.has_edge(X, Y):
                            cpdag.remove_edge(X, Y)
                        if cpdag.has_edge(Z, X) and cpdag.has_edge(X, Z):
                            cpdag.remove_edge(X, Z)
    return cpdag
