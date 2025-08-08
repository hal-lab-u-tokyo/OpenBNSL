from pgmpy.base import PDAG as PgmpyPDAG
from openbnsllib.base import PDAG as OpenBNSLPDAG


def to_pgmpy(openbnsl_pdag: OpenBNSLPDAG, node_labels: list[str]) -> PgmpyPDAG:
    """
    Convert an OpenBNSL PDAG to a pgmpy PDAG.

    Args:
        openbnsl_pdag (OpenBNSLPDAG): The OpenBNSL PDAG to convert.

    Returns:
        PgmpyPDAG: The converted pgmpy PDAG.
    """

    num_vars = openbnsl_pdag.num_vars
    if len(node_labels) != num_vars:
        raise ValueError(
            f"Number of node labels ({len(node_labels)}) does not match the number of variables ({num_vars})."
        )

    pgmpy_pdag = PgmpyPDAG()
    pgmpy_pdag.add_nodes_from(node_labels)
    for i in range(num_vars):
        for j in range(num_vars):
            if openbnsl_pdag.has_edge(i, j):
                pgmpy_pdag.add_edge(node_labels[i], node_labels[j])
    return pgmpy_pdag


def to_openbnsl(pgmpy_pdag: PgmpyPDAG, var_str2idx: dict[str, int]) -> OpenBNSLPDAG:
    """
    Convert a pgmpy PDAG to an OpenBNSL PDAG.

    Args:
        pgmpy_pdag (PgmpyPDAG): The pgmpy PDAG to convert.

    Returns:
        OpenBNSLPDAG: The converted OpenBNSL PDAG.
    """
    n = pgmpy_pdag.number_of_nodes()
    if len(var_str2idx) != n:
        raise ValueError(
            f"Number of variable labels ({len(var_str2idx)}) does not match the number of nodes ({n})."
        )
    obnsl_pdag = OpenBNSLPDAG(n)
    for u, v in pgmpy_pdag.edges():
        u_idx = var_str2idx[u]
        v_idx = var_str2idx[v]
        obnsl_pdag.add_edge(u_idx, v_idx)
    return obnsl_pdag
