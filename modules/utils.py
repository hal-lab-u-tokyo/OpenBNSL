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
