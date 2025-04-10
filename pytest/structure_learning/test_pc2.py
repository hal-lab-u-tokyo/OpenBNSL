import sys
import pytest
import numpy as np
from pgmpy.base import PDAG
from pgmpy.utils import get_example_model

import openbnsllib

from modules.structural_distance import structural_errors, PDAG2CPDAG
from modules.util import to_pgmpy


@pytest.mark.parametrize(
    "model_name, citest_type, sample_size, seed",
    [
        ("cancer", openbnsllib.citest.ChiSquare(0.05), int(1e5), 0),  # 5 nodes
        ("asia", openbnsllib.citest.ChiSquare(0.05), int(1e5), 0),  # 8 nodes
        ("child", openbnsllib.citest.ChiSquare(0.05), int(1e5), 0),  # 20 nodes
    ],
)
def test_pc2(model_name, citest_type, sample_size, seed):

    model_original = get_example_model(model_name)
    samples = model_original.simulate(sample_size, seed=seed)
    samples = samples[sorted(samples.columns)]

    df_wrapper = openbnsllib.base.DataframeWrapper(samples)
    _pdag = openbnsllib.structure_learning.PC2(df_wrapper, citest_type)

    # Convert _pdag to pgmpy PDAG
    model_estimated = PDAG()
    node_labels = list(samples.columns)
    model_estimated.add_nodes_from(node_labels)
    for i in range(len(node_labels)):
        for j in range(len(node_labels)):
            if _pdag.has_edge(i, j):
                model_estimated.add_edge(node_labels[i], node_labels[j])

    model_original = PDAG2CPDAG(model_original)
    model_estimated = PDAG2CPDAG(model_estimated)
    errors = structural_errors(model_original, model_estimated)

    if errors[0] > 0:
        print(f"Model: {model_name}, [SHD, ME, EE, DE, ED, MD, RD]: {errors}")
    assert errors[0] == 0, f"SHD too high for {model_name}: {errors[0]}"
