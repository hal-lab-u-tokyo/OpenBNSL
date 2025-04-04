import sys
import pytest
import numpy as np
from pgmpy.base import PDAG
from pgmpy.utils import get_example_model

sys.path.append("/workspace/build")
import openbnsllib

sys.path.append("/workspace")
from modules.structural_distance import structural_errors, PDAG2CPDAG


@pytest.mark.parametrize(
    "model_name, score_type, sample_size, seed",
    [
        ("cancer", openbnsllib.score.BDeu(1.0), int(1e5), 0),  # 5 nodes
        ("asia", openbnsllib.score.BDeu(1.0), int(1e5), 0),  # 8 nodes
        ("child", openbnsllib.score.BDeu(1.0), int(1e5), 0),  # 20 nodes
    ],
)
def test_exhaustive_search(model_name, score_type, sample_size, seed):

    model_original = get_example_model(model_name)
    samples = model_original.simulate(sample_size, seed=seed)
    samples = samples[sorted(samples.columns)]

    df_wrapper = openbnsllib.base.DataframeWrapper(samples)
    ndarray = openbnsllib.exhaustive_search(df_wrapper, score_type, max_parents=3)

    # Convert ndarray to PDAG
    model_estimated = PDAG()
    model_estimated.add_nodes_from(range(samples.shape[1]))
    for i in range(ndarray.shape[0]):
        for j in range(ndarray.shape[1]):
            if ndarray[i, j] > 0:
                model_estimated.add_edge(samples.columns[i], samples.columns[j])

    model_original = PDAG2CPDAG(model_original)
    model_estimated = PDAG2CPDAG(model_estimated)
    errors = structural_errors(model_original, model_estimated)

    if errors[0] > 0:
        print(f"Model: {model_name}, [SHD, ME, EE, DE, ED, MD, RD]: {errors}")
    assert errors[0] == 0, f"SHD too high for {model_name}: {errors[0]}"
