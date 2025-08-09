import sys
import pytest
import numpy as np
from pgmpy.base import PDAG
from pgmpy.utils import get_example_model

import openbnsllib
from modules.utils import to_pgmpy
from modules.structural_distance import structural_errors


# @pytest.mark.parametrize("model_name", ["cancer", "asia", "child", "alarm"])
@pytest.mark.parametrize("model_name", ["cancer"])
@pytest.mark.parametrize("score_type", [openbnsllib.score.BDeu(1.0)])
@pytest.mark.parametrize("sample_size", [int(1e5)])
@pytest.mark.parametrize("seed", [0])
def test_simulated_annealing(model_name, score_type, sample_size, seed):
    model_original = get_example_model(model_name)
    samples = model_original.simulate(sample_size, seed=seed)
    samples = samples[sorted(samples.columns)]
    df_wrapper = openbnsllib.base.DataframeWrapper(samples)
    _pdag = openbnsllib.structure_learning.simulated_annealing(
        df_wrapper,
        score_type,
        max_parents=3,
        max_iters=1000,
        init_temp=1.0,
        cooling_rate=0.9995,
        is_deterministic=True,
        seed=seed,
        num_chains=1,
    )
    model_estimated = to_pgmpy(_pdag, list(samples.columns))
    errors = structural_errors(model_original, model_estimated)

    # TODO: Implementation
    # if errors["SHD"] > 0:
    #     print(f"Model: {model_name}, [SHD, ME, EE, DE, ED, MD, RD]: {errors}")
    # assert errors["SHD"] == 0, f"Failed for model {model_name} with Errors: {errors}"
