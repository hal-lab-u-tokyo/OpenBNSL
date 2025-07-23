import pytest
import random
from pgmpy.utils import get_example_model
from pgmpy.estimators import PC

import openbnsllib

from modules.utils import to_pgmpy
from modules.structural_distance import structural_errors


@pytest.mark.parametrize(
    "model_name",
    [
        "cancer",  # 5 nodes
        "asia",  # 8 nodes
        # "child",   # 20 nodes
        # "alarm",   # 37 nodes
    ],
)
@pytest.mark.parametrize(
    "max_cond_vars",
    [
        5,
    ],
)
@pytest.mark.parametrize(
    "citest_type_str",
    [
        "chi_square",
        "g_sq",
    ],
)
@pytest.mark.parametrize(
    "level",
    [
        0.01,
        0.05,
    ],
)
@pytest.mark.parametrize(
    "sample_size",
    [
        int(1e5),
    ],
)
@pytest.mark.parametrize(
    "seed",
    [
        0,
        # 0, 1, 2, 3, 4, 5, 6, 7, 8, 9,
    ],
)
def test_rai(model_name, max_cond_vars, citest_type_str, level, sample_size, seed):
    random.seed(seed)

    model_original = get_example_model(model_name)
    samples = model_original.simulate(sample_size, seed=seed)
    samples = samples[sorted(samples.columns)]

    # pgmpy implementation
    estimator = PC(samples)
    expected_pgmpy = estimator.estimate(
        variant="orig",
        ci_test=citest_type_str,
        significance_level=level,
        max_cond_vars=max_cond_vars,
    )

    # our implementation
    if citest_type_str == "chi_square":
        citest_type = openbnsllib.citest.ChiSquare(level)
    elif citest_type_str == "g_sq":
        citest_type = openbnsllib.citest.GSquare(level)
    else:
        raise ValueError(f"Unsupported citest_type_str: {citest_type_str}")
    df_wrapper = openbnsllib.base.DataframeWrapper(samples)
    _pdag = openbnsllib.structure_learning.RAI(
        df_wrapper,
        citest_type,
        max_cond_vars,
    )

    expected_obnsl = to_pgmpy(_pdag, list(samples.columns))
    error_pgmpy = structural_errors(expected_pgmpy, model_original)
    error_obnsl = structural_errors(expected_obnsl, model_original)
    print(f"[pgmpy] Model: {model_name}, Seed: {seed}, {error_pgmpy}")
    print(f"[obnsl] Model: {model_name}, Seed: {seed}, {error_obnsl}")
