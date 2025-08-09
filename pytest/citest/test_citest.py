from typing import Literal
import pytest
import random
import numpy as np
from pgmpy.utils import get_example_model
from pgmpy.estimators.CITests import chi_square, g_sq

import openbnsllib

@pytest.mark.parametrize("model_name", ["cancer", "asia", "child", "alarm"])
@pytest.mark.parametrize("sample_size", [int(1e5)])
@pytest.mark.parametrize("citest_type_str", ["chi_square", "g_sq"])
@pytest.mark.parametrize("level", [0.01])
@pytest.mark.parametrize("seed", [0])
def test_citest(
    model_name: Literal["asia"],
    citest_type_str: Literal["chi_square"] | Literal["g_sq"],
    level: float,
    sample_size: int,
    seed: int,
):
    random.seed(seed)

    model_original = get_example_model(model_name)
    samples = model_original.simulate(sample_size, seed=seed)
    samples = samples[sorted(samples.columns)]

    df_wrapper = openbnsllib.base.DataframeWrapper(samples)
    n = df_wrapper.num_of_vars

    remaining_indices = list(range(n))
    x_idx = random.choice(remaining_indices)
    remaining_indices.remove(x_idx)
    y_idx = random.choice(remaining_indices)
    remaining_indices.remove(y_idx)
    num_sepset_candidates = random.randint(0, min(5, len(remaining_indices)))
    selected_sepset_indices = random.sample(remaining_indices, num_sepset_candidates)
    sorted_sepset_indices = sorted(selected_sepset_indices)

    # pgmpy implementation
    if citest_type_str == "chi_square":
        pgmpy_citest_func = chi_square
    elif citest_type_str == "g_sq":
        pgmpy_citest_func = g_sq
    else:
        raise ValueError(f"Unsupported citest_type_str: {citest_type_str}")
    x_name = samples.columns[x_idx]
    y_name = samples.columns[y_idx]
    z_names = [samples.columns[i] for i in sorted_sepset_indices]
    (chi2, p_value, dof) = pgmpy_citest_func(
        X=x_name,
        Y=y_name,
        Z=z_names,
        data=samples,
        boolean=False,
        significance_level=level,
    )
    if np.isnan(p_value):
        pytest.skip("pgmpy produced NaN p‑value (dof == 0); test not informative.")

    expected_result = p_value >= level
    # print(
    #     f"[pgmpy] stat: {chi2}, "
    #     f"p_value: {p_value}, "
    #     f"dof: {dof}, "
    #     f"level: {level}, "
    #     f"result: {expected_result}"
    # )

    # our implementation
    if citest_type_str == "chi_square":
        citest_type = openbnsllib.citest.ChiSquare(level)
    elif citest_type_str == "g_sq":
        citest_type = openbnsllib.citest.GSquare(level)
    else:
        raise ValueError(f"Unsupported citest_type_str: {citest_type_str}")
    var_indices = sorted(sorted_sepset_indices + [x_idx, y_idx])
    ct = openbnsllib.base.ContingencyTable(var_indices, df_wrapper)
    computed_result = openbnsllib.citest.citest(
        x_idx, y_idx, sorted_sepset_indices, ct, citest_type
    )

    # Check the result
    msg = (
        f"CItest for {x_idx}, {y_idx} | {sorted_sepset_indices} of {model_name}: "
        f"pgmpy: {expected_result}, ours: {computed_result}"
    )
    assert expected_result == computed_result, msg
