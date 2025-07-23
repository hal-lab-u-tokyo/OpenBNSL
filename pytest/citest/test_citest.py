from typing import Literal
import pytest
import random
from pgmpy.utils import get_example_model
from pgmpy.estimators.CITests import chi_square, g_sq

import openbnsllib


@pytest.mark.parametrize(
    "model_name",
    [
        "cancer",  # e.g. 5 nodes
        "asia",  # e.g. 8 nodes
        # "child",   # e.g. 20 nodes
        # "alarm",   # e.g. 37 nodes
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
        1,
        2,
        # 0, 1, 2, 3, 4, 5, 6, 7, 8, 9,
    ],
)
def test_citest(
    model_name: Literal["asia"],
    citest_type_str: Literal["chi_square"] | Literal["g_sq"],
    level: float,
    sample_size: int,
    seed: (
        Literal[0]
        | Literal[1]
        | Literal[2]
        | Literal[3]
        | Literal[4]
        | Literal[5]
        | Literal[6]
        | Literal[7]
        | Literal[8]
        | Literal[9]
    ),
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
        # boolean=True,
        boolean=False,
        significance_level=level,
    )
    expected_result = p_value >= level
    print(
        f"[pgmpy] stat: {chi2}, p_value: {p_value}, dof: {dof}, level: {level}, "
        f"expected_result: {expected_result}"
    )

    # our implementation
    if citest_type_str == "chi_square":
        citest_type = openbnsllib.citest.ChiSquare(level)
    elif citest_type_str == "g_sq":
        citest_type = openbnsllib.citest.GSquare(level)
    else:
        raise ValueError(f"Unsupported citest_type_str: {citest_type_str}")
    var_indices = sorted(sorted_sepset_indices + [x_idx, y_idx])
    ct = openbnsllib.base.buildContingencyTable(var_indices, df_wrapper)
    computed_result = openbnsllib.citest.citest(
        x_idx, y_idx, sorted_sepset_indices, ct, citest_type
    )

    msg = (
        f"CItest for {x_idx}, {y_idx} | {sorted_sepset_indices} of {model_name}: "
        f"pgmpy: {expected_result}, ours: {computed_result}"
    )
    assert expected_result == computed_result, msg
