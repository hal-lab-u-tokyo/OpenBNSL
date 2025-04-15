import pytest
import random
from pgmpy.utils import get_example_model
from pgmpy.estimators.CITests import chi_square, g_sq

import openbnsllib

@pytest.mark.parametrize(
    "model_name, sample_size, seed",
    [
        ("cancer", int(1e4), 0),  # e.g. 5 nodes
        ("asia", int(1e4), 0),  # e.g. 8 nodes
        ("child", int(1e4), 0),  # e.g. 20 nodes
        ("alarm", int(1e4), 0),  # e.g. 37 nodes
    ],
)
@pytest.mark.parametrize(
    "citest_name, citest_type, pgmpy_citest_func",
    [
        (
            "ChiSquare",
            openbnsllib.citest.ChiSquare(0.05),
            lambda X, Y, Z, data: chi_square(
                X=X, Y=Y, Z=Z, data=data, boolean=True, significance_level=0.05
            ),
        ),
        (
            "GSquare",
            openbnsllib.citest.GSquare(0.05),
            lambda X, Y, Z, data: g_sq(
                X=X, Y=Y, Z=Z, data=data, boolean=True, significance_level=0.05
            ),
        ),
    ]
)
def test_citest(
    model_name, sample_size, seed,
    citest_name, citest_type, pgmpy_citest_func):
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
    num_sepset_candidates = random.randint(0, min(3, len(remaining_indices)))
    selected_sepset_indices = random.sample(remaining_indices, num_sepset_candidates)
    sorted_sepset_indices = sorted(selected_sepset_indices)

    # pgmpy implementation
    x_name = samples.columns[x_idx]
    y_name = samples.columns[y_idx]
    z_names = [samples.columns[i] for i in sorted_sepset_indices]
    expected_result = pgmpy_citest_func(
        X=x_name,
        Y=y_name,
        Z=z_names,
        data=samples,
    )

    # our implementation
    var_indices = sorted(selected_sepset_indices + [x_idx, y_idx])
    ct = openbnsllib.base.buildContingencyTable(var_indices, df_wrapper)
    computed_result = openbnsllib.citest.citest(
        x_idx, y_idx, selected_sepset_indices, ct, citest_type
    )

    msg = (
        f"CItest for {x_idx}, {y_idx} | {selected_sepset_indices} of {model_name}: "
        f"pgmpy: {expected_result}, ours: {computed_result}"
    )
    assert expected_result == computed_result, msg