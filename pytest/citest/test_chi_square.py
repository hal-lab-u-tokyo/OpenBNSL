import pytest
import random
from pgmpy.utils import get_example_model
from pgmpy.estimators.CITests import chi_square

import openbnsllib


@pytest.mark.parametrize(
    "model_name, sample_size, seed",
    [
        ("cancer", int(1e5), 0),  # e.g. 5 nodes
        ("cancer", int(1e5), 1),  # e.g. 5 nodes
        ("cancer", int(1e5), 2),  # e.g. 5 nodes
        ("asia", int(1e5), 0),  # e.g. 8 nodes
        ("asia", int(1e5), 1),  # e.g. 8 nodes
        ("asia", int(1e5), 2),  # e.g. 8 nodes
        ("child", int(1e5), 0),  # e.g. 20 nodes
        ("child", int(1e5), 1),  # e.g. 20 nodes
        ("child", int(1e5), 2),  # e.g. 20 nodes
        ("alarm", int(1e5), 0),  # e.g. 37 nodes
        ("alarm", int(1e5), 1),  # e.g. 37 nodes
        ("alarm", int(1e5), 2),  # e.g. 37 nodes
    ],
)
def test_chi_square(model_name, sample_size, seed):
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
    expected_result = chi_square(
        X=samples.columns[x_idx],
        Y=samples.columns[y_idx],
        Z=[samples.columns[i] for i in selected_sepset_indices],
        data=samples,
        boolean=True,
        significance_level=0.05,
    )

    # our implementation
    var_indices = sorted(selected_sepset_indices + [x_idx, y_idx])
    ct = openbnsllib.base.buildContingencyTable(var_indices, df_wrapper)
    citest_type = openbnsllib.citest.ChiSquare(0.05)
    computed_result = openbnsllib.citest.citest(
        x_idx, y_idx, selected_sepset_indices, ct, citest_type
    )

    assert (
        expected_result == computed_result
    ), f"Expected result: {expected_result}, Computed result: {computed_result}"
