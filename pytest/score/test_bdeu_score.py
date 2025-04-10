import pytest
import random
from pgmpy.utils import get_example_model
from pgmpy.estimators import BDeuScore
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
def test_local_score(model_name, sample_size, seed):
    random.seed(seed)

    model_original = get_example_model(model_name)
    samples = model_original.simulate(sample_size, seed=seed)
    samples = samples[sorted(samples.columns)]
    df_wrapper = openbnsllib.base.DataframeWrapper(samples)

    n = df_wrapper.num_of_vars

    child_idx = random.randint(0, n - 1)
    candidate_parent_indices = [i for i in range(n) if i != child_idx]
    num_parent_candidates = random.randint(1, min(3, len(candidate_parent_indices)))
    selected_parent_indices = random.sample(
        candidate_parent_indices, num_parent_candidates
    )
    sorted_parent_indices = sorted(selected_parent_indices)

    # pgmpy implementation
    pgmpy_scoring = BDeuScore(samples, 1.0)
    columns = sorted(samples.columns)
    child_variable_name = columns[child_idx]
    parent_variable_names = [columns[i] for i in selected_parent_indices]
    expected_score = pgmpy_scoring.local_score(
        child_variable_name, parent_variable_names
    )

    # our implementation
    var_indices = sorted(selected_parent_indices + [child_idx])
    ct = openbnsllib.base.buildContingencyTable(var_indices, df_wrapper)
    score_type = openbnsllib.score.BDeu(1.0)
    computed_score = openbnsllib.score.calculate_local_score(
        child_idx, selected_parent_indices, ct, score_type
    )

    assert (
        abs(expected_score - computed_score) < 1e-6
    ), f"Expected score: {expected_score}, Computed score: {computed_score}"
