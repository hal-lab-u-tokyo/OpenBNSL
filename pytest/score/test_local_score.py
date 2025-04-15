import pytest
import random
from pgmpy.utils import get_example_model
from pgmpy.estimators import K2, BDeu
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
    "score_name, score_type, pgmpy_score_func",
    [
        ("K2", openbnsllib.score.K2(), lambda s: K2(s)),
        ("BDeu", openbnsllib.score.BDeu(1.0), lambda s: BDeu(s, 1.0)),
    ]
)
def test_local_score(
    model_name, sample_size, seed,
    score_name, score_type, pgmpy_score_func
):
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
    pgmpy_scoring = pgmpy_score_func(samples)
    columns = sorted(samples.columns)
    child_variable_name = columns[child_idx]
    parent_variable_names = [columns[i] for i in selected_parent_indices]
    expected_score = pgmpy_scoring.local_score(
        child_variable_name, parent_variable_names
    )

    # our implementation
    var_indices = sorted(selected_parent_indices + [child_idx])
    ct = openbnsllib.base.buildContingencyTable(var_indices, df_wrapper)
    computed_score = openbnsllib.score.calculate_local_score(
        child_idx, selected_parent_indices, ct, score_type
    )

    msg = (
        f"{score_name} score for {child_idx} <= {selected_parent_indices} of {model_name}: "
        f"pgmpy: {expected_score}, ours: {computed_score}"
    )
    assert abs(expected_score - computed_score) < 1e-6, msg
