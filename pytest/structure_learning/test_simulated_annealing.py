import pytest
import time
from pgmpy.utils import get_example_model

import openbnsllib
from helpers.pgmpy_bridge import to_openbnsl


@pytest.mark.parametrize("model_name", ["cancer", "asia", "child", "alarm"])
@pytest.mark.parametrize("score_type", [openbnsllib.score.BDeu(1.0)])
@pytest.mark.parametrize("sample_size", [int(1e5)])
@pytest.mark.parametrize("seed", [0])
def test_simulated_annealing(model_name, score_type, sample_size, seed):

    original_pdag_pgmpy = get_example_model(model_name)
    samples = original_pdag_pgmpy.simulate(sample_size, seed=seed)
    samples = samples[sorted(samples.columns)]

    df_wrapper = openbnsllib.base.DataframeWrapper(samples)
    original_pdag_obnsl = to_openbnsl(original_pdag_pgmpy, df_wrapper.col_str2idx)

    start = time.perf_counter()
    learned_pdag_obnsl = openbnsllib.structure_learning.simulated_annealing(
        df_wrapper, score_type, max_parents=3
    )
    elapsed = time.perf_counter() - start
    print(f"[OpenBNSL] model={model_name}, seed={seed}, time={elapsed:.2f}s")

    original_score = original_pdag_obnsl.score(df_wrapper, openbnsllib.score.BDeu(1.0))
    learned_score = learned_pdag_obnsl.score(df_wrapper, openbnsllib.score.BDeu(1.0))
    score_ratio = (
        learned_score / original_score if original_score != 0 else float("inf")
    )

    msg = f"score: {learned_score:.2f}/{original_score:.2f} ({score_ratio:.2f}) in {elapsed:.2f}s"
    # print(msg)
    assert score_ratio >= 0.99, msg
