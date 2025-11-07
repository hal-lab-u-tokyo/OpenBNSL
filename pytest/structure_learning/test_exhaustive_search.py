import pytest
import time
from pgmpy.utils import get_example_model

import openbnsllib
from helpers.pgmpy_bridge import to_pgmpy
from helpers.structural_distance import structural_errors

@pytest.mark.parametrize("model_name", ["cancer", "asia", "child"])
@pytest.mark.parametrize("score_type", [openbnsllib.score.BDeu(1.0)])
@pytest.mark.parametrize("sample_size", [int(1e5)])
@pytest.mark.parametrize("seed", [0])
def test_exhaustive_search(model_name, score_type, sample_size, seed):

    original_pdag_pgmpy = get_example_model(model_name)
    samples = original_pdag_pgmpy.simulate(sample_size, seed=seed)
    samples = samples[sorted(samples.columns)]

    df_wrapper = openbnsllib.base.DataframeWrapper(samples)

    start = time.perf_counter()
    learned_pdag_obnsl = openbnsllib.structure_learning.exhaustive_search(
        df_wrapper, score_type, max_parents=3
    )
    elapsed = time.perf_counter() - start
    print(f"[OpenBNSL] model={model_name}, seed={seed}, time={elapsed:.2f}s")

    learned_pdag_pgmpy = to_pgmpy(learned_pdag_obnsl, list(samples.columns))
    errors = structural_errors(original_pdag_pgmpy, learned_pdag_pgmpy)

    if errors["SHD"] > 0:
        print(f"Model: {model_name}, [SHD, ME, EE, DE, ED, MD, RD]: {errors}")
    assert errors["SHD"] == 0, f"Failed for model {model_name} with Errors: {errors}"
