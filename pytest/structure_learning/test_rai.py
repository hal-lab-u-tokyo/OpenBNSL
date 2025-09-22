import pytest
import random
import time
from pgmpy.utils import get_example_model

import openbnsllib

from helpers.pgmpy_bridge import to_pgmpy, to_openbnsl
from helpers.structural_distance import structural_errors


@pytest.mark.parametrize(
    "model_name",
    [
        # "asia",
        "cancer",
        "earthquake",
        "sachs",
        # "survey",  # Small networks
        # "alarm",
        # "barley",
        # "child",
        # "insurance",
        # "mildew",
        # "water",  # Medium networks
        # "hailfinder", "hepar2", "win95pts", # Large networks
        # "andes", "diabetes", "link", "munin1", "pathfinder", "pigs", # Extra large networks
        # "munin", "munin2", "munin3", "munin4", # Very extra large networks
    ],
)
@pytest.mark.parametrize("seed", [0])
def test_rai(model_name, seed):
    random.seed(seed)

    model_original = get_example_model(model_name)
    samples = model_original.simulate(int(1e6), seed=seed)
    samples = samples[sorted(samples.columns)]

    df_wrapper = openbnsllib.base.DataframeWrapper(samples)
    citest_type = openbnsllib.citest.ChiSquare(level=0.01)

    t0 = time.perf_counter()
    _pdag = openbnsllib.structure_learning.rai(
        df_wrapper, citest_type, max_cond_vars=len(samples.columns)
    )
    t1 = time.perf_counter()
    elapsed = t1 - t0
    print(f"[TIME] model={model_name}, seed={seed}, time={elapsed:.2f}s")

    expected_obnsl = to_pgmpy(_pdag, list(samples.columns))
    error_obnsl = structural_errors(model_original, expected_obnsl)

    msg = f"Structural errors for {model_name}): {error_obnsl} in {elapsed:.2f}s"
    assert error_obnsl["SHD"] == 0, msg
