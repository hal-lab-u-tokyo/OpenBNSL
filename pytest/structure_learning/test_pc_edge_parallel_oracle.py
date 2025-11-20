import pytest
import random
from pgmpy.utils import get_example_model

import openbnsllib

from helpers.pgmpy_bridge import to_pgmpy, to_openbnsl
from helpers.structural_distance import structural_errors


@pytest.mark.parametrize(
    "model_name",
    [
        "asia",
        "cancer",
        "earthquake",
        "sachs",
        "survey",  # Small networks
        "alarm",
        "barley",
        "child",
        "insurance",
        "mildew",
        "water",  # Medium networks
        # "hailfinder", "hepar2", "win95pts", # Large networks
        # "andes", "diabetes", "link", "munin1", "pathfinder", "pigs", # Extra large networks
        # "munin", "munin2", "munin3", "munin4", # Very extra large networks
    ],
)
@pytest.mark.parametrize("seed", [0])
def test_pc(model_name, seed):
    random.seed(seed)

    model_original = get_example_model(model_name)
    samples = model_original.simulate(int(1e3), seed=seed)  # dummy
    samples = samples[sorted(samples.columns)]

    df_wrapper = openbnsllib.base.DataframeWrapper(samples)
    oracle_graph = to_openbnsl(model_original, df_wrapper.col_str2idx)
    citest_type = openbnsllib.citest.OracleGraph(oracle_graph)

    _pdag = openbnsllib.structure_learning.pc_edge_parallel(
        df_wrapper, citest_type, max_cond_vars=len(samples.columns)
    )
    expected_obnsl = to_pgmpy(_pdag, list(samples.columns))
    error_obnsl = structural_errors(model_original, expected_obnsl)

    msg = f"Structural errors for {model_name}: {error_obnsl}"
    # print(msg)
    assert error_obnsl["SHD"] == 0, msg
