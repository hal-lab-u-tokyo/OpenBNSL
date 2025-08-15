import pytest
import random
from pgmpy.utils import get_example_model

import openbnsllib

from helpers.pgmpy_bridge import to_pgmpy, to_openbnsl
from helpers.structural_distance import structural_errors


@pytest.mark.parametrize("model_name", ["cancer", "asia", "child", "alarm"])
@pytest.mark.parametrize("stable", [True, False])
@pytest.mark.parametrize("seed", [0])
def test_pc(model_name, stable, seed):
    random.seed(seed)

    model_original = get_example_model(model_name)
    samples = model_original.simulate(int(1e3), seed=seed)  # dummy
    samples = samples[sorted(samples.columns)]

    df_wrapper = openbnsllib.base.DataframeWrapper(samples)
    oracle_graph = to_openbnsl(model_original, df_wrapper.col_str2idx)
    citest_type = openbnsllib.citest.OracleGraph(oracle_graph)

    _pdag = openbnsllib.structure_learning.pc(
        df_wrapper, citest_type, max_cond_vars=len(samples.columns), stable=stable
    )
    expected_obnsl = to_pgmpy(_pdag, list(samples.columns))
    error_obnsl = structural_errors(model_original, expected_obnsl)

    msg = f"Structural errors for {model_name} (stable={stable}): {error_obnsl}"
    # print(msg)
    assert error_obnsl["SHD"] == 0, msg
