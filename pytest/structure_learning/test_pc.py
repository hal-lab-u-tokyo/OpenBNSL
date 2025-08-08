import pytest
import random
from pgmpy.utils import get_example_model
from pgmpy.base import PDAG
from pgmpy.estimators import PC

import openbnsllib

from modules.utils import to_pgmpy, to_openbnsl
from modules.structural_distance import structural_errors


@pytest.mark.parametrize("model_name", ["cancer", "asia", "child", "alarm"])
# @pytest.mark.parametrize("model_name", ["cancer", "asia"])
@pytest.mark.parametrize("max_cond_vars", [5])  # pgmpy default is 5
@pytest.mark.parametrize("citest_type_str", ["oracle"])
@pytest.mark.parametrize("level", [0.05])
@pytest.mark.parametrize("sample_size", [int(1e5)])
@pytest.mark.parametrize("seed", [0])
def test_pc(model_name, max_cond_vars, citest_type_str, level, sample_size, seed):
    random.seed(seed)

    model_original = get_example_model(model_name)
    samples = model_original.simulate(sample_size, seed=seed)
    samples = samples[sorted(samples.columns)]

    # # pgmpy implementation
    # estimator = PC(samples)
    # expected_pgmpy = estimator.estimate(
    #     variant="orig",
    #     ci_test=citest_type_str,
    #     significance_level=level,
    #     max_cond_vars=max_cond_vars,
    # )

    # our implementation
    df_wrapper = openbnsllib.base.DataframeWrapper(samples)
    if citest_type_str == "chi_square":
        citest_type = openbnsllib.citest.ChiSquare(level)
    elif citest_type_str == "g_sq":
        citest_type = openbnsllib.citest.GSquare(level)
    elif citest_type_str == "oracle":
        oracle_graph = to_openbnsl(model_original, df_wrapper.col_str2idx)
        citest_type = openbnsllib.citest.OracleGraph(oracle_graph)
    else:
        raise ValueError(f"Unsupported citest_type_str: {citest_type_str}")
    _pdag = openbnsllib.structure_learning.PC(
        df_wrapper,
        citest_type,
        max_cond_vars,
    )
    expected_obnsl = to_pgmpy(_pdag, list(samples.columns))

    # error_pgmpy = structural_errors(expected_pgmpy, model_original)
    error_obnsl = structural_errors(expected_obnsl, model_original)
    print(f"Structural errors for {model_name}: {error_obnsl}")

    # assert error_pgmpy == error_obnsl, (
    #     f"Structural errors do not match: "
    #     f"pgmpy: {error_pgmpy}, "
    #     f"openbnsllib: {error_obnsl}"
    # )

    model_original.to_graphviz().draw(
        f"tmp_images/{model_name}_original.png", prog="dot", format="png"
    )
    expected_obnsl.to_graphviz().draw(
        f"tmp_images/{model_name}_output.png", prog="dot", format="png"
    )
