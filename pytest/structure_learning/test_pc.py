import pytest
from pgmpy.base import PDAG
from pgmpy.utils import get_example_model

import openbnsllib

from modules.utils import to_pgmpy
from modules.structural_distance import structural_errors


@pytest.mark.parametrize(
    "model_name, citest_type, sample_size, seed",
    [
        ("cancer", openbnsllib.citest.ChiSquare(0.05), int(1e5), 0),  # 5 nodes
        ("asia", openbnsllib.citest.ChiSquare(0.05), int(1e5), 0),  # 8 nodes
        ("child", openbnsllib.citest.ChiSquare(0.05), int(1e5), 0),  # 20 nodes
    ],
)
def test_pc2(model_name, citest_type, sample_size, seed):
    model_original = get_example_model(model_name)
    samples = model_original.simulate(sample_size, seed=seed)
    samples = samples[sorted(samples.columns)]
    df_wrapper = openbnsllib.base.DataframeWrapper(samples)
    _pdag = openbnsllib.structure_learning.PC(df_wrapper, citest_type)
    model_estimated = _pdag.to_pgmpy()
    errors = structural_errors(model_original, model_estimated)

    if errors["SHD"] > 0:
        print(f"Model: {model_name}, [SHD, ME, EE, DE, ED, MD, RD]: {errors}")
    assert errors["SHD"] == 0, f"Failed for model {model_name} with Errors: {errors}"
