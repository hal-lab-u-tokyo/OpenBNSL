import pytest
from pgmpy.utils import get_example_model
import openbnsllib


@pytest.mark.parametrize("model_name", ["cancer", "asia", "child", "alarm"])
@pytest.mark.parametrize("sample_size", [int(1e5)])
@pytest.mark.parametrize("seed", [0])
def test_contingency_table(model_name, sample_size, seed):

    model_original = get_example_model(model_name)
    samples = model_original.simulate(sample_size, seed=seed)
    samples = samples[sorted(samples.columns)]
    df_wrapper = openbnsllib.base.DataframeWrapper(samples)

    ct1 = openbnsllib.base.ContingencyTable([0], df_wrapper)

    n_unique = df_wrapper.num_of_values[0]
    assert len(ct1.counts) <= n_unique
    assert sum(ct1.counts.values()) == df_wrapper.num_of_datapoints

    if df_wrapper.num_of_vars >= 2:
        ct2 = openbnsllib.base.ContingencyTable([0, 1], df_wrapper)

        n_cells = df_wrapper.num_of_values[0] * df_wrapper.num_of_values[1]
        assert len(ct2.counts) <= n_cells
        assert sum(ct2.counts.values()) == df_wrapper.num_of_datapoints
