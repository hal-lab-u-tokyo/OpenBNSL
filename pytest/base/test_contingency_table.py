import pytest
from pgmpy.utils import get_example_model
import openbnsllib


@pytest.mark.parametrize(
    "model_name, sample_size, seed",
    [
        ("cancer", int(1e5), 0),  # e.g. 5 nodes
        ("asia", int(1e5), 0),  # e.g. 8 nodes
        ("child", int(1e5), 0),  # e.g. 20 nodes
    ],
)
def test_contingency_table(model_name, sample_size, seed):
    # Simulate samples from the example model.
    model_original = get_example_model(model_name)
    samples = model_original.simulate(sample_size, seed=seed)

    # Ensure columns are sorted (lexicographical order)
    samples = samples[sorted(samples.columns)]

    # Create a DataframeWrapper instance.
    df_wrapper = openbnsllib.base.DataframeWrapper(samples)

    # --- Test 1: Single-column contingency table ---
    # Use first column (index 0)
    ct_single = openbnsllib.base.buildContingencyTable([0], df_wrapper)

    # The table should have as many cells as the number of unique values in column 0.
    expected_cells = df_wrapper.num_of_values[0]
    assert (
        len(ct_single.counts) == expected_cells
    ), f"Expected {expected_cells} cells, but got {len(ct_single.counts)}"

    # Sum of counts must equal the total number of datapoints.
    total_counts = sum(ct_single.counts)
    assert (
        total_counts == df_wrapper.num_of_datapoints
    ), f"Total counts {total_counts} does not equal number of datapoints {df_wrapper.num_of_datapoints}"

    # --- Test 2: Two-column contingency table (if available) ---
    if df_wrapper.num_of_vars >= 2:
        ct_double = openbnsllib.base.buildContingencyTable([0, 1], df_wrapper)

        # Expect the number of cells to be the product of unique values in columns 0 and 1.
        expected_cells = df_wrapper.num_of_values[0] * df_wrapper.num_of_values[1]
        assert (
            len(ct_double.counts) == expected_cells
        ), f"Expected {expected_cells} cells, but got {len(ct_double.counts)}"

        # Sum of counts must equal the total number of datapoints.
        total_counts = sum(ct_double.counts)
        assert (
            total_counts == df_wrapper.num_of_datapoints
        ), f"Total counts {total_counts} does not equal number of datapoints {df_wrapper.num_of_datapoints}"


def test_contingency_table_unsorted_input():
    model_name = "cancer"
    sample_size = int(1e5)
    seed = 0

    model_original = get_example_model(model_name)
    samples = model_original.simulate(sample_size, seed=seed)
    samples = samples[sorted(samples.columns)]
    df_wrapper = openbnsllib.base.DataframeWrapper(samples)

    unsorted_vars = [1, 0]

    with pytest.raises(Exception) as exc_info:
        openbnsllib.base.buildContingencyTable(unsorted_vars, df_wrapper)

    assert "sorted" in str(exc_info.value)
