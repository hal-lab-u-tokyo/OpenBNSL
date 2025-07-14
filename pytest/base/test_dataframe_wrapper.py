import pytest
from pgmpy.utils import get_example_model

import openbnsllib


@pytest.mark.parametrize(
    "model_name, sample_size, seed",
    [
        ("cancer", int(1e3), 0),  # 5 nodes
        ("asia", int(1e3), 0),  # 8 nodes
        ("child", int(1e3), 0),  # 20 nodes
    ],
)
def test_dataframe_wrapper_basic(model_name, sample_size, seed):
    # Get example model and simulate samples.
    model_original = get_example_model(model_name)
    samples = model_original.simulate(sample_size, seed=seed)

    # Force the columns to be in lexicographical order.
    samples = samples[sorted(samples.columns)]

    # Create the DataframeWrapper instance.
    df_wrapper = openbnsllib.base.DataframeWrapper(samples)

    # Check that the number of variables and datapoints match the DataFrame.
    assert df_wrapper.num_of_vars == len(
        samples.columns
    ), f"Expected {len(samples.columns)} variables, got {df_wrapper.num_of_vars}"
    assert (
        df_wrapper.num_of_datapoints == samples.shape[0]
    ), f"Expected {samples.shape[0]} datapoints, got {df_wrapper.num_of_datapoints}"

    # Check that column names are stored in lexicographical order.
    sorted_cols = sorted(samples.columns)
    assert (
        df_wrapper.col_idx2str == sorted_cols
    ), f"Expected col_idx2str {sorted_cols}, got {df_wrapper.col_idx2str}"

    # Check that each column name maps to the correct index.
    for idx, col in enumerate(sorted_cols):
        assert col in df_wrapper.col_str2idx, f"Column {col} missing in col_str2idx"
        assert (
            df_wrapper.col_str2idx[col] == idx
        ), f"Expected {col} to map to {idx}, got {df_wrapper.col_str2idx[col]}"

    # Check that for each column, the unique values (val_idx2str)
    # are the lexicographically sorted unique values in the DataFrame.
    for i, col in enumerate(sorted_cols):
        unique_vals = sorted(samples[col].unique())
        assert (
            df_wrapper.val_idx2str[i] == unique_vals
        ), f"For column {col}: expected unique values {unique_vals}, got {df_wrapper.val_idx2str[i]}"
        # Also, number of unique values must match.
        assert df_wrapper.num_of_values[i] == len(
            unique_vals
        ), f"For column {col}: expected num_of_values {len(unique_vals)}, got {df_wrapper.num_of_values[i]}"

    # Too Slow to check in Python for loop
    # Check that the row-major representation is the transpose of the column-major one.
    # num_rows, num_cols = samples.shape[0], len(samples.columns)
    # for j in range(num_cols):
    #     col_data = df_wrapper.data_column_major[j]
    #     for i in range(num_rows):
    #         assert (
    #             df_wrapper.data_row_major[i][j] == col_data[i]
    #         ), f"Mismatch at row {i}, col {j}: expected {col_data[i]}, got {df_wrapper.data_row_major[i][j]}"


def test_dataframe_wrapper_invalid_input():
    # Test with a non-DataFrame input.
    with pytest.raises(Exception):
        DataframeWrapper(42)
