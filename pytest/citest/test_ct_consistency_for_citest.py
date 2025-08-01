# pytest/citest/test_ct_consistency_for_citest.py
# -------------------------------------------------------------------
# Purpose
# -------------------------------------------------------------------
# Ensure that the sparse ContingencyTable used by openbnsllib stores
# exactly the same *frequencies* as the dense contingency table that
# pgmpy builds inside its CI tests.  pgmpy includes zero–frequency
# cells whereas openbnsllib drops them, so we pad the sparse vector
# with zeros to equalise the lengths before multiset comparison.
# -------------------------------------------------------------------
from typing import Literal
import random
import numpy as np
import pandas as pd
import pytest
from pgmpy.utils import get_example_model

import openbnsllib


@pytest.mark.parametrize("model_name", ["cancer", "asia", "child", "alarm"])
@pytest.mark.parametrize("sample_size", [int(1e5)])
@pytest.mark.parametrize("seed", [0, 1])
def test_contingency_table_consistency(
    model_name: str,
    sample_size: int,
    seed: int,
):
    random.seed(seed)

    # ------------------------------------------------------------ #
    # 1) Generate synthetic data from pgmpy’s Bayesian‑network model
    # ------------------------------------------------------------ #
    model = get_example_model(model_name)
    samples = model.simulate(sample_size, seed=seed)
    samples = samples[sorted(samples.columns)]          # deterministic column order
    df_wrapper = openbnsllib.base.DataframeWrapper(samples)

    # ------------------------------------------------------------ #
    # 2) Randomly pick variables:   X, Y | Z
    # ------------------------------------------------------------ #
    n_vars    = df_wrapper.num_of_vars
    remaining = list(range(n_vars))

    x_idx = random.choice(remaining)
    remaining.remove(x_idx)

    y_idx = random.choice(remaining)
    remaining.remove(y_idx)

    z_idx = sorted(random.sample(remaining, k=random.randint(0, min(5, len(remaining)))))

    x_name = samples.columns[x_idx]
    y_name = samples.columns[y_idx]
    z_names = [samples.columns[i] for i in z_idx]

    var_indices = sorted(z_idx + [x_idx, y_idx])        # order expected by C++ CT

    # ------------------------------------------------------------ #
    # 3) Reproduce *pgmpy‑style* dense count vector
    # ------------------------------------------------------------ #
    counts_list = []
    if not z_names:
        ct_dense = pd.crosstab(samples[x_name], samples[y_name])
        counts_list.append(ct_dense.values.flatten())
    else:
        for _, grp in samples.groupby(z_names, observed=True, sort=False):
            # Even empty groups produce a zero‑only table in pgmpy → include them
            ct_dense = pd.crosstab(grp[x_name], grp[y_name])
            counts_list.append(ct_dense.values.flatten())

    pgmpy_counts = (
        np.concatenate(counts_list) if counts_list else np.array([], dtype=int)
    )
    # pgmpy_counts **contains zeros** by design → そのまま使用

    # ------------------------------------------------------------ #
    # 4) Build openbnsllib’s sparse contingency table
    # ------------------------------------------------------------ #
    ct = openbnsllib.base.ContingencyTable(var_indices, df_wrapper)
    obnsl_counts = np.fromiter(ct.counts.values(), dtype=int, count=len(ct.counts))
    # **no zeros** – sparse representation

    # ------------------------------------------------------------ #
    # 5) Pad sparse vector with zeros so that lengths match
    # ------------------------------------------------------------ #
    if len(obnsl_counts) < len(pgmpy_counts):
        pad_len = len(pgmpy_counts) - len(obnsl_counts)
        obnsl_counts = np.pad(obnsl_counts, (0, pad_len))
    elif len(obnsl_counts) > len(pgmpy_counts):
        raise AssertionError(
            "Sparse table has more non‑zero entries than pgmpy’s dense table!"
        )

    # ------------------------------------------------------------ #
    # 6) Compare totals and multiset equality (order‑insensitive)
    # ------------------------------------------------------------ #
    assert (
        pgmpy_counts.sum() == obnsl_counts.sum()
    ), "Sample‑size check failed: totals differ!"

    pg_sorted    = np.sort(pgmpy_counts)
    obnsl_sorted = np.sort(obnsl_counts)

    if not np.array_equal(pg_sorted, obnsl_sorted):
        # Detailed debug info to help pinpoint the mismatch
        max_len        = len(pg_sorted)
        diff           = pg_sorted - obnsl_sorted
        mismatch_idx   = np.where(diff != 0)[0]

        print(
            f"\n🔍  Count vectors differ at {len(mismatch_idx)} position(s). "
            f"Showing first few mismatches:"
        )
        for i in mismatch_idx[:10]:
            print(f"  idx {i}: pgmpy={pg_sorted[i]}  obnsl={obnsl_sorted[i]}")

        print("\n[openbnsllib raw] (key -> count)")
        for k, v in ct.counts.items():
            print(f"  {k}: {v}")
        print()
        raise AssertionError("Contingency counts do not match.")
