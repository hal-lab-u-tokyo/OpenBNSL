import time
import pytest
from pgmpy.utils import get_example_model
import openbnsllib


@pytest.mark.parametrize("model_name", ["cancer", "asia", "child", "alarm"])
@pytest.mark.parametrize("sample_size", [int(1e5)])
@pytest.mark.parametrize("seed", [0])
def test_marginalize_correctness_and_speed(model_name, sample_size, seed, capsys):

    model = get_example_model(model_name)
    samples = model.simulate(sample_size, seed=seed)
    samples = samples[sorted(samples.columns)]
    dfw = openbnsllib.base.DataframeWrapper(samples)

    num_vars = dfw.num_vars

    test_pairs = []
    if num_vars >= 2:
        S = [0, 1]
        for T in [[0], [1], [0, 1]]:
            test_pairs.append((S, T))
    if num_vars >= 3:
        S = [0, 1, 2]
        for T in [[0, 2], [1, 2], [0], [2]]:
            test_pairs.append((S, T))

    report_lines = []
    seen_S = set()
    for S, _ in test_pairs:
        key = tuple(S)
        if key in seen_S:
            continue
        seen_S.add(key)

        ct_S = openbnsllib.base.ContingencyTable(S, dfw)

        for _, T in filter(lambda p: p[0] == S, test_pairs):
            t0 = time.perf_counter()
            ct_T_from_marg = ct_S.marginalize_to(T)
            t1 = time.perf_counter()
            t_marg = t1 - t0

            t2 = time.perf_counter()
            ct_T_direct = openbnsllib.base.ContingencyTable(T, dfw)
            t3 = time.perf_counter()
            t_direct = t3 - t2

            assert ct_T_from_marg.counts == ct_T_direct.counts
            assert sum(ct_T_from_marg.counts.values()) == sum(
                ct_T_direct.counts.values()
            )

            speedup = (t_direct / t_marg) if t_marg > 0 else float("inf")
            report_lines.append(
                f"{model_name}: n={dfw.num_datapoints}, S={S} -> T={T} | "
                f"marg={t_marg:.6f}s, direct={t_direct:.6f}s, speedup≈{speedup:.2f}x"
            )

    print("\n".join(report_lines))
