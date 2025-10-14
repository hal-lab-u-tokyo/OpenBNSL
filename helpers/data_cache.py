import os
import tempfile
import pandas as pd
from typing import Callable


def _cache_dir(model_name: str, seed: int) -> str:
    return os.path.join("data_cache", model_name, f"seed_{seed}")


def _cache_path(model_name: str, seed: int, compress: bool) -> str:
    fname = "samples.pkl.gz" if compress else "samples.pkl"
    return os.path.join(_cache_dir(model_name, seed), fname)


def _atomic_write_pickle(df: pd.DataFrame, path: str, compress: bool) -> None:
    os.makedirs(os.path.dirname(path), exist_ok=True)
    suffix = ".pkl.gz" if compress else ".pkl"
    with tempfile.NamedTemporaryFile(
        delete=False, dir=os.path.dirname(path), suffix=suffix
    ) as tmp:
        tmp_path = tmp.name
    try:
        df.to_pickle(tmp_path, compression=("gzip" if compress else None))
        os.replace(tmp_path, path)
    finally:
        if os.path.exists(tmp_path):
            try:
                os.remove(tmp_path)
            except:
                pass


def get_samples_head(
    *,
    model_name: str,
    seed: int,
    num_samples: int,
    generator: Callable[[int, int], pd.DataFrame],
    compress: bool = True,
) -> pd.DataFrame:
    """
    Get a head of samples for a specific model and seed, using cached samples if available.
    """
    cpath = _cache_path(model_name, seed, compress)
    if os.path.exists(cpath):
        df = pd.read_pickle(cpath, compression=("gzip" if compress else None))
        if len(df) >= num_samples:
            print(f"Using cached samples for model '{model_name}' with seed {seed}.")
            df = df[sorted(df.columns)]
            return df.head(num_samples).copy()

    print(f"Generating samples for model '{model_name}' with seed {seed}.")
    df_new = generator(num_samples, seed)
    df_new = df_new[sorted(df_new.columns)]
    _atomic_write_pickle(df_new, cpath, compress)
    return df_new.head(num_samples).copy()
