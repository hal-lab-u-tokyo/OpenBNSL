# Helpers

- [`./pgmpy_bridge.py`](./pgmpy_bridge.py) 
  Conversion utilities between OpenBNSL’s C++ PDAG (`openbnsllib.base.PDAG`) and pgmpy’s Python PDAG (`pgmpy.base.PDAG`).  
  Provides `to_pgmpy()` and `to_openbnsl()` for bidirectional conversion.

- [`./structural_distance.py`](./structural_distance.py) 
  Structural comparison metrics for PDAGs. Includes functions to:
  - Convert PDAGs to CPDAGs (`PDAG2CPDAG`)
  - Compute Structural Hamming Distance (SHD) and its components (ME, EE, DE, ED, MD, RD) via `structural_errors()`.

These helpers are pure Python and are used in tests and benchmarks for model evaluation and interoperability.
