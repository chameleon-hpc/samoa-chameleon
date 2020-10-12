# Python library for reading Samoa XDMF files

## Dependencies

- `h5py`
- `pyqtree`
- `numba`
- `lxml`
- `numpy`
- `pandas`

## Precompilation

Run `./sampler.py` to compile the DG sampling functions to a native library using Numba/LLVM. This dramatically improves performance.