# Find Peaks Python Bindings

This directory contains Python bindings for the Find Peaks library, allowing you to use both the C and C++ implementations from Python and compare them with the reference SciPy implementation.

## Overview

The Python bindings provide:
- Direct access to the C implementation (`find_peaks_c_py` module)
- Direct access to the C++ implementation (`find_peaks_cpp_py` module)
- Identical interfaces to SciPy's `find_peaks` function for easy comparison
- Performance benchmarks between all three implementations
- Validation tests to ensure compatibility

## Building the Python Bindings

The bindings are built automatically when you build the main project with Python bindings enabled:

```bash
# From the main project directory
mkdir build && cd build
cmake .. -DBUILD_PYTHON_BINDINGS=ON
cmake --build .
```

Or build just the Python modules:

```bash
cmake --build . --target find_peaks_cpp_py find_peaks_c_py
```

## Usage

After building, you can import and use the bindings directly:

```python
import sys
import numpy as np
from scipy.signal import find_peaks as scipy_find_peaks

# Add the path to your built modules
sys.path.insert(0, '/path/to/built/modules')

# Import the C and C++ implementations
from find_peaks_cpp_py import find_peaks as find_peaks_cpp
from find_peaks_c_py import find_peaks as find_peaks_c

# Create a test signal
signal = np.array([0, 1, 0, 2, 0, 3, 0, 2, 0, 1, 0])

# Run with SciPy implementation
scipy_peaks, scipy_props = scipy_find_peaks(signal, height=1.0)

# Run with C++ implementation
cpp_peaks, cpp_props = find_peaks_cpp(signal, height=1.0)

# Run with C implementation
c_peaks, c_props = find_peaks_c(signal, height=1.0)

# Compare results
print(f"SciPy found {len(scipy_peaks)} peaks at positions: {scipy_peaks}")
print(f"C++ found {len(cpp_peaks)} peaks at positions: {cpp_peaks}")
print(f"C found {len(c_peaks)} peaks at positions: {c_peaks}")
```

## Parameter Compatibility

The Python bindings accept the same parameters as SciPy's `find_peaks` function:

- `height`: Required height of peaks. Either a number or a 2-element tuple/array
- `threshold`: Required threshold of peaks. Either a number or a 2-element tuple/array
- `distance`: Required minimal horizontal distance between neighboring peaks
- `prominence`: Required prominence of peaks. Either a number or a 2-element tuple/array
- `width`: Required width of peaks. Either a number or a 2-element tuple/array
- `wlen`: Window length for prominence calculation
- `rel_height`: Relative height for width calculation (from 0 to 1)
- `plateau_size`: Required size of flat top plateau. Either a number or a 2-element tuple/array

## Return Value Format

Like SciPy, both implementations return a tuple with two elements:
1. A numpy array containing the indices of detected peaks
2. A dictionary containing properties of the detected peaks:
   - `prominences`: Prominence values for each peak
   - `left_bases`, `right_bases`: Indices of bases for each peak
   - `widths`: Width values for each peak
   - `width_heights`: Height levels at which widths were measured
   - `left_ips`, `right_ips`: Interpolated positions of width measurement
   - `plateau_sizes`: Size of the flat peak plateau for each peak
   - `left_edges`, `right_edges`: Indices of plateau edges

## Validation and Testing

Run the included validation tests to ensure the C and C++ implementations match SciPy's results:

```bash
cd build
ctest -R test_find_peaks_scipy
```

Or run the tests directly:

```bash
cd bindings

# Python needs to know where to load find_peaks_c_py find_peaks_cpp_py form.
# So before running pytest, the path to the built binaries (.so or .pyd) of bindings
# need  to be added as an environment variable:

# export FP_BINDING_PATH=<path_to_module_dir>
# or 
# set FP_BINDING_PATH=<path_to_module_dir>

python -m pytest test_find_peaks_scipy.py -v
```

## Performance Benchmarking

Use the included Jupyter notebook to run performance benchmarks comparing the C, C++, and SciPy implementations:

```bash
cd bindings
jupyter-notebook
# Open the notebook file and run the benchmarks
```

<!-- ## Implementation Details

The bindings use pybind11 to expose the C and C++ functionality to Python. The binding code converts between Python/NumPy data types and native C/C++ types, ensuring that the interface matches SciPy's while leveraging the performance benefits of the native implementations. -->
