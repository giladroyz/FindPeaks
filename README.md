![GitHub Actions Workflow Status](https://github.com/username/repository/actions/workflows/main.yml/badge.svg)

# Find Peaks Library

C/C++ library for peak detection and analysis in signal data.

The find_peaks algorithm is a direct C & C++ conversion of SciPy's [find_peaks](https://docs.scipy.org/doc/scipy/reference/generated/scipy.signal.find_peaks.html) algorithm.

## Features

The library offers comprehensive peak detection with customizable criteria and detailed peak property analysis.

- **Flexible peak detection** with customizable filtering criteria
- **Detailed peak analysis** including position, height, prominence, width, and plateau characteristics
- **Both C and C++ interfaces** for easy integration into any project
- **Documented examples** and parameter descriptions

## Key Peak Detection Parameters

- **Height**: Filter peaks based on their absolute height
- **Threshold**: Require peaks to exceed neighbors by a specific amount
- **Distance**: Ensure minimum separation between detected peaks
- **Prominence**: Filter peaks based on how much they stand out from baseline
- **Width**: Filter peaks based on their width at a configurable height level
- **Plateau Size**: Filter peaks based on the number of samples at their maximum value

## C++ Interface Usage

```cpp
#include "find_peaks.hpp"
#include <vector>
#include <iostream>

int main() {
    // Sample signal data
    std::vector<double> signal = {0, 1, 3, 1, 0, 2, 4, 2, 0, 1, 2, 1, 0};
    
    // Create custom peak detection conditions
    PeakConditions conditions;
    conditions.set_height(2.0);           // Minimum height of 2.0
    conditions.set_prominence(1.0);       // Minimum prominence of 1.0
    conditions.set_distance(2);           // At least 2 samples between peaks
    conditions.set_width(1.0, 3.0);       // Width between 1.0 and 3.0 samples
    conditions.set_rel_height(0.7);       // Measure width at 70% of peak height
    
    // Find peaks with custom conditions
    std::vector<peak_result_t> peaks = find_peaks(signal, conditions);
    
    // Access detailed peak information
    for (const auto& peak : peaks) {
        std::cout << "Peak at position " << peak.peak << ":\n";
        std::cout << "  Height: " << peak.peak_height << "\n";
        std::cout << "  Prominence: " << peak.prominence.prominence << "\n";
        std::cout << "  Width: " << peak.width.width << " at height " << peak.width.width_height << "\n";
        std::cout << "  Left threshold: " << peak.threshold.left_threshold << "\n";
        std::cout << "  Right threshold: " << peak.threshold.right_threshold << "\n";
        
        // Plateau information
        if (peak.plateau.plateau_size > 1) {
            std::cout << "  Plateau size: " << peak.plateau.plateau_size << "\n";
        }
    }
    
    return 0;
}
```

## C Interface Usage

```c
#include "find_peaks.h"
#include <stdio.h>
#include <stdlib.h>

int main() {
    // Sample signal data
    double signal[] = {0, 1, 3, 1, 0, 2, 4, 2, 0, 1, 2, 1, 0};
    size_t signal_size = sizeof(signal) / sizeof(signal[0]);
    
    // Create conditions and customize
    fp_conditions_t cond;
    fp_init_conditions(&cond);
    
    fp_cond_set_height_mn(&cond, 2.0);          // Minimum height of 2.0
    fp_cond_set_prominence_mn(&cond, 1.0);      // Minimum prominence of 1.0
    fp_cond_set_distance_mn(&cond, 2);          // At least 2 samples between peaks
    fp_cond_set_width_mn_mx(&cond, 1.0, 3.0);   // Width between 1.0 and 3.0 samples
    fp_cond_set_rel_height_mn(&cond, 0.7);      // Measure width at 70% of peak height
    
    // Variables to store results
    peak_result_t* peaks = NULL;
    size_t num_peaks = 0;
    
    // Find peaks
    if (find_peaks(signal, signal_size, cond, &peaks, &num_peaks) == FP_STATUS_OK) {
        printf("Detected %zu peaks\n", num_peaks);
        
        // Access detailed peak information
        for (size_t i = 0; i < num_peaks; i++) {
            printf("Peak at position %zu:\n", peaks[i].peak);
            printf("  Height: %f\n", peaks[i].peak_height);
            printf("  Prominence: %f\n", peaks[i].prominence.prominence);
            printf("  Width: %f at height %f\n", 
                   peaks[i].width.width, peaks[i].width.width_height);
            printf("  Left threshold: %f\n", peaks[i].threshold.left_threshold);
            printf("  Right threshold: %f\n", peaks[i].threshold.right_threshold);
            
            // Plateau information
            if (peaks[i].plateau.plateau_size > 1) {
                printf("  Plateau size: %zu\n", peaks[i].plateau.plateau_size);
            }
        }
        
        // Free allocated memory
        fp_free(peaks);
    }
    
    return 0;
}
```

## Peak Detection Algorithm

The algorithm have the following steps:
1. Identifies all local maxima in the signal
2. Analyzes peak properties including height, prominence, width, etc.
3. Filters peaks based on the specified conditions
4. Returns detailed information about each qualifying peak

### Reference Implementation

This library's algorithm is based on SciPy's `scipy.signal.find_peaks` function:

- **SciPy Reference**: [scipy.signal.find_peaks](https://docs.scipy.org/doc/scipy/reference/generated/scipy.signal.find_peaks.html)
- **SciPy GitHub**: [scipy/signal/_peak_finding.py](https://github.com/scipy/scipy/blob/main/scipy/signal/_peak_finding.py)

The C/C++ implementation follows the same functionality and produces results compatible with the Python version.

## Peak Result Structure

Each detected peak includes comprehensive information:

```cpp
typedef struct{
    size_t peak;                     // Index of the peak in the input data
    double peak_height;              // Height (value) of the peak

    lpr_peak_plateau_t plateau;      // Information about the peak's plateau
    lr_peak_threshold_t threshold;   // Threshold values for the peak
    lpr_peak_prominence_t prominence; // Prominence information for the peak
    whlr_peak_width_t width;         // Width information for the peak
} peak_result_t;
```

## Performance Considerations

- The algorithm has O(nlogn) avarage time complexity (quick sort + O(n) operations)
- For large datasets, use the distance parameter to improve performance by reducing the number of peaks analyzed


## Building and Integration

This library can be used as a header-only library or built with CMake for all components.

### Prerequisites

- CMake 3.14 or higher
- C++11 compatible compiler
- Python 3.x with numpy & scipy installed (if building Python bindings)

### Integration Options

1. **Code-only**: Simply include the appropriate files (`find_peaks.hpp` and `find_peaks.cpp` for C++ or `find_peaks.h` and `find_peaks.c` for C) in your project.

2. **Static Library**: Link against the built libraries:
    - For C++: `find_peaks_cpp`
    - For C: `find_peaks_c`

### CMake Integration

To use find_peaks in your CMake project:

```cmake
# Add find_peaks as a subdirectory
add_subdirectory(path/to/find_peaks)

# For C++ usage
target_link_libraries(your_target PRIVATE find_peaks_cpp)
target_include_directories(your_target PRIVATE path/to/find_peaks/src/cpp)

# For C usage
target_link_libraries(your_target PRIVATE find_peaks_c)
target_include_directories(your_target PRIVATE path/to/find_peaks/src/c)
```

## Testing and Validation

The library includes comprehensive test suites to ensure correctness and compatibility with the original SciPy implementation:

### C++ and C Tests

Unit tests are implemented using the Google Test framework and cover:
- Basic peak detection functionality
- Parameter filtering (height, prominence, width, etc.)
- Edge cases (empty signals, flat signals, etc.)
- Compatibility with SciPy reference results

Run the tests with:
```bash
cd build
cmake --build . --target test_find_peaks_c test_find_peaks_cpp

ctest -V
```

Or run specific test suites:
```bash
./test_find_peaks_cpp  # C++ implementation tests
./test_find_peaks_c    # C implementation tests
```

### Python Bindings and Validation

Python bindings allow easy comparison with the reference SciPy implementation. The bindings folder includes:
- Direct Python interfaces to both C and C++ implementations
- Validation tests comparing results against SciPy
- Performance benchmarks comparing execution speed

For details on using the Python bindings, see the [Python Bindings README](bindings/README.md).

## License

This project contains code derived from scipy.signal.find_peaks, which is licensed under the BSD 3-clause license.
The original code is:

Copyright (c) 2001-2002 Enthought, Inc. 2003, SciPy Developers.
Available at: https://github.com/scipy/scipy

This derivative work is distributed under the same BSD 3-clause license. See the LICENSE file for the complete license text.

## Contributing

Contributions are welcome! Please feel free to submit a Pull Request.
