"""
Tests for comparing the scipy.signal.find_peaks function to the C/C++ implementations.
This file contains the Python equivalent of the C++ tests in test_find_peaks_cpp.cpp.
"""
import os

import pytest
import numpy as np
from scipy.signal import find_peaks
import sys

try:
    # Add path to find the Python bindings
    sys.path.insert(0, os.environ['FP_BINDING_PATH'])
    from find_peaks_wrapper import find_peaks_cpp, find_peaks_c
    HAS_CPP_BINDINGS = True
    print('Test Binding!')
except ImportError:
    HAS_CPP_BINDINGS = False
    print('Doesn\'t Test Binding!')

# Test fixtures - sample signals matching those in C++ tests
@pytest.fixture
def simple_signal():
    return np.array([0, 1, 0, 2, 0, 3, 0, 2, 0, 1, 0], dtype=np.float64)

@pytest.fixture
def plateau_signal():
    return np.array([0, 1, 1, 1, 0, 2, 2, 0, 3, 3, 3, 0], dtype=np.float64)

@pytest.fixture
def noisy_signal():
    return np.array([1.2, 0.8, 1.9, 1.5, 2.7, 1.8, 1.1, 2.5, 3.2, 2.8, 1.6, 0.9], dtype=np.float64)

@pytest.fixture
def flat_signal():
    return np.array([1, 1, 1, 1, 1], dtype=np.float64)

@pytest.fixture
def empty_signal():
    return np.array([], dtype=np.float64)

# Utility function to compare peak results
def compare_peaks(scipy_result, cpp_result=None, c_result=None):
    """Compare peak indices and properties between different implementations."""
    if cpp_result is not None:
        assert len(scipy_result[0]) == len(cpp_result[0])
        assert np.all(scipy_result[0] == cpp_result[0])
        cpp_result_dict = {key: cpp_result[1].get(key, []) for key in scipy_result[1]}
        key_comp = scipy_result[1].keys() == cpp_result_dict.keys()
        val_comp = all([np.all(scipy_result[1][key] == cpp_result_dict[key]) for key in scipy_result[1]])
        assert key_comp and val_comp

    if c_result is not None:
        assert len(scipy_result[0]) == len(c_result[0])
        assert np.all(scipy_result[0] == c_result[0])
        c_result_dict = {key: c_result[1].get(key, []) for key in scipy_result[1]}
        key_comp = scipy_result[1].keys() == c_result_dict.keys()
        val_comp = all([np.all(scipy_result[1][key] == c_result_dict[key]) for key in scipy_result[1]])
        assert key_comp and val_comp

# Test basic peak detection with default conditions
def test_basic_peak_detection(simple_signal):
    # With SciPy
    peaks, _ = find_peaks(simple_signal)

    # Check if we have expected peaks
    assert len(peaks) == 5
    assert peaks[0] == 1  # Peak at index 1 with height 1
    assert peaks[1] == 3  # Peak at index 3 with height 2
    assert peaks[2] == 5  # Peak at index 5 with height 3
    assert peaks[3] == 7  # Peak at index 7 with height 3
    assert peaks[4] == 9  # Peak at index 9 with height 3

    # Compare with C/C++ implementations if available
    if HAS_CPP_BINDINGS:
        cpp_peaks = find_peaks_cpp(simple_signal)
        c_peaks = find_peaks_c(simple_signal)
        compare_peaks((peaks, _), cpp_peaks, c_peaks)

# Test detection with height filtering
def test_height_filtering(simple_signal):
    # With SciPy
    peaks, properties = find_peaks(simple_signal, height=2.5)

    # Check if only peaks higher than 2.5 are found
    assert len(peaks) == 1
    assert peaks[0] == 5  # Only the peak with height 3
    assert simple_signal[peaks[0]] == 3.0

    # Compare with C/C++ implementations if available
    if HAS_CPP_BINDINGS:
        cpp_peaks = find_peaks_cpp(simple_signal, height=2.5)
        c_peaks = find_peaks_c(simple_signal, height=2.5)
        compare_peaks((peaks, properties), cpp_peaks, c_peaks)

# Test detection with minimum distance between peaks
def test_distance_filtering(simple_signal):
    # With SciPy
    peaks, properties = find_peaks(simple_signal, distance=3)

    # Check if peaks are at least 3 samples apart
    assert len(peaks) == 3
    assert peaks[0] == 1  # First peak at index 1
    assert peaks[1] == 5  # Second peak at index 5
    assert peaks[2] == 9  # Second peak at index 5

    # Compare with C/C++ implementations if available
    if HAS_CPP_BINDINGS:
        cpp_peaks = find_peaks_cpp(simple_signal, distance=3)
        c_peaks = find_peaks_c(simple_signal, distance=3)
        compare_peaks((peaks, properties), cpp_peaks, c_peaks)

# Test plateau detection
def test_plateau_detection(plateau_signal):
    # With SciPy
    peaks, properties = find_peaks(plateau_signal, plateau_size=1)

    # Check if we found all plateaus
    assert len(peaks) == 3

    # First plateau (middle is at index 2)
    assert peaks[0] == 2
    # Second plateau (middle is at index 5)
    assert peaks[1] == 5
    # Third plateau (middle is at index 9)
    assert peaks[2] == 9

    # Check plateau sizes if available in properties
    if 'plateau_sizes' in properties:
        assert properties['plateau_sizes'][0] == 3  # First plateau size
        assert properties['plateau_sizes'][1] == 2  # Second plateau size
        assert properties['plateau_sizes'][2] == 3  # Third plateau size

    # Compare with C/C++ implementations if available
    if HAS_CPP_BINDINGS:
        cpp_peaks = find_peaks_cpp(plateau_signal, plateau_size=1)
        c_peaks = find_peaks_c(plateau_signal, plateau_size=1)
        compare_peaks((peaks, properties), cpp_peaks, c_peaks)

# Test prominence calculation and filtering
def test_prominence_filtering(simple_signal):
    # With SciPy
    peaks, properties = find_peaks(simple_signal, prominence=1.5)
    # Only the peak with height 3 has enough prominence
    assert len(peaks) > 0
    assert 5 in peaks  # The peak at index 5 with height 3 should be included

    # Check all prominences
    for i, peak in enumerate(peaks):
        assert properties['prominences'][i] >= 1.5

    # Compare with C/C++ implementations if available
    if HAS_CPP_BINDINGS:
        cpp_peaks = find_peaks_cpp(simple_signal, prominence=1.5)
        c_peaks = find_peaks_c(simple_signal, prominence=1.5)
        compare_peaks((peaks, properties), cpp_peaks, c_peaks)

# Test width calculation and filtering
def test_width_filtering(simple_signal):
    # With SciPy
    peaks, properties = find_peaks(simple_signal, width=2.0, rel_height=0.5)

    # Check that all peaks have width >= 2.0
    for i in range(len(peaks)):
        assert properties['widths'][i] >= 2.0

    # Compare with C/C++ implementations if available
    if HAS_CPP_BINDINGS:
        cpp_peaks = find_peaks_cpp(simple_signal, width=2.0, rel_height=0.5)
        c_peaks = find_peaks_c(simple_signal, width=2.0, rel_height=0.5)
        compare_peaks((peaks, properties), cpp_peaks, c_peaks)

# Test threshold filtering
def test_threshold_filtering(simple_signal):
    # With SciPy
    peaks, properties = find_peaks(simple_signal, threshold=1.5)

    # Check that all peaks exceed neighbors by at least 1.5
    if 'left_thresholds' in properties and 'right_thresholds' in properties:
        for i in range(len(peaks)):
            assert properties['left_thresholds'][i] >= 1.5 or properties['right_thresholds'][i] >= 1.5

    # Compare with C/C++ implementations if available
    if HAS_CPP_BINDINGS:
        cpp_peaks = find_peaks_cpp(simple_signal, threshold=1.5)
        c_peaks = find_peaks_c(simple_signal, threshold=1.5)
        compare_peaks((peaks, properties), cpp_peaks, c_peaks)

# Test with noisy data
def test_noisy_signal(noisy_signal):
    # With SciPy
    peaks, _ = find_peaks(noisy_signal)

    # Verify that peaks are correctly identified in noisy signal
    assert len(peaks) > 0

    # Check if the highest peak is detected (index 8, value 3.2)
    assert 8 in peaks
    assert noisy_signal[8] == 3.2

    # Compare with C/C++ implementations if available
    if HAS_CPP_BINDINGS:
        cpp_peaks = find_peaks_cpp(noisy_signal)
        c_peaks = find_peaks_c(noisy_signal)
        compare_peaks((peaks, _), cpp_peaks, c_peaks)

# Test with flat signal (no peaks)
def test_flat_signal(flat_signal):
    # With SciPy
    peaks, _ = find_peaks(flat_signal)
    assert len(peaks) == 0

    # Compare with C/C++ implementations if available
    if HAS_CPP_BINDINGS:
        cpp_peaks = find_peaks_cpp(flat_signal)
        c_peaks = find_peaks_c(flat_signal)
        compare_peaks((peaks, _), cpp_peaks, c_peaks)

# Test with empty signal
def test_empty_signal(empty_signal):
    # With SciPy
    peaks, _ = find_peaks(empty_signal)
    assert len(peaks) == 0

    # Compare with C/C++ implementations if available
    if HAS_CPP_BINDINGS:
        cpp_peaks = find_peaks_cpp(empty_signal)
        c_peaks = find_peaks_c(empty_signal)
        compare_peaks((peaks, _), cpp_peaks, c_peaks)

# Test combined filtering criteria
def test_combined_filters(simple_signal):
    # With SciPy
    peaks, properties = find_peaks(
        simple_signal,
        height=1.5,          # Height >= 1.5
        prominence=1.0,      # Prominence >= 1.0
        distance=2,          # At least 2 samples between peaks
        width=(1.0, 4.0)     # Width between 1.0 and 4.0
    )

    # Check that all peaks satisfy all conditions
    for i, peak in enumerate(peaks):
        assert simple_signal[peak] >= 1.5  # Height condition
        assert properties['prominences'][i] >= 1.0  # Prominence condition
        assert properties['widths'][i] >= 1.0  # Width lower bound
        assert properties['widths'][i] <= 4.0  # Width upper bound

    # Make sure peaks are at least 2 samples apart
    for i in range(1, len(peaks)):
        assert peaks[i] - peaks[i-1] >= 2

    # Compare with C/C++ implementations if available
    if HAS_CPP_BINDINGS:
        cpp_peaks = find_peaks_cpp(
            simple_signal,
            height=1.5,
            prominence=1.0,
            distance=2,
            width=(1.0, 4.0)
        )
        c_peaks = find_peaks_c(
            simple_signal,
            height=1.5,
            prominence=1.0,
            distance=2,
            width=(1.0, 4.0)
        )
        compare_peaks((peaks, properties), cpp_peaks, c_peaks)

# Test against specific known signal for compatibility validation
def test_scipy_compatibility():
    # This is a known test signal with peaks that match SciPy's find_peaks results
    scipy_signal = np.array([0, 1, 0, 2, 0, 3, 0, 2, 0, 1, 0])

    # These are the expected peak indices from SciPy's find_peaks with height=0.5
    expected_peaks = np.array([1, 3, 5, 7, 9])

    # With SciPy
    peaks, _ = find_peaks(scipy_signal, height=0.5)

    # Check if we got the expected peaks
    assert len(peaks) == len(expected_peaks)
    assert np.all(peaks == expected_peaks)

    # Compare with C/C++ implementations if available
    if HAS_CPP_BINDINGS:
        cpp_peaks = find_peaks_cpp(scipy_signal, height=0.5)
        c_peaks = find_peaks_c(scipy_signal, height=0.5)
        compare_peaks((peaks, _), cpp_peaks, c_peaks)

if __name__ == "__main__":
    pytest.main()
