"""
Python wrapper for C++ and C implementations of peak finding algorithms.

This module provides a Python interface to optimized peak finding implementations
in C++ and C. The functions find local maxima (peaks) in sequential data and
provide optional filtering capabilities based on various properties like height,
prominence, width, etc.

Functions:
----------
find_peaks_cpp : Find peaks in data using the C++ implementation
find_peaks_c : Find peaks in data using the C implementation
_unpack_condition_args : Helper function to parse condition arguments

Notes:
------
The C++ and C implementations provide significant performance improvements
over pure Python implementations while maintaining the same interface.
"""

import sys

import find_peaks_cpp_py as fp_cpp  # C++ implementation bindings
import find_peaks_c_py as fp_c      # C implementation bindings
import numpy as np

def _unpack_condition_args(interval, isint = False):
    """
    Parse condition arguments for peak finding functions.

    Parameters
    ----------
    interval : number or ndarray or sequence
        Either a number, ndarray, or a 2-element sequence. The first value is 
        interpreted as `imin` and the second, if supplied, as `imax`.
    isint : bool, optional
        If True, the default range will use integer limits, otherwise
        float limits will be used. Default is False.

    Returns
    -------
    imin, imax : number or ndarray or None
        Minimal and maximal values derived from the interval argument.

    Notes
    -----
    This is a helper function used internally by the peak finding functions
    to standardize how interval conditions are specified.
    """

    if interval is None:
        # If no interval is provided, set default min/max bounds
        if isint:
            interval = (0, sys.maxsize)  # Use full range for integer values
        else:
            interval = (-sys.float_info.max, sys.float_info.max)  # Use full range for float values

    try:
        # Try to unpack interval as a pair of values
        imin, imax = interval
    except (TypeError, ValueError):
        # If interval is a single value, use it as the minimum and set a default maximum
        if isint:
            imin, imax = (interval, sys.maxsize)  # For integer parameters
        else:
            # Using a large finite value instead of infinity to avoid floating point issues
            imin, imax = (interval, 1000000)  # For float parameters

    return imin, imax


def find_peaks_cpp(x, height=None, threshold=None, distance=None,
               prominence=None, width=None, wlen=None, rel_height=0.5,
               plateau_size=None):
    """
    Find peaks in data using the C++ implementation.

    Finds local maxima in `x` whose properties match the specified conditions.

    Parameters
    ----------
    x : ndarray
        A signal with peaks.
    height : number or ndarray or sequence, optional
        Required height of peaks. Either a number, None, an array matching `x` or
        a 2-element sequence of the former. The first element is always interpreted
        as the minimal and the second, if supplied, as the maximal required height.
    threshold : number or ndarray or sequence, optional
        Required threshold of peaks, the vertical distance to its neighboring
        samples. Either a number, None, an array matching `x` or a 2-element
        sequence of the former. The first element is always interpreted as the
        minimal and the second, if supplied, as the maximal required threshold.
    distance : int, optional
        Required minimal horizontal distance (>= 1) in samples between
        neighbouring peaks. The smaller peaks are removed first until the condition
        is fulfilled for all remaining peaks. Default is 0 (no distance constraint).
    prominence : number or ndarray or sequence, optional
        Required prominence of peaks. Either a number, None, an array matching `x`
        or a 2-element sequence of the former. The first element is always
        interpreted as the minimal and the second, if supplied, as the maximal
        required prominence.
    width : number or ndarray or sequence, optional
        Required width of peaks in samples. Either a number, None, an array
        matching `x` or a 2-element sequence of the former. The first element is
        always interpreted as the minimal and the second, if supplied, as the
        maximal required width.
    wlen : int, optional
        Used for calculation of the peaks prominences, thus it is only used if one
        of the arguments `prominence` or `width` is given. Default is 0 (uses entire
        signal).
    rel_height : float, optional
        Used for calculation of the peaks width, thus it is only used if `width`
        is given. Default is 0.5.
    plateau_size : number or ndarray or sequence, optional
        Required size of the flat top of peaks in samples. Either a number, None,
        an array matching `x` or a 2-element sequence of the former. The first
        element is always interpreted as the minimal and the second, if supplied,
        as the maximal required plateau size.

    Returns
    -------
    peaks : ndarray
        Indices of peaks in `x` that satisfy all given conditions.
    properties : dict
        A dictionary containing properties of the detected peaks:
        - 'peak_heights': Height of each peak.
        - 'left_thresholds', 'right_thresholds': Vertical distance to the
          neighbouring samples.
        - 'prominences': Prominence of each peak.
        - 'left_bases', 'right_bases': Indices of bases of prominence.
        - 'width_heights': Height at which the peak width is measured.
        - 'left_ips', 'right_ips': x-coordinates of left and right intersection
          points at the respective peak width height.
        - 'plateau_sizes': Size of the flat top of each peak.

    Examples
    --------
    >>> import numpy as np
    >>> from find_peaks_wrapper import find_peaks_cpp
    >>> x = np.array([0, 1, 0, 2, 0, 3, 0, 2, 0, 1, 0])
    >>> peaks, properties = find_peaks_cpp(x, height=1.5)
    >>> peaks
    array([3, 5, 7])
    """

    # Convert all filtering condition arguments to min/max pairs
    height = _unpack_condition_args(height)          # Peak height constraints
    threshold = _unpack_condition_args(threshold)    # Threshold constraints
    prominence = _unpack_condition_args(prominence)  # Prominence constraints
    width = _unpack_condition_args(width)           # Width constraints
    plateau_size = _unpack_condition_args(plateau_size, True)  # Plateau size (integer) constraints

    # Set defaults for null parameters
    if distance is None:
        distance = 0  # No distance filtering by default

    if wlen is None:
        wlen = 0  # Use the entire signal for prominence calculation by default

    # Call the C++ implementation with the processed parameters
    result = fp_cpp._find_peaks(x, height, threshold, distance, prominence, width, wlen, rel_height, plateau_size)

    # Reorganize the results: result[0] is peak indices, result[1] is a list of property dicts per peak
    d = {}
    for r in result[1]:  # Iterate through each peak's properties
        for key in r:    # Iterate through each property type
            if key not in d:
                d[key] = []  # Initialize a list for this property if not already present
            d[key] += [r[key]]  # Append this peak's property value

    # Convert property lists to numpy arrays for consistency
    for key in d:
        d[key] = np.array(d[key])

    return result[0], d  # Return peak indices and property dictionary

def find_peaks_c(x, height=None, threshold=None, distance=None,
                   prominence=None, width=None, wlen=None, rel_height=0.5,
                   plateau_size=None):
    """
    Find peaks in data using the C implementation.

    Finds local maxima in `x` whose properties match the specified conditions.
    This is functionally identical to find_peaks_cpp but uses the C implementation.

    Parameters
    ----------
    x : ndarray
        A signal with peaks.
    height : number or ndarray or sequence, optional
        Required height of peaks. Either a number, None, an array matching `x` or
        a 2-element sequence of the former. The first element is always interpreted
        as the minimal and the second, if supplied, as the maximal required height.
    threshold : number or ndarray or sequence, optional
        Required threshold of peaks, the vertical distance to its neighboring
        samples. Either a number, None, an array matching `x` or a 2-element
        sequence of the former. The first element is always interpreted as the
        minimal and the second, if supplied, as the maximal required threshold.
    distance : int, optional
        Required minimal horizontal distance (>= 1) in samples between
        neighbouring peaks. The smaller peaks are removed first until the condition
        is fulfilled for all remaining peaks. Default is 0 (no distance constraint).
    prominence : number or ndarray or sequence, optional
        Required prominence of peaks. Either a number, None, an array matching `x`
        or a 2-element sequence of the former. The first element is always
        interpreted as the minimal and the second, if supplied, as the maximal
        required prominence.
    width : number or ndarray or sequence, optional
        Required width of peaks in samples. Either a number, None, an array
        matching `x` or a 2-element sequence of the former. The first element is
        always interpreted as the minimal and the second, if supplied, as the
        maximal required width.
    wlen : int, optional
        Used for calculation of the peaks prominences, thus it is only used if one
        of the arguments `prominence` or `width` is given. Default is 0 (uses entire
        signal).
    rel_height : float, optional
        Used for calculation of the peaks width, thus it is only used if `width`
        is given. Default is 0.5.
    plateau_size : number or ndarray or sequence, optional
        Required size of the flat top of peaks in samples. Either a number, None,
        an array matching `x` or a 2-element sequence of the former. The first
        element is always interpreted as the minimal and the second, if supplied,
        as the maximal required plateau size.

    Returns
    -------
    peaks : ndarray
        Indices of peaks in `x` that satisfy all given conditions.
    properties : dict
        A dictionary containing properties of the detected peaks:
        - 'peak_heights': Height of each peak.
        - 'left_thresholds', 'right_thresholds': Vertical distance to the
          neighbouring samples.
        - 'prominences': Prominence of each peak.
        - 'left_bases', 'right_bases': Indices of bases of prominence.
        - 'width_heights': Height at which the peak width is measured.
        - 'left_ips', 'right_ips': x-coordinates of left and right intersection
          points at the respective peak width height.
        - 'plateau_sizes': Size of the flat top of each peak.

    Examples
    --------
    >>> import numpy as np
    >>> from find_peaks_wrapper import find_peaks_c
    >>> x = np.array([0, 1, 0, 2, 0, 3, 0, 2, 0, 1, 0])
    >>> peaks, properties = find_peaks_c(x, height=1.5)
    >>> peaks
    array([3, 5, 7])
    """

    # Convert all filtering condition arguments to min/max pairs
    height = _unpack_condition_args(height)          # Peak height constraints
    threshold = _unpack_condition_args(threshold)    # Threshold constraints
    prominence = _unpack_condition_args(prominence)  # Prominence constraints
    width = _unpack_condition_args(width)           # Width constraints
    plateau_size = _unpack_condition_args(plateau_size, True)  # Plateau size (integer) constraints

    # Set defaults for null parameters
    if distance is None:
        distance = 0  # No distance filtering by default

    if wlen is None:
        wlen = 0  # Use the entire signal for prominence calculation by default

    # Call the C implementation with the processed parameters
    result = fp_c._find_peaks(x, height, threshold, distance, prominence, width, wlen, rel_height, plateau_size)

    # Reorganize the results: result[0] is peak indices, result[1] is a list of property dicts per peak
    d = {}
    for r in result[1]:  # Iterate through each peak's properties
        for key in r:    # Iterate through each property type
            if key not in d:
                d[key] = []  # Initialize a list for this property if not already present
            d[key] += [r[key]]  # Append this peak's property value

    # Convert property lists to numpy arrays for consistency
    for key in d:
        d[key] = np.array(d[key])

    return result[0], d  # Return peak indices and property dictionary

