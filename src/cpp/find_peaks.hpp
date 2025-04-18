#ifndef FIND_PEAKS_FIND_PEAKS_HPP
#define FIND_PEAKS_FIND_PEAKS_HPP

#include <string>
#include <sstream>
#include <limits>
#include <vector>

/**
 * @file find_peaks.hpp
 * @brief Header file for peak detection and analysis in signal data
 *
 * This file defines data structures and functions for finding and analyzing
 * peaks in signal data with customizable detection criteria and detailed
 * peak property analysis.
 */

namespace findPeaks {

    /**
     * @struct lmr_peak_index_t
     * @brief Stores the indices defining a peak's position
     *
     * This structure contains the left edge, middle point, and right edge indices
     * of a detected peak, defining its position and extent in the input data.
     */
    typedef struct {
        size_t left_edge;   /**< Index of the leftmost sample belonging to the peak */
        size_t mid_point;   /**< Index of the peak's highest point (or middle of plateau) */
        size_t right_edge;  /**< Index of the rightmost sample belonging to the peak */
    } lmr_peak_index_t;

    /**
     * @struct lpr_peak_prominence_t
     * @brief Stores peak prominence information
     *
     * Prominence quantifies how much a peak stands out from the surrounding baseline.
     * This structure contains the prominence value and the indices of the left and right
     * base points used to calculate it.
     */
    typedef struct {
        size_t left_base;   /**< Index of the left base point used for prominence calculation */
        double prominence;  /**< Calculated prominence value of the peak */
        size_t right_base;  /**< Index of the right base point used for prominence calculation */
    } lpr_peak_prominence_t;

    /**
     * @struct whlr_peak_width_t
     * @brief Stores peak width information
     *
     * Width characterizes the horizontal extent of a peak at a specific height.
     * This structure contains the calculated width and the interpolated positions
     * where the peak crosses the width height level.
     */
    typedef struct {
        double width;        /**< Calculated width of the peak */
        double width_height; /**< Height level at which width is measured */
        double left_ip;      /**< Interpolated position of the left width crossing point */
        double right_ip;     /**< Interpolated position of the right width crossing point */
    } whlr_peak_width_t;

    /**
     * @struct lr_peak_threshold_t
     * @brief Stores peak threshold information
     *
     * Threshold represents how much a peak rises above its neighboring valleys.
     * This structure contains the threshold values for the left and right sides of the peak.
     */
    typedef struct {
        double left_threshold;  /**< Height difference between peak and left neighbor valley */
        double right_threshold; /**< Height difference between peak and right neighbor valley */
    } lr_peak_threshold_t;

    /**
     * @struct lpr_peak_plateau_t
     * @brief Stores information about flat peak plateaus
     *
     * Some peaks may have a flat top (plateau) rather than a single highest point.
     * This structure contains information about the plateau size and its edges.
     */
    typedef struct {
        size_t plateau_size; /**< Number of samples in the peak's plateau */
        size_t left_edge;    /**< Index of the leftmost sample in the plateau */
        size_t right_edge;   /**< Index of the rightmost sample in the plateau */
    } lpr_peak_plateau_t;

    /**
     * @struct peak_result_t
     * @brief Complete peak information structure
     *
     * This structure aggregates all information about a detected peak,
     * including its position, height, and various characteristics like
     * prominence, width, threshold, and plateau information.
     */
    typedef struct {
        size_t peak;         /**< Index of the peak in the input data */
        double peak_height;  /**< Height (value) of the peak */

        lpr_peak_plateau_t plateau;      /**< Information about the peak's plateau */
        lr_peak_threshold_t threshold;    /**< Threshold values for the peak */
        lpr_peak_prominence_t prominence; /**< Prominence information for the peak */
        whlr_peak_width_t width;         /**< Width information for the peak */

    } peak_result_t;

    /**
     * @class RangeFloat
     * @brief Represents a range of floating-point values with lower and upper bounds
     *
     * This class allows creating flexible numeric ranges for floating-point values,
     * supporting various initialization methods and providing string conversion.
     */
    class RangeFloat {
    public:
        /**
         * @brief Indicates whether the range was explicitly initialized
         *
         * Note: This field is typically used for tracking range configuration status
         */
        bool initialized;

        /**
         * @brief Lower bound of the numeric range
         */
        double lower;

        /**
         * @brief Upper bound of the numeric range
         */
        double upper;

        /**
         * @brief Default constructor
         * @details Creates a range spanning the entire possible double value range
         * Sets lower bound to minimum double and upper bound to maximum double
         */
        RangeFloat() {
            this->lower = -std::numeric_limits<double>::max();//std::numeric_limits<double>::min();
            this->upper = std::numeric_limits<double>::max();
        }

        /**
         * @brief Copy constructor
         * @param other Another RangeFloat to copy from
         * @details Creates a new RangeFloat with the same lower and upper bounds
         */
        RangeFloat(const RangeFloat &other) {
            this->lower = other.lower;
            this->upper = other.upper;
        }

        /**
         * @brief Constructor with only lower bound
         * @param lower The lower bound of the range
         * @details Sets the lower bound to the specified value and upper bound to maximum double
         */
        RangeFloat(double lower) {
            this->lower = lower;
            this->upper = std::numeric_limits<double>::max();
        }

        /**
         * @brief Constructor with both lower and upper bounds
         * @param lower The lower bound of the range
         * @param upper The upper bound of the range
         * @details Creates a range with specified lower and upper bounds
         */
        RangeFloat(double lower, double upper) {
            this->lower = lower;
            this->upper = upper;
        }

        /**
         * @brief Converts the range to a string representation
         * @return A string in the format "(lower, upper)"
         *
         * This operator allows easy string conversion for debugging or display purposes
         */
        explicit operator std::string() const {
            std::ostringstream ss;
            ss << "(" << lower << ", " << upper << ")";
            return ss.str();
        }
    };

    /**
     * @class RangeInt
     * @brief Represents a range of integer values (size_t) with lower and upper bounds
     *
     * Similar to RangeFloat, this class provides flexible numeric ranges for integer values,
     * supporting various initialization methods and string conversion.
     */
    class RangeInt {
    public:
        /**
         * @brief Lower bound of the numeric range
         */
        size_t lower;

        /**
         * @brief Upper bound of the numeric range
         */
        size_t upper;

        /**
         * @brief Default constructor
         * @details Creates a range from 0 to maximum possible size_t value
         */
        RangeInt() {
            this->lower = 0;
            this->upper = std::numeric_limits<size_t>::max();
        }

        /**
         * @brief Copy constructor
         * @param other Another RangeInt to copy from
         * @details Creates a new RangeInt with the same lower and upper bounds
         */
        RangeInt(const RangeInt &other) {
            this->lower = other.lower;
            this->upper = other.upper;
        }

        /**
         * @brief Constructor with only lower bound
         * @param lower The lower bound of the range
         * @details Sets the lower bound to the specified value and upper bound to maximum size_t
         */
        RangeInt(size_t lower) {
            this->lower = lower;
            this->upper = std::numeric_limits<size_t>::max();
        }

        /**
         * @brief Constructor with both lower and upper bounds
         * @param lower The lower bound of the range
         * @param upper The upper bound of the range
         * @details Creates a range with specified lower and upper bounds
         */
        RangeInt(size_t lower, size_t upper) {
            this->lower = lower;
            this->upper = upper;
        }

        /**
         * @brief Converts the range to a string representation
         * @return A string in the format "(lower, upper)"
         *
         * This operator allows easy string conversion for debugging or display purposes
         */
        explicit operator std::string() const {
            std::ostringstream ss;
            ss << "(" << lower << ", " << upper << ")";
            return ss.str();
        }
    };

    /**
     * @class PeakConditions
     * @brief Configures comprehensive criteria for peak detection in signal processing
     *
     * This class allows detailed customization of peak detection parameters,
     * supporting flexible configuration through various setter methods.
     */
    class PeakConditions {
    public:
        /**
         * @brief Range of acceptable peak heights
         */
        RangeFloat height;

        /**
         * @brief Range of acceptable peak detection thresholds
         */
        RangeFloat threshold;

        /**
         * @brief Minimum horizontal distance between neighboring peaks
         */
        size_t distance;

        /**
         * @brief Range of acceptable peak prominences
         */
        RangeFloat prominence;

        /**
         * @brief Range of acceptable peak widths
         */
        RangeFloat width;

        /**
         * @brief Window length for prominence calculation
         */
        size_t wlen;

        /**
         * @brief Relative height at which peak width is measured (0.0-1.0)
         */
        double rel_height;

        /**
         * @brief Range of acceptable plateau sizes
         */
        RangeInt plateau_size;

        /**
         * @brief Comprehensive constructor for peak detection conditions
         * @param height Range of peak heights (default: entire range)
         * @param threshold Range of peak detection thresholds (default: entire range)
         * @param distance Minimum horizontal distance between peaks (default: 1)
         * @param prominence Range of peak prominences (default: entire range)
         * @param width Range of peak widths (default: entire range)
         * @param wlen Window length for prominence calculation (default: 0)
         * @param rel_height Relative height for width measurement (default: 0.5)
         * @param plateau_size Range of acceptable plateau sizes (default: entire range)
         *
         * Allows flexible initialization of peak detection parameters with sensible defaults
         */
        PeakConditions(const RangeFloat &height = RangeFloat(),
                       const RangeFloat &threshold = RangeFloat(),
                       size_t distance = 1,
                       const RangeFloat &prominence = RangeFloat(),
                       const RangeFloat &width = RangeFloat(),
                       size_t wlen = 0,
                       double rel_height = 0.5,
                       const RangeInt &plateau_size = RangeInt()) :
                height(height), threshold(threshold),
                distance(distance), prominence(prominence),
                width(width), wlen(wlen),
                rel_height(rel_height), plateau_size(plateau_size) {}

        /**
         * @brief Sets peak height range using variadic template
         * @tparam Args Argument types for RangeFloat construction
         * @param args Arguments to construct RangeFloat
         *
         * Supports flexible height range configuration
         */
        template<class... Args>
        void set_height(Args &&... args) { this->height = RangeFloat(std::forward<Args>(args)...); }

        /**
         * @brief Sets peak threshold range using variadic template
         * @tparam Args Argument types for RangeFloat construction
         * @param args Arguments to construct RangeFloat
         *
         * Supports flexible threshold range configuration
         */
        template<class... Args>
        void set_threshold(Args &&... args) { this->threshold = RangeFloat(std::forward<Args>(args)...); }

        /**
         * @brief Sets minimum horizontal distance between peaks
         * @param distance Number of samples between peaks
         */
        void set_distance(size_t distance) { this->distance = distance; }

        /**
         * @brief Sets peak prominence range using variadic template
         * @tparam Args Argument types for RangeFloat construction
         * @param args Arguments to construct RangeFloat
         *
         * Supports flexible prominence range configuration
         */
        template<class... Args>
        void set_prominence(Args &&... args) { this->prominence = RangeFloat(std::forward<Args>(args)...); }

        /**
         * @brief Sets peak width range using variadic template
         * @tparam Args Argument types for RangeFloat construction
         * @param args Arguments to construct RangeFloat
         *
         * Supports flexible width range configuration
         */
        template<class... Args>
        void set_width(Args &&... args) { this->width = RangeFloat(std::forward<Args>(args)...); }

        /**
         * @brief Sets window length for prominence calculation
         * @param wlen Window length in samples
         */
        void set_wlen(size_t wlen) { this->wlen = wlen; }

        /**
         * @brief Sets relative height for width measurement
         * @param rel_height Relative height as a fraction between 0 and 1
         */
        void set_rel_height(double rel_height) { this->rel_height = rel_height; }

        /**
         * @brief Sets plateau size range using variadic template
         * @tparam Args Argument types for RangeInt construction
         * @param args Arguments to construct RangeInt
         *
         * Supports flexible plateau size range configuration
         */
        template<class... Args>
        void set_plateau_size(Args &&... args) { this->plateau_size = RangeInt(std::forward<Args>(args)...); }
    };

    /**
     * @brief Finds and analyzes peaks in a numerical signal
     *
     * This function detects peaks in a given signal based on configurable criteria,
     * providing comprehensive peak analysis including position, height, prominence,
     * width, and plateau characteristics.
     *
     * @param x Input signal represented as a vector of double-precision values
     * @param conditions Peak detection criteria (optional)
     *
     * @return Vector of peak_result_t structures containing detailed peak information
     *
     * @details
     * The function performs sophisticated peak detection with the following key features:
     * - Flexible peak detection criteria through PeakConditions
     * - Comprehensive peak property analysis
     * - Efficient O(nlog(n)) time complexity algorithm
     *
     * Detection Criteria:
     * - Height range
     * - Threshold values
     * - Minimum distance between peaks
     * - Prominence range
     * - Width range
     * - Plateau size
     *
     * @note
     * - Returns an empty vector if no peaks meet the specified conditions
     * - Peaks are returned in order of appearance in the input signal
     *
     * @complexity O(nlog(n)), where n is the length of the input signal
     *
     * @throws No explicit exceptions, but may have undefined behavior
     *         with extremely large or malformed input signals
     *
     * @see PeakConditions
     * @see peak_result_t
     *
     * @par Example
     * @code{.cpp}
     * // Basic peak detection
     * std::vector<double> signal = {1.0, 2.5, 1.7, 3.2, 2.1, 1.8};
     * auto peaks = find_peaks(signal);
     *
     * // Customized peak detection
     * PeakConditions conditions;
     * conditions.set_height(2.0, 4.0);    // Peaks between 2.0 and 4.0
     * conditions.set_distance(2);          // At least 2 samples between peaks
     * conditions.set_prominence(0.5);      // Minimum prominence of 0.5
     *
     * auto filtered_peaks = find_peaks(signal, conditions);
     * @endcode
     *
     * @warning
     * - Always validate peak detection results visually or statistically
     * - Adjust detection parameters based on signal characteristics
     *
     * @par Potential Applications
     * - Signal processing
     * - Sensor data analysis
     * - Biomedical signal interpretation
     * - Financial time series analysis
     * - Vibration and oscillation studies
     */
    std::vector<peak_result_t> find_peaks(
            const std::vector<double> &x,                      ///< Input signal vector
            PeakConditions conditions = PeakConditions()       ///< Peak detection conditions (optional)
    );

}

#endif //FIND_PEAKS_FIND_PEAKS_HPP
