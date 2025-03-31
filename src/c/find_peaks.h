#ifndef FIND_PEAKS_FIND_PEAKS_H
#define FIND_PEAKS_FIND_PEAKS_H

/** ---------------------------------------------------------------------------------------
 * This section contains utility functions to configure conditions used by the find_peaks
 * algorithm. These functions manipulate range structures and condition parameters that
 * control peak detection behavior.
 */

/**
 * @struct fp_double_range_t
 * @brief Represent a range of double values
 */
typedef struct{
    double min;
    double max;
} fp_double_range_t;

/**
 * @struct fp_int_range_t
 * @brief Represent a range of integer values
 */
typedef struct{
    size_t min;
    size_t max;
} fp_int_range_t;


/**
 * @struct fp_conditions_t
 * @brief Configuration parameters for peak detection
 *
 * This structure contains all the configurable parameters that control which
 * features in the data are identified as peaks and which are filtered out.
 * Each field represents a different criterion for peak detection and filtering.
 */
typedef struct{

    /**
     * @brief Height range for peak filtering
     *
     * Only peaks with heights within this range will be detected.
     * Default range is [-DBL_MAX, DBL_MAX] (all peaks regardless of height).
     */
    fp_double_range_t height;

    /**
     * @brief Threshold range for peak filtering
     *
     * Only peaks that rise above their neighboring valleys by an amount within
     * this range will be detected. This measures how distinct a peak is from its
     * immediate surroundings.
     * Default range is [-DBL_MAX, DBL_MAX] (all peaks regardless of threshold).
     */
    fp_double_range_t threshold;

    /**
     * @brief Minimum distance between peaks
     *
     * Ensures peaks are separated by at least this many samples. When multiple
     * peaks are found within this distance, only the highest one is kept.
     * Default value is 1 (adjacent peaks allowed).
     */
    size_t distance;

    /**
     * @brief Prominence range for peak filtering
     *
     * Only peaks with prominence values within this range will be detected.
     * Prominence measures how much a peak stands out from its surrounding baseline.
     * Default range is [-DBL_MAX, DBL_MAX] (all peaks regardless of prominence).
     */
    fp_double_range_t prominence;

    /**
     * @brief Width range for peak filtering
     *
     * Only peaks with widths within this range will be detected. Width is
     * measured at a height determined by rel_height.
     * Default range is [-DBL_MAX, DBL_MAX] (all peaks regardless of width).
     */
    fp_double_range_t width;

    /**
     * @brief Window length for prominence and width calculations
     *
     * Used to limit the evaluated area for prominence and width calculations.
     * Default value is 0 (use the full data extent).
     */
    size_t wlen;

    /**
     * @brief Relative height for width calculation
     *
     * Determines the height level at which peak width is measured, as a proportion
     * of the peak height. For example, 0.5 means width is measured at half the
     * peak's height above its base.
     * Default value is 0.5 (half height).
     */
    double rel_height;

    /**
     * @brief Plateau size range for peak filtering
     *
     * Only peaks with plateau sizes within this range will be detected.
     * Plateau size is the number of consecutive samples at the peak's maximum value.
     * Default range is [0, SIZE_MAX] (all peaks regardless of plateau size).
     */
    fp_int_range_t plateau_size;
} fp_conditions_t;


/**
 * @brief Set minimum and maximum values for a double range
 *
 * @param range Pointer to a fp_double_range_t structure to be modified
 * @param min Minimum value for the range
 * @param max Maximum value for the range
 */
void fp_set_double_range_mn_mx(fp_double_range_t* range, double min, double max);

/**
 * @brief Set minimum value for a double range (maximum set to DBL_MAX)
 *
 * @param range Pointer to a fp_double_range_t structure to be modified
 * @param min Minimum value for the range
 */
void fp_set_double_range_mn(fp_double_range_t* range, double min);

/**
 * @brief Initialize a double range with default values (-DBL_MAX to DBL_MAX)
 *
 * @param range Pointer to a fp_double_range_t structure to be initialized
 */
void fp_init_double_range(fp_double_range_t* range);


/**
 * @brief Set minimum and maximum values for an integer range
 *
 * @param range Pointer to a fp_int_range_t structure to be modified
 * @param min Minimum value for the range
 * @param max Maximum value for the range
 */
void fp_set_int_range_mn_mx(fp_int_range_t* range, size_t min, size_t max);

/**
 * @brief Set minimum value for an integer range (maximum set to SIZE_MAX)
 *
 * @param range Pointer to a fp_int_range_t structure to be modified
 * @param min Minimum value for the range
 */
void fp_set_int_range_mn(fp_int_range_t* range, size_t min);

/**
 * @brief Initialize an integer range with default values (0 to SIZE_MAX)
 *
 * @param range Pointer to a fp_int_range_t structure to be initialized
 */
void fp_init_int_range(fp_int_range_t* range);


/**
 * @brief Set minimum and maximum height condition for peak detection
 *
 * Height condition defines the absolute height range a peak must have to be detected.
 *
 * @param cond Pointer to a fp_conditions_t structure to be modified
 * @param min Minimum peak height to be considered
 * @param max Maximum peak height to be considered
 */
void fp_cond_set_height_mn_mx(fp_conditions_t *cond, double min, double max);

/**
 * @brief Set minimum height condition for peak detection (maximum set to DBL_MAX)
 *
 * Height condition defines the absolute height range a peak must have to be detected.
 *
 * @param cond Pointer to a fp_conditions_t structure to be modified
 * @param min Minimum peak height to be considered
 */
void fp_cond_set_height_mn(fp_conditions_t *cond, double min);

/**
 * @brief Initialize height condition for peak detection to default range (-DBL_MAX to DBL_MAX)
 *
 * Height condition defines the absolute height range a peak must have to be detected.
 *
 * @param cond Pointer to a fp_conditions_t structure to be modified
 */
void fp_cond_init_height(fp_conditions_t *cond);


/**
 * @brief Set minimum and maximum threshold condition for peak detection
 *
 * Threshold condition defines how much a data point needs to exceed its neighbors to be considered a peak.
 *
 * @param cond Pointer to a fp_conditions_t structure to be modified
 * @param min Minimum threshold value
 * @param max Maximum threshold value
 */
void fp_cond_set_threshold_mn_mx(fp_conditions_t *cond, double min, double max);


/**
 * @brief Set minimum threshold condition for peak detection (maximum set to DBL_MAX)
 *
 * Threshold condition defines how much a data point needs to exceed its neighbors to be considered a peak.
 *
 * @param cond Pointer to a fp_conditions_t structure to be modified
 * @param min Minimum threshold value
 */
void fp_cond_set_threshold_mn(fp_conditions_t *cond, double min);

/**
 * @brief Initialize threshold condition to default range (-DBL_MAX to DBL_MAX)
 *
 * Threshold condition defines how much a data point needs to exceed its neighbors to be considered a peak.
 *
 * @param cond Pointer to a fp_conditions_t structure to be modified
 */
void fp_cond_init_threshold(fp_conditions_t *cond);


/**
 * @brief Set minimum distance between detected peaks
 *
 * Distance condition ensures peaks are separated by at least the specified number of samples.
 *
 * @param cond Pointer to a fp_conditions_t structure to be modified
 * @param distance Minimum samples between peaks
 */
void fp_cond_set_distance_mn(fp_conditions_t *cond, size_t distance);

/**
 * @brief Initialize distance condition to default value (1)
 *
 * Distance condition ensures peaks are separated by at least the specified number of samples.
 *
 * @param cond Pointer to a fp_conditions_t structure to be modified
 */
void fp_cond_init_distance(fp_conditions_t *cond);

/**
 * @brief Set minimum and maximum prominence condition for peak detection
 *
 * Prominence quantifies how much a peak stands out from surrounding baseline.
 *
 * @param cond Pointer to a fp_conditions_t structure to be modified
 * @param min Minimum peak prominence
 * @param max Maximum peak prominence
 */
void fp_cond_set_prominence_mn_mx(fp_conditions_t *cond, double min, double max);

/**
 * @brief Set minimum prominence condition for peak detection (maximum set to DBL_MAX)
 *
 * Prominence quantifies how much a peak stands out from surrounding baseline.
 *
 * @param cond Pointer to a fp_conditions_t structure to be modified
 * @param min Minimum peak prominence
 */
void fp_cond_set_prominence_mn(fp_conditions_t *cond, double min);

/**
 * @brief Initialize prominence condition to default range (-DBL_MAX to DBL_MAX)
 *
 * Prominence quantifies how much a peak stands out from surrounding baseline.
 *
 * @param cond Pointer to a fp_conditions_t structure to be modified
 */
void fp_cond_init_prominence(fp_conditions_t *cond);

/**
 * @brief Set minimum and maximum width condition for peak detection
 *
 * Defines the Required width condition of peaks in samples.
 *
 * @param cond Pointer to a fp_conditions_t structure to be modified
 * @param min Minimum peak width
 * @param max Maximum peak width
 */
void fp_cond_set_width_mn_mx(fp_conditions_t *cond, double min, double max);

/**
 * @brief Set minimum width condition for peak detection (maximum set to DBL_MAX)
 *
 * Defines the Required width condition of peaks in samples.
 *
 * @param cond Pointer to a fp_conditions_t structure to be modified
 * @param min Minimum peak width
 */
void fp_cond_set_width_mn(fp_conditions_t *cond, double min);


/**
 * @brief Initialize width condition to default range (-DBL_MAX to DBL_MAX)
 *
 * Defines the Required width condition of peaks in samples.
 *
 * @param cond Pointer to a fp_conditions_t structure to be modified
 */
void fp_cond_init_width(fp_conditions_t *cond);

/**
 * @brief Set a window length in samples that optionally limits the evaluated area for each peak
 *
 * Defines the size of the window used for calculating peak widths.
 * The peak is always placed in the middle of the window therefore the given length is rounded up to the next odd integer
 *
 * @param cond Pointer to a fp_conditions_t structure to be modified
 * @param wlen Window length in samples
 */
void fp_cond_set_wlen_mn(fp_conditions_t *cond, size_t wlen);

/**
 * @brief Initialize window length (wlen) that limits the evaluated area for each peak to default value (0)
 *
 * Defines the size of the window used for calculating peak widths.
 * The peak is always placed in the middle of the window therefore the given length is rounded up to the next odd integer
 *
 * @param cond Pointer to a fp_conditions_t structure to be modified
 */
void fp_cond_init_wlen(fp_conditions_t *cond);

/**
 * @brief Set relative height for peak width calculation
 *
 * Defines the relative height at which the peak width is measured as a percentage of its prominence.
 * 1.0 calculates the width of the peak at its lowest contour line while 0.5 evaluates at half the prominence height.
 * Must be at least 0
 *
 * @param cond Pointer to a fp_conditions_t structure to be modified
 * @param rel_height Relative height (0.0 to 1.0, where 1.0 is the peak height)
 */
void fp_cond_set_rel_height_mn(fp_conditions_t *cond, double rel_height);

/**
 * @brief Initialize relative height for peak width calculation to default value (0.5)
 *
 * Defines the relative height at which the peak width is measured as a percentage of its prominence.
 * 1.0 calculates the width of the peak at its lowest contour line while 0.5 evaluates at half the prominence height.
 * Must be at least 0
 *
 * @param cond Pointer to a fp_conditions_t structure to be modified
 */
void fp_cond_init_rel_height(fp_conditions_t *cond);

/**
 * @brief Set minimum and maximum plateau size condition for peak detection
 *
 * Plateau size defines the acceptable range for the number of samples at a peak's maximum value.
 *
 * @param cond Pointer to a fp_conditions_t structure to be modified
 * @param min Minimum number of samples at peak level
 * @param max Maximum number of samples at peak level
 */
void fp_cond_set_plateau_size_mn_mx(fp_conditions_t *cond, size_t min, size_t max);

/**
 * @brief Set minimum plateau size condition for peak detection (maximum set to SIZE_MAX)
 *
 * Plateau size defines the acceptable range for the number of samples at a peak's maximum value.
 *
 * @param cond Pointer to a fp_conditions_t structure to be modified
 * @param min Minimum number of samples at peak level
 */
void fp_cond_set_plateau_size_mn(fp_conditions_t *cond, size_t min);

/**
 * @brief Initialize plateau size condition to default range (0 to SIZE_MAX)
 *
 * Plateau size defines the acceptable range for the number of samples at a peak's maximum value.
 *
 * @param cond Pointer to a fp_conditions_t structure to be modified
 */
void fp_cond_init_plateau_size(fp_conditions_t *cond);


/**
 * @brief Initialize all conditions in a conditions structure to their default values
 *
 * This function initializes all peak detection parameters to their default values,
 * providing a clean starting point for condition configuration.
 *
 * @param cond Pointer to a fp_conditions_t structure to be initialized
 */
void fp_init_conditions(fp_conditions_t* cond);

/**
 * @brief Get a new conditions structure with all default values
 *
 * Creates and returns a new fp_conditions_t structure with all conditions
 * initialized to their default values.
 *
 * @return A new fp_conditions_t structure with default values
 */
fp_conditions_t fp_get_default_conditions();



/** ---------------------------------------------------------------------------
 * This section defines the data structures used by the find_peaks algorithm
 * to represent peak information and results.
 */


/**
 * @struct lmr_peak_index_t
 * @brief Stores the indices defining a peak's position
 *
 * This structure contains the left edge, middle point, and right edge indices
 * of a detected peak, defining its position and extent in the input data.
 */
typedef struct{
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
typedef struct{
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
typedef struct{
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
typedef struct{
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
typedef struct{
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
typedef struct{
    size_t peak;         /**< Index of the peak in the input data */
    double peak_height;  /**< Height (value) of the peak */

    lpr_peak_plateau_t plateau;      /**< Information about the peak's plateau */
    lr_peak_threshold_t threshold;    /**< Threshold values for the peak */
    lpr_peak_prominence_t prominence; /**< Prominence information for the peak */
    whlr_peak_width_t width;         /**< Width information for the peak */

} peak_result_t;


typedef enum {
    FP_STATUS_OK = 0,
    FP_STATUS_MEMORY_ALLOCATION_ERROR = 1,
    FP_STATUS_UNKNOWN_ERROR = 2
} fp_status;

/**
 * @brief Function for allocating memory objects for the find_peaks functions
 *
 * @param size Size of the memory allocated in bytes
 */
void* fp_malloc(size_t size);

/**
 * @brief Function for reallocating memory objects for the find_peaks functions
 *
 * @param struct_ptr Pointer to a memory allocated struct
 * @param size New size of the memory allocated in bytes
 */
void* fp_realloc(void* struct_ptr, size_t size);

/**
 * @brief Function for releasing the memory allocated objects allocated inside the find_peaks function
 *
 * @param struct_ptr Pointer to a memory allocated struct
 */
void fp_free(void* struct_ptr);

/**
 * @brief Find peaks in a one-dimensional array according to specified conditions
 *
 * This function identifies local maxima (peaks) in an array of double values and filters them
 * based on various criteria including height, prominence, width, threshold, and distance between peaks.
 * The function first identifies all local maxima, then applies the filtering conditions to select
 * peaks that meet all criteria.
 *
 * @param x Pointer to an array of double values to search for peaks
 * @param size The size of the input array
 * @param cond A structure containing all filtering conditions:
 *        - height: min/max allowed peak heights
 *        - prominence: min/max allowed peak prominences
 *        - width: min/max allowed peak widths
 *        - threshold: min/max thresholds for peak detection
 *        - plateau_size: min/max size of peak plateaus
 *        - distance: minimum distance between peaks
 *        - wlen: window length for prominence calculation
 *        - rel_height: relative height for width calculation
 * @param results Pointer to a pointer where the array of peak results will be stored
 *               (will be allocated by this function)
 * @param results_size Pointer where the number of detected peaks will be stored
 *
 * @return 0 on success, non-zero on error (typically memory allocation failure)
 *
 * @note The caller is responsible for freeing the memory allocated for *results
 * @note If no peaks match the criteria, *results will be NULL and *results_size will be 0
 */
fp_status find_peaks(const double* x, size_t size, fp_conditions_t cond, peak_result_t** results, size_t* results_size);


#endif //FIND_PEAKS_FIND_PEAKS_H
