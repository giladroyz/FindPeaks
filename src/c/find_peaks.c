/*
 * This implementation is derived from scipy.signal.find_peaks
 * Original code copyright:
 * Copyright (c) 2001-2002 Enthought, Inc. 2003, SciPy Developers.
 * Original code available under the BSD 3-clause license.
 * See LICENSE file for the complete license text.
 */

#include <stdlib.h>
#include <float.h>

#include "find_peaks.h"


void fp_set_double_range_mn_mx(fp_double_range_t* range, double min, double max)
{
    range->min = min;
    range->max = max;
}

void fp_set_double_range_mn(fp_double_range_t* range, double min)
{
    fp_set_double_range_mn_mx(range, min, DBL_MAX);
}

void fp_init_double_range(fp_double_range_t* range)
{
    fp_set_double_range_mn_mx(range, -DBL_MAX, DBL_MAX);
}

void fp_set_int_range_mn_mx(fp_int_range_t* range, size_t min, size_t max)
{
    range->min = min;
    range->max = max;
}

void fp_set_int_range_mn(fp_int_range_t* range, size_t min)
{
    fp_set_int_range_mn_mx(range, min, SIZE_MAX);
}

void fp_init_int_range(fp_int_range_t* range)
{
    fp_set_int_range_mn_mx(range, 0, SIZE_MAX);
}

void fp_cond_set_height_mn_mx(fp_conditions_t *cond, double min, double max)
{
    fp_set_double_range_mn_mx((&cond->height), min, max);
}

void fp_cond_set_height_mn(fp_conditions_t *cond, double min)
{
    fp_set_double_range_mn((&cond->height), min);
}

void fp_cond_init_height(fp_conditions_t *cond)
{
    fp_init_double_range((&cond->height));
}

void fp_cond_set_threshold_mn_mx(fp_conditions_t *cond, double min, double max)
{
    fp_set_double_range_mn_mx(&(cond->threshold), min, max);
}

void fp_cond_set_threshold_mn(fp_conditions_t *cond, double min)
{
    fp_set_double_range_mn(&(cond->threshold), min);
}

void fp_cond_init_threshold(fp_conditions_t *cond)
{
    fp_init_double_range(&(cond->threshold));
}

void fp_cond_set_distance_mn(fp_conditions_t *cond, size_t distance)
{
    cond->distance = distance;
}

void fp_cond_init_distance(fp_conditions_t *cond)
{
    cond->distance = 1;
}

void fp_cond_set_prominence_mn_mx(fp_conditions_t *cond, double min, double max)
{
    fp_set_double_range_mn_mx(&(cond->prominence), min, max);
}

void fp_cond_set_prominence_mn(fp_conditions_t *cond, double min)
{
    fp_set_double_range_mn(&(cond->prominence), min);
}

void fp_cond_init_prominence(fp_conditions_t *cond)
{
    fp_init_double_range(&(cond->prominence));
//    fp_cond_set_prominence_mn_mx(cond, -DBL_MAX, DBL_MAX);
}

void fp_cond_set_width_mn_mx(fp_conditions_t *cond, double min, double max)
{
    fp_set_double_range_mn_mx(&(cond->width), min, max);
}

void fp_cond_set_width_mn(fp_conditions_t *cond, double min)
{
    fp_set_double_range_mn(&(cond->width), min);
}

void fp_cond_init_width(fp_conditions_t *cond)
{
    fp_init_double_range(&(cond->width));
}

void fp_cond_set_wlen_mn(fp_conditions_t *cond, size_t wlen)
{
    cond->wlen = wlen;
}

void fp_cond_init_wlen(fp_conditions_t *cond)
{
    cond->wlen = 0;
}

void fp_cond_set_rel_height_mn(fp_conditions_t *cond, double rel_height)
{
    cond->rel_height = rel_height;
}

void fp_cond_init_rel_height(fp_conditions_t *cond)
{
    cond->rel_height = 0.5;
}

void fp_cond_set_plateau_size_mn_mx(fp_conditions_t *cond, size_t min, size_t max)
{
    fp_set_int_range_mn_mx(&(cond->plateau_size), min, max);
}

void fp_cond_set_plateau_size_mn(fp_conditions_t *cond, size_t min)
{
    fp_set_int_range_mn(&(cond->plateau_size), min);
}

void fp_cond_init_plateau_size(fp_conditions_t *cond)
{
    fp_init_int_range(&(cond->plateau_size));
}

void fp_init_conditions(fp_conditions_t* cond)
{
    fp_cond_init_height(cond);
    fp_cond_init_threshold(cond);
    fp_cond_init_distance(cond);
    fp_cond_init_prominence(cond);
    fp_cond_init_width(cond);
    fp_cond_init_wlen(cond);
    fp_cond_init_rel_height(cond);
    fp_cond_init_plateau_size(cond);
}

fp_conditions_t fp_get_default_conditions()
{
    fp_conditions_t cond;
    fp_init_conditions(&cond);

    return cond;
}

typedef struct{
    double value;
    size_t index;
} value_and_index_t;

int argsort_compare(const void* a, const void* b)
{
    value_and_index_t i = *(value_and_index_t*)a;
    value_and_index_t j = *(value_and_index_t*)b;
    return (i.value > j.value) - ((i.value < j.value));
//    return ((value_and_index_t*)a)->value > ((value_and_index_t*)b)->value;
}

fp_status argsort(const double* x, size_t size, size_t* idx)
{
    value_and_index_t* comp_array = (value_and_index_t*) fp_malloc(size * sizeof(value_and_index_t));
    if(comp_array == NULL)
        return FP_STATUS_MEMORY_ALLOCATION_ERROR;

    for(size_t i = 0; i < size; i++)
    {
        comp_array[i].value = x[i];
        comp_array[i].index = i;
    }

    qsort(comp_array, size, sizeof(value_and_index_t),argsort_compare);

    for(size_t i = 0; i < size; i++)
        idx[i] = comp_array[i].index;

    fp_free(comp_array);

    return FP_STATUS_OK;
}

size_t max_i(size_t x, size_t y)
{
    return x >= y ? x : y;
}

size_t min_i(size_t x, size_t y)
{
    return x <= y ? x : y;
}

double max_f(double x, double y)
{
    return x >= y ? x : y;
}

double min_f(double x, double y)
{
    return x <= y ? x : y;
}

fp_status local_maxima_1d(const double* x, size_t size, lmr_peak_index_t** peaks, size_t* peaks_size)
{
    if( x == NULL || size == 0)
    {
        *peaks = NULL;
        *peaks_size = 0;
        return FP_STATUS_OK;
    }

    *peaks = (lmr_peak_index_t*) fp_malloc(size * sizeof(lmr_peak_index_t));

    size_t i_ahead;
    size_t i_max = size - 1;

    size_t peak_idx = 0;

    if(*peaks == NULL)
        return FP_STATUS_MEMORY_ALLOCATION_ERROR;

    for(size_t i = 1; i < i_max; i++)
    {
        // `i` is the Pointer to current sample, first one can't be maxima

        //Test if previous sample is smaller
        if(x[i - 1] < x[i])
        {
            //Index to look ahead of current sample
            i_ahead = i + 1;

            //Find next sample that is unequal to x[i]
            while(i_ahead < i_max && x[i_ahead] == x[i])
                i_ahead++;

            //Maxima is found if next unequal sample is smaller than x[i]
            if(x[i_ahead] < x[i])
            {
                lmr_peak_index_t peak;
                peak.left_edge = i;
                peak.right_edge = i_ahead - 1;
                peak.mid_point = (peak.left_edge + peak.right_edge) / 2;

                (*peaks)[peak_idx] = peak;

                peak_idx++;

                //Skip samples that can't be maximum
                i = i_ahead;
            }
        }
    }

    if(peak_idx == 0)
    {
        fp_free(*peaks);
        *peaks = NULL;
        *peaks_size = peak_idx;
        return FP_STATUS_OK;
    }

    *peaks = fp_realloc(*peaks, peak_idx * sizeof(lmr_peak_index_t));

    if(*peaks == NULL)
        return FP_STATUS_MEMORY_ALLOCATION_ERROR;

    *peaks_size = peak_idx;

    return FP_STATUS_OK;
}

fp_status select_by_peak_distance(const size_t* peaks, size_t size, const double* priority, size_t distance, int* keep)
{
    size_t i, k, j;

    //Create map from `i` (index for `peaks` sorted by `priority`) to `j` (index
    //for `peaks` sorted by position). This allows to iterate `peaks` and `keep`
    //with `j` by order of `priority` while still maintaining the ability to
    //step to neighbouring peaks with (`j` + 1) or (`j` - 1).
    size_t* priority_to_position = (size_t*) fp_malloc(size * sizeof(size_t));

    if(priority_to_position == NULL)
        return FP_STATUS_MEMORY_ALLOCATION_ERROR;

    argsort(priority, size, priority_to_position);

    //    //Round up because actual peak distance can only be natural number
    //    size_t distance_ = distance;
    //    distance_ = distance_ == distance ? distance_ : distance_ + 1;

    //Highest priority first -> iterate in reverse order (decreasing)
    for(i = 0; i < size; i++)
    {
        //"Translate" `i` to `j` which points to current peak whose
        //neighbours are to be evaluated
        j = priority_to_position[size - 1 - i];

        //Skip evaluation for peak already marked as "don't keep"
        if(keep[j] == 0)
            continue;

        k = 1;
        //Flag "earlier" peaks for removal until minimal distance is exceeded
        while(k <= j && peaks[j] - peaks[j - k] < distance)
        {
            keep[j - k] = 0;
            k++;
        }

        k = j + 1;
        //Flag "later" peaks for removal until minimal distance is exceeded
        while(k < size && peaks[k] - peaks[j] < distance)
        {
            keep[k] = 0;
            k++;
        }
    }

    fp_free(priority_to_position);

    return FP_STATUS_OK;
}

void peak_prominences(const double* x,
                      size_t size_x,
                      const size_t* peaks,
                      size_t size_peaks,
                      size_t wlen,
                      const int *mask,
                      lpr_peak_prominence_t* prominences)
{
    size_t i;
    double left_min, right_min;

    size_t peak, i_min, i_max;
    size_t half_wlen = wlen / 2;

//    printf("wlen: %d, %d", wlen, half_wlen);
//    std::cout << "wlen: " << wlen << ", " << half_wlen << std::endl;

    for(size_t peak_nr = 0; peak_nr < size_peaks; peak_nr++)
    {
        if(!mask[peak_nr])
            continue;

        lpr_peak_prominence_t prominence;

        peak = peaks[peak_nr];
        i_min = 0;
        i_max = size_x - 1;

        if (wlen >= 2)
        {
            //Adjust window around the evaluated peak (within bounds);
            //if wlen is even the resulting window length is is implicitly
            //rounded to next odd integer
            i_min = max_i(peak - half_wlen, i_min);
            i_max = min_i(peak + half_wlen, i_max);
        }

        //Find the left base in interval [i_min, peak]
        i = peak;
        prominence.left_base = peak;
        left_min = x[peak];

        while(i_min <= i && x[i] <= x[peak])
        {
            if(x[i] < left_min)
            {
                left_min = x[i];
                prominence.left_base = i;
            }

            if(i == 0 && i_min == 0)
                break;

            i--;
        }

        //Find the right base in interval [peak, i_max]
        i = peak;
        prominence.right_base = peak;
        right_min = x[peak];

        while(i <= i_max && x[i] <= x[peak])
        {
            if(x[i] < right_min)
            {
                right_min = x[i];
                prominence.right_base = i;
            }
            i++;
        }

        prominence.prominence = x[peak] - max_f(left_min, right_min);

        prominences[peak_nr] = prominence;
    }
}

void peak_widths(const double* x,
                 const size_t* peaks,
                 size_t size_peaks,
                 double rel_height,
                 const lpr_peak_prominence_t* prominences,
                 const int *mask,
                 whlr_peak_width_t* widths)
{

    size_t peak, i, i_max, i_min;
    double height, left_ip, right_ip;

    for(size_t p = 0; p < size_peaks; p++)
    {
        if(!mask[p])
            continue;

        whlr_peak_width_t width_data;

        i_min = prominences[p].left_base;
        i_max = prominences[p].right_base;
        peak = peaks[p];

        height = x[peak] - prominences[p].prominence * rel_height;
        width_data.width_height = x[peak] - prominences[p].prominence * rel_height;

        //Find intersection point on left side
        i = peak;
        while(i_min < i && height < x[i])
            i -= 1;
        left_ip = (double)i;
        //Interpolate if true intersection height is between samples
        if(x[i] < height)
            left_ip += (height - x[i]) / (x[i + 1] - x[i]);

        //Find intersection point on right side
        i = peak;
        while(i < i_max && height < x[i])
            i += 1;
        right_ip = (double )i;
        //Interpolate if true intersection height is between samples
        if(x[i] < height)
            right_ip -= (height - x[i]) / (x[i - 1] - x[i]);

        width_data.width = right_ip - left_ip;
        width_data.left_ip = left_ip;
        width_data.right_ip = right_ip;

        widths[p] = width_data;
    }
}

void peak_thresholds(const double* x,
                     const size_t* peaks,
                     size_t size_peaks,
                     lr_peak_threshold_t* thresholds)
{

    for(size_t peak_idx = 0; peak_idx < size_peaks; peak_idx++)
    {
        lr_peak_threshold_t thr;
        thr.left_threshold = x[peaks[peak_idx]] - x[peaks[peak_idx]-1];
        thr.right_threshold = x[peaks[peak_idx]] - x[peaks[peak_idx]+1];

        thresholds[peak_idx] = thr;
    }
}

void peak_plateaus(const lmr_peak_index_t* peaks, size_t size_peaks, lpr_peak_plateau_t* plateaus)
{

    for(size_t p = 0; p < size_peaks; p++)
    {
        lpr_peak_plateau_t plateau;
        plateau.right_edge = peaks[p].right_edge;
        plateau.left_edge = peaks[p].left_edge;
        plateau.plateau_size = plateau.right_edge - plateau.left_edge + 1;

        plateaus[p] = plateau;
    }
}

void peak_heights(const double* x, const lmr_peak_index_t* peaks, size_t size_peaks, double* heights)
{
    for(size_t p = 0; p < size_peaks; p++)
        heights[p] = x[peaks[p].mid_point];
}

void peak_indices(const lmr_peak_index_t* peaks, size_t size_peaks, size_t* peak_indices)
{
    for(size_t p = 0; p < size_peaks; p++)
        peak_indices[p] = peaks[p].mid_point;
}

typedef struct
{
    size_t max_peaks_count;

    size_t peaks_size;
    peak_result_t* peaks_data;

//    size_t* peak_idx;
    int* mask;
    size_t* priority_to_position;

} fp_spec;

void* fp_malloc(size_t size)
{
    return malloc(size);
}

void* fp_realloc(void* struct_ptr, size_t size)
{
    if(struct_ptr == NULL)
        return NULL;

    return realloc(struct_ptr, size);
}

void fp_free(void* struct_ptr)
{
    if(struct_ptr != NULL)
        free(struct_ptr);
}

fp_status find_peaks(const double* x, size_t size, fp_conditions_t cond, peak_result_t** results, size_t* results_size)
{
    size_t i, counter;
    peak_result_t r;

    /* Variables for storing peak information */
    lmr_peak_index_t* peaks = NULL;       /* Raw peak data */
    size_t peaks_size = 0;                /* Number of peaks found */

    size_t* peak_idx = NULL;              /* Array of peak indices */
    double* heights = NULL;               /* Array of peak heights */

    /* Variables for storing peak characteristics */
    lpr_peak_plateau_t* plateaus = NULL;      /* Information about peak plateaus */
    lr_peak_threshold_t* thresholds = NULL;   /* Threshold information for each peak */
    lpr_peak_prominence_t* prominences = NULL; /* Prominence information for each peak */
    whlr_peak_width_t* widths = NULL;         /* Width information for each peak */

    int* mask = NULL;                     /* Mask array to track which peaks pass filters */

    fp_status ret;                              /* Return value */

    /* Step 1: Find all local maxima in the input array */
    ret = local_maxima_1d(x, size, &peaks, &peaks_size);
    if(ret)
        goto cleanup;

    /* Step 2: Allocate memory for peak properties */
    peak_idx = (size_t*) fp_malloc(peaks_size * sizeof(size_t));
    heights = (double*) fp_malloc(peaks_size * sizeof(double));
    mask = (int*) fp_malloc(peaks_size * sizeof(int));

    /* Check memory allocation success */
    if(peak_idx == NULL || heights == NULL || mask == NULL)
    {
        ret = FP_STATUS_MEMORY_ALLOCATION_ERROR;
        goto cleanup;
    }

    /* Allocate memory for more detailed peak characteristics */
    plateaus = (lpr_peak_plateau_t*) fp_malloc(peaks_size * sizeof(lpr_peak_plateau_t));
    thresholds = (lr_peak_threshold_t*) fp_malloc(peaks_size * sizeof(lr_peak_threshold_t));
    prominences = (lpr_peak_prominence_t*) fp_malloc(peaks_size * sizeof(lpr_peak_prominence_t));
    widths = (whlr_peak_width_t*) fp_malloc(peaks_size * sizeof(whlr_peak_width_t));

    /* Check memory allocation success */
    if(plateaus == NULL || thresholds == NULL || prominences == NULL || widths == NULL)
    {
        ret = FP_STATUS_MEMORY_ALLOCATION_ERROR;
        goto cleanup;
    }

    /* Initialize all peaks as valid (passing filters) */
    for(i = 0; i < peaks_size; i++)
        mask[i] = 1;

    /* Step 3: Extract peak indices from the raw peak data */
    peak_indices(peaks, peaks_size, peak_idx);

    /* Step 4: Calculate height of each peak */
    peak_heights(x, peaks, peaks_size, heights);

    /* Step 5: Identify plateau characteristics for each peak */
    peak_plateaus(peaks, peaks_size, plateaus);

    /* Step 6: Calculate threshold information for each peak */
    peak_thresholds(x, peak_idx, peaks_size, thresholds);

    /* Step 7: Filter peaks based on height, plateau size, and threshold criteria */
    for(i = 0; i < peaks_size; i++)
    {
        if(!mask[i])
            continue;  /* Skip peaks that have already been filtered out */

        /* Apply height filter */
        if(heights[i] > cond.height.max || heights[i] < cond.height.min)
            mask[i] = 0;

        /* Apply plateau size filter */
        if(plateaus[i].plateau_size > cond.plateau_size.max || plateaus[i].plateau_size < cond.plateau_size.min)
            mask[i] = 0;

        /* Apply threshold filter */
        if(min_f(thresholds[i].right_threshold, thresholds[i].left_threshold) < cond.threshold.min ||
           max_f(thresholds[i].right_threshold, thresholds[i].left_threshold) > cond.threshold.max)
            mask[i] = 0;
    }

    /* Step 8: Filter peaks based on minimum distance between peaks */
    ret = select_by_peak_distance(peak_idx, peaks_size, heights, cond.distance, mask);
    if(ret)
        goto cleanup;

    /* Step 9: Calculate peak prominences */
    peak_prominences(x, size, peak_idx, peaks_size, cond.wlen, mask, prominences);

    /* Step 10: Calculate peak widths */
    peak_widths(x, peak_idx, peaks_size, cond.rel_height, prominences, mask, widths);

    /* Step 11: Filter peaks based on prominence and width criteria */
    for(i = 0; i < peaks_size; i++)
    {
        if(!mask[i])
            continue;  /* Skip peaks that have already been filtered out */

        /* Apply prominence filter */
        if(prominences[i].prominence > cond.prominence.max || prominences[i].prominence < cond.prominence.min)
            mask[i] = 0;

        /* Apply width filter */
        if(widths[i].width > cond.width.max || widths[i].width < cond.width.min)
            mask[i] = 0;
    }

    /* Step 12: Count how many peaks passed all filters */
    counter = 0;
    for(i = 0; i < peaks_size; i++)
        if(mask[i])
            counter++;

    *results_size = counter;

    /* If no peaks pass the filters, return early */
    if(*results_size == 0)
    {
        ret = FP_STATUS_OK;
        goto cleanup;
    }

    /* Step 13: Allocate memory for the final results array */
    *results = (peak_result_t*) fp_malloc(*results_size * sizeof(peak_result_t));

    if(*results == NULL)
    {
        ret = FP_STATUS_MEMORY_ALLOCATION_ERROR;
        goto cleanup;
    }

    /* Step 14: Fill the results array with peaks that passed all filters */
    for(i = 0, counter = 0; i < peaks_size; i++)
    {
        if(!mask[i])
            continue;  /* Skip filtered peaks */

        /* Store all peak properties in the result structure */
        r.peak = peak_idx[i];
        r.peak_height = heights[i];

        r.plateau = plateaus[i];
        r.threshold = thresholds[i];
        r.prominence = prominences[i];
        r.width = widths[i];

        (*results)[counter] = r;
        counter++;
    }

    cleanup:
    /* Free all allocated memory */
    fp_free(peaks);
    fp_free(peak_idx);
    fp_free(heights);
    fp_free(plateaus);
    fp_free(thresholds);
    fp_free(prominences);
    fp_free(widths);
    fp_free(mask);

    return ret;  /* Return FP_STATUS_OK (0) on success, non-zero on error */
}