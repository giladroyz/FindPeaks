/*
 * This implementation is derived from scipy.signal.find_peaks
 * Original code copyright:
 * Copyright (c) 2001-2002 Enthought, Inc. 2003, SciPy Developers.
 * Original code available under the BSD 3-clause license.
 * See LICENSE file for the complete license text.
 */

#include <vector>
#include <numeric>
#include <algorithm>
#include <iostream>

#include "find_peaks.hpp"

namespace findPeaks {

    size_t max_int(size_t x, size_t y) {
        return x >= y ? x : y;
    }

    size_t min_int(size_t x, size_t y) {
        return x <= y ? x : y;
    }

    double max_double(double x, double y) {
        return x >= y ? x : y;
    }

    double min_double(double x, double y) {
        return x <= y ? x : y;
    }

    std::vector<size_t> argsort(const std::vector<double> &x) {
        std::vector<size_t> idx(x.size());
        std::iota(idx.begin(), idx.end(), 0);

        stable_sort(idx.begin(), idx.end(), [&x](size_t i1, size_t i2) { return x[i1] < x[i2]; });

        return idx;
    }

    std::vector<lmr_peak_index_t> _local_maxima_1d(const std::vector<double> &x) {
        std::vector<lmr_peak_index_t> peaks;

        if (x.empty())
            return peaks;

        size_t i_ahead;
        size_t i_max = x.size() - 1;

        for (size_t i = 1; i < i_max; i++) {
            // `i` is the Pointer to current sample, first one can't be maxima

            //Test if previous sample is smaller
            if (x[i - 1] < x[i]) {
                //Index to look ahead of current sample
                i_ahead = i + 1;

                //Find next sample that is unequal to x[i]
                while (i_ahead < i_max && x[i_ahead] == x[i])
                    i_ahead++;

                //Maxima is found if next unequal sample is smaller than x[i]
                if (x[i_ahead] < x[i]) {
                    lmr_peak_index_t peak;
                    peak.left_edge = i;
                    peak.right_edge = i_ahead - 1;
                    peak.mid_point = (peak.left_edge + peak.right_edge) / 2;

                    peaks.push_back(peak);

                    //Skip samples that can't be maximum
                    i = i_ahead;
                }
            }
        }

        return peaks;
    }

    std::vector<bool> _select_by_peak_distance(const std::vector<size_t> &peaks,
                                               const std::vector<double> &priority,
                                               size_t distance) {
        size_t i, k, j;
        std::vector<bool> keep(peaks.size(), true);

        //Create map from `i` (index for `peaks` sorted by `priority`) to `j` (index
        //for `peaks` sorted by position). This allows to iterate `peaks` and `keep`
        //with `j` by order of `priority` while still maintaining the ability to
        //step to neighbouring peaks with (`j` + 1) or (`j` - 1).
        std::vector<size_t> priority_to_position = argsort(priority);

        //    //Round up because actual peak distance can only be natural number
        //    size_t distance_ = distance;
        //    distance_ = distance_ == distance ? distance_ : distance_ + 1;

        //Highest priority first -> iterate in reverse order (decreasing)
        for (i = 0; i < peaks.size(); i++) {
            //"Translate" `i` to `j` which points to current peak whose
            //neighbours are to be evaluated
            j = priority_to_position[peaks.size() - 1 - i];

            //Skip evaluation for peak already marked as "don't keep"
            if (keep[j] == 0)
                continue;

            k = 1;
            //Flag "earlier" peaks for removal until minimal distance is exceeded
            while (k <= j && peaks[j] - peaks[j - k] < distance) {
                keep[j - k] = false;
                k++;
            }

            k = j + 1;
            //Flag "later" peaks for removal until minimal distance is exceeded
            while (k < peaks.size() && peaks[k] - peaks[j] < distance) {
                keep[k] = false;
                k++;
            }
        }

        return keep;
    }

    std::vector<lpr_peak_prominence_t> _peak_prominences(const std::vector<double> &x,
                                                         const std::vector<size_t> &peaks,
                                                         size_t wlen) {
        std::vector<lpr_peak_prominence_t> prominences;

        size_t i;
        double left_min, right_min;

        size_t peak, i_min, i_max;
        size_t half_wlen = wlen / 2;

//    std::cout << "wlen: " << wlen << ", " << half_wlen << std::endl;

        for (size_t peak_nr = 0; peak_nr < peaks.size(); peak_nr++) {
            lpr_peak_prominence_t prominence;

            peak = peaks[peak_nr];
            i_min = 0;
            i_max = x.size() - 1;

            if (wlen >= 2) {
                //Adjust window around the evaluated peak (within bounds);
                //if wlen is even the resulting window length is is implicitly
                //rounded to next odd integer
                i_min = max_int(peak - half_wlen, i_min);
                i_max = min_int(peak + half_wlen, i_max);
            }

            //Find the left base in interval [i_min, peak]
            i = peak;
            prominence.left_base = peak;
            left_min = x[peak];

            while (i_min <= i && x[i] <= x[peak]) {
                if (x[i] < left_min) {
                    left_min = x[i];
                    prominence.left_base = i;
                }

                if (i == 0 && i_min == 0)
                    break;

                i--;
            }

            //Find the right base in interval [peak, i_max]
            i = peak;
            prominence.right_base = peak;
            right_min = x[peak];

            while (i <= i_max && x[i] <= x[peak]) {
                if (x[i] < right_min) {
                    right_min = x[i];
                    prominence.right_base = i;
                }
                i++;
            }

            prominence.prominence = x[peak] - max_double(left_min, right_min);

            prominences.push_back(prominence);
        }

        return prominences;
    }

    std::vector<whlr_peak_width_t> _peak_widths(const std::vector<double> &x,
                                                const std::vector<size_t> &peaks,
                                                double rel_height,
                                                const std::vector<lpr_peak_prominence_t> &prominences) {
        std::vector<whlr_peak_width_t> widths;

        size_t peak, i, i_max, i_min;
        double height, left_ip, right_ip;

        for (size_t p = 0; p < peaks.size(); p++) {
            whlr_peak_width_t width_data;

            i_min = prominences[p].left_base;
            i_max = prominences[p].right_base;
            peak = peaks[p];

            height = x[peak] - prominences[p].prominence * rel_height;
            width_data.width_height = x[peak] - prominences[p].prominence * rel_height;

            //Find intersection point on left side
            i = peak;
            while (i_min < i && height < x[i])
                i -= 1;
            left_ip = (double) i;
            //Interpolate if true intersection height is between samples
            if (x[i] < height)
                left_ip += (height - x[i]) / (x[i + 1] - x[i]);

            //Find intersection point on right side
            i = peak;
            while (i < i_max && height < x[i])
                i += 1;
            right_ip = (double) i;
            //Interpolate if true intersection height is between samples
            if (x[i] < height)
                right_ip -= (height - x[i]) / (x[i - 1] - x[i]);

            width_data.width = right_ip - left_ip;
            width_data.left_ip = left_ip;
            width_data.right_ip = right_ip;

            widths.push_back(width_data);
        }

        return widths;
    }

    std::vector<lr_peak_threshold_t> _peak_thresholds(const std::vector<double> &x,
                                                      const std::vector<size_t> &peaks) {
        std::vector<lr_peak_threshold_t> thresholds;

        for (size_t peak: peaks) {
            lr_peak_threshold_t thr;
            thr.left_threshold = x[peak] - x[peak - 1];
            thr.right_threshold = x[peak] - x[peak + 1];

            thresholds.push_back(thr);
        }

        return thresholds;
    }

    std::vector<lpr_peak_plateau_t> _peak_plateaus(const std::vector<lmr_peak_index_t> &peaks) {
        std::vector<lpr_peak_plateau_t> plateaus;

        for (size_t p = 0; p < peaks.size(); p++) {
            lpr_peak_plateau_t plateau;
            plateau.right_edge = peaks[p].right_edge;
            plateau.left_edge = peaks[p].left_edge;
            plateau.plateau_size = plateau.right_edge - plateau.left_edge + 1;

            plateaus.push_back(plateau);
        }

        return plateaus;
    }

    std::vector<size_t> _peak_indices(const std::vector<lmr_peak_index_t> &peaks) {
        std::vector<size_t> peak_indices;

        for (size_t p = 0; p < peaks.size(); p++)
            peak_indices.push_back(peaks[p].mid_point);

        return peak_indices;
    }

    std::vector<double> _peak_heights(const std::vector<double> &x,
                                      const std::vector<lmr_peak_index_t> &peaks) {
        std::vector<double> heights;

        for (size_t p = 0; p < peaks.size(); p++)
            heights.push_back(x[peaks[p].mid_point]);

        return heights;
    }

    void operator&=(std::vector<bool> &a, const std::vector<bool> &b) {
        for (size_t i = 0; i < a.size(); i++)
            a[i] = a[i] & b[i];
    }

    template<typename T>
    std::vector<T> apply_mask(const std::vector<T> &vec, const std::vector<bool> &mask) {
        std::vector<T> result;

        for (size_t p = 0; p < vec.size(); p++)
            if (mask[p]) result.push_back(vec[p]);

        return result;
    }

    std::vector<peak_result_t> find_peaks(const std::vector<double> &x, PeakConditions conditions) {
        std::vector<peak_result_t> results;

        std::vector<lmr_peak_index_t> peaks = _local_maxima_1d(x);
        std::vector<size_t> peak_indices = _peak_indices(peaks);
        std::vector<double> heights = _peak_heights(x, peaks);

        std::vector<lpr_peak_plateau_t> plateaus = _peak_plateaus(peaks);
        std::vector<lr_peak_threshold_t> thresholds = _peak_thresholds(x, peak_indices);

        std::vector<bool> good_peak_mask(peaks.size(), false);
        for (size_t p = 0; p < peak_indices.size(); p++) {
            peak_result_t r;
            r.peak = peak_indices[p];
            r.peak_height = heights[p];

            r.plateau = plateaus[p];
            r.threshold = thresholds[p];

            if (r.peak_height > conditions.height.upper || r.peak_height < conditions.height.lower)
                continue;
            if (r.plateau.plateau_size > conditions.plateau_size.upper ||
                r.plateau.plateau_size < conditions.plateau_size.lower)
                continue;
            if (min_double(r.threshold.right_threshold, r.threshold.left_threshold) < conditions.threshold.lower ||
                max_double(r.threshold.right_threshold, r.threshold.left_threshold) > conditions.threshold.upper)
                continue;

            good_peak_mask[p] = true;
        }

        peaks = apply_mask(peaks, good_peak_mask);
        peak_indices = _peak_indices(peaks);
        heights = apply_mask(heights, good_peak_mask);
        plateaus = apply_mask(plateaus, good_peak_mask);
        thresholds = apply_mask(thresholds, good_peak_mask);

        std::vector<bool> distance_mask = _select_by_peak_distance(peak_indices, heights, conditions.distance);

        std::vector<lpr_peak_prominence_t> prominences = _peak_prominences(x, peak_indices, conditions.wlen);
        std::vector<whlr_peak_width_t> widths = _peak_widths(x, peak_indices, conditions.rel_height, prominences);

        for (size_t p = 0; p < peak_indices.size(); p++) {
            peak_result_t r;
            r.peak = peak_indices[p];
            r.peak_height = heights[p];

            r.plateau = plateaus[p];
            r.threshold = thresholds[p];
            r.prominence = prominences[p];
            r.width = widths[p];

            if (!distance_mask[p])
                continue;
            if (r.prominence.prominence > conditions.prominence.upper ||
                r.prominence.prominence < conditions.prominence.lower)
                continue;
            if (r.width.width > conditions.width.upper || r.width.width < conditions.width.lower)
                continue;

            results.push_back(r);
        }

        return results;
    }

}