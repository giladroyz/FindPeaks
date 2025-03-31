#include <stdio.h>
#include "find_peaks.h"

int basic_example() {
    // Sample signal data
    double signal[] = {0, 1, 3, 1, 0, 2, 4, 2, 0, 1, 2, 1, 0};
    size_t signal_size = sizeof(signal) / sizeof(signal[0]);

    // Get default conditions
    fp_conditions_t cond = fp_get_default_conditions();

    // Variables to store results
    peak_result_t* peaks = NULL;
    size_t num_peaks = 0;

    printf("Original array: ");
    for (int i = 0; i < signal_size; i++) {
        printf("%.1f ", signal[i]);
    }
    printf("\n");

    // Find peaks
    if (find_peaks(signal, signal_size, cond, &peaks, &num_peaks) == FP_STATUS_OK) {
        printf("Detected %zu peaks\n", num_peaks);

        // Print peak positions and heights
        for (size_t i = 0; i < num_peaks; i++) {
            printf("Peak at position %zu with height %0.1f\n",
                   peaks[i].peak, peaks[i].peak_height);
        }

        // Free allocated memory
        fp_free(peaks);
    }

    return 0;
}

int medium_example()
{
    // Example usage
    double signal[] = {5, 3, 8, 4, 9, 1, 6, 5, 6};
    int signal_size = sizeof(signal) / sizeof(signal[0]);

    fp_conditions_t cond = fp_get_default_conditions();
    fp_cond_set_threshold_mn_mx(&cond, 2, 10); // Minimum relative height of 2.0 and maximum of 10.0

    // Variables to store results
    peak_result_t* peaks = NULL;
    size_t num_peaks = 0;

    int ret = find_peaks(signal, signal_size, cond, &peaks, &num_peaks);

    if(ret != FP_STATUS_OK)
    {
        printf("find_peaks had an error!");
        return 1;
    }

    printf("Original array: ");
    for (int i = 0; i < signal_size; i++) {
        printf("%.1f ", signal[i]);
    }
    printf("\n");

    printf("Detected %zu peaks\n", num_peaks);

    // Access detailed peak information
    for (size_t i = 0; i < num_peaks; i++) {
        printf("Peak at position %zu:\n", peaks[i].peak);
        printf("  Height: %0.1f\n", peaks[i].peak_height);
        printf("  Left threshold: %0.1f\n", peaks[i].threshold.left_threshold);
        printf("  Right threshold: %0.1f\n", peaks[i].threshold.right_threshold);
    }

    fp_free(peaks);

    return 0;
}

int advance_example() {
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

    printf("Original array: ");
    for (int i = 0; i < signal_size; i++) {
        printf("%.1f ", signal[i]);
    }
    printf("\n");

    // Find peaks
    if (find_peaks(signal, signal_size, cond, &peaks, &num_peaks) == FP_STATUS_OK) {
        printf("Detected %zu peaks\n", num_peaks);

        // Access detailed peak information
        for (size_t i = 0; i < num_peaks; i++) {
            printf("Peak at position %zu:\n", peaks[i].peak);
            printf("  Height: %0.1f\n", peaks[i].peak_height);
            printf("  Prominence: %0.1f\n", peaks[i].prominence.prominence);
            printf("  Width: %f at height %0.1f\n",
                   peaks[i].width.width, peaks[i].width.width_height);
            printf("  Left threshold: %0.1f\n", peaks[i].threshold.left_threshold);
            printf("  Right threshold: %0.1f\n", peaks[i].threshold.right_threshold);

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

int main()
{
    printf("===== Basic Example =====\n\n");
    basic_example();
    printf("=========================\n\n");

    printf("===== Medium Example =====\n\n");
    medium_example();
    printf("==========================\n\n");

    printf("===== Advance Example =====\n\n");
    advance_example();
    printf("===========================\n\n");
}