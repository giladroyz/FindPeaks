#include <gtest/gtest.h>
extern "C" {
    #include "find_peaks.h"
}
#include <vector>
#include <cmath>

// Test fixture for find_peaks C tests
class FindPeaksCTest : public ::testing::Test {
protected:
    // Clean up allocated memory after each test
    void TearDown() override {
        if (peaks != nullptr) {
            fp_free(peaks);
            peaks = nullptr;
        }
    }

    // Sample signals for testing
    double simple_signal[11] = {0, 1, 0, 2, 0, 3, 0, 2, 0, 1, 0};
    double plateau_signal[12] = {0, 1, 1, 1, 0, 2, 2, 0, 3, 3, 3, 0};
    double noisy_signal[12] = {1.2, 0.8, 1.9, 1.5, 2.7, 1.8, 1.1, 2.5, 3.2, 2.8, 1.6, 0.9};
    double flat_signal[5] = {1, 1, 1, 1, 1};
    
    // Variables for storing peak detection results
    peak_result_t* peaks = nullptr;
    size_t num_peaks = 0;
};

// Test basic peak detection with default conditions
TEST_F(FindPeaksCTest, BasicPeakDetection) {
    fp_conditions_t cond = fp_get_default_conditions();
    
    ASSERT_EQ(FP_STATUS_OK, find_peaks(simple_signal, 11, cond, &peaks, &num_peaks));
    ASSERT_EQ(5, num_peaks);
    EXPECT_EQ(1, peaks[0].peak); // Peak at index 1 with height 1
    EXPECT_EQ(3, peaks[1].peak); // Peak at index 3 with height 2
    EXPECT_EQ(5, peaks[2].peak); // Peak at index 5 with height 3
    EXPECT_EQ(7, peaks[3].peak); // Peak at index 7 with height 2
    EXPECT_EQ(9, peaks[4].peak); // Peak at index 9 with height 1

}

// Test detection with height filtering
TEST_F(FindPeaksCTest, HeightFiltering) {
    fp_conditions_t cond = fp_get_default_conditions();
    fp_cond_set_height_mn(&cond, 2.5); // Only peaks >= 2.5
    
    ASSERT_EQ(FP_STATUS_OK, find_peaks(simple_signal, 11, cond, &peaks, &num_peaks));
    ASSERT_EQ(1, num_peaks);
    EXPECT_EQ(5, peaks[0].peak); // Only the peak with height 3
    EXPECT_DOUBLE_EQ(3.0, peaks[0].peak_height);
}

// Test detection with minimum distance between peaks
TEST_F(FindPeaksCTest, DistanceFiltering) {
    fp_conditions_t cond = fp_get_default_conditions();
    fp_cond_set_distance_mn(&cond, 3); // Peaks must be at least 3 samples apart
    
    ASSERT_EQ(FP_STATUS_OK, find_peaks(simple_signal, 11, cond, &peaks, &num_peaks));
    ASSERT_EQ(3, num_peaks);
    EXPECT_EQ(1, peaks[0].peak); // First peak at index 1
    EXPECT_EQ(5, peaks[1].peak); // Second peak at index 5, the peak in index 3 is skipped
    EXPECT_EQ(9, peaks[2].peak); // Second peak at index 9, the peak in index 7 is skipped
}

// Test plateau detection
TEST_F(FindPeaksCTest, PlateauDetection) {
    fp_conditions_t cond = fp_get_default_conditions();
    
    ASSERT_EQ(FP_STATUS_OK, find_peaks(plateau_signal, 12, cond, &peaks, &num_peaks));
    ASSERT_EQ(3, num_peaks);
    
    // Check first plateau peak
    EXPECT_EQ(2, peaks[0].peak); // The middle of the plateau
    EXPECT_EQ(3, peaks[0].plateau.plateau_size); // 3 samples wide
    EXPECT_EQ(1, peaks[0].plateau.left_edge);
    EXPECT_EQ(3, peaks[0].plateau.right_edge);
    
    // Check second plateau peak
    EXPECT_EQ(5, peaks[1].peak); // The middle of the plateau
    EXPECT_EQ(2, peaks[1].plateau.plateau_size); // 2 samples wide
    
    // Check third plateau peak
    EXPECT_EQ(9, peaks[2].peak); // The middle of the plateau
    EXPECT_EQ(3, peaks[2].plateau.plateau_size); // 3 samples wide
}

// Test prominence calculation and filtering
TEST_F(FindPeaksCTest, ProminenceFiltering) {
    fp_conditions_t cond = fp_get_default_conditions();
    fp_cond_set_prominence_mn(&cond, 1.5); // Only peaks with prominence >= 1.5
    
    ASSERT_EQ(FP_STATUS_OK, find_peaks(simple_signal, 11, cond, &peaks, &num_peaks));
    ASSERT_EQ(3, num_peaks);

    EXPECT_EQ(3, peaks[0].peak); // Peak with height 2 has enough prominence
    EXPECT_GT(peaks[0].prominence.prominence, 1.5);

    EXPECT_EQ(5, peaks[1].peak); // Peak with height 3 has enough prominence
    EXPECT_GT(peaks[1].prominence.prominence, 1.5);

    EXPECT_EQ(7, peaks[2].peak); // Peak with height 2 has enough prominence
    EXPECT_GT(peaks[2].prominence.prominence, 1.5);
}

// Test width calculation and filtering
TEST_F(FindPeaksCTest, WidthFiltering) {
    fp_conditions_t cond = fp_get_default_conditions();
    fp_cond_set_width_mn(&cond, 2.0); // Only peaks with width >= 2.0
    fp_cond_set_rel_height_mn(&cond, 0.5); // Measure width at 50% of peak height
    
    ASSERT_EQ(FP_STATUS_OK, find_peaks(simple_signal, 11, cond, &peaks, &num_peaks));
    
    // We expect only peaks that are wide enough at half height
    for (size_t i = 0; i < num_peaks; i++) {
        EXPECT_GE(peaks[i].width.width, 2.0);
    }
}

// Test threshold filtering
TEST_F(FindPeaksCTest, ThresholdFiltering) {
    fp_conditions_t cond = fp_get_default_conditions();
    fp_cond_set_threshold_mn(&cond, 1.5); // Peak must exceed neighbors by at least 1.5
    
    ASSERT_EQ(FP_STATUS_OK, find_peaks(simple_signal, 11, cond, &peaks, &num_peaks));
    
    for (size_t i = 0; i < num_peaks; i++) {
        EXPECT_GE(peaks[i].threshold.left_threshold, 1.5);
        EXPECT_GE(peaks[i].threshold.right_threshold, 1.5);
    }
}

// Test with noisy data
TEST_F(FindPeaksCTest, NoisySignal) {
    fp_conditions_t cond = fp_get_default_conditions();
    
    ASSERT_EQ(FP_STATUS_OK, find_peaks(noisy_signal, 12, cond, &peaks, &num_peaks));
    
    // Verify that peaks are correctly identified in noisy signal
    ASSERT_GT(num_peaks, 0);
    
    // Check if the highest peak is detected (index 8, value 3.2)
    bool found_max_peak = false;
    for (size_t i = 0; i < num_peaks; i++) {
        if (peaks[i].peak == 8) {
            found_max_peak = true;
            EXPECT_DOUBLE_EQ(3.2, peaks[i].peak_height);
            break;
        }
    }
    EXPECT_TRUE(found_max_peak);
}

// Test with flat signal (no peaks)
TEST_F(FindPeaksCTest, FlatSignal) {
    fp_conditions_t cond = fp_get_default_conditions();
    
    ASSERT_EQ(FP_STATUS_OK, find_peaks(flat_signal, 5, cond, &peaks, &num_peaks));
    EXPECT_EQ(0, num_peaks);
}

// Test combined filtering criteria
TEST_F(FindPeaksCTest, CombinedFilters) {
    fp_conditions_t cond = fp_get_default_conditions();
    fp_cond_set_height_mn(&cond, 1.5);         // Height >= 1.5
    fp_cond_set_prominence_mn(&cond, 1.0);     // Prominence >= 1.0
    fp_cond_set_distance_mn(&cond, 2);         // At least 2 samples between peaks
    fp_cond_set_width_mn_mx(&cond, 1.0, 4.0);  // Width between 1.0 and 4.0
    
    ASSERT_EQ(FP_STATUS_OK, find_peaks(simple_signal, 11, cond, &peaks, &num_peaks));
    
    // Check that all peaks satisfy all conditions
    for (size_t i = 0; i < num_peaks; i++) {
        EXPECT_GE(peaks[i].peak_height, 1.5);
        EXPECT_GE(peaks[i].prominence.prominence, 1.0);
        EXPECT_GE(peaks[i].width.width, 1.0);
        EXPECT_LE(peaks[i].width.width, 4.0);
    }
    
    // Make sure peaks are at least 2 samples apart
    for (size_t i = 1; i < num_peaks; i++) {
        EXPECT_GE(peaks[i].peak - peaks[i-1].peak, 2);
    }
}

// Test against SciPy-like sample for compatibility
TEST_F(FindPeaksCTest, SciPyCompatibility) {
    // This is a known test signal with peaks that match SciPy's find_peaks results
    double scipy_signal[11] = {0, 1, 0, 2, 0, 3, 0, 2, 0, 1, 0};
    
    // These are the expected peak indices from SciPy's find_peaks
    size_t expected_peaks[5] = {1, 3, 5, 7, 9};
    
    fp_conditions_t cond = fp_get_default_conditions();
    fp_cond_set_height_mn(&cond, 0.5);  // Same as SciPy test
    
    ASSERT_EQ(FP_STATUS_OK, find_peaks(scipy_signal, 11, cond, &peaks, &num_peaks));
    
    ASSERT_EQ(5, num_peaks);
    
    for (size_t i = 0; i < num_peaks; i++) {
        EXPECT_EQ(expected_peaks[i], peaks[i].peak);
    }
}
