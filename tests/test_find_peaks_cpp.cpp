#include <gtest/gtest.h>
#include "find_peaks.hpp"
#include <vector>
#include <cmath>

using namespace findPeaks;

// Test fixture for find_peaks C++ tests
class FindPeaksCppTest : public ::testing::Test {
protected:
    // Sample signals for testing
    std::vector<double> simple_signal = {0, 1, 0, 2, 0, 3, 0, 2, 0, 1, 0};
    std::vector<double> plateau_signal = {0, 1, 1, 1, 0, 2, 2, 0, 3, 3, 3, 0};
    std::vector<double> noisy_signal = {1.2, 0.8, 1.9, 1.5, 2.7, 1.8, 1.1, 2.5, 3.2, 2.8, 1.6, 0.9};
    std::vector<double> flat_signal = {1, 1, 1, 1, 1};
    std::vector<double> empty_signal = {};
};

// Test basic peak detection with default conditions
TEST_F(FindPeaksCppTest, BasicPeakDetection) {
    auto peaks = find_peaks(simple_signal);
    ASSERT_EQ(5, peaks.size());
    EXPECT_EQ(1, peaks[0].peak); // Peak at index 1 with height 1
    EXPECT_EQ(3, peaks[1].peak); // Peak at index 3 with height 2
    EXPECT_EQ(5, peaks[2].peak); // Peak at index 5 with height 3
    EXPECT_EQ(7, peaks[3].peak); // Peak at index 7 with height 2
    EXPECT_EQ(9, peaks[4].peak); // Peak at index 9 with height 1
}

// Test detection with height filtering
TEST_F(FindPeaksCppTest, HeightFiltering) {
    PeakConditions conditions;
    conditions.set_height(2.5); // Only peaks >= 2.5

    auto peaks = find_peaks(simple_signal, conditions);
    ASSERT_EQ(1, peaks.size());
    EXPECT_EQ(5, peaks[0].peak); // Only the peak with height 3
    EXPECT_DOUBLE_EQ(3.0, peaks[0].peak_height);
}

// Test detection with minimum distance between peaks
TEST_F(FindPeaksCppTest, DistanceFiltering) {
    PeakConditions conditions;
    conditions.set_distance(3); // Peaks must be at least 3 samples apart

    auto peaks = find_peaks(simple_signal, conditions);
    ASSERT_EQ(3, peaks.size());
    EXPECT_EQ(1, peaks[0].peak); // First peak at index 1
    EXPECT_EQ(5, peaks[1].peak); // Second peak at index 5, the peak in index 3 is skipped
    EXPECT_EQ(9, peaks[2].peak); // Second peak at index 9, the peak in index 7 is skipped
}

// Test plateau detection
TEST_F(FindPeaksCppTest, PlateauDetection) {
    auto peaks = find_peaks(plateau_signal);
    ASSERT_EQ(3, peaks.size());

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
TEST_F(FindPeaksCppTest, ProminenceFiltering) {
    PeakConditions conditions;
    conditions.set_prominence(1.5); // Only peaks with prominence >= 1.5

    auto peaks = find_peaks(simple_signal, conditions);
    ASSERT_EQ(3, peaks.size());

    EXPECT_EQ(3, peaks[0].peak); // Peak with height 2 has enough prominence
    EXPECT_GT(peaks[0].prominence.prominence, 1.5);

    EXPECT_EQ(5, peaks[1].peak); // Peak with height 3 has enough prominence
    EXPECT_GT(peaks[1].prominence.prominence, 1.5);

    EXPECT_EQ(7, peaks[2].peak); // Peak with height 2 has enough prominence
    EXPECT_GT(peaks[2].prominence.prominence, 1.5);
}

// Test width calculation and filtering
TEST_F(FindPeaksCppTest, WidthFiltering) {
    PeakConditions conditions;
    conditions.set_width(2.0); // Only peaks with width >= 2.0
    conditions.set_rel_height(0.5); // Measure width at 50% of peak height

    auto peaks = find_peaks(simple_signal, conditions);

    // We expect only peaks that are wide enough at half height
    for (const auto& peak : peaks) {
        EXPECT_GE(peak.width.width, 2.0);
    }
}

// Test threshold filtering
TEST_F(FindPeaksCppTest, ThresholdFiltering) {
    PeakConditions conditions;
    conditions.set_threshold(1.5); // Peak must exceed neighbors by at least 1.5

    auto peaks = find_peaks(simple_signal, conditions);

    for (const auto& peak : peaks) {
        EXPECT_GE(peak.threshold.left_threshold, 1.5);
        EXPECT_GE(peak.threshold.right_threshold, 1.5);
    }
}

// Test with noisy data
TEST_F(FindPeaksCppTest, NoisySignal) {
    auto peaks = find_peaks(noisy_signal);

    // Verify that peaks are correctly identified in noisy signal
    ASSERT_GT(peaks.size(), 0);

    // Check if the highest peak is detected (index 8, value 3.2)
    bool found_max_peak = false;
    for (const auto& peak : peaks) {
        if (peak.peak == 8) {
            found_max_peak = true;
            EXPECT_DOUBLE_EQ(3.2, peak.peak_height);
            break;
        }
    }
    EXPECT_TRUE(found_max_peak);
}

// Test with flat signal (no peaks)
TEST_F(FindPeaksCppTest, FlatSignal) {
    auto peaks = find_peaks(flat_signal);
    EXPECT_EQ(0, peaks.size());
}

// Test with empty signal
TEST_F(FindPeaksCppTest, EmptySignal) {
    auto peaks = find_peaks(empty_signal);
    EXPECT_EQ(0, peaks.size());
}

// Test combined filtering criteria
TEST_F(FindPeaksCppTest, CombinedFilters) {
    PeakConditions conditions;
    conditions.set_height(1.5);         // Height >= 1.5
    conditions.set_prominence(1.0);     // Prominence >= 1.0
    conditions.set_distance(2);         // At least 2 samples between peaks
    conditions.set_width(1.0, 4.0);     // Width between 1.0 and 4.0

    auto peaks = find_peaks(simple_signal, conditions);

    // Check that all peaks satisfy all conditions
    for (const auto& peak : peaks) {
        EXPECT_GE(peak.peak_height, 1.5);
        EXPECT_GE(peak.prominence.prominence, 1.0);
        EXPECT_GE(peak.width.width, 1.0);
        EXPECT_LE(peak.width.width, 4.0);
    }

    // Make sure peaks are at least 2 samples apart
    for (size_t i = 1; i < peaks.size(); i++) {
        EXPECT_GE(peaks[i].peak - peaks[i-1].peak, 2);
    }
}

// Test against SciPy-like sample for compatibility
TEST_F(FindPeaksCppTest, SciPyCompatibility) {
    // This is a known test signal with peaks that match SciPy's find_peaks results
    std::vector<double> scipy_signal = {0, 1, 0, 2, 0, 3, 0, 2, 0, 1, 0};

    // These are the expected peak indices from SciPy's find_peaks
    std::vector<size_t> expected_peaks = {1, 3, 5, 7, 9};

    PeakConditions conditions;
    conditions.set_height(0.5);         // Same as SciPy test

    auto peaks = find_peaks(scipy_signal, conditions);

    ASSERT_EQ(expected_peaks.size(), peaks.size());

    for (size_t i = 0; i < peaks.size(); i++) {
        EXPECT_EQ(expected_peaks[i], peaks[i].peak);
    }
}
