#include <iostream>
#include <vector>

#include "find_peaks.hpp"

using namespace findPeaks;

int basic_example() {
    // Sample signal data
    std::vector<double> signal = {0, 1, 3, 1, 0, 2, 4, 2, 0, 1, 2, 1, 0};

    // Find peaks with default conditions
    std::vector<peak_result_t> peaks = find_peaks(signal);

    std::cout << "Original array: ";
    for (const auto& value : signal) {
        std::cout << value << " ";
    }
    std::cout << std::endl;

    // Print peak positions and heights
    std::cout << "Detected " << peaks.size() << " peaks\n";
    for (const auto& peak : peaks) {
        std::cout << "Peak at position " << peak.peak
                  << " with height " << peak.peak_height << "\n";
    }

    return 0;
}

int medium_example()
{
    // Sample signal data
    std::vector<double> data = {5, 3, 8, 4, 9, 1, 5, 4, 6};

    PeakConditions cond;
    cond.set_threshold(2, 10);      // Minimum relative height of 2.0 and maximum of 10.0

    std::vector<peak_result_t> peaks = find_peaks(data, cond);

    std::cout << "Original array: ";
    for (const auto& value : data) {
        std::cout << value << " ";
    }
    std::cout << std::endl;

    // Access detailed peak information
    for (const auto& peak : peaks) {
        std::cout << "Peak at position " << peak.peak << ":\n";
        std::cout << "  Height: " << peak.peak_height << "\n";
        std::cout << "  Left threshold: " << peak.threshold.left_threshold << "\n";
        std::cout << "  Right threshold: " << peak.threshold.right_threshold << "\n";
    }

    return 0;
}

int advance_example() {
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

int main()
{
    std::cout << "===== Basic Example =====\n\n";
    basic_example();
    std::cout << "=========================\n\n";

    std::cout << "===== Medium Example =====\n\n";
    medium_example();
    std::cout << "==========================\n\n";

    std::cout << "===== Advance Example =====\n\n";
    advance_example();
    std::cout << "===========================\n\n";
}