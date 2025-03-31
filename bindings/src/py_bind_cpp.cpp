#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <pybind11/iostream.h>
#include <iostream>
#include <cmath>

#include "find_peaks.hpp"

#define STRINGIFY(x) #x
#define MACRO_STRINGIFY(x) STRINGIFY(x)

namespace py = pybind11;
using namespace findPeaks;

/**
 * Find peaks in a 1D signal with C++ implementation
 * 
 * @param x Input signal array
 * @param height Height range constraint for peaks
 * @param threshold Threshold range for peak detection
 * @param distance Minimum distance between neighboring peaks
 * @param prominence Prominence range constraint for peaks
 * @param width Width range constraint for peaks
 * @param wlen Window length for calculating prominence
 * @param rel_height Relative height for width calculation
 * @param plateau_size Plateau size range constraint
 * @return Tuple containing peak indices and properties
 */
py::tuple find_peaks_cpp(const py::array_t<double>& x,
                                  const py::tuple& height,
                                  const py::tuple& threshold,
                                  size_t distance,
                                  const py::tuple& prominence,
                                  const py::tuple& width,
                                  size_t wlen,
                                  double rel_height,
                                  const py::tuple& plateau_size)
{
    // Initialize peak detection conditions
    PeakConditions cond;

    // Parse optional parameters and set conditions if provided
    if(!height.is_none())
        cond.set_height(RangeFloat(py::cast<double>(height[0]), py::cast<double>(height[1])));
    if(!threshold.is_none())
        cond.set_threshold(RangeFloat(py::cast<double>(threshold[0]), py::cast<double>(threshold[1])));
    if(distance != 0)
        cond.set_distance(distance);
    if(!prominence.is_none())
        cond.set_prominence(RangeFloat(py::cast<double>(prominence[0]), py::cast<double>(prominence[1])));
    if(!width.is_none())
        cond.set_width(RangeFloat(py::cast<double>(width[0]), py::cast<double>(width[1])));
    if(wlen >= 1)
        cond.set_wlen(wlen);
    if(rel_height != 0.5)
        cond.set_rel_height(rel_height);
    if(!plateau_size.is_none())
        cond.set_plateau_size(RangeInt(py::cast<size_t>(plateau_size[0]), py::cast<size_t>(plateau_size[1])));

    // Get input array size and convert to vector
    size_t size = x.size();
    auto x_vec = std::vector<double>(x.data(), x.data() + size);

    // Call the C++ implementation to find peaks
    std::vector<peak_result_t> result = find_peaks(x_vec, cond);

    // Create NumPy array for peak indices
    py::array_t<size_t> result_peaks(result.size());
    for(size_t i = 0; i < result.size(); i++)
        result_peaks.mutable_at(i) = result[i].peak;

    // Create list of dictionaries containing peak properties
    py::list result_params;
    for(auto r: result)
    {
        py::dict d;

        // Add peak height
        d["peak_heights"] = r.peak_height;

        // Add threshold properties
        d["left_thresholds"] = r.threshold.left_threshold;
        d["right_thresholds"] = r.threshold.right_threshold;

        // Add prominence properties
        d["prominences"] = r.prominence.prominence;
        d["left_bases"] = r.prominence.left_base;
        d["right_bases"] = r.prominence.right_base;

        // Add width properties
        d["widths"] = r.width.width;
        d["width_heights"] = r.width.width_height;
        d["left_ips"] = r.width.left_ip;
        d["right_ips"] = r.width.right_ip;

        // Add plateau properties
        d["plateau_sizes"] = r.plateau.plateau_size;
        d["left_edges"] = r.plateau.left_edge;
        d["right_edges"] = r.plateau.right_edge;

        result_params.append(d);
    }

    // Return tuple of peak indices and properties
    py::tuple total_res(2);
    total_res[0] = result_peaks;
    total_res[1] = result_params;

    return total_res;
}


PYBIND11_MODULE(find_peaks_cpp_py, m) {
m.doc() = R"pbdoc(
        Pybind11 example plugin
        -----------------------

        .. currentmodule:: find_peaks_cpp_py

        .. autosummary::
           :toctree: _generate

           _find_peaks

    )pbdoc";

// Define Python module interface with all parameters and documentation
m.def("_find_peaks", &find_peaks_cpp, py::call_guard<py::scoped_ostream_redirect>(),
        py::arg("x"),
        py::arg("height") = py::none(),
        py::arg("threshold") = py::none(),
        py::arg("distance") = 0,
        py::arg("prominence") = py::none(),
        py::arg("width") = py::none(),
        py::arg("wlen") = -1,
        py::arg("rel_height") = 0.5,
        py::arg("plateau_size") = py::none(),
        R"pbdoc(
            Find peaks in a 1D signal x.
            This wrapper function wraps the C++ function in a way that make the interface identical
            to the scipy.signal.find_peaks function in order to compare them and test the functionality of the C function.
            See documentation of scipy.signal.find_peaks for parameters and interface.
        )pbdoc");

// Set version information
#ifdef VERSION_INFO
    m.attr("__version__") = MACRO_STRINGIFY(VERSION_INFO);
#else
    m.attr("__version__") = "dev";
#endif
}