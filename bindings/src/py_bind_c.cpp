#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <pybind11/iostream.h>
#include <iostream>
#include <cmath>

extern "C" {
    #include "find_peaks.h"
}

#define STRINGIFY(x) #x
#define MACRO_STRINGIFY(x) STRINGIFY(x)

namespace py = pybind11;

/**
 * Find peaks in a 1D signal using the C implementation
 * 
 * Wraps the C find_peaks function to create the same interface as scipy.signal.find_peaks.
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
py::tuple find_peaks_c(const py::array_t<double>& x,
                         const py::tuple& height,
                         const py::tuple& threshold,
                         size_t distance,
                         const py::tuple& prominence,
                         const py::tuple& width,
                         size_t wlen,
                         double rel_height,
                         const py::tuple& plateau_size)
{
    // Initialize peak detection conditions with default values
    fp_conditions_t cond = fp_get_default_conditions();

    // Parse optional parameters and set conditions if provided
    if(!height.is_none())
        fp_cond_set_height_mn_mx(&cond, py::cast<double>(height[0]), py::cast<double>(height[1]));
    if(!threshold.is_none())
        fp_cond_set_threshold_mn_mx(&cond, py::cast<double>(threshold[0]), py::cast<double>(threshold[1]));
    if(!prominence.is_none())
        fp_cond_set_prominence_mn_mx(&cond, py::cast<double>(prominence[0]), py::cast<double>(prominence[1]));
    if(!width.is_none())
        fp_cond_set_width_mn_mx(&cond, py::cast<double>(width[0]), py::cast<double>(width[1]));
    if(!plateau_size.is_none())
        fp_cond_set_plateau_size_mn_mx(&cond, py::cast<size_t>(plateau_size[0]), py::cast<size_t>(plateau_size[1]));

    // Set scalar parameters
    fp_cond_set_distance_mn(&cond, distance);
    fp_cond_set_wlen_mn(&cond, wlen);
    fp_cond_set_rel_height_mn(&cond, rel_height);

    // Get input array size
    size_t size = x.size();
    auto x_vec = std::vector<double>(x.data(), x.data() + size);

    // Call C implementation to find peaks
    peak_result_t *result = NULL;
    size_t result_size = 0;
    find_peaks(x.data(), size, cond, &result, &result_size);

    // Create NumPy array for peak indices
    py::array_t<size_t> result_peaks(result_size);
    for(size_t i = 0; i < result_size; i++)
        result_peaks.mutable_at(i) = result[i].peak;

    // Create list of dictionaries containing peak properties
    py::list result_params;
    for(size_t i = 0; i < result_size; i++)
    {
        peak_result_t r = result[i];
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

    // Free memory allocated by C function
    fp_free(result);

    return total_res;
}


PYBIND11_MODULE(find_peaks_c_py, m) {
    m.doc() = R"pbdoc(
        Pybind11 wrapper for find_peaks C function.
        -----------------------

        .. currentmodule:: find_peaks_c_py

        .. autosummary::
           :toctree: _generate

           _find_peaks

    )pbdoc";

    // Define Python module interface with all parameters and documentation
    m.def("_find_peaks", &find_peaks_c, py::call_guard<py::scoped_ostream_redirect>(),
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
            This wrapper function wraps the C function in a way that make the interface identical
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
