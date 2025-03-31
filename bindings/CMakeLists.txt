# CMake configuration for find_peaks Python bindings

# Python Bindings Configuration
# Includes the pybind11 library for creating Python bindings

# Let pybind11 find the python interpreter for our pytest
set(PYBIND11_FINDPYTHON ON)

include(FetchContent)
FetchContent_Declare(
        pybind11
        GIT_REPOSITORY https://github.com/pybind/pybind11
        GIT_TAG        v2.13.6
        EXCLUDE_FROM_ALL
)

FetchContent_GetProperties(pybind11)
if(NOT pybind11_POPULATED)
    FetchContent_Populate(pybind11)
    add_subdirectory(${pybind11_SOURCE_DIR} ${pybind11_BINARY_DIR})
endif()

# Python Bindings for C++ Library
# Creates a Python module named find_peaks_cpp_py from the source file
pybind11_add_module(find_peaks_cpp_py ${CMAKE_CURRENT_SOURCE_DIR}/src/py_bind_cpp.cpp
        ${FIND_PEAKS_SRC_DIR}/cpp/find_peaks.cpp)
# Adds the src directory to the include path for this module
target_include_directories(find_peaks_cpp_py PUBLIC ${FIND_PEAKS_SRC_DIR}/cpp)

# Python Bindings for C Library
# Creates a Python module named find_peaks_c_py from the source file
pybind11_add_module(find_peaks_c_py ${CMAKE_CURRENT_SOURCE_DIR}/src/py_bind_c.cpp
        ${FIND_PEAKS_SRC_DIR}/c/find_peaks.c)
# Adds the src directory to the include path for this module
target_include_directories(find_peaks_c_py PUBLIC ${FIND_PEAKS_SRC_DIR}/c)

# Add Python tests
if(BUILD_TESTS)
    add_test(
            NAME test_find_peaks_scipy
            COMMAND ${CMAKE_COMMAND} -E env FP_BINDING_PATH=$<TARGET_FILE_DIR:find_peaks_c_py> ${Python_EXECUTABLE} -m pytest -s ${CMAKE_CURRENT_SOURCE_DIR}/test_find_peaks_scipy.py -v
    )
    set_tests_properties(test_find_peaks_scipy PROPERTIES DEPENDS "find_peaks_c_py;find_peaks_cpp_py")
endif()
