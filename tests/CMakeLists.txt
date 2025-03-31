# CMake configuration for find_peaks tests
# Uses Google Test framework for C and C++ implementations

# fix gtest compiling errors
set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)

# Get GTest package
include(FetchContent)
FetchContent_Declare(
        googletest
        GIT_REPOSITORY https://github.com/google/googletest
        GIT_TAG        v1.16.0
        EXCLUDE_FROM_ALL
)

FetchContent_GetProperties(googletest)
if(NOT googletest_POPULATED)
    FetchContent_Populate(googletest)
    add_subdirectory(${googletest_SOURCE_DIR} ${googletest_BINARY_DIR})
endif()

include_directories(${GTEST_INCLUDE_DIRS})

#set(CMAKE_CXX_FLAGS "/MT")

# C++ tests
add_executable(test_find_peaks_cpp
    gtest_runner.cpp
    test_find_peaks_cpp.cpp
)
target_include_directories(test_find_peaks_cpp PUBLIC ${FIND_PEAKS_SRC_DIR}/cpp)
target_link_libraries(test_find_peaks_cpp PUBLIC
    GTest::gtest
#    GTest::gtest_main
    find_peaks_cpp
)

# C tests
add_executable(test_find_peaks_c
    gtest_runner.cpp
    test_find_peaks_c.cpp
)
target_include_directories(test_find_peaks_c PUBLIC ${FIND_PEAKS_SRC_DIR}/c)
target_link_libraries(test_find_peaks_c PUBLIC
    GTest::gtest
#    GTest::gtest_main
    find_peaks_c
)

# Register tests with CTest
add_test(NAME test_cpp COMMAND test_find_peaks_cpp)
add_test(NAME test_c COMMAND test_find_peaks_c)
