cmake_minimum_required(VERSION 3.28)
project(HAlign4Tests LANGUAGES C CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
option(HALIGN4_ENABLE_OFFSET_COMPACT "Enable WFA offset compaction in tests" OFF)


include(CTest)  # Provide ctest infrastructure and add_test support

# If user does not specify CMAKE_BUILD_TYPE, default to Release so optimizations are enabled by default
if(NOT DEFINED CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE "Release" CACHE STRING "Choose the type of build (default: Release)" FORCE)
    message(STATUS "CMAKE_BUILD_TYPE not set -> defaulting to: ${CMAKE_BUILD_TYPE}")
else()
    message(STATUS "CMAKE_BUILD_TYPE is: ${CMAKE_BUILD_TYPE}")
endif()
get_filename_component(HALIGN4_ROOT "${CMAKE_CURRENT_LIST_DIR}/.." ABSOLUTE)


find_package(OpenMP REQUIRED)
find_package(CURL)
find_package(ZLIB)

# --- Reuse your implementation: build a support static library ---
add_library(halign4_test_support STATIC
        "${HALIGN4_ROOT}/src/utils/preprocess.cpp"
        "${HALIGN4_ROOT}/src/utils/seq_io.cpp"
        "${HALIGN4_ROOT}/src/utils/cmd.cpp"
        "${HALIGN4_ROOT}/src/utils/file_io.cpp"
        "${HALIGN4_ROOT}/src/mash/xxhash.c"
        "${HALIGN4_ROOT}/src/mash/hash.cpp"
        "${HALIGN4_ROOT}/src/consensus/selector.cpp"
        "${HALIGN4_ROOT}/src/consensus/consensus.cpp"
        "${HALIGN4_ROOT}/src/seed/minimizer.cpp"
        "${HALIGN4_ROOT}/src/seed/anchor.cpp"
        "${HALIGN4_ROOT}/src/mash/mash.cpp"
        "${HALIGN4_ROOT}/src/align/ref_aligner.cpp"
        "${HALIGN4_ROOT}/src/align/ksw2_extz2_sse.c"
        "${HALIGN4_ROOT}/src/align/cigar.cpp"
        "${HALIGN4_ROOT}/src/align/align.cpp"
)

target_include_directories(halign4_test_support PUBLIC
        "${HALIGN4_ROOT}/include"
        "${HALIGN4_ROOT}/extern/CLI11/include"
        "${HALIGN4_ROOT}/extern/spdlog/include"
        "${HALIGN4_ROOT}/extern/cereal/include"
        "${HALIGN4_ROOT}/extern/doctest/include"
        "${HALIGN4_ROOT}/extern/WFA2-lib"
)

# Note: Explicitly link WFA2 static library here to avoid exposing only include directories without the corresponding library,
#       otherwise tests using WFA2 interfaces may get undefined symbol errors at link time.
target_link_libraries(halign4_test_support PUBLIC
        OpenMP::OpenMP_CXX
        wfa2::wfa2cpp_static
)

# Enable the same aggressive optimizations for the test support library as the main project in Release builds
if (MSVC)
    target_compile_options(halign4_test_support PRIVATE
        $<$<CONFIG:Release>:/O2>
        $<$<CONFIG:Release>:/arch:AVX2>
        $<$<CONFIG:Release>:/GL>
    )
    target_link_options(halign4_test_support PRIVATE
        $<$<CONFIG:Release>:/LTCG>
    )
else()
    target_compile_options(halign4_test_support PRIVATE
        $<$<CONFIG:Release>:-O3>
        $<$<CONFIG:Release>:-march=native>
    )
endif()

if (CURL_FOUND)
    target_link_libraries(halign4_test_support PUBLIC CURL::libcurl)
    target_compile_definitions(halign4_test_support PUBLIC HALIGN4_HAVE_LIBCURL=1)
endif()

if (ZLIB_FOUND)
    target_link_libraries(halign4_test_support PUBLIC ZLIB::ZLIB)
    target_compile_definitions(halign4_test_support PUBLIC HALIGN4_HAVE_ZLIB=1)
endif()

target_compile_definitions(halign4_test_support PUBLIC
        HALIGN4_ENABLE_OFFSET_COMPACT=$<BOOL:${HALIGN4_ENABLE_OFFSET_COMPACT}>
        $<$<CONFIG:Debug>:_DEBUG>
)

# ---------------------------------------------------------------------------
# WFA2-lib (Wavefront Alignment) subproject trimming: compile only the core library, not tests/benchmark
# ---------------------------------------------------------------------------
set(BUILD_TESTING OFF CACHE BOOL "Disable tests for subprojects" FORCE)
set(CMAKE_DISABLE_FIND_PACKAGE_PkgConfig ON CACHE BOOL "Avoid PkgConfig dependency for WFA2-lib" FORCE)

# The OPENMP option for WFA2-lib is OFF by default. Our main project enables OpenMP,
# so enable it here to ensure WFA2-lib correctly turns on WFA_PARALLEL.
set(OPENMP ON CACHE BOOL "Enable OpenMP for WFA2-lib" FORCE)

# Key: Only build the core library, disable WFA2-lib's built-in benchmark/tests.
# (These two switches are added by us in extern/WFA2-lib/CMakeLists.txt, default ON; here we force them off in the main project.)
set(WFA2LIB_BUILD_BENCHMARK OFF CACHE BOOL "Do not build WFA2-lib benchmark" FORCE)
set(WFA2LIB_BUILD_TESTS OFF CACHE BOOL "Do not register WFA2-lib tests" FORCE)

# Note: This CMakeLists is under test/, so we need to go back to the project root before entering extern/WFA2-lib.
# CMake requires: when add_subdirectory() includes a "non-current-directory subdirectory", the binary dir must be explicitly specified,
# otherwise it will error: "not a subdirectory ... must be explicitly specified".
#
# Here we fix the WFA2-lib build directory inside the test build tree to avoid polluting the main project build directory,
# and to avoid multiple add_subdirectory calls conflicting with each other.
add_subdirectory(${HALIGN4_ROOT}/extern/WFA2-lib ${CMAKE_CURRENT_BINARY_DIR}/wfa2lib-build)

# Note: We cannot directly link to the main project target `halign4` (not defined in this CMakeLists),
# otherwise configuration will fail. We link OpenMP + WFA2 static libraries into the test support library,
# so halign4_tests obtains these dependencies indirectly through halign4_test_support, keeping test logic unchanged.
target_link_libraries(halign4_test_support PRIVATE OpenMP::OpenMP_CXX wfa2::wfa2cpp_static)


add_executable(halign4_tests
        "${CMAKE_CURRENT_LIST_DIR}/doctest_main.cpp"
        "${CMAKE_CURRENT_LIST_DIR}/test_consensus.cpp"
        "${CMAKE_CURRENT_LIST_DIR}/test_read_fasta.cpp"
        "${CMAKE_CURRENT_LIST_DIR}/test_minimizer.cpp"
        "${CMAKE_CURRENT_LIST_DIR}/test_mash.cpp"
        "${CMAKE_CURRENT_LIST_DIR}/test_jaccard.cpp"
        "${CMAKE_CURRENT_LIST_DIR}/test_fasta_writer.cpp"
        "${CMAKE_CURRENT_LIST_DIR}/test_align.cpp"
        "${CMAKE_CURRENT_LIST_DIR}/test_anchor.cpp"
)

target_link_libraries(halign4_tests PRIVATE halign4_test_support)

# The Release optimization option is also applied to the test executable.
if (MSVC)
    target_compile_options(halign4_tests PRIVATE
        $<$<CONFIG:Release>:/O2>
        $<$<CONFIG:Release>:/arch:AVX2>
        $<$<CONFIG:Release>:/GL>
    )
    target_link_options(halign4_tests PRIVATE
        $<$<CONFIG:Release>:/LTCG>
    )
else()
    target_compile_options(halign4_tests PRIVATE
        $<$<CONFIG:Release>:-O3>
        $<$<CONFIG:Release>:-march=native>
        $<$<CONFIG:Release>:-flto>
        $<$<CONFIG:Release>:-pipe>
    )
    target_link_options(halign4_tests PRIVATE
        $<$<CONFIG:Release>:-flto>
    )
endif()

# It is recommended to enable it explicitly (even if include(CTest) has already been done, there will be no conflict).
enable_testing()

# Register a test: Run this executable file
# add_test(NAME halign4_tests COMMAND halign4_tests)

add_test(NAME consensus
        COMMAND halign4_tests --test-suite=consensus)

add_test(NAME read_fasta
        COMMAND halign4_tests --test-suite=read_fasta)


add_test(NAME minimizer
        COMMAND halign4_tests --test-suite=minimizer)

add_test(NAME mash
        COMMAND halign4_tests --test-suite=mash)

add_test(NAME jaccard
        COMMAND halign4_tests --test-suite=jaccard)

add_test(NAME write_fasta
        COMMAND halign4_tests --test-suite=write_fasta)

# Align Correctness Test (Quick)
add_test(NAME align
        COMMAND halign4_tests --test-suite=align)

# Align performance test (time-consuming, includes high similarity sequence testing)
add_test(NAME align_perf
        COMMAND halign4_tests --test-suite=align_perf)

# Anchor correctness test (quick)
add_test(NAME anchor
        COMMAND halign4_tests --test-suite=anchor)



