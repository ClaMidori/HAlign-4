cmake_minimum_required(VERSION 3.18)
project(HAlign4)

set(CMAKE_CXX_STANDARD 20)
option(HALIGN4_NATIVE_ARCH "Enable -march=native" OFF)
set(HALIGN4_WFA_OFFSET_BITS "32" CACHE STRING "WFA offset width in bits (16 or 32)")
set_property(CACHE HALIGN4_WFA_OFFSET_BITS PROPERTY STRINGS 16 32)

# If the user did not specify CMAKE_BUILD_TYPE, default to Release so optimizations are enabled by default
if(NOT DEFINED CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE "Release" CACHE STRING "Choose the type of build (default: Release)" FORCE)
    message(STATUS "CMAKE_BUILD_TYPE not set -> defaulting to: ${CMAKE_BUILD_TYPE}")
else()
    message(STATUS "CMAKE_BUILD_TYPE is: ${CMAKE_BUILD_TYPE}")
endif()


include_directories(
        ${PROJECT_SOURCE_DIR}/include
        ${PROJECT_SOURCE_DIR}/extern/CLI11/include
        ${PROJECT_SOURCE_DIR}/extern/spdlog/include
        ${PROJECT_SOURCE_DIR}/extern/cereal/include
        ${PROJECT_SOURCE_DIR}/extern/WFA2-lib
)

add_executable(halign4
        src/halign4.cpp
        src/utils/preprocess.cpp
        src/utils/seq_io.cpp
        src/utils/cmd.cpp
        src/utils/file_io.cpp
        src/mash/hash.cpp
        src/mash/xxhash.c
        src/mash/mash.cpp
        src/consensus/selector.cpp
        src/consensus/consensus.cpp
        src/align/ref_aligner.cpp
        src/align/ksw2_extz2_sse.c
        src/align/ksw2_gg2_sse.c
        src/align/psw_gg3_sse.c
        src/align/cigar.cpp
        src/align/align.cpp
        src/seed/minimizer.cpp
        src/seed/anchor.cpp

)

find_package(CURL)
if (CURL_FOUND)
    target_link_libraries(halign4 PRIVATE CURL::libcurl)
    target_compile_definitions(halign4 PRIVATE HALIGN4_HAVE_LIBCURL=1)
else()
    message(STATUS "libcurl not found; downloadFile will fall back to calling curl/wget at runtime")
endif()


find_package(ZLIB)
if (ZLIB_FOUND)
    target_link_libraries(halign4 PRIVATE ZLIB::ZLIB)
    target_compile_definitions(halign4 PRIVATE HALIGN4_HAVE_ZLIB=1)
else()
    message(STATUS "zlib not found; compression support will be disabled or limited")
endif()

find_package(OpenMP REQUIRED)

# ---------------------------------------------------------------------------
# Trim down the WFA2-lib (Wavefront Alignment) subproject: only build the core library, not tests/benchmarks
# ---------------------------------------------------------------------------
set(BUILD_TESTING OFF CACHE BOOL "Disable tests for subprojects" FORCE)
set(CMAKE_DISABLE_FIND_PACKAGE_PkgConfig ON CACHE BOOL "Avoid PkgConfig dependency for WFA2-lib" FORCE)

# WFA2-lib's OpenMP option is OFF by default. Since our top-level project enables OpenMP,
# explicitly turn it on here so WFA_PARALLEL is enabled inside WFA2-lib.
set(OPENMP ON CACHE BOOL "Enable OpenMP for WFA2-lib" FORCE)
set(WFA2_OFFSET_BITS ${HALIGN4_WFA_OFFSET_BITS} CACHE STRING "WFA2 internal offset width (16 or 32)" FORCE)

# Critical: only build the core library, turn off WFA2-lib's bundled benchmark/tests.
# (These two options are added in extern/WFA2-lib/CMakeLists.txt and default to ON; we force them off here.)
set(WFA2LIB_BUILD_BENCHMARK OFF CACHE BOOL "Do not build WFA2-lib benchmark" FORCE)
set(WFA2LIB_BUILD_TESTS OFF CACHE BOOL "Do not register WFA2-lib tests" FORCE)

add_subdirectory(extern/WFA2-lib)

target_link_libraries(halign4 PRIVATE OpenMP::OpenMP_CXX wfa2::wfa2cpp_static)
target_compile_definitions(halign4 PRIVATE HALIGN4_WFA_OFFSET_BITS=${HALIGN4_WFA_OFFSET_BITS})

# Add aggressive optimization flags in Release builds (handled separately for GCC/Clang vs MSVC)
if (MSVC)
    # MSVC optimizations: /O2 enables optimizations, /arch:AVX2 requests AVX2 (change to AVX/AVX512 if desired), /GL enables whole program optimization
    target_compile_options(halign4 PRIVATE
        $<$<CONFIG:Release>:/O2>
        $<$<CONFIG:Release>:/arch:AVX2>
        $<$<CONFIG:Release>:/GL>
    )
    # Enable link-time code generation (LTCG)
    target_link_options(halign4 PRIVATE
        $<$<CONFIG:Release>:/LTCG>
    )
else()
    # GCC / Clang: -O3, -march=native, -flto, -pipe
    target_compile_options(halign4 PRIVATE
        $<$<CONFIG:Release>:-O3>
            $<$<AND:$<CONFIG:Release>,$<BOOL:${HALIGN4_NATIVE_ARCH}>>:-march=native>

    )
endif()


target_compile_definitions(halign4 PRIVATE
        $<$<CONFIG:Debug>:_DEBUG>
)
