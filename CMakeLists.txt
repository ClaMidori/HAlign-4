cmake_minimum_required(VERSION 3.18)
project(HAlign4)

set(CMAKE_CXX_STANDARD 20)
option(HALIGN4_NATIVE_ARCH "Enable -march=native" OFF)

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
        src/utils/msa_validation.cpp
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
        src/align/cigar.cpp
        src/align/align.cpp
        src/align/native_msa.cpp
        src/align/strategy.cpp
        src/align/homology_plan.cpp
        src/align/lossless_extractor.cpp
        src/align/block_local_msa.cpp
        src/ultralong/native_pipeline.cpp
        src/align/profile_projection.cpp
        src/align/ordered_sam.cpp
        src/align/anchored_insertions.cpp
        src/seed/minimizer.cpp
        src/seed/anchor.cpp

)

# Diagnostic extractor can also be built independently without abPOA/WFA2.
add_executable(halign4_pr3d3a_extract
    tools/pr3d3a_native_cli.cpp
    src/align/lossless_extractor.cpp)
target_include_directories(halign4_pr3d3a_extract PRIVATE "${PROJECT_SOURCE_DIR}/include")

# Static libcurl distributions can have missing transitive dependencies.
# Optional network download support stays available through curl/wget executables.
option(HALIGN4_WITH_CURL "Link HAlign4 to libcurl (OFF uses curl/wget fallback)" OFF)
if (HALIGN4_WITH_CURL)
    find_package(CURL REQUIRED)
endif()
if (HALIGN4_WITH_CURL AND CURL_FOUND)
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

# Critical: only build the core library, turn off WFA2-lib's bundled benchmark/tests.
# (These two options are added in extern/WFA2-lib/CMakeLists.txt and default to ON; we force them off here.)
set(WFA2LIB_BUILD_BENCHMARK OFF CACHE BOOL "Do not build WFA2-lib benchmark" FORCE)
set(WFA2LIB_BUILD_TESTS OFF CACHE BOOL "Do not register WFA2-lib tests" FORCE)

add_subdirectory(extern/WFA2-lib)

# Integrate the version-pinned abPOA submodule directly. Do not build its CLI
# during a normal HAlign-4 build. This is a library-only integration at PR 2.
if(NOT EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/extern/abPOA/include/abpoa.h")
    message(FATAL_ERROR "abPOA submodule missing: run git submodule update --init --recursive")
endif()
add_subdirectory(extern/abPOA ${CMAKE_CURRENT_BINARY_DIR}/abpoa-build EXCLUDE_FROM_ALL)
target_include_directories(halign4 PRIVATE "${CMAKE_CURRENT_SOURCE_DIR}/extern/abPOA/include")
target_link_libraries(halign4 PRIVATE OpenMP::OpenMP_CXX wfa2::wfa2cpp_static abpoa)

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

# PR3D.4D.1: 64-bit sparse out-of-core FASTA access, a separate indexing tool.
# No change to the proven --ultralong aligner or legacy commands.
add_library(halign4_sparse_fasta STATIC src/ultralong/sparse_fasta.cpp)
target_include_directories(halign4_sparse_fasta PUBLIC "${PROJECT_SOURCE_DIR}/include")
add_executable(halign4_fasta64_index tools/fasta64_cli.cpp)
target_link_libraries(halign4_fasta64_index PRIVATE halign4_sparse_fasta)

# PR3D.4D.2: independent, 64-bit index-backed homology WINDOW PLANNER only.
# The proven --ultralong MSA backend remains unchanged in this PR.
add_library(halign4_window64 STATIC
    src/ultralong/window64_plan.cpp
    src/align/homology_plan.cpp)
target_include_directories(halign4_window64 PUBLIC "${PROJECT_SOURCE_DIR}/include")
target_link_libraries(halign4_window64 PUBLIC halign4_sparse_fasta)
add_executable(halign4_window64_plan tools/window64_cli.cpp)
target_link_libraries(halign4_window64_plan PRIVATE halign4_window64)

# PR3D.4D.3: real out-of-core MSA (index-backed windows -> local abPOA -> streamed FASTA).
# The original memory-resident --ultralong and legacy behavior are not altered.
target_sources(halign4 PRIVATE src/ultralong/window64_msa.cpp)
target_link_libraries(halign4 PRIVATE halign4_window64)
add_executable(halign4_window64_msa
  tools/window64_msa_cli.cpp
  src/ultralong/window64_msa.cpp
  src/align/native_msa.cpp)
target_include_directories(halign4_window64_msa PRIVATE
  "${PROJECT_SOURCE_DIR}/include" "${PROJECT_SOURCE_DIR}/extern/abPOA/include")
target_link_libraries(halign4_window64_msa PRIVATE
  halign4_window64 abpoa OpenMP::OpenMP_CXX)
if(ZLIB_FOUND)
  target_link_libraries(halign4_window64_msa PRIVATE ZLIB::ZLIB)
endif()
