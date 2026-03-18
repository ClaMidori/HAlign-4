#!/usr/bin/env bash
set -euxo pipefail

# Do not use -march=native in conda-build environment (portability disaster)
# If you added HALIGN4_NATIVE_ARCH option in CMakeLists as suggested, keep it OFF here
# Also avoid injecting native flags into sub-projects
export CFLAGS="${CFLAGS} -O3"
export CXXFLAGS="${CXXFLAGS} -O3"

cmake -S . -B build \
  -DCMAKE_BUILD_TYPE=Release \
  -DCMAKE_INSTALL_PREFIX="${PREFIX}" \
  -DCMAKE_PREFIX_PATH="${PREFIX}" \
  -DCMAKE_INSTALL_LIBDIR=lib \
  -DBUILD_TESTING=OFF \
  -DWFA2LIB_BUILD_BENCHMARK=OFF \
  -DWFA2LIB_BUILD_TESTS=OFF \
  -DHALIGN4_NATIVE_ARCH=OFF

cmake --build build -j "${CPU_COUNT}"

# Your project currently has no install() rules, so manually install executables here
install -d "${PREFIX}/bin"
install -m 0755 build/halign4 "${PREFIX}/bin/halign4"
