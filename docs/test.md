# Tests

This document explains how to run HAlign-4 tests under the `test/` directory.

HAlign-4 uses **doctest** for unit tests and CTest for registration.

---

## 1. Quick start (recommended)

Use the helper script:

```bash
cd test
./run_tests.sh -t Release -j 8
```

What it does:

- configures the test build into `build-test/` (or your chosen `-B` directory)
- builds the test binaries
- runs `ctest` (verbose by default)

The test suite validates shared project behavior. To test both runtime modes, build and run it once for each configuration:

```bash
cd ..
cmake -S test -B build-test-original \
  -DCMAKE_BUILD_TYPE=Release \
  -DHALIGN4_ENABLE_OFFSET_COMPACT=OFF
cmake --build build-test-original -j
ctest --test-dir build-test-original -V

cmake -S test -B build-test-compact \
  -DCMAKE_BUILD_TYPE=Release \
  -DHALIGN4_ENABLE_OFFSET_COMPACT=ON
cmake --build build-test-compact -j
ctest --test-dir build-test-compact -V
```

The helper script uses its own test build directory and does not replace `build-original` or `build-compact`.

Performance cases can be slow and include hardware-dependent timing thresholds. Run the correctness tests without performance cases when validating a build on a constrained or shared machine:

```bash
ctest --test-dir build-test-compact -V -E 'align_perf'
```

---

## 2. Common workflows

### 2.1 Run a single CTest test by name

Example: run the `align` test only:

```bash
cd test
./run_tests.sh -t Release -- -- -R align
```

(Everything after `--` is passed to `ctest`.)

### 2.2 Enable perf tests

Some tests are heavier and are guarded by an environment variable.

```bash
cd test
./run_tests.sh -t Release --perf
```

This sets:

- `HALIGN4_RUN_PERF=1`

### 2.3 Build directory management

To build tests into a custom directory:

```bash
cd test
./run_tests.sh -B ../build-test-release -t Release -j 16
```

To clean rebuild:

```bash
cd test
./run_tests.sh --clean -t Release
```

The helper script does not expose project-specific CMake options. Configure compact-mode tests directly:

```bash
cmake -S test -B build-test-compact \
  -DCMAKE_BUILD_TYPE=Release \
  -DHALIGN4_ENABLE_OFFSET_COMPACT=ON
cmake --build build-test-compact -j
ctest --test-dir build-test-compact -V
```

For normal mode, use `-DHALIGN4_ENABLE_OFFSET_COMPACT=OFF` or omit the option.

---

## 3. Run CTest directly (advanced)

If you already have a configured build directory:

```bash
ctest --test-dir build-test -V
```

---

## 4. Using `test/data` as example inputs

The repository includes small datasets that are useful for smoke tests and documentation examples:

- `test/data/mt1x.fasta.gz` (minimal dataset)
- `test/data/covid-ref.fasta.gz` / `test/data/covid-test.fasta.gz`

For runnable CLI examples, see [`docs/usage.md`](usage.md).
