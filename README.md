# HAlign-4

[![Downloads](https://anaconda.org/malab/halign4/badges/downloads.svg)](https://anaconda.org/malab/halign4)
[![License](https://anaconda.org/malab/halign4/badges/license.svg)](https://anaconda.org/malab/halign4)
[![Platforms](https://anaconda.org/malab/halign4/badges/platforms.svg)](https://anaconda.org/malab/halign4)

[HAlign 4: A New Strategy for Rapidly Aligning Millions of Sequences.](https://doi.org/10.1093/bioinformatics/btae718)

Documentation:

- Detailed usage & examples: [`docs/usage.md`](docs/usage.md)
- Source install & dependencies: [`docs/install.md`](docs/install.md)
- Tests: [`docs/test.md`](docs/test.md)

---

## Install (Conda)

Conda is the recommended installation method for end users.

```bash
conda install -c malab halign4
```

Verify:

```bash
halign4 --version
halign4 -h
```

Source installation: see [`docs/install.md`](docs/install.md).

---

## Install from source

Requirements: CMake >= 3.18, C++20 compiler, OpenMP, and Make or Ninja. `zlib` is recommended for gzipped FASTA input. `libcurl` is optional for URL input.

```bash
git clone https://github.com/ClaMidori/HAlign-4.git
cd HAlign-4
```

Build normal mode:

```bash
cmake -S . -B build-original \
  -DCMAKE_BUILD_TYPE=Release \
  -DHALIGN4_ENABLE_OFFSET_COMPACT=OFF
cmake --build build-original -j
```

Build compact mode:

```bash
cmake -S . -B build-compact \
  -DCMAKE_BUILD_TYPE=Release \
  -DHALIGN4_ENABLE_OFFSET_COMPACT=ON
cmake --build build-compact -j
```

Compact mode enables experimental WFA offset compaction. Normal mode uses standard WFA2 behavior. Keep separate build directories when switching modes.

Verify either binary:

```bash
./build-original/halign4 --version
./build-compact/halign4 --version
```

See [`docs/install.md`](docs/install.md) for dependencies and external MSA tools.

---

## Quick start

The repository includes small datasets under `test/data/` which are perfect for a first run.

Normal mode:

```bash
./build-original/halign4 \
  -i test/data/mt1x.fasta.gz \
  -o mt1x-normal.fasta
```

Compact mode:

```bash
./build-compact/halign4 \
  -i test/data/mt1x.fasta.gz \
  -o mt1x-compact.fasta
```

---

## Parameters (overview)

The most important parameters are:

- `-i/--input`: input FASTA (required)
- `-o/--output`: output aligned FASTA (required)
- `-w/--workdir`: working directory (optional; default: `./tmp-<random>`)
- `-p/--msa-cmd`: MSA method (keyword: `minipoa`/`mafft`/`clustalo`, or a custom template)
- `-c/--center-path`: provide a reference/center FASTA (optional)
- `--keep-length`: keep reference length coordinate rules

For the full parameter list and detailed examples, see [`docs/usage.md`](docs/usage.md).

---

## Tests

Run the standard test suite:

```bash
cd test
./run_tests.sh -t Release -j 8
```

Run performance tests:

```bash
./run_tests.sh -t Release --perf
```

See [`docs/test.md`](docs/test.md) for filters, clean builds, and advanced CTest usage.

---

## Citation

If you use HAlign-4 in academic work, please cite:

HAlign 4: a new strategy for rapidly aligning millions of sequences. Bioinformatics, 2024, 40(12): btae718. https://doi.org/10.1093/bioinformatics/btae718
