[![CMake](https://github.com/zhucongx/LatticeMonteCarlo/actions/workflows/cmake.yml/badge.svg?branch=main)](https://github.com/zhucongx/LatticeMonteCarlo/actions/workflows/cmake.yml)


# LatticeMC

LatticeMC is a lattice Monte Carlo framework for multicomponent alloys. It supports Kinetic Monte Carlo (KMC), Canonical Monte Carlo (CMC), Simulated Annealing (SA), and post-processing tools (Analysis/Reformat) for logs and configurations.

Core idea: total energy and migration barriers are predicted from JSON coefficient files (similar to cluster-expansion or ML models). Local updates and LRU caches are used to accelerate event evaluation. The project also provides extended cfg/xyz I/O (with optional gzip) and rich analysis output.

### Features

- Multiple simulation methods: `KineticMc*` (KMC with OMP/MPI variants), `CanonicalMc*` (CMC), and `SimulatedAnnealing`.
- Predictors:
  - Total energy: `EnergyPredictor`
  - Vacancy migration barriers: `VacancyMigrationPredictor*` (LRU-backed variant included)
  - Local energy changes: `EnergyChangePredictorPairSite`
  - Time–temperature interpolation
- Parallelization: OpenMP and MPI (some algorithms mix OMP + MPI)
- I/O: `.cfg`/`.xyz` (gz supported), `lattice/element/map` triplet, logging, and analysis outputs

### Requirements

- C++ compiler with C++20 support (GCC/Clang/ICC/ICX/MPICXX)
- CMake ≥ 3.14
- Required libraries:
  - Boost (filesystem, iostreams)
  - OpenMP (enabled by most methods)
  - MPI (required for `*Mpi/*Ompi` methods only)
- Auto-fetched (via FetchContent):
  - nlohmann_json
  - Eigen3 (if not found on system)
- Linear algebra (optional):
  - MKL or generic BLAS/LAPACK. The build prefers MKL; if not found it falls back to generic BLAS/LAPACK, and otherwise to Eigen-only.

Note (macOS): with Apple Clang you typically need `libomp` for OpenMP (`brew install libomp`). Make sure CMake detects `OpenMP::OpenMP_CXX` successfully.

### Build

Using CMake (recommended):

- `cmake -S . -B build -DCMAKE_BUILD_TYPE=Release`
- `cmake --build build -j`
- Executable is placed at `bin/lmc.exe`

Using the script:

- `./build.sh r d`  # Release + default compiler
- Examples: `./build.sh d g` (Debug + gcc/g++), `./build.sh r m` (Release + mpicc/mpicxx)

Tips:

- On macOS, installing `libomp` helps CMake locate OpenMP.
- oneAPI + MKL: source your oneAPI environment (e.g., `source /opt/intel/oneapi/setvars.sh`) so that `MKLROOT` is visible to CMake.

### Tests

The unit tests are optional and need no input files. They download GoogleTest at configure time:

- `cmake -S . -B build -DCMAKE_BUILD_TYPE=Release -DLMC_BUILD_TESTS=ON`
- `cmake --build build -j`
- `ctest --test-dir build --output-on-failure` (or run `bin/lmc_test.exe` directly)

Test sources live next to the code they cover, in `lmc/<module>/test/`.

### Run

The executable is `bin/lmc.exe`. Provide a parameter file via `-p <file>`.

Sample parameter files are under `script/`:

- KMC: `script/kmc_param.txt`
- CMC: `script/cmc_param.txt`
- Simulated annealing: `script/sa_param.txt`
- Analysis: `script/ansys_param.txt`
- Reformat (batch cfg ↔ map conversion): `script/format_param.txt`

The number of MPI processes is fixed by the method, and the program stops with a message if it is started with another number:

| `simulation_method` | MPI processes | OpenMP threads |
|---|---|---|
| `KineticMcFirstSerial` | 1 | at start-up only |
| `KineticMcFirstMpi` | 12 (one per jump of the vacancy) | no |
| `KineticMcChainOmpi` | 12 (one per jump of the vacancy) | yes |
| `CanonicalMcSerial`, `SimulatedAnnealing` | 1 | no |
| `CanonicalMcOmp`, `Ansys`, `Reformat` | 1 | yes |

A step costs a few microseconds, which is less than it costs to start a group of threads or to exchange messages between processes. The single-process methods are therefore the fastest ones, and spare cores are better spent on independent runs. Steps per second for 10976 sites (Al with 5 % Mg and Zn) at 450 K on an Apple M3 Max:

| `simulation_method` | processes × threads | steps per second |
|---|---|---|
| `KineticMcFirstSerial` | 1 × any | 270 000 |
| `KineticMcFirstMpi` | 12 × 1 | 115 000 to 145 000 |
| `KineticMcChainOmpi` | 12 × 1 | 13 000 to 18 000 |
| `CanonicalMcSerial` | 1 × 1 | 320 000 |
| `CanonicalMcOmp` | 1 × 12 | 130 000 |

Memory grows with the number of lattice sites: the tables of the barrier model, which the KMC methods use, take about 7 kB per site, and those of the swap model of CMC and simulated annealing about 2.4 kB per site.

Examples:

- KMC with the sample file, which selects `KineticMcChainOmpi`:
  - `mpirun -np 12 bin/lmc.exe -p script/kmc_param.txt`

- Single-process KMC (set `simulation_method KineticMcFirstSerial` in the parameter file):
  - `bin/lmc.exe -p script/kmc_param.txt`

- CMC:
  - `bin/lmc.exe -p script/cmc_param.txt`

- Simulated annealing:
  - `bin/lmc.exe -p script/sa_param.txt`

- Analysis:
  - `bin/lmc.exe -p script/ansys_param.txt`

- Reformat:
  - `bin/lmc.exe -p script/format_param.txt`

The program exits with status 0 only if the run completed. Any problem with the input is reported as one `Error: ...` line and a non-zero exit status.

### Parameter File Cheatsheet

One `key value` pair per line, separated by blanks or tabs. Text after `#` is a comment. A key that is not known, a value of the wrong type (integers are written with digits only, booleans as `true`/`false`) and a missing key that the chosen method needs are errors.

Common keys (only relevant ones are used by each method):

- `simulation_method`: one of `KineticMcFirstSerial`, `KineticMcFirstMpi`, `KineticMcChainOmpi`, `CanonicalMcSerial`, `CanonicalMcOmp`, `SimulatedAnnealing`, `Ansys`, `Reformat`
- `json_coefficients_filename`: JSON coefficients for energy/barrier models
- `config_filename` / `map_filename`: input configuration (.cfg) or the triplet (lattice/element/map). Choose one.
- `element_set`: element list, e.g., `Al Mg Zn`
- `temperature` / `initial_temperature`: in Kelvin. `temperature` may be omitted for KMC when a time-temperature table is given.
- `log_dump_steps` / `config_dump_steps` / `maximum_steps` / `thermodynamic_averaging_steps`
- `restart_steps` / `restart_energy` / `restart_time`
- `seed`: seed of the random number generator (optional). Two runs with the same seed, executable and thread/process count follow the same trajectory. If the key is omitted a time-based seed is used and printed at startup, so any run can be repeated by adding `seed <value>`. A restarted run continues with a fresh random stream derived from `seed` and `restart_steps`.

KMC-specific:

- `time_temperature_filename`: time-temperature table (optional). Every line that starts with two numbers is a `time temperature` point, other lines are ignored; the temperature is interpolated linearly and held constant outside the table.
- `early_stop`: `true|false` (stop when vacancy neighborhood becomes pure solvent)
- `solute_disp`: `true|false` (record solute center-of-mass displacement)

SA-specific:

- `factor`: FCC supercell factor (x=y=z=factor)
- `solvent_element` / `solute_element_set` / `solute_number_set`

Analysis/Reformat-specific:

- `log_type`: `kinetic_mc` | `canonical_mc` | `simulated_annealing`
- `config_type`: `config` | `map`
- `initial_steps` / `increment_steps`
- `smallest_cluster_criteria` / `solvent_bond_criteria` / `escape_temperature`

See `script/` for ready-to-use samples.

### Outputs

- Logs:
  - KMC: `kmc_log.txt`
  - CMC: `cmc_log.txt`
  - SA: `sa_log.txt`
- Configurations: `N.cfg.gz` at dump intervals and `end.cfg.gz` at completion. Reformat supports writing/reading `lattice.txt` / `element.txt` / `mapN.txt`.
- Analysis: `ansys_frame_log.txt`, `ansys_cluster_log.txt`, and `ansys/*.xyz.gz` files for visualization/statistics.

### FAQ / Tips

- MKL/BLAS: If MKL isn’t found, the build falls back to generic BLAS/LAPACK, and then to Eigen-only (lowest performance). Setting `MKLROOT` can help detection.
- OpenMP (macOS): install `libomp` and ensure CMake finds `OpenMP::OpenMP_CXX`.
- MPI: only required for `*Mpi/*Ompi` methods.
- Executable name: the output is `bin/lmc.exe` for consistency across platforms.
- Large model files: JSON coefficient files can be large; consider hosting externally or use Git LFS.

### Repository Layout (high level)

- `lmc/`: sources (`cfg`/`pred`/`mc`/`api`/`ansys`)
- `script/`: sample parameter files and helper scripts
- `bin/`: executable output directory
- `lmc/<module>/test/`: unit tests
- `test/`: local example runs and large model files (not tracked, not unit tests)
