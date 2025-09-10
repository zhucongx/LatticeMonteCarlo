[![CMake](https://github.com/zhucongx/LatticeMC/actions/workflows/cmake.yml/badge.svg?branch=main)](https://github.com/zhucongx/LatticeMC/actions/workflows/cmake.yml)


# LatticeMC

LatticeMC is a lattice Monte Carlo framework for multicomponent alloys. It supports Kinetic Monte Carlo (KMC), Canonical Monte Carlo (CMC), Simulated Annealing (SA), and post-processing tools (Analysis/Reformat) for logs and configurations.

Core idea: total energy and migration barriers are predicted from JSON coefficient files (similar to cluster-expansion or ML models). Local updates and LRU caches are used to accelerate event evaluation. The project also provides extended cfg/xyz I/O (with optional gzip) and rich analysis output.

### Features

- Multiple simulation methods: `KineticMc*` (KMC with OMP/MPI variants), `CanonicalMc*` (CMC), and `SimulatedAnnealing`.
- Predictors:
  - Total energy: `EnergyPredictor`
  - Vacancy migration barriers: `VacancyMigrationPredictor*` (LRU-backed variant included)
  - Local energy changes: `EnergyChangePredictorPairSite`
  - Time–temperature interpolation and a rate corrector
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

### Run

The executable is `bin/lmc.exe`. Provide a parameter file via `-p <file>`.

Sample parameter files are under `script/`:

- KMC: `script/kmc_param.txt`
- CMC: `script/cmc_param.txt`
- Simulated annealing: `script/sa_param.txt`
- Analysis: `script/ansys_param.txt`
- Reformat (batch cfg ↔ map conversion): `script/format_param.txt`

Examples:

- Single-process KMC:
  - `bin/lmc.exe -p script/kmc_param.txt`

- MPI-parallel KMC (4 processes):
  - `mpirun -np 4 bin/lmc.exe -p script/kmc_param.txt`

- CMC:
  - `bin/lmc.exe -p script/cmc_param.txt`

- Simulated annealing:
  - `bin/lmc.exe -p script/sa_param.txt`

- Analysis:
  - `bin/lmc.exe -p script/ansys_param.txt`

- Reformat:
  - `bin/lmc.exe -p script/format_param.txt`

### Parameter File Cheatsheet

Common keys (only relevant ones are read by each method):

- `simulation_method`: one of `KineticMcFirstOmp`, `KineticMcFirstMpi`, `KineticMcChainOmpi`, `CanonicalMcSerial`, `CanonicalMcOmp`, `SimulatedAnnealing`, `Ansys`, `Reformat`
- `json_coefficients_filename`: JSON coefficients for energy/barrier models
- `config_filename` / `map_filename`: input configuration (.cfg) or the triplet (lattice/element/map). Choose one.
- `element_set`: element list, e.g., `Al Mg Zn`
- `temperature` / `initial_temperature`: in Kelvin
- `log_dump_steps` / `config_dump_steps` / `maximum_steps` / `thermodynamic_averaging_steps`
- `restart_steps` / `restart_energy` / `restart_time`

KMC-specific:

- `time_temperature_filename`: time-temperature table (optional)
- `rate_corrector`: `true|false`
- `vacancy_trajectory`: initial vacancy displacement (3 numbers)
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

See `script/` for ready-to-use samples. The `script/kmc_param.txt` file has been verified to contain proper line endings so that `solute_disp` doesn’t merge with the next key.

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
- `test/`: example/large model files (not unit tests)
