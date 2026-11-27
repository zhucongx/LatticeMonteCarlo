# Parameters

The program reads one parameter file, given with `-p <file>`. The directory `script/` holds a sample for every kind of run.

## Rules of the file

- One key per line, followed by its value. Blanks or tabs separate them, as many as you like.
- Text from a word that begins with `#` to the end of the line is a comment. Blank lines are allowed.
- File names are relative to the directory the program is started in. They cannot contain blanks.
- Whole numbers are written with digits only: `1000000`, not `1e6`. Other numbers may have a decimal point and an exponent, but no unit.
- `true`, `false`, `1` and `0` are the values of a switch. Upper and lower case do not matter.
- The entries of a list stand on one line, separated by blanks.
- A key that is not known, a value that does not fit its key, and a missing key that the method needs stop the program before the run begins. A key that the method does not use has no effect.
- A key that is given twice takes its last value, and the program prints a warning.

## Keys of every method

| Key | Value | Default | Meaning |
|---|---|---|---|
| `simulation_method` | name | required | What to do, see [Methods](methods.md) |
| `json_coefficients_filename` | file | required | The coefficients of the model, see [Input files](input-files.md#coefficient-file) |

## Start configuration

For KMC and CMC. One of the two file names is required; if both are given, `map_filename` is used.

| Key | Value | Default | Meaning |
|---|---|---|---|
| `config_filename` | file | | A [configuration](input-files.md#configuration), `.cfg` or `.cfg.gz` |
| `map_filename` | file | | A [map](input-files.md#lattice-element-and-map). The files `lattice.txt` and `element.txt` have to be in the directory the program is started in |
| `element_set` | list | required | The elements of the model, without the vacancy. An element of the configuration that is missing here is an error |

## Steps and output

| Key | Value | Default | Methods | Meaning |
|---|---|---|---|---|
| `maximum_steps` | integer | required | KMC, CMC, SA | Number of the last step. A run makes the steps 0 to `maximum_steps` |
| `log_dump_steps` | integer > 0 | required | KMC, CMC, SA | Steps between two rows of the log. In KMC and CMC the first rows are closer together, see [Output files](output-files.md#when-a-row-is-written) |
| `config_dump_steps` | integer > 0 | required | KMC, CMC | Steps between two configurations. Use a power of ten for `log_dump_steps` and a multiple of it here: a restart and `Ansys` need a [row](output-files.md#when-a-row-is-written) of the log for the step of a configuration |
| `thermodynamic_averaging_steps` | integer | 0 | CMC | Number of steps over which `average_energy` is taken. With 0 the column stays 0 |
| `seed` | integer | from the clock | KMC, CMC, SA | Seed of the random numbers, see [Restarting](restart.md#repeating-a-run) |

## Temperature

| Key | Value | Default | Methods | Meaning |
|---|---|---|---|---|
| `temperature` | K, > 0 | required | KMC, CMC | Temperature of the run. May be left out for KMC with a time-temperature table |
| `time_temperature_filename` | file | none | KMC | [Table](input-files.md#time-temperature-table) of the temperature against the time. It replaces `temperature` from the first step on |
| `initial_temperature` | K, > 0 | required | SA | Temperature at the start of the [schedule](methods.md#simulatedannealing) |

## Restart

| Key | Value | Default | Methods | Meaning |
|---|---|---|---|---|
| `restart_steps` | integer | 0 | KMC, CMC | Step at which the run goes on. A value above 0 makes it a [restart](restart.md): the log is continued instead of overwritten |
| `restart_energy` | eV | 0 | KMC, CMC | Value of the `energy` column at that step |
| `restart_time` | s | 0 | KMC | Value of the `time` column at that step |

## KMC

| Key | Value | Default | Meaning |
|---|---|---|---|
| `early_stop` | switch | `false` | Stop when all neighbors of the vacancy up to the third shell are solvent atoms |
| `solute_disp` | switch | `false` | Write the displacement of the center of mass of the solutes to the log |

## SA

| Key | Value | Default | Meaning |
|---|---|---|---|
| `factor` | integer > 2 | required | The box has `factor` × `factor` × `factor` unit cells, which is 4 `factor`³ sites |
| `solvent_element` | name | required | Element on the sites that hold no solute |
| `solute_element_set` | list | required | The solute elements: all elements of the coefficient file but the solvent. The vacancy `X` can be added, with the number 1 |
| `solute_number_set` | list | required | Number of atoms of each of them, in the same order. 0 is allowed |

## Ansys and Reformat

| Key | Value | Default | Methods | Meaning |
|---|---|---|---|---|
| `log_type` | name | required | both | Which log the run wrote: `kinetic_mc`, `canonical_mc` or `simulated_annealing` |
| `config_type` | name | required | both | `config` or `map`: the kind of file to read |
| `element_set` | list | required | both | The elements of the model |
| `solvent_element` | name | required | both | The solvent |
| `initial_steps` | integer | 0 | both | First step to read |
| `increment_steps` | integer > 0 | required | both | Steps between two frames |
| `smallest_cluster_criteria` | integer | required | `Ansys` | Smallest number of connected solute atoms and vacancies that is a cluster |
| `solvent_bond_criteria` | integer | required | `Ansys` | A solvent atom with more neighbors than this among them belongs to the cluster |
| `escape_temperature` | K, > 0 | required | `Ansys` | Temperature of the escape times |

[Analysis](analysis.md) explains the last five.

## Keys without effect

These keys are accepted, so that older parameter files can still be read, but nothing uses their value: `decrement_temperature`, `early_stop_steps`, and, with a warning, `vacancy_trajectory` and `rate_corrector`. `rate_corrector true` is an error, because the time is not rescaled any more, see [Methods](methods.md#time-and-vacancy-concentration). `config_dump_steps` has no effect on SA. In a KMC run, `thermodynamic_averaging_steps` only costs memory, up to 8 bytes times its value in every process.

## Required keys at a glance

| Method | Required besides `simulation_method` and `json_coefficients_filename` |
|---|---|
| `KineticMcFirstSerial`, `KineticMcFirstMpi`, `KineticMcSecondSerial`, `KineticMcSecondMpi`, `KineticMcSecondOmpi` | `config_filename` or `map_filename`, `element_set`, `log_dump_steps`, `config_dump_steps`, `maximum_steps`, and `temperature` or `time_temperature_filename` |
| `CanonicalMcSerial`, `CanonicalMcOmp` | `config_filename` or `map_filename`, `element_set`, `log_dump_steps`, `config_dump_steps`, `maximum_steps`, `temperature` |
| `SimulatedAnnealing` | `factor`, `solvent_element`, `solute_element_set`, `solute_number_set`, `log_dump_steps`, `maximum_steps`, `initial_temperature` |
| `Ansys` | `log_type`, `config_type`, `solvent_element`, `element_set`, `increment_steps`, `smallest_cluster_criteria`, `solvent_bond_criteria`, `escape_temperature` |
| `Reformat` | `log_type`, `config_type`, `solvent_element`, `element_set`, `increment_steps` |
