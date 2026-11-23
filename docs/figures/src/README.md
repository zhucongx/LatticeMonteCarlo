# Sources of the figures

Every figure of the documentation can be made again from the files in this directory. The comment at the top of a file says which figure it belongs to and gives the commands. Intermediate files go to a scratch directory, not here.

| Figures | Made with |
|---|---|
| `files` | `fig_files.tex`, then `pdf_to_png.py` |
| `log_rows`, `log_plot`, `vacancy_concentration` | `plots.py` |
| `kmc_step`, `kmc_second_order`, `second_order`, `metropolis`, `cluster_definition`, `barrier_groups`, `energy_change`, `restart` | `animate_<name>.py`, which draws every frame as a line drawing |
| `clustering`, `quench`, `vacancy_jumps`, `trapped_vacancy`, `annealing`, `jump_path`, `environment`, `lattice`, `cluster_shapes`, `frame_fields` | renderings of the `render_*.py` scripts, put together by `animate_<name>.py` |

The figures of the first two rows are stills. All others are animations: a WebP file for the pages and an MP4 file of the same content, which is not part of the repository.

Shared files:

- `figstyle.tex`, `style.py` and `scene.py` hold the look of the line drawings, the plots and the renderings.
- `tikz_film.py` and `tikz_frames.py` make an animation of line drawings, and `turn_film.py` one of renderings that turn.
- `solute_box.py` draws the box of a KMC run with its solute atoms.
- `alloy_time.py` gives the factor between the time of a KMC run and the [time in an alloy](../../methods.md#time-and-vacancy-concentration), for the animations and for `plots.py`.

Needed: LaTeX, Ghostscript, Python with matplotlib, seaborn, pandas, SciPy and Pillow, the matplotlib style `xzc_sans` (not part of the repository), Blender, the WebP tools (`img2webp`, `cwebp`, `webpmux`), and `ffmpeg`.

The data come from the runs of the [first example](../../README.md#a-first-run), from `Ansys` runs on its KMC part, from the run with the [time-temperature table](../../input-files.md#time-temperature-table), and from the simulated annealing run on the page [Methods](../../methods.md#simulatedannealing). For the animations these runs were repeated with more output, or continued; the comments at the top of the scripts say how. The start configuration and the coefficient file are not part of the repository.
