# Sources of the figures

Every figure of the documentation can be made again from the files in this directory. The comment at the top of a file says which figure it belongs to and gives the commands. Intermediate files go to a scratch directory, not here.

| Figures | Made with |
|---|---|
| `files`, `kmc_step`, `kmc_second_order`, `metropolis`, `cluster_definition`, `barrier_groups`, `energy_change`, `path`, `restart` | `fig_<name>.tex` |
| `lattice`, `cluster_shapes`, `environment`, `frame_fields`, `annealing` | renderings (`render_*.py`), put together in `fig_<name>.tex` |
| `time_temperature`, `log_rows`, `log_plot`, `vacancy_concentration` | `plots.py` |
| `clustering`, `vacancy_jumps`, `second_order`, `annealing`, `jump_path`, `environment` (animations) | `animate_<name>.py`, most of them from renderings |

`figstyle.tex`, `style.py` and `scene.py` hold the look that the line drawings, the plots and the renderings share. `pdf_to_png.py` and `pdf_to_png_dithered.py` make the PNG of a still from its PDF. `alloy_time.py` gives the factor between the time of a KMC run and the [time in an alloy](../../methods.md#time-and-vacancy-concentration), for the animations and for `plots.py`.

Needed: LaTeX, Ghostscript, Python with matplotlib, seaborn, pandas, SciPy, scikit-learn and Pillow, the matplotlib style `xzc_sans` (not part of the repository), Blender, the WebP tools (`img2webp`, `cwebp`, `webpmux`), and `ffmpeg` for the videos.

The data come from the runs of the [first example](../../README.md#a-first-run), from an `Ansys` run on its KMC part, and from the simulated annealing run on the page [Methods](../../methods.md#simulatedannealing). For the animations these runs were repeated with more output, or continued; the comments at the top of the scripts say how. The start configuration and the coefficient file are not part of the repository.
