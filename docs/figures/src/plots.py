# Data plots of the documentation: time_temperature.png, log_rows.png, log_plot.png and vacancy_concentration.png.
#   python3 plots.py <run folder> <output folder>
# The run folder holds cmc/cmc_log.txt and kmc/kmc_log.txt, the logs of the two runs on the page README.md, and in ansys/ansys the files
# 0.xyz.gz and 10000000.xyz.gz of Ansys for the KMC run. Rows of the KMC log after step 10 000 000 are not used. Needs LaTeX, Ghostscript,
# pandas and Pillow besides the modules and the matplotlib style of style.py.
import gzip
import os
import re
import subprocess
import sys
import tempfile
from bisect import bisect_left

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import pandas as pd
from matplotlib.collections import LineCollection
from PIL import Image
from style import *

runs, out = sys.argv[1:3]
PX = 1600 / 10.5         # pixels of the PNG per inch: 14 pt text is 14 px high in a picture of 1600 px that is shown 760 px wide


def canvas(width, height, *panels):
    """ width of the PNG in pixels, height in inches; panels: left, bottom, width, height of the axes in inches """
    w = width / PX
    fig = plt.figure(figsize=(w, height))
    return fig, [fig.add_axes((l / w, b / height, aw / w, ah / height)) for l, b, aw, ah in panels]


def panel_label(ax, letter):
    ax.annotate(r"\textbf{(%s)}" % letter, (0, 1), (-65, 8), xycoords="axes fraction", textcoords="offset points",
                fontsize=18, ha="left", va="bottom")


def note(ax, x, y, text, size=14, **kw):
    """ text on a white box, so that the grid does not run through it """
    ax.text(x, y, text, fontsize=size, bbox=dict(facecolor="white", edgecolor="none", pad=1.5), **kw)


def save_png(fig, name, width):
    """ PNG on white with a palette of 256 colors """
    with tempfile.TemporaryDirectory() as tmp:
        pdf, raw = os.path.join(tmp, "figure.pdf"), os.path.join(tmp, "figure.png")
        fig.savefig(pdf, bbox_inches=fig.bbox_inches)
        subprocess.run(["gs", "-q", "-dNOPAUSE", "-dBATCH", "-sDEVICE=png16m", "-r450", "-dTextAlphaBits=4", "-dGraphicsAlphaBits=4",
                        "-sOutputFile=" + raw, pdf], check=True)
        im = Image.open(raw).convert("RGB")
        im = im.resize((width, round(im.height * width / im.width)), Image.LANCZOS)
    plt.close(fig)
    path = os.path.join(out, name)
    im.quantize(colors=256, method=Image.Quantize.MAXCOVERAGE, kmeans=4).save(path, optimize=True)
    print(name, im.size, os.path.getsize(path) // 1024, "kB")


# --------------------------------------------------------------------------------- time-temperature table (input-files.md)
TABLE = ((0, 750), (10, 300), (40, 300), (60, 400))


def temperature(time, table=TABLE):
    """ TimeTemperatureInterpolator::GetTemperature """
    i = bisect_left([point[0] for point in table], time)
    if i == len(table):
        return table[-1][1]
    if i == 0:
        return table[0][1]
    (t0, T0), (t1, T1) = table[i - 1], table[i]
    return T0 + (time - t0) / (t1 - t0) * (T1 - T0)


def color(T):
    return TEMP(TNORM(T))


fig, (ax,) = canvas(1100, 3.86, (0.85, 0.71, 6.25, 3.02))
t = np.unique(np.concatenate([np.linspace(t0, t1, 101 if T0 != T1 else 2) for (t0, T0), (t1, T1) in zip(TABLE, TABLE[1:])]))
p = np.array([(x, temperature(x)) for x in t])                                               # short pieces, colored by the temperature
ax.add_collection(LineCollection(np.stack([p[:-1], p[1:]], axis=1), colors=color(0.5 * (p[1:, 1] + p[:-1, 1])), lw=2.2, capstyle="round",
                                 zorder=2))
p = np.array([(x, temperature(x)) for x in (TABLE[-1][0], 90)])
ax.plot(*p.T, color=color(p[-1, 1]), lw=2.2, ls=(0, (3.5, 1.6)), zorder=2)
ax.scatter(*zip(*TABLE), c=color([point[1] for point in TABLE]), edgecolor="k", linewidth=.5, s=75, zorder=3)
PLACES = ((9, 0, "left", "center"), (-4, -8, "right", "top"), (4, -8, "left", "top"), (-4, 8, "right", "bottom"))
for (x, y), (dx, dy, ha, va) in zip(TABLE, PLACES):                                          # the lines of the file
    ax.annotate(r"\texttt{%d %d}" % (x, y), (x, y), (dx, dy), textcoords="offset points", fontsize=15, ha=ha, va=va)
note(ax, 7, 540, "linear between\ntwo points", ha="left", va="center")
note(ax, 71.5, 418, "constant after\nthe last point", ha="center", va="bottom")
ax.set_xlim(-3, 83); ax.set_ylim(230, 800); ax.set_xticks(range(0, 81, 10)); ax.set_yticks(range(300, 801, 100))
ax.set_xlabel("Time (s)"); ax.set_ylabel("Temperature (K)"); grid(ax)
save_png(fig, "time_temperature.png", 1100)


# --------------------------------------------------------------------------------- rows of a log (output-files.md)
def log_dump_interval(steps, log_dump_steps):
    """ GetLogDumpInterval of McAbstract.cpp """
    if steps > 10 * log_dump_steps:
        return log_dump_steps
    interval, remaining = 1, (steps + 1) // 100
    while remaining > 0:
        interval, remaining = interval * 10, remaining // 10
    return min(interval, log_dump_steps)


steps = np.array([1] + [10 ** k + d for k in range(2, 10) for d in (-2, -1)] + [10 ** 10])   # the interval changes at 10^k - 1
curves = {lds: np.array([log_dump_interval(int(s), lds) for s in steps]) for lds in (10 ** 2, 10 ** 4, 10 ** 6)}
common = curves[10 ** 4] == curves[10 ** 6]
LINE = dict(where="post", lw=2.4, solid_joinstyle="miter", solid_capstyle="butt")
fig, (ax,) = canvas(1100, 4.24, (0.85, 0.71, 6.13, 3.38))
ax.step(steps[common], curves[10 ** 6][common], color=DARK, **LINE)                          # dark: at least two of the curves coincide
for (lds, y), c in zip(curves.items(), (BLUE, YELLOW, RED)):
    alone = np.all([y != other for key, other in curves.items() if key != lds], axis=0)
    k = np.argmax(alone) - 1
    ax.step(steps[k:], y[k:], color=c, **LINE)
    note(ax, 7e9, 1.35 * lds, r"\texttt{log\_dump\_steps %d}" % lds, ha="right", va="bottom")
for k, s in enumerate(("every step", "every 10th step", "every 100th step")):
    note(ax, 10 ** (k + 1.9), 1.35 * 10 ** k, s, ha="right", va="bottom")
log_x(ax, 12); log_y(ax, 10); ax.minorticks_off(); ax.set_xlim(1, 1e10); ax.set_ylim(0.5, 1e7)
ax.set_xlabel("Number of the step"); ax.set_ylabel("Steps between two rows"); grid(ax)
save_png(fig, "log_rows.png", 1100)


# --------------------------------------------------------------------------------- energies in the logs (output-files.md)
PANELS = ((0.97, 0.72, 4.00, 3.1), (6.22, 0.72, 4.00, 3.1))
cmc, kmc = (pd.read_csv(os.path.join(runs, m, m + "_log.txt"), sep=r"\s+") for m in ("cmc", "kmc"))
cmc, kmc = cmc[cmc["steps"] > 0], kmc[(kmc["time"] > 0) & (kmc["steps"] <= 10 ** 7)]         # no zero on a log axis
fig, axes = canvas(1600, 4.24, *PANELS)
for ax, letter, x, y, lw, xlabel, label in ((axes[0], "a", cmc["steps"], cmc["energy"], 1.0, r"\texttt{steps}", "CMC, 800 K"),
                                            (axes[1], "b", kmc["time"], kmc["energy"], 1.3, r"\texttt{time} (s)", "KMC, 450 K")):
    ax.plot(x, y, color=BLUE, lw=lw)
    log_x(ax, 9); ax.xaxis.minorticks_off(); ax.margins(y=0.1)
    ax.set_xlim(10 ** np.floor(np.log10(x.min())), 10 ** np.ceil(np.log10(x.max())))
    note(ax, 0.045, 0.06, label, 15, transform=ax.transAxes)
    ax.set_xlabel(xlabel); ax.set_ylabel(r"\texttt{energy} (eV)"); grid(ax); panel_label(ax, letter)
save_png(fig, "log_plot.png", 1600)


# --------------------------------------------------------------------------------- vacancies of an alloy and of the box (methods.md)
SITES, T = 10976, np.linspace(250, 900, 651)


def profile_energy(step):
    """ the field profile_energy of all sites in a file of Ansys """
    with gzip.open(os.path.join(runs, "ansys", "ansys", step + ".xyz.gz"), "rt") as f:
        f.readline()
        fields = re.search(r"Properties=(\S+)", f.readline()).group(1).split(":")
        column = sum(int(width) for name, width in zip(fields[:fields.index("profile_energy"):3], fields[2::3]))
        return np.array([float(line.split()[column]) for line in f])


def c_equilibrium(energy, temp):
    """ all vacancies in equilibrium: the concentration in pure Al times the mean of exp(-E/kT) over the sites """
    return 1.67 * np.exp(0.7 - 0.66 / (kB * temp)) * np.exp(-energy[:, None] / (kB * np.atleast_1d(temp))).mean(axis=0)


solution, clusters = profile_energy("0"), profile_energy("10000000")                         # first and last frame of the KMC run
c, low, c_quench = 1 / SITES, T[T <= 600], c_equilibrium(solution, 800.0)[0]
fig, (a, b) = canvas(1600, 4.24, *PANELS)
for ax, y in ((a, lambda v: v), (b, lambda v: c / v)):
    ax.plot(T, y(c_equilibrium(solution, T)), color=RED, lw=2.2, zorder=3, label="solid solution")
    ax.plot(low, y(c_equilibrium(clusters, low)), color=RED, lw=2.2, ls=(0, (3.5, 1.6)), zorder=3, label="with clusters")
    ax.plot([250, 800], [y(c_quench)] * 2, color=BLUE, lw=2.2, zorder=3)
    ax.plot([800], [y(c_quench)], "o", color=BLUE, zorder=4, **MARK)
    ax.axhline(y(c), color="k", lw=1.2, ls="--", zorder=2)
a.annotate("", xy=(285, 1e-8), xytext=(285, 0.75 * c_quench), zorder=2,
           arrowprops=dict(arrowstyle="-|>", color=BLUE, lw=1.6, shrinkA=0, shrinkB=0))
note(a, 525, 7e-4, "after a quench from 800 K", 15, color=BLUE, ha="center", va="bottom")
note(a, 300, 4.8e-5, "one in %d sites" % SITES, 15, va="top")
note(a, 300, 1.7e-6, "mobile part", 15, color=BLUE, va="center")
legend(a, loc="lower right", title="equilibrium", fontsize=14, title_fontsize=14, labelcolor=RED)
a.get_legend().get_title().set_color(RED)
note(b, 520, 3e4, "equilibrium", 15, color=RED, va="center")
note(b, 300, 0.16, "after a quench", 15, color=BLUE, va="top")
for ax, letter, label, lim in ((a, "a", "Vacancies per site", (1e-14, 1e-2)), (b, "b", r"$c/c_{\text{alloy}}$", (1e-2, 1e10))):
    log_y(ax, 8); ax.yaxis.set_minor_locator(LogLocator(base=10.0, subs=(1,), numticks=99))  # labels at every second power of ten
    ax.set_ylim(*lim); ax.set_xlim(250, 900); ax.set_xticks(range(300, 901, 100))
    ax.set_xlabel(r"$T$ (K)"); ax.set_ylabel(label); grid(ax); panel_label(ax, letter)
save_png(fig, "vacancy_concentration.png", 1600)
