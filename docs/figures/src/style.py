# Common settings of the data plots: the matplotlib style xzc_sans, which is not part of the repository, with the palette of the figures.
import os

import matplotlib
matplotlib.use("pdf")
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
from cycler import cycler
from matplotlib.colors import LinearSegmentedColormap, Normalize
from matplotlib.ticker import LogLocator, NullFormatter

sns.set_theme(style="ticks")
plt.style.use(["xzc_sans"])
BLUE, RED, YELLOW, GREEN, GRAY, DARK = "#3F78B5", "#C8423B", "#E0A21B", "#3E9B5C", "#7E858F", "#2B2B2B"
AL, MG, ZN = "#9BA1A9", "#DE8D3A", "#3F78B5"                       # species: Al gray, Mg orange, Zn blue
SPECIES = {"Al": AL, "Mg": MG, "Zn": ZN}
mpl.rcParams["axes.prop_cycle"] = cycler(color=[BLUE, RED, YELLOW, GREEN, GRAY, DARK])
TEMP = LinearSegmentedColormap.from_list("temperature", [BLUE, "#4FA08B", YELLOW, RED])      # 300 K to 800 K
TNORM = Normalize(vmin=300.0, vmax=800.0)
BIND = LinearSegmentedColormap.from_list("binding", [BLUE, "#F4F4F4", RED])                  # bound: blue, repelled: red
BARRIER = LinearSegmentedColormap.from_list("barrier", ["#FFF1CC", YELLOW, "#7A4A00"])
MAP = LinearSegmentedColormap.from_list("map", ["#FFF6E5", "#F2C879", YELLOW, RED, "#5E1F24"][::-1])
MARK = dict(markeredgecolor="k", markeredgewidth=.5, markersize=7)
kB = 8.617333262145e-5
OUT = "."


def set_size(w, h, ax):
    """ w, h: width, height of the axes in inches """
    p = ax.figure.subplotpars
    ax.figure.set_size_inches(float(w) / (p.right - p.left), float(h) / (p.top - p.bottom))


def save(fig, name):
    fig.savefig(os.path.join(OUT, name), bbox_inches="tight", dpi=400)
    plt.close(fig)
    print(name)


def log_x(ax, n):
    ax.set_xscale("log")
    ax.xaxis.set_major_locator(LogLocator(base=10.0, numticks=n))
    ax.xaxis.set_minor_formatter(NullFormatter())


def log_y(ax, n):
    ax.set_yscale("log")
    ax.yaxis.set_major_locator(LogLocator(base=10.0, numticks=n))
    ax.yaxis.set_minor_formatter(NullFormatter())


def grid(ax):
    ax.grid(visible=True, which="major", axis="both")


def legend(ax, **kw):
    ax.legend(**{**dict(fontsize=12, borderpad=0.25, labelspacing=0.25, handlelength=1.5, handletextpad=0.45, columnspacing=0.8), **kw})
