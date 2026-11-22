# The two plots of annealing.png, from the log of a simulated annealing run.
#   python3 plot_annealing.py <folder with sa_log.txt> <output folder>
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import style
from style import *
import pandas as pd

data, style.OUT = sys.argv[1], sys.argv[2]
log = pd.read_csv(os.path.join(data, "sa_log.txt"), sep=r"\s+")
last = log["steps"].iloc[-1]                              # maximum_steps of the run
x, T = log["steps"] / 1e6, log["temperature"]
bare = T.iloc[0] * np.exp(-3 * log["steps"] / last)       # the factor exp(-3 / maximum_steps) of every step alone

# The log shows where the temperature was changed in another way: the ratio of two rows then differs from that of the bare schedule.
ratio = np.exp(np.diff(np.log(T / bare)))
down, up = np.where(np.abs(ratio - 0.99) < 1e-6)[0] + 1, np.where(np.abs(ratio - 1.10) < 1e-6)[0] + 1
other = np.sum(np.abs(ratio - 1) > 1e-6) - len(down) - len(up)
low = log["energy"].idxmin()
print(f"x0.99: {len(down)} times, rows {log['steps'][down[0]]} to {log['steps'][down[-1]]}; x1.10: {len(up)} times, rows "
      f"{list(log['steps'][up])}; other changes: {other}; lowest energy {log['energy'][low]:.6f} eV at step {log['steps'][low]}, {T[low]:.3f} K")


def finish(ax, name, values, offsets):
    """axes, and the two states that are drawn below the plots"""
    for i, label, shift in zip((0, low), "cd", offsets):
        ax.plot(x[i], values[i], "o", color="w", zorder=5, clip_on=False, **{**MARK, "markeredgewidth": 1.3})
        ax.annotate(rf"\textbf{{({label})}}", (x[i], values[i]), shift, textcoords="offset points", ha="center", va="center", fontsize=14)
    ax.set_xlabel(r"\texttt{steps} ($10^6$)"); ax.set_xlim(0, 5); grid(ax); set_size(4.3, 2.9, ax); save(fig, name)


fig, ax = plt.subplots()
ax.plot(x, T, color=RED, lw=1.8, label="run")
ax.axvspan(0, x[down[-1]], color=YELLOW, alpha=0.30, lw=0, label=rf"$\times\,0.99$ ({len(down)} times)")
ax.plot(x[up], T[up - 1] - 22, "^", color=RED, label=rf"$\times\,1.10$ ({len(up)} times)", **{**MARK, "markersize": 6})
ax.plot(x, bare, color="k", lw=1.2, ls="--")
ax.text(1.25, 430, rf"${T.iloc[0]:.0f}\ \text{{K}}\,\exp(-3\,\texttt{{steps}}/\texttt{{maximum\_steps}})$", fontsize=14, va="center")
ax.set_ylabel(r"\texttt{temperature} (K)"); ax.set_ylim(0, 900)
legend(ax, loc="upper right", fontsize=14)
finish(ax, "plot_annealing_temperature.pdf", T, ((16, 6), (10, 13)))

fig, ax = plt.subplots()
ax.plot(x, log["energy"], color=BLUE, lw=1, label=r"\texttt{energy}")
ax.plot(x, log["lowest_energy"], color=RED, lw=1.4, label=r"\texttt{lowest\_energy}")
ax.set_ylabel(r"Energy (eV)"); ax.set_ylim(-9.3, 1.0)
legend(ax, loc="upper right", fontsize=14)
finish(ax, "plot_annealing_energy.pdf", log["energy"], ((16, 8), (10, 13)))
