# Plot of the acceptance probability, panel (b) of metropolis.png.   python3 plot_metropolis.py <output folder>
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import style
from style import *

style.OUT = sys.argv[1]
x = np.linspace(-3, 5.2, 821)
fig, ax = plt.subplots()
ax.plot(x, np.minimum(1, np.exp(-x)), color=RED, lw=2.2)
ax.axvline(0, color="k", lw=1, ls="--")
ax.text(-1.5, 0.915, "always", ha="center", va="baseline", fontsize=14)
ax.text(-1.5, 0.83, "accepted", ha="center", va="baseline", fontsize=14)
ax.text(1.15, 0.5, r"$\exp(-\Delta E/k_{\text{B}}T)$", fontsize=14)
ax.set_xlim(-3, 5.2)
ax.set_ylim(0, 1.15)
ax.set_xlabel(r"$\Delta E/k_{\text{B}}T$")
ax.set_ylabel("Acceptance probability")
grid(ax)
set_size(4.8, 2.7, ax)              # fig_metropolis.tex uses the height to align the two panels
save(fig, "plot_metropolis.pdf")
