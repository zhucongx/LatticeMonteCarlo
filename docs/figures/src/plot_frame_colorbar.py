# The color bar of frame_fields.png, with the range and the colors of render_frame.py.   python3 plot_frame_colorbar.py <output folder>
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import style
from style import *

style.OUT = sys.argv[1]
fig, ax = plt.subplots(figsize=(3.0, 0.32))
cb = fig.colorbar(mpl.cm.ScalarMappable(norm=Normalize(-0.15, 0.15), cmap=BIND), cax=ax, orientation="horizontal")
cb.set_ticks([-0.15, 0, 0.15], labels=["$-0.15$", "$0$", "$0.15$"])
cb.set_ticks([-0.1, -0.05, 0.05, 0.1], minor=True)
save(fig, "plot_frame_colorbar.pdf")
