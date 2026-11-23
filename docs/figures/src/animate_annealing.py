# annealing.webp and annealing.mp4: a simulated annealing run, with the temperature and the energy of its log. For the folder R of the
# run and any scratch folder T:
#   blender -b -P render_annealing.py -- T/pictures 1040 64 R/lowest_energy.cfg.gz $(python3 animate_annealing.py list R)
#   python3 animate_annealing.py make R T/pictures T/frames <folder for annealing.webp> <folder for annealing.mp4>
# The first command makes a picture of every configuration that is listed, the second one puts the plots next to them. Pictures and
# frames that exist already are kept. The program writes lowest_energy.cfg.gz only; the numbered configurations, one every 1000 steps,
# are from a build in which SimulatedAnnealing::Dump writes them as well. Needs LaTeX, img2webp, ffmpeg, pandas, Pillow and SciPy.
import io
import os
import subprocess
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from style import *
import pandas as pd
from matplotlib.lines import Line2D
from matplotlib.patches import Patch, Rectangle
from PIL import Image
from scipy.interpolate import PchipInterpolator

FPS = 12
FIRST, END, LOWEST = 0.5, 0.6, 2.5                         # time in s of the start, the end of the run and the lowest energy
# frames and steps: slowly while the cluster forms, faster later
KNOTS = np.array(((0, 0), (10, 250_000), (44, 420_000), (62, 800_000), (78, 1_444_000), (92, 5_000_000)))
steps = [1000 * round(s / 1000) for s in PchipInterpolator(*KNOTS.T)(range(KNOTS[-1, 0] + 1)).tolist()]
names = [str(s) for s in steps] + ["lowest_energy"]

mode, run = sys.argv[1:3]
if mode == "list":
    print("\n".join(os.path.join(run, n + ".cfg.gz") for n in names))
    sys.exit()

pictures, frames, webp, mp4 = sys.argv[3:7]
os.makedirs(frames, exist_ok=True)
log = pd.read_csv(os.path.join(run, "sa_log.txt"), sep=r"\s+")
low = log["energy"].idxmin()
rows = [log.index[log["steps"] == s][0] for s in steps] + [low]
x, T = log["steps"] / 1e6, log["temperature"]
last = log["steps"].iloc[-1]                               # maximum_steps of the run
bare = T.iloc[0] * np.exp(-3 * log["steps"] / last)        # the factor exp(-3 / maximum_steps) of every step alone
# The log shows where the temperature was changed in another way: the ratio of two rows then differs from that of the bare schedule.
ratio = np.exp(np.diff(np.log(T / bare)))
down, up = np.where(np.abs(ratio - 0.99) < 1e-6)[0] + 1, np.where(np.abs(ratio - 1.10) < 1e-6)[0] + 1
assert np.sum(np.abs(ratio - 1) > 1e-6) == len(down) + len(up)
print(f"x0.99: {len(down)} times, rows {log['steps'][down[0]]} to {log['steps'][down[-1]]}; x1.10: {len(up)} times, rows "
      f"{list(log['steps'][up])}; lowest energy {log['energy'][low]:.6f} eV at step {log['steps'][low]}, {T[low]:.3f} K")

# ------------------------------------------------------------------ what does not change: the axes, the legends, and all of the run,
# pale. On top of it, in full color, comes what the run has reached in a frame.
W, H, DPI = 1920, 1080, 160
PALE = 0.3
mpl.rcParams.update({"savefig.bbox": None, "font.size": 16})
TEXT = dict(fontsize=19, color=DARK)
SMALL = 15                                                 # font size of the legends: 13 pixels in the WebP
TRIANGLE = {**MARK, "markersize": 7}
fig = plt.figure(figsize=(W / DPI, H / DPI))
ax_t, ax_e = (fig.add_axes([0.593, bottom, 0.389, 0.412]) for bottom in (0.562, 0.109))
ax_t.set_ylim(0, 870); ax_e.set_ylim(-9.2, 0.8)
reached = []                                               # what is drawn in full color up to the step of the frame
for ax, key, color, width in ((ax_t, "temperature", RED, 2.0), (ax_e, "energy", BLUE, 1.2), (ax_e, "lowest_energy", RED, 1.5)):
    ax.plot(x, log[key], color=color, lw=1.2, alpha=PALE)
    reached += ax.plot(x, log[key], color=color, lw=width)
# the schedule without the two other factors
ax_t.plot(x, bare, color="k", lw=1.2, ls="--", alpha=PALE)
reached += ax_t.plot(x, bare, color="k", lw=1.2, ls="--")
ax_t.text(1.33, 392, rf"${T.iloc[0]:.0f}\ \text{{K}}\,\exp(-3\,\texttt{{steps}}/\texttt{{maximum\_steps}})$", fontsize=SMALL, va="center")
# the steps in which the temperature was multiplied by 0.99, and those in which it was multiplied by 1.10
end_down = x[down[-1]]
shade_pale, shade = (ax_t.add_patch(Rectangle((0, 0), end_down, 1, transform=ax_t.get_xaxis_transform(), color=YELLOW, alpha=alpha, lw=0))
                     for alpha in (0.30 * PALE, 0.30))
ax_t.plot(x[up], T[up - 1] - 24, "^", color=RED, alpha=PALE, **TRIANGLE)
triangles = [ax_t.plot(x[u], T[u - 1] - 24, "^", color=RED, **TRIANGLE)[0] for u in up]
legend(ax_t, loc="upper right", fontsize=SMALL, handles=[
    Line2D([], [], color=RED, lw=2.0, label="run"),
    Patch(color=YELLOW, alpha=0.30, lw=0, label=rf"$\times\,0.99$ ({len(down)} times)"),
    Line2D([], [], color=RED, marker="^", ls="", label=rf"$\times\,1.10$ ({len(up)} times)", **TRIANGLE)])
legend(ax_e, loc="upper right", fontsize=SMALL, handles=[
    Line2D([], [], color=BLUE, lw=1.2, label=r"\texttt{energy}"), Line2D([], [], color=RED, lw=1.5, label=r"\texttt{lowest\_energy}")])
plots = []
for ax, key, color, unit in ((ax_t, "temperature", RED, "K"), (ax_e, "energy", BLUE, "eV")):
    dot, = ax.plot([], [], "o", color=color, zorder=5, clip_on=False, **{**MARK, "markersize": 10, "markeredgewidth": 0.8})
    text = ax.text(0.5, 0.91, "", transform=ax.transAxes, ha="right", va="top", **{**TEXT, "color": color})
    ax.set_xlim(0, 5); grid(ax)
    plots.append((log[key], dot, text, unit))
ax_t.set_ylabel(r"\texttt{temperature} (K)"); ax_e.set_ylabel("Energy (eV)")           # the lower plot has two columns of the log
ax_t.tick_params(labelbottom=False)
ax_e.set_xlabel(r"\texttt{steps} ($10^6$)")
fig.align_ylabels()
caption = fig.text(495 / W, 34 / H, "", ha="center", va="center", **TEXT)


def ball(color, size):
    """the shaded ball of figstyle.tex: a disk of radius 25, and the center of the shading to the upper left"""
    u = (np.arange(4 * size) + 0.5) / (4 * size) * 50 - 25
    X, Y = np.meshgrid(u, -u)
    c = np.array(mpl.colors.to_rgb(color))
    shades = [0.22 * c + 0.78, 0.62 * c + 0.38, 0.96 * c + 0.04, 0.92 * c, 0.74 * c, 0.55 * c]
    rgb = [np.interp(np.hypot(X + 9, Y - 9), (0, 6, 13, 20, 27, 50), [s[k] for s in shades]) for k in range(3)]
    im = np.dstack([*rgb, np.hypot(X, Y) <= 25]) * 255
    return Image.fromarray(im.round().astype("uint8")).resize((size, size), Image.LANCZOS)


base = Image.new("RGBA", (W, H), "white")
for left, color, name in ((46, "#E08A2E", "Mg"), (208, "#3F78B5", "Zn"), (362, None, "vacancy")):
    if color:
        base.alpha_composite(ball(color, 42), (left - 21, 17))
    else:
        fig.add_artist(mpl.lines.Line2D([left / W], [1 - 38 / H], marker="o", markersize=38 * 72 / DPI, mfc="w", mec=RED, mew=2.4))
    fig.text((left + 34) / W, 1 - 38 / H, name, ha="left", va="center_baseline", **TEXT)

# ------------------------------------------------------------------ the frames
for i, (name, row) in enumerate(zip(names, rows)):
    path = os.path.join(frames, f"{i:03d}.png")
    if os.path.exists(path):
        continue
    # The last frame is the whole figure: the run is over, and the point is where its lowest energy was found.
    now = x.iloc[-1] if name == "lowest_energy" else x[row]
    for artist in reached:                                 # the same line in every frame, cut at the step of the frame
        artist.set_clip_path(Rectangle((0, 0), now, 1, transform=artist.axes.get_xaxis_transform()))
    done = min(now, end_down)
    shade.set_width(done); shade_pale.set_x(done); shade_pale.set_width(end_down - done)
    for u, triangle in zip(up, triangles):
        triangle.set_visible(x[u] <= now)
    for v, dot, text, unit in plots:
        dot.set_data([x[row]], [v[row]])
        text.set_text(f"${round(v[row], 2 if unit == 'eV' else None) + 0:{'.2f' if unit == 'eV' else 'd'}}$ {unit}")
    s = f"${log['steps'][row]:,}$".replace(",", r"\,")
    caption.set_text(rf"\texttt{{lowest\_energy.cfg.gz}} (step {s})" if name == "lowest_energy" else
                     f"step {s}" + ", start of the run" * (row == 0) + ", end of the run" * (row == rows[-2]))
    buffer = io.BytesIO()
    fig.savefig(buffer, format="png", dpi=DPI, transparent=True)
    im = base.copy()
    im.alpha_composite(Image.open(os.path.join(pictures, name + ".png")).convert("RGBA"), (0, 0), (10, 10))
    im.alpha_composite(Image.open(buffer).convert("RGBA"))
    im.convert("RGB").save(path + ".part", format="png")   # an interrupted frame does not count
    os.replace(path + ".part", path)
    print(path, name, flush=True)

# ------------------------------------------------------------------ the two files
times = [FIRST] + [1 / FPS] * (len(names) - 3) + [END, LOWEST]
args = []
for i, t in enumerate(times):
    small = os.path.join(frames, f"small_{i:03d}.png")
    Image.open(os.path.join(frames, f"{i:03d}.png")).resize((760, 428), Image.LANCZOS).save(small)
    args += ["-d", str(round(1000 * t)), small]
webp, mp4 = os.path.join(webp, "annealing.webp"), os.path.join(mp4, "annealing.mp4")
# A key frame every 9 to 17 frames: between them, pixels that change little keep the errors of the frames before. Without -exact the
# colors of the pixels that stay the same are thrown away, and faint colored lines appear on the white next to them.
subprocess.run(["img2webp", "-loop", "0", "-sharp_yuv", "-kmin", "9", "-kmax", "17", "-lossy", "-exact", "-q", "80", "-m", "6", *args,
                "-o", webp], check=True)
# the colors of the frames, in the form that players expect at this size
colors = "scale=out_color_matrix=bt709:out_range=tv,setparams=color_primaries=bt709:color_trc=iec61966-2-1"
movie = subprocess.Popen(["ffmpeg", "-loglevel", "error", "-y", "-f", "rawvideo", "-pix_fmt", "rgb24", "-s", f"{W}x{H}", "-r", str(2 * FPS),
                          "-i", "-", "-vf", colors, "-c:v", "libx264", "-pix_fmt", "yuv420p", "-crf", "18", "-movflags", "+faststart", mp4],
                         stdin=subprocess.PIPE)
for i, t in enumerate(times):
    movie.stdin.write(Image.open(os.path.join(frames, f"{i:03d}.png")).tobytes() * round(2 * FPS * t))
movie.stdin.close()
movie.wait()
print(f"{len(names)} frames, {sum(times):.1f} s, " + ", ".join(f"{f}: {os.path.getsize(f) / 1000:.0f} kB" for f in (webp, mp4)))
