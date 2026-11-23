# trapped_vacancy.webp and trapped_vacancy.mp4: the vacancy of a KMC run among the clusters. It makes most of its steps at the clusters,
# in about half of the time, and the few steps in the matrix take the other half. All times are times in the alloy, see alloy_time.py.
# For the folder R of the run, the file F of Ansys for the start of the KMC example, and any scratch folder T:
#   python3 trapped_vacancy_data.py R F T/run.npz
#   blender -b -P render_trapped_vacancy.py -- T/run.npz T/pictures 1040 64
#   python3 animate_trapped_vacancy.py T/run.npz T/pictures T/frames <output folder>
# Frames that exist already are kept. Needs LaTeX, cwebp and webpmux, ffmpeg, Pillow and SciPy.
# The run continues the KMC example of docs/README.md: its last configuration, 10000000.cfg.gz, is the start.cfg.gz of a run with the same
# parameter file, but maximum_steps 3000000 and log_dump_steps 1. The folder R holds start.cfg.gz and kmc_log.txt of this run.
import io
import json
import os
import subprocess
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from style import *
import trapped_vacancy_data as data
from matplotlib.collections import LineCollection
from PIL import Image
from scipy.ndimage import binary_dilation, distance_transform_edt
from trapped_vacancy_data import MATRIX, NAMED, OTHER

FPS, START, END = 12, 0.6, 2.0                             # pictures per second, time in s of the first and of the last picture
CLUSTER = ["#8A63A9", "#3E9B5C", "#D487AC", "#E6C343"]     # the colors of cluster_id 1 to 4 in the figure frame_fields, see frame_fields_data.py
INK = ["white", "white", "black", "black"]                 # and of their numbers
COLOR = {MATRIX: DARK, OTHER: GRAY, **{k + 1: c for k, c in enumerate(CLUSTER)}}
ROW = {MATRIX: 0, OTHER: 1, **{k: NAMED + 2 - k for k in range(1, NAMED + 1)}}      # rows of the upper panel, from the bottom

source, pictures, frames, out = sys.argv[1:5]
os.makedirs(frames, exist_ok=True)
d = data.load(source)
where, time, shown = d["where"], d["alloy"] * 1e3, d["shown"]          # the time in the alloy, in ms
rows = len(time)
labels = []                                                # what render_trapped_vacancy.py has written for every picture
for i in range(len(shown)):
    with open(os.path.join(pictures, f"{i:03d}.json")) as f:
        labels.append(json.load(f))
    assert labels[i]["step"] == shown[i]
# a visit: the rows first ... last - 1, during which the vacancy stays at one cluster or in the matrix
first = np.r_[0, np.flatnonzero(np.diff(where)) + 1]
last = np.r_[first[1:], rows - 1]
steps, spent, at_clusters = last - first, time[last] - time[first], where[first] != MATRIX
print(f"{rows - 1} steps, {d['time'][-1] * 1e6:.3f} us of the run, {time[-1]:.3f} ms in the alloy, factor {float(d['factor']):.1f}")
print(f"at clusters: {steps[at_clusters].sum()} steps ({steps[at_clusters].sum() / (rows - 1):.2%}) in {spent[at_clusters].sum():.3f} ms "
      f"({spent[at_clusters].sum() / time[-1]:.2%} of the time in the alloy)")

# ------------------------------------------------------------------ what does not change: the axes, the pale curves, the legend
W, H, DPI = 1920, 1080, 160
BOX = 1040                                                 # size of the pictures, which are at the left
mpl.rcParams.update({"savefig.bbox": None, "font.size": 16})
TEXT = dict(fontsize=19, color=DARK)
fig = plt.figure(figsize=(W / DPI, H / DPI))
top = fig.add_axes([0.590, 0.562, 0.392, 0.412])
low = fig.add_axes([0.590, 0.109, 0.392, 0.412])
bars = [[(time[a], ROW[where[a]]), (time[b], ROW[where[a]])] for a, b in zip(first, last)]
# the curve of the lower panel: a point at both ends of every visit, and every 2000 steps
knots = np.unique(np.r_[first, np.arange(0, rows, 2000), rows - 1])
curve = [[(time[a], a / 1e6), (time[b], b / 1e6)] for a, b in zip(knots[:-1], knots[1:])]
top.add_collection(LineCollection(bars, colors=[COLOR[s] for s in where[first]], lw=11, alpha=0.25, capstyle="butt"))
low.add_collection(LineCollection(curve, colors=[COLOR[s] for s in where[knots[:-1]]], lw=1.2, alpha=0.30))
past = [top.add_collection(LineCollection([], lw=11, capstyle="butt")), low.add_collection(LineCollection([], lw=3.0, capstyle="round"))]
dots = [ax.plot([], [], "o", zorder=5, clip_on=False, **{**MARK, "markersize": 12, "markeredgewidth": 0.8})[0] for ax in (top, low)]
top.set_ylim(-0.6, NAMED + 1.6); top.tick_params(labelbottom=False); top.tick_params(axis="y", which="minor", left=False)
top.set_yticks([ROW[MATRIX], ROW[OTHER], *range(2, NAMED + 2)], ["matrix", "other", *[""] * NAMED])
low.set_ylabel(r"\texttt{steps} ($10^6$)"); low.set_ylim(-0.12, 3.12); low.set_xlabel("time in the alloy (ms)")
for ax in (top, low):
    ax.set_xlim(-0.026 * time[-1], 1.054 * time[-1]); ax.set_xticks(np.arange(0, 1.054 * time[-1], 5)); grid(ax)
caption = [fig.text(left / W, 34 / H, "", ha="left", va="center_baseline", **TEXT) for left in (140, 345, 620)]
ring = mpl.patches.Circle((0, 0), 1, transform=None, fc="w", ec=RED, lw=2.4)          # the vacancy, on top of the picture as in the legend
fig.add_artist(ring)


def ball(color, size):
    """the shaded ball of figstyle.tex: a disk of radius 25, and the center of the shading to the upper left"""
    u = (np.arange(4 * size) + 0.5) / (4 * size) * 50 - 25
    x, y = np.meshgrid(u, -u)
    c = np.array(mpl.colors.to_rgb(color))
    shades = [0.22 * c + 0.78, 0.62 * c + 0.38, 0.96 * c + 0.04, 0.92 * c, 0.74 * c, 0.55 * c]
    rgb = [np.interp(np.hypot(x + 9, y - 9), (0, 6, 13, 20, 27, 50), [s[k] for s in shades]) for k in range(3)]
    im = np.dstack([*rgb, np.hypot(x, y) <= 25]) * 255
    return Image.fromarray(im.round().astype("uint8")).resize((size, size), Image.LANCZOS)


base = Image.new("RGBA", (W, H), "white")
for left, color, size, name in ((46, "#E08A2E", 42, "Mg"), (208, "#3F78B5", 42, "Zn"), (362, None, 38, "vacancy"),
                                (590, "#F3D0AB", 22, ""), (626, "#B2C9E1", 22, "in no cluster")):
    if color:
        base.alpha_composite(ball(color, size), (left - size // 2, 38 - size // 2))
    else:
        fig.add_artist(mpl.lines.Line2D([left / W], [1 - 38 / H], marker="o", markersize=size * 72 / DPI, mfc="w", mec=RED, mew=2.4))
    fig.text((left + size / 2 + 13) / W, 1 - 38 / H, name, ha="left", va="center_baseline", **TEXT)
# The numbers of the clusters: a ball in the color of the cluster, at its row and next to the cluster in the box. The place in the box is
# the same in all pictures: the one nearest to the middle of the cluster that is CLEAR pixels away from every large ball of every picture.
CLEAR = 52
taken = np.zeros((BOX, BOX), bool)
for label in labels:
    x, y = (np.array(label["large"]) * BOX).astype(int).T
    taken[BOX - 1 - y, x] = True
free = distance_transform_edt(~taken) >= CLEAR
free[:80] = free[-80:] = free[:, :40] = free[:, -40:] = False
ys, xs = np.nonzero(free)
marks = Image.new("RGBA", (W, H), (0, 0, 0, 0))
for k in range(1, NAMED + 1):
    x, y = fig.transFigure.inverted().transform(top.transData.transform((top.get_xlim()[0], ROW[k]))) * (W, H)    # the left end of its row
    middle = np.median([label[f"cluster {k}"] for label in labels], axis=0) * BOX
    nearest = np.hypot(xs - middle[0], ys - (BOX - middle[1])).argmin()
    for x, y in ((x - 40, H - y), (xs[nearest], ys[nearest])):
        marks.alpha_composite(ball(CLUSTER[k - 1], 46), (round(x) - 23, round(y) - 23))
        fig.text(x / W, 1 - y / H, rf"\textbf{{{k}}}", ha="center", va="center_baseline", fontsize=17, color=INK[k - 1])

# ------------------------------------------------------------------ the frames
for i, (row, label) in enumerate(zip(shown, labels)):
    path = os.path.join(frames, f"{i:03d}.png")
    if os.path.exists(path):
        continue
    seen = first <= row                                    # the visits and the pieces of the curve so far, the last one up to this row
    past[0].set_segments([[(time[a], ROW[where[a]]), (time[min(b, row)], ROW[where[a]])] for a, b in zip(first[seen], last[seen])])
    past[0].set_color([COLOR[s] for s in where[first[seen]]])
    seen = knots[:-1] < row
    past[1].set_segments([[(time[a], a / 1e6), (time[min(b, row)], min(b, row) / 1e6)] for a, b in zip(knots[:-1][seen], knots[1:][seen])])
    past[1].set_color([COLOR[s] for s in where[knots[:-1][seen]]] or ["k"])
    for dot, y in zip(dots, (ROW[where[row]], row / 1e6)):
        dot.set_data([time[row]], [y])
        dot.set_markerfacecolor(COLOR[where[row]])
    for part, text in zip(caption, (f"{d['temperature']:.0f} K", rf"$t$ = ${time[row]:.2f}$ ms", f"step ${row:,}$".replace(",", r"\,"))):
        part.set_text(text)
    x, y = label["vacancy"]
    ring.set_center((x * BOX, H - BOX + y * BOX))
    ring.set_radius((label["rim"][0] - x) * BOX)
    buffer = io.BytesIO()
    fig.savefig(buffer, format="png", dpi=DPI, transparent=True)
    im = base.copy()
    im.alpha_composite(Image.open(os.path.join(pictures, f"{i:03d}.png")).convert("RGBA").resize((BOX, BOX), Image.LANCZOS))
    im.alpha_composite(marks)
    im.alpha_composite(Image.open(buffer).convert("RGBA"))
    im.convert("RGB").save(path + ".part", format="png")   # an interrupted frame does not count
    os.replace(path + ".part", path)
    print(path, row, flush=True)

# ------------------------------------------------------------------ the files
times = [START] + [1 / FPS] * (len(shown) - 2) + [END]
webp, movie = (os.path.join(out, "trapped_vacancy." + e) for e in ("webp", "mp4"))
# The WebP file. Every frame draws the blocks of 8 x 8 pixels in which something has changed, and nothing else. The pixels next to these
# blocks keep their colors, although they are transparent: the edges of the blocks would get their color otherwise.
args = []
for i, v in enumerate(times):
    new = np.array(Image.open(os.path.join(frames, f"{i:03d}.png")).resize((760, 428), Image.LANCZOS)).astype(int)
    h, w = new.shape[:2]
    changed = np.abs(new - drawn).max(axis=2) > 4 if i else np.ones((h, w), bool)
    blocks = np.pad(changed, ((0, -h % 8), (0, -w % 8))).reshape(-(-h // 8), 8, -(-w // 8), 8).any(axis=(1, 3))
    draw, near = (np.kron(b, np.ones((8, 8), bool))[:h, :w] for b in (blocks, binary_dilation(blocks, np.ones((3, 3)))))
    drawn = np.where(draw[..., None], new, drawn) if i else new
    small = os.path.join(frames, f"small_{i:03d}")
    Image.fromarray(np.dstack([np.where(near[..., None], new, 255), 255 * draw]).astype("uint8")).save(small + ".png")
    subprocess.run(["cwebp", "-quiet", "-q", "80", "-m", "6", "-sharp_yuv", "-exact", small + ".png", "-o", small + ".webp"], check=True)
    args += ["-frame", small + ".webp", f"+{round(1000 * v)}+0+0+0+b"]
subprocess.run(["webpmux", *args, "-loop", "0", "-o", webp], check=True)
# The MP4 file, with the colors of the frames in the form that players expect at this size
colors = "scale=out_color_matrix=bt709:out_range=tv,setparams=color_primaries=bt709:color_trc=iec61966-2-1"
mp4 = subprocess.Popen(["ffmpeg", "-loglevel", "error", "-y", "-f", "rawvideo", "-pix_fmt", "rgb24", "-s", f"{W}x{H}", "-r", str(2 * FPS),
                        "-i", "-", "-vf", colors, "-c:v", "libx264", "-pix_fmt", "yuv420p", "-crf", "18", "-movflags", "+faststart", movie],
                       stdin=subprocess.PIPE)
for i, v in enumerate(times):
    content = Image.open(os.path.join(frames, f"{i:03d}.png")).tobytes()
    for _ in range(round(2 * FPS * v)):
        mp4.stdin.write(content)
mp4.stdin.close()
mp4.wait()
print(f"{len(shown)} pictures, {sum(times):.1f} s, {os.path.getsize(webp) / 1000:.0f} kB WebP, {os.path.getsize(movie) / 1000:.0f} kB MP4")
