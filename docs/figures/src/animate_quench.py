# quench.webp and quench.mp4: a KMC run that follows a time-temperature table, with the temperature of the table and the energy of the
# run against the time. For the folder R of the run, the file F of Ansys for its start and any scratch folder T:
#   blender -b -P render_quench.py -- T/pictures 1040 64 0 119 R
#   python3 animate_quench.py R F T/pictures T/frames quench.webp
#   python3 animate_quench.py R F T/pictures T/frames quench.mp4
# Pictures and frames that exist already are kept. Needs LaTeX, cwebp and webpmux, ffmpeg, pandas, Pillow and SciPy.
# The run is the KMC example of docs/README.md with the sample table of docs/input-files.md as time_temperature.txt, and with these
# lines in the parameter file in place of temperature, config_dump_steps and maximum_steps:
#   time_temperature_filename time_temperature.txt
#   config_dump_steps 10000
#   maximum_steps 90000000
# A configuration every 10000 steps puts one within 0.03 ms of the time of every picture. The time in the alloy is that of an alloy which
# keeps the vacancies that it has in equilibrium at the temperature at which the table begins, see alloy_time.py.
import glob
import io
import json
import os
import subprocess
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from style import *
import pandas as pd
from alloy_time import factor as alloy_factor
from matplotlib.collections import LineCollection
from PIL import Image
from scipy.ndimage import binary_dilation

FPS, START, END = 12, 0.5, 2.5                             # pictures per second, time in s of the first and of the last picture

run, frame, pictures, frames, output = sys.argv[1:6]
os.makedirs(frames, exist_ok=True)
log = pd.read_csv(os.path.join(run, "kmc_log.txt"), sep=r"\s+", usecols=["steps", "time", "temperature", "energy"]).set_index("steps")
t, T, E = 1e3 * log["time"].values, log["temperature"].values, log["energy"].values           # the time in ms
lines = []                                                 # the table as the program reads it: the lines that begin with two numbers
with open(os.path.join(run, "time_temperature.txt")) as f:
    for line in f:
        try:
            lines.append((float(line.split()[0]), float(line.split()[1]), " ".join(line.split()[:2])))
        except (IndexError, ValueError):
            pass
tt, TT, labels = zip(*sorted(lines))
tt = 1e3 * np.array(tt)
# TimeTemperatureInterpolator::GetTemperature: linear between two points, constant before the first and after the last one
if abs(T - np.interp(t, tt, TT)).max() > 1e-6:
    sys.exit("the temperatures of the log are not those of the table")
# The alloy keeps the vacancies that it has in equilibrium at the temperature at which the table begins, and the box holds one in all
# its sites: a time of the run times this factor is the time in the alloy.
factor = alloy_factor(frame, TT[0])
shown = []
for f in sorted(glob.glob(os.path.join(pictures, "[0-9][0-9][0-9].png"))):
    with open(f[:-4] + ".json") as g:
        shown.append((f, json.load(g)))

# ------------------------------------------------------------------ what does not change: the axes, the complete curves, the labels
W, H, DPI = 1920, 1080, 160
BOX = 1040                                                 # size of the pictures, which are at the left
LEFT, RIGHT = -0.25, 6.45                                  # the time axis in ms
mpl.rcParams.update({"savefig.bbox": None, "font.size": 16})
TEXT = dict(fontsize=19, color=DARK)
LINE = dict(lw=2.2, capstyle="round", zorder=2)
DASHES = dict(lw=2.2, ls=(0, (3.5, 1.6)), zorder=2)
fig = plt.figure(figsize=(W / DPI, H / DPI))
top, low = (fig.add_axes([0.593, bottom, 0.389, 0.3725]) for bottom in (0.5225, 0.109))
color = lambda temperature: TEMP(TNORM(temperature))
pale = lambda c: 1 - 0.3 * (1 - np.asarray(c))             # as with alpha 0.3 on white


def pieces(x):
    """short pieces of the curve of the table through the times x, and the color of each"""
    p = np.column_stack([x, np.interp(x, tt, TT)])
    return np.stack([p[:-1], p[1:]], axis=1), color(0.5 * (p[:-1, 1] + p[1:, 1]))


def note(x, y, text, **kw):
    """text on a white box, so that the grid does not run through it"""
    top.annotate(text, (x, y), textcoords="offset points", fontsize=17, bbox=dict(facecolor="white", edgecolor="none", pad=1.5), **kw)


knots = np.unique(np.concatenate([np.linspace(a, b, 101) for a, b in zip(tt, tt[1:])]))
segments, colors = pieces(knots)
top.add_collection(LineCollection(segments, colors=pale(colors), **LINE))
top.plot([tt[-1], RIGHT], [TT[-1]] * 2, color=pale(color(TT[-1])), **DASHES)
run_solid = top.add_collection(LineCollection([], **LINE))
run_dashed, = top.plot([], [], color=color(TT[-1]), **DASHES)
top.scatter(tt, TT, c=color(TT), edgecolor="k", linewidth=.5, s=75, zorder=3)
PLACES = ((9, 0, "left", "center"), (4, -9, "left", "top"), (4, -9, "left", "top"), (-4, 8, "right", "bottom"))
for x, y, label, (dx, dy, ha, va) in zip(tt, TT, labels, PLACES):                              # the lines of the file
    top.annotate(rf"\texttt{{{label}}}", (x, y), (dx, dy), textcoords="offset points", fontsize=17, ha=ha, va=va)
note(0.4, 600, "linear between\ntwo points", xytext=(14, 0), ha="left", va="center")
note(5.25, TT[-1], "constant after\nthe last point", xytext=(0, 10), ha="center", va="bottom")
low.plot(t, E, color=BLUE, lw=1.2, alpha=0.3)
past, = low.plot([], [], color=BLUE, lw=1.2)
value = low.text(0.97, 0.91, "", transform=low.transAxes, ha="right", va="top", **{**TEXT, "color": BLUE})
dots = [ax.plot([], [], "o", color=BLUE, zorder=5, clip_on=False, **{**MARK, "markersize": 11, "markeredgewidth": 0.8})[0] for ax in (top, low)]
top.set_ylabel(r"\texttt{temperature} (K)"); top.set_ylim(185, 870); top.set_yticks(range(300, 801, 100)); top.tick_params(labelbottom=False)
low.set_ylabel(r"\texttt{energy} (eV)"); low.set_ylim(-19.5, 3.5); low.set_xlabel(r"\texttt{time} (ms)")
for ax in (top, low):
    ax.set_xlim(LEFT, RIGHT); ax.set_xticks(range(7)); grid(ax)
alloy = top.secondary_xaxis("top", functions=(lambda x: x * factor, lambda x: x / factor))
alloy.set_xlabel("time in the alloy (ms)"); alloy.set_xticks(np.arange(0, factor * RIGHT, 0.5))
fig.align_ylabels()
caption = [fig.text(left / W, 34 / H, "", ha="left", va="center_baseline", **{**TEXT, "fontsize": size})
           for left, size in ((60, 19), (215, 19), (500, 16), (765, 19))]
ring = mpl.patches.Circle((0, 0), 1, transform=None, fill=False, ec=RED, lw=2.2)
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
for left, c, size, name in ((46, "#E08A2E", 42, "Mg"), (208, "#3F78B5", 42, "Zn"), (362, None, 38, "vacancy"),
                            (590, "#F3D0AB", 22, ""), (626, "#B2C9E1", 22, "in no cluster")):
    if c:
        base.alpha_composite(ball(c, size), (left - size // 2, 38 - size // 2))
    else:
        fig.add_artist(mpl.lines.Line2D([left / W], [1 - 38 / H], marker="o", markersize=size * 72 / DPI, mfc="w", mec=RED, mew=2.4))
    fig.text((left + size / 2 + 13) / W, 1 - 38 / H, name, ha="left", va="center_baseline", **TEXT)

# ------------------------------------------------------------------ the frames
for i, (picture, d) in enumerate(shown):
    path = os.path.join(frames, f"{i:03d}.png")
    if os.path.exists(path):
        continue
    row = log.index.get_loc(d["step"])
    segments, colors = pieces(np.append(knots[knots < min(t[row], tt[-1])], min(t[row], tt[-1])))     # the table up to this time
    run_solid.set_segments(segments); run_solid.set_color(colors if len(colors) else "k")
    run_dashed.set_data([tt[-1], max(t[row], tt[-1])], [TT[-1]] * 2)
    past.set_data(t[:row + 1], E[:row + 1])
    for dot, v in zip(dots, (T, E)):
        dot.set_data([t[row]], [v[row]])
    dots[0].set_markerfacecolor(color(T[row]))
    value.set_text(f"${round(E[row], 1) + 0:.1f}$ eV")
    for part, text in zip(caption, (f"{T[row]:.0f} K", f"$t$ = {factor * t[row]:.3f} ms", rf"\texttt{{time}} {t[row]:.2f} ms",
                                    f"step ${d['step']:,}$".replace(",", r"\,"))):
        part.set_text(text)
    x, y = d["vacancy"]
    ring.set_center((x * BOX, H - BOX + y * BOX))
    ring.set_radius((d["rim"][0] - x) * BOX)
    buffer = io.BytesIO()
    fig.savefig(buffer, format="png", dpi=DPI, transparent=True)
    im = base.copy()
    im.alpha_composite(Image.open(picture).convert("RGBA").resize((BOX, BOX), Image.LANCZOS))
    im.alpha_composite(Image.open(buffer).convert("RGBA"))
    im.convert("RGB").save(path + ".part", format="png")   # an interrupted frame does not count
    os.replace(path + ".part", path)
    print(path, d["step"], flush=True)

# ------------------------------------------------------------------ the file
times = [START] + [1 / FPS] * (len(shown) - 2) + [END]
if output.endswith(".webp"):
    # As in animate_clustering.py: every frame draws the blocks of 8 x 8 pixels in which something has changed, and nothing else. The
    # pixels next to these blocks keep their colors, although they are transparent: the edges of the blocks would get their color
    # otherwise.
    args = []
    for i, v in enumerate(times):
        new = np.array(Image.open(os.path.join(frames, f"{i:03d}.png")).resize((760, 428), Image.LANCZOS)).astype(int)
        h, w = new.shape[:2]
        changed = np.abs(new - last).max(axis=2) > 4 if i else np.ones((h, w), bool)
        blocks = np.pad(changed, ((0, -h % 8), (0, -w % 8))).reshape(-(-h // 8), 8, -(-w // 8), 8).any(axis=(1, 3))
        draw, near = (np.kron(b, np.ones((8, 8), bool))[:h, :w] for b in (blocks, binary_dilation(blocks, np.ones((3, 3)))))
        last = np.where(draw[..., None], new, last) if i else new
        small = os.path.join(frames, f"small_{i:03d}")
        Image.fromarray(np.dstack([np.where(near[..., None], new, 255), 255 * draw]).astype("uint8")).save(small + ".png")
        subprocess.run(["cwebp", "-quiet", "-q", "80", "-m", "6", "-sharp_yuv", "-exact", small + ".png", "-o", small + ".webp"], check=True)
        args += ["-frame", small + ".webp", f"+{round(1000 * v)}+0+0+0+b"]
    subprocess.run(["webpmux", *args, "-loop", "0", "-o", output], check=True)
else:
    # the colors of the frames, in the form that players expect at this size
    colors = "scale=out_color_matrix=bt709:out_range=tv,setparams=color_primaries=bt709:color_trc=iec61966-2-1"
    movie = subprocess.Popen(["ffmpeg", "-loglevel", "error", "-y", "-f", "rawvideo", "-pix_fmt", "rgb24", "-s", f"{W}x{H}", "-r", str(2 * FPS),
                              "-i", "-", "-vf", colors, "-c:v", "libx264", "-pix_fmt", "yuv420p", "-crf", "18", "-threads", "4",
                              "-movflags", "+faststart", output], stdin=subprocess.PIPE)
    for i, v in enumerate(times):
        movie.stdin.write(Image.open(os.path.join(frames, f"{i:03d}.png")).tobytes() * round(2 * FPS * v))
    movie.stdin.close()
    movie.wait()
print(f"{len(shown)} pictures, {sum(times):.1f} s, {output}: {os.path.getsize(output) / 1000:.0f} kB")
