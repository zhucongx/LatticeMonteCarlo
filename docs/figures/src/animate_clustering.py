# clustering.webp and clustering.mp4: clusters of solute atoms form in a KMC run. For the folders R of the run and any scratch folder T:
#   blender -b -P render_clustering.py -- T/rest 1040 64 24 24 0 99 R ...
#   blender -b -P render_clustering.py -- T/turn 1040 64 54 24 0 99 R ...
#   python3 animate_clustering.py T/rest T/rest_frames clustering.webp F R ...
#   python3 animate_clustering.py T/turn T/turn_frames clustering.mp4 F R ...
# F is the file of Ansys for the start of the run. The time shown is that of an alloy which keeps the vacancies that this solid solution
# has in equilibrium at the temperature of the run: the time of the run times the factor of alloy_time.py.
# The camera turns in the MP4 file only: a moving camera makes the WebP file too large. Pictures and frames that exist already are kept.
# Needs LaTeX, cwebp and webpmux, ffmpeg, pandas, Pillow and SciPy.
# The run is the KMC example of docs/README.md with maximum_steps 130000000 and config_dump_steps 200000. Its early configurations are
# wanted more often, so the run is made twice more, in folders of their own: 1200000 steps with config_dump_steps 2000 and
# log_dump_steps 1000, and 6000000 steps with config_dump_steps 20000. The seed is the same, so the trajectory is the same.
import glob
import io
import json
import os
import subprocess
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from style import *
import pandas as pd
from alloy_time import factor
from PIL import Image
from scipy.ndimage import binary_dilation

FPS, START, END = 12, 0.6, 2.0                             # pictures per second, time in s of the first and of the last picture

pictures, frames, output, frame = sys.argv[1:5]
os.makedirs(frames, exist_ok=True)
log = pd.concat([pd.read_csv(os.path.join(f, "kmc_log.txt"), sep=r"\s+", usecols=["steps", "time", "temperature", "energy"])
                 for f in sys.argv[5:]])
log = log.drop_duplicates("steps").sort_values("steps").set_index("steps")
shown = []
for f in sorted(glob.glob(os.path.join(pictures, "[0-9][0-9][0-9].png"))):
    with open(f[:-4] + ".json") as g:
        shown.append((f, json.load(g)))
log = log.loc[:shown[-1][1]["step"]]
KELVIN = log["temperature"].values[0]                      # the run is at one temperature
alloy = log["time"].values * factor(frame, KELVIN)         # the time in the alloy
t, E = np.maximum(alloy, 1e-3), log["energy"].values       # the axis begins at 1e-3 s

# ------------------------------------------------------------------ what does not change: the axes, the complete curve, the legend
W, H, DPI = 1920, 1080, 160
BOX = 1040                                                 # size of the pictures, which are at the left
mpl.rcParams.update({"savefig.bbox": None, "font.size": 16})
TEXT = dict(fontsize=19, color=DARK)
fig = plt.figure(figsize=(W / DPI, H / DPI))
ax = fig.add_axes([0.609, 0.19, 0.362, 0.65])
ax.plot(t, E, color=BLUE, lw=1.2, alpha=0.3)
past, = ax.plot([], [], color=BLUE, lw=1.2)
dot, = ax.plot([], [], "o", color=BLUE, zorder=5, clip_on=False, **{**MARK, "markersize": 10, "markeredgewidth": 0.8})
value = ax.text(0.895, 0.86, "", transform=ax.transAxes, ha="right", va="center", **{**TEXT, "color": BLUE})   # between grid lines
log_x(ax, 6); ax.set_xlim(1e-3, 2e0); ax.set_ylim(-30, 2); grid(ax)
ax.set_xlabel("time in the alloy (s)"); ax.set_ylabel(r"\texttt{energy} (eV)")
caption = [fig.text(left / W, 34 / H, "", ha="left", va="center_baseline", **TEXT) for left in (140, 345, 620)]
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
for left, color, size, name in ((46, "#E08A2E", 42, "Mg"), (208, "#3F78B5", 42, "Zn"), (362, None, 38, "vacancy"),
                                (590, "#F3D0AB", 22, ""), (626, "#B2C9E1", 22, "in no cluster")):
    if color:
        base.alpha_composite(ball(color, size), (left - size // 2, 38 - size // 2))
    else:
        fig.add_artist(mpl.lines.Line2D([left / W], [1 - 38 / H], marker="o", markersize=size * 72 / DPI, mfc="w", mec=RED, mew=2.4))
    fig.text((left + size / 2 + 13) / W, 1 - 38 / H, name, ha="left", va="center_baseline", **TEXT)


def time_text(v):
    """time with three digits"""
    v, unit = (v, "s") if v >= 0.9995 else (v * 1e3, "ms") if v >= 0.9995e-3 else (v * 1e6, r"\textmu s")
    return f"{v:.{2 if v < 9.995 else 1 if v < 99.95 else 0}f} {unit}"


# ------------------------------------------------------------------ the frames
for i, (picture, d) in enumerate(shown):
    path = os.path.join(frames, f"{i:03d}.png")
    if os.path.exists(path):
        continue
    row = log.index.get_loc(d["step"])
    past.set_data(t[:row + 1], E[:row + 1])
    dot.set_data([t[row]], [E[row]])
    value.set_text(f"${round(E[row], 1) + 0:.1f}$ eV")
    for part, text in zip(caption, (f"{KELVIN:.0f} K", "$t$ = " + time_text(alloy[row]), f"step ${d['step']:,}$".replace(",", r"\,"))):
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
    # Every frame draws the blocks of 8 x 8 pixels in which something has changed, and nothing else. The pixels next to these blocks keep
    # their colors, although they are transparent: the edges of the blocks would get their color otherwise. (img2webp decides for single
    # pixels. That makes a larger file, and what does not change flickers with the key frames.)
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
                              "-i", "-", "-vf", colors, "-c:v", "libx264", "-pix_fmt", "yuv420p", "-crf", "18", "-movflags", "+faststart",
                              output], stdin=subprocess.PIPE)
    for i, v in enumerate(times):
        movie.stdin.write(Image.open(os.path.join(frames, f"{i:03d}.png")).tobytes() * round(2 * FPS * v))
    movie.stdin.close()
    movie.wait()
print(f"{len(shown)} pictures, {sum(times):.1f} s, {output}: {os.path.getsize(output) / 1000:.0f} kB")
