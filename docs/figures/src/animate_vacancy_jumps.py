# vacancy_jumps.webp and vacancy_jumps.mp4: some single steps of a KMC run around the vacancy, with the barrier and the time of each.
#   python3 animate_vacancy_jumps.py <run folder> <folder of the pictures of render_vacancy_jumps.py> <folder for the frames> <output folder>
# Frames that exist already are kept. Needs LaTeX, img2webp, ffmpeg and Pillow.
import io
import os
import subprocess
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from style import *
import vacancy_jumps_data as data
from matplotlib.collections import LineCollection
from PIL import Image

run, pictures, frames, out = sys.argv[1:5]
os.makedirs(frames, exist_ok=True)
FIRST, STEPS, SUB = data.FIRST, data.STEPS, data.SUB
LAST = STEPS * SUB                                         # number of the last picture
d = data.load(run)
step, time, barrier = d["step"], d["time"] * 1e9, d["barrier"]
element = [d["species"][i] for i in d["selected"][:STEPS]]  # what jumps in every step
COLOR = {**SPECIES, "Al": GRAY}                            # the gray of the marked Al atom in the pictures

# ------------------------------------------------------------------ text and panels
W, H, DPI = 1920, 1080, 160
mpl.rcParams["savefig.bbox"] = None                        # frames of exactly this size
fig = plt.figure(figsize=(W / DPI, H / DPI))
top = fig.add_axes([0.590, 0.562, 0.392, 0.412])
low = fig.add_axes([0.590, 0.109, 0.392, 0.412])
rise = [[(step[n], time[n]), (step[n + 1], time[n + 1])] for n in range(STEPS)]
low.add_collection(LineCollection(rise, colors=[COLOR[e] for e in element], lw=1.2, alpha=0.30))
past = low.add_collection(LineCollection([], lw=3.0, capstyle="round"))
points = {}
for e in ("Al", "Mg", "Zn"):
    k = [n for n in range(STEPS) if element[n] == e]
    top.plot(step[k], barrier[k], "o", color=COLOR[e], alpha=0.25, markeredgewidth=0, markersize=7)
    points[e] = top.plot([], [], "o", color=COLOR[e], **MARK)[0]
dots = [ax.plot([], [], "o", zorder=5, clip_on=False, **{**MARK, "markersize": 12, "markeredgewidth": 0.8})[0] for ax in (top, low)]
jump = top.text(0.97, 0.91, "", transform=top.transAxes, ha="right", va="top", fontsize=18, color=DARK)
now = low.text(0.03, 0.91, "", transform=low.transAxes, ha="left", va="top", fontsize=18, color=DARK)
top.set_ylabel(r"\texttt{Ea} (eV)"); top.set_ylim(0.27, 0.67); top.tick_params(labelbottom=False)
low.set_ylabel(r"\texttt{time} (ns)"); low.set_ylim(1750, 1905); low.set_xlabel(r"\texttt{steps}")
for ax in (top, low):
    ax.set_xlim(FIRST - 1, FIRST + STEPS + 1); grid(ax)
fig.align_ylabels([top, low])
caption = fig.text(508 / W, 33 / H, "", ha="center", va="center", fontsize=18)


def ball(color, size):
    """the shaded ball of figstyle.tex: a disk of radius 25 with the center of the shading to the upper left"""
    v = (np.arange(4 * size) + 0.5) / (4 * size) * 50 - 25
    X, Y = np.meshgrid(v, -v)
    c = np.array(mpl.colors.to_rgb(color))
    shades = [0.22 * c + 0.78, 0.62 * c + 0.38, 0.96 * c + 0.04, 0.92 * c, 0.74 * c, 0.55 * c]
    rgb = [np.interp(np.hypot(X + 9, Y - 9), (0, 6, 13, 20, 27, 50), [q[i] for q in shades]) for i in range(3)]
    im = np.dstack([*rgb, np.hypot(X, Y) <= 25]) * 255
    return Image.fromarray(im.round().astype("uint8")).resize((size, size), Image.LANCZOS)


base = Image.new("RGBA", (W, H), "white")
for x0, color, name in ((46, "#B8BCC2", "Al"), (170, "#E08A2E", "Mg"), (318, "#3F78B5", "Zn"), (458, None, "vacancy")):
    if color:
        base.alpha_composite(ball(color, 40), (x0 - 20, 18))
    else:
        fig.add_artist(mpl.lines.Line2D([x0 / W], [1 - 38 / H], marker="o", markersize=36 * 72 / DPI, mfc="w", mec=RED, mew=2.4))
    fig.text((x0 + 34) / W, 1 - 38 / H, name, ha="left", va="center_baseline", fontsize=18)


def layer(n):
    """text and panels during step n of the window; all values are those of row FIRST + n of the log"""
    for e, marks in points.items():
        k = [m for m in range(n) if element[m] == e]
        marks.set_data(step[k], barrier[k])
    past.set_segments(rise[:n]); past.set_color([COLOR[e] for e in element[:n]] or ["k"])
    caption.set_text(rf"step ${step[n]}$")
    now.set_text(rf"${time[n]:.2f}$ ns")
    dots[1].set_data([step[n]], [time[n]])
    color = COLOR[element[n]] if n < STEPS else DARK
    dots[0].set_data([step[n]], [barrier[n]]) if n < STEPS else dots[0].set_data([], [])       # none after the last jump
    jump.set_text(rf"{element[n]} jumps, ${barrier[n]:.3f}$ eV" if n < STEPS else "")
    for dot in dots:
        dot.set_markerfacecolor(color)
    with io.BytesIO() as f:
        fig.savefig(f, format="png", dpi=DPI, transparent=True)
        return Image.open(f).convert("RGBA")


# ------------------------------------------------------------------ the frames
missing, layers = 0, {}
for number in range(LAST + 1):
    n, k = divmod(number, SUB)
    path, picture = (os.path.join(folder, f"{number:03d}.png") for folder in (frames, pictures))
    if os.path.exists(path):
        continue
    if not os.path.exists(picture):
        missing += 1
        continue
    if n not in layers:
        layers = {n: layer(n)}
    im = base.copy()
    im.alpha_composite(Image.open(picture).convert("RGBA"), (0, 0), (22, 20))
    im.alpha_composite(layers[n])
    im.convert("RGB").save(path + ".part", format="png")
    os.replace(path + ".part", path)
if missing:
    sys.exit(f"{missing} pictures are missing")

# ------------------------------------------------------------------ the WebP, with every second picture
# It is nearly lossless, because lossy partial pictures leave tinted lines at their edges. A pixel that changes by at most 3 of 255 levels
# keeps its color, so that a picture holds only what moves; white has to be exactly white, or the faded trail would leave a trace.
args, total, shown = [], 0, None
for number in range(0, LAST + 1, 2):
    im = np.asarray(Image.open(os.path.join(frames, f"{number:03d}.png")).resize((760, 428), Image.LANCZOS)).astype(int)
    if shown is not None:
        keep = (abs(im - shown).max(2) <= 3) & (im < 255).any(2)
        im = np.where(keep[..., None], shown, im)
    shown = im
    small = os.path.join(frames, f"small_{number:03d}.png")
    Image.fromarray(im.astype("uint8")).save(small)
    ms = {0: 450, LAST: 2000}.get(number, 150 if number % SUB == 0 else 75)
    total += ms
    args += ["-d", str(ms), small]
webp, movie = (os.path.join(out, "vacancy_jumps." + e) for e in ("webp", "mp4"))
subprocess.run(["img2webp", "-loop", "0", "-kmax", "0", "-near_lossless", "60", "-q", "75", "-m", "4", *args, "-o", webp], check=True)

# ------------------------------------------------------------------ the MP4, with all pictures
colors = "scale=out_color_matrix=bt709:out_range=tv,setparams=color_primaries=bt709:color_trc=iec61966-2-1"
mp4 = subprocess.Popen(["ffmpeg", "-loglevel", "error", "-y", "-f", "rawvideo", "-pix_fmt", "rgb24", "-s", f"{W}x{H}", "-r", "30", "-i", "-",
                        "-vf", colors, "-c:v", "libx264", "-pix_fmt", "yuv420p", "-crf", "18", "-movflags", "+faststart", movie],
                       stdin=subprocess.PIPE)
for number in range(LAST + 1):
    content = Image.open(os.path.join(frames, f"{number:03d}.png")).convert("RGB").tobytes()
    for _ in range({0: 15, LAST: 60}.get(number, 4 if number % SUB == 0 else 1)):
        mp4.stdin.write(content)
mp4.stdin.close()
mp4.wait()
print(f"{LAST // 2 + 1} pictures, {total / 1000:.1f} s, {os.path.getsize(webp) // 1024} kB WebP, {os.path.getsize(movie) // 1024} kB MP4")
