# frame_fields.webp and frame_fields.mp4: the clusters of one frame of Ansys in a full turn, colored by the fields cluster_id and
# binding_energy of its file. Make them in this folder with, for the folder R of the Ansys run and any scratch folder T:
#   blender -b -P render_frame_turn.py -- T/pictures R/ansys/10000000.xyz.gz
#   python3 animate_frame_fields.py T/pictures T/frames R/ansys/10000000.xyz.gz <folder for the WebP> <folder for the MP4>
# Frames that exist already are kept. Needs LaTeX, cwebp and webpmux, ffmpeg, Pillow and SciPy.
# The frame is step 10 000 000 of the KMC run of docs/README.md, analyzed with smallest_cluster_criteria 10 and solvent_bond_criteria 3.
import io
import json
import os
import subprocess
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from style import *
import frame_fields_data as data
from PIL import Image
from scipy.ndimage import binary_dilation

FPS, COUNT = 24, 180                                       # of the MP4; the WebP has every second picture
pictures, frames, frame, webp, mp4 = sys.argv[1:6]
os.makedirs(frames, exist_ok=True)
members = data.load(frame)[0]
sizes = {a["cluster_id"]: a["cluster_size"] for a in members}
for c, n in sizes.items():                                 # the vacancy is a member, but it does not count
    assert n == sum(a["cluster_id"] == c and a["species"] != "X" for a in members)
assert all(abs(a["binding_energy"]) <= data.VMAX for a in members)

# ------------------------------------------------------------------ what does not change: labels and legends
W, H, DPI = 1920, 1080, 160
BOX = 960, 900                                             # size of a picture; the two are side by side at the top
ROWS = 908, 975, 1035                                      # labels of the panels, and the two rows of the legends
mpl.rcParams.update({"savefig.bbox": None, "font.size": 19})
TEXT = dict(fontsize=19, color=DARK, va="center_baseline")
fig = plt.figure(figsize=(W / DPI, H / DPI))
at = lambda x, y: (x / W, 1 - y / H)                       # pixels from the upper left corner
for left, letter, name in ((0, "a", "clusters"), (BOX[0], "b", "the same atoms")):
    fig.text(*at(left + 12, 10), rf"\textbf{{({letter})}}", ha="left", va="top", fontsize=21, color="k")
    fig.text(*at(left + BOX[0] / 2, ROWS[0]), name, ha="center", **TEXT)


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
LEFT, PITCH, SIZE = 372, 60, 50                            # (a) the clusters: where the keys end, and the balls
for row, key in zip(ROWS[1:], ("cluster_id", "cluster_size")):
    fig.text(*at(LEFT, row), rf"\texttt{{{key}}}".replace("_", r"\_"), ha="right", **TEXT)
for c, n in sorted(sizes.items()):                         # dark numbers on the two light colors
    x = LEFT + 20 + SIZE // 2 + PITCH * (c - 1)
    base.alpha_composite(ball("#" + data.CLUSTER[c - 1], SIZE), (x - SIZE // 2, ROWS[1] - SIZE // 2))
    fig.text(*at(x, ROWS[1]), rf"\textbf{{{c}}}", ha="center", va="center_baseline", fontsize=16, color="k" if c in (3, 4) else "w")
    fig.text(*at(x, ROWS[2]), str(n), ha="center", **TEXT)
x = 882                                                    # the vacancy, between the two legends
fig.add_artist(mpl.lines.Line2D([at(x, 0)[0]], [at(0, ROWS[1])[1]], marker="o", markersize=38 * 72 / DPI, mfc="w", mec=RED, mew=2.4))
fig.text(*at(x + 32, ROWS[1]), "vacancy", ha="left", **TEXT)
LEFT, LENGTH, THICK = 1522, 300, 30                        # (b) the color bar
fig.text(*at(LEFT - 22, ROWS[1]), r"\texttt{binding\_energy} (eV)", ha="right", **TEXT)
bar = fig.colorbar(mpl.cm.ScalarMappable(norm=Normalize(-data.VMAX, data.VMAX), cmap=BIND), orientation="horizontal",
                   cax=fig.add_axes([LEFT / W, 1 - (ROWS[1] + THICK / 2) / H, LENGTH / W, THICK / H]))
bar.set_ticks([-data.VMAX, 0, data.VMAX], labels=[f"${v:g}$" for v in (-data.VMAX, 0, data.VMAX)])
bar.set_ticks([v * data.VMAX / 3 for v in (-2, -1, 1, 2)], minor=True)
bar.ax.tick_params(pad=9)
rings = [fig.add_artist(mpl.patches.Circle((0, 0), 1, transform=None, fill=False, ec=RED, lw=2.0)) for _ in range(2)]

# ------------------------------------------------------------------ the frames
for k in range(COUNT):
    path = os.path.join(frames, f"{k:03d}.png")
    if os.path.exists(path):
        continue
    with open(os.path.join(pictures, f"id_{k:03d}.json")) as f:
        d = json.load(f)
    im = base.copy()
    for ring, left, name in zip(rings, (0, BOX[0]), ("id", "binding")):                # the same camera in both pictures
        im.alpha_composite(Image.open(os.path.join(pictures, f"{name}_{k:03d}.png")).convert("RGBA"), (left, 0))
        ring.set_center((left + d["vacancy"][0] * BOX[0], H - BOX[1] + d["vacancy"][1] * BOX[1]))
        ring.set_radius((d["rim"][0] - d["vacancy"][0]) * BOX[0])
    buffer = io.BytesIO()
    fig.savefig(buffer, format="png", dpi=DPI, transparent=True)
    im.alpha_composite(Image.open(buffer).convert("RGBA"))
    im.convert("RGB").save(path + ".part", format="png")   # an interrupted frame does not count
    os.replace(path + ".part", path)
    print(path, flush=True)

# ------------------------------------------------------------------ the WebP: every second picture, a turn in 7.5 s without a pause
# The first frame is complete. Every later frame draws all blocks of 8 x 8 pixels in which something changes during the turn, and
# nothing else: labels and legends stay as they are, and no trace is left where an atom has been. The pixels next to these blocks keep
# their colors, although they are transparent: the edges of the blocks would get their color otherwise. The encoder uses the same
# step for the whole picture (-segments 1) and a finer one for the color (-sns 100): with its defaults it saves on the clusters, and
# small pale blue atoms between red ones turn green.
def small(k):
    """frame k at the size of the WebP"""
    return np.array(Image.open(os.path.join(frames, f"{k:03d}.png")).resize((760, 428), Image.LANCZOS)).astype(int)


first = small(0)
h, w = first.shape[:2]
moving = np.zeros((h, w), bool)
for k in range(2, COUNT, 2):
    moving |= np.abs(small(k) - first).max(axis=2) > 4
tiles = lambda mask: np.pad(mask, ((0, -h % 8), (0, -w % 8))).reshape(-(-h // 8), 8, -(-w // 8), 8).any(axis=(1, 3))
blocks = tiles(moving)
assert not (blocks & tiles((first < 250).any(axis=2) & ~moving)).any(), "a label or a legend is too close to what moves"
draw, near = (np.kron(b, np.ones((8, 8), bool))[:h, :w] for b in (blocks, binary_dilation(blocks, np.ones((3, 3)))))
args = []
for k in range(0, COUNT, 2):
    new, name = small(k), os.path.join(frames, f"small_{k:03d}")
    part = np.dstack([np.where(near[..., None], new, 255), 255 * draw]) if k else np.dstack([new, np.full((h, w), 255)])
    Image.fromarray(part.astype("uint8")).save(name + ".png")
    subprocess.run(["cwebp", "-quiet", "-q", "37", "-m", "6", "-segments", "1", "-sns", "100", "-sharp_yuv", "-exact", name + ".png",
                    "-o", name + ".webp"], check=True)
    args += ["-frame", name + ".webp", f"+{round(1000 * (k + 2) / FPS) - round(1000 * k / FPS)}+0+0+0+b"]
webp, mp4 = os.path.join(webp, "frame_fields.webp"), os.path.join(mp4, "frame_fields.mp4")
subprocess.run(["webpmux", *args, "-loop", "0", "-o", webp], check=True)

# ------------------------------------------------------------------ the MP4: all pictures
# the colors of the frames, in the form that players expect at this size
colors = "scale=out_color_matrix=bt709:out_range=tv,setparams=color_primaries=bt709:color_trc=iec61966-2-1"
subprocess.run(["ffmpeg", "-loglevel", "error", "-y", "-framerate", str(FPS), "-i", os.path.join(frames, "%03d.png"), "-vf", colors,
                "-c:v", "libx264", "-threads", "4", "-pix_fmt", "yuv420p", "-crf", "18", "-movflags", "+faststart", mp4], check=True)
print(f"{COUNT} pictures, {COUNT / FPS:.1f} s, " + ", ".join(f"{f}: {os.path.getsize(f) / 1000:.0f} kB" for f in (webp, mp4)))
