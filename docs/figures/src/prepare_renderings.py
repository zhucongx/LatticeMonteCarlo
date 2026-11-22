# Puts the pictures of render_lattice.py, render_clusters.py, render_environment.py and render_frame.py on white and crops them; the
# pictures of one group get the same crop. Writes the positions of labels, relative to the cropped picture, for the TikZ files.
#   python3 prepare_renderings.py <folder of the renderings> <output folder>
import json
import os
import sys

from PIL import Image

src, out = sys.argv[1:3]
os.makedirs(out, exist_ok=True)
GROUPS = [["lattice_cell"], ["shell1", "shell2", "shell3"],
          ["cluster_site", "cluster_pair1", "cluster_pair2", "cluster_pair3", "cluster_triplet111", "cluster_triplet112",
           "cluster_triplet113", "cluster_triplet123"],
          ["env_all"], ["env_mmm", "env_mm2"], ["frame_solutes", "frame_clusters"], ["frame_binding"]]
LABELS = {"env_mmm": [(f"\\envnum{{{k}}}{{%.4f}}{{%.4f}}", str(k)) for k in range(1, 19)]                    # the numbers of the 18 sites
          + [(f"\\coordinate ({k}) at (%.4f,%.4f);", k) for k in "oxyz"],                               # the axes of the local frame
          "frame_solutes": [("\\coordinate (vacancy) at (%.4f,%.4f);", "vacancy")]}
LABELS["env_mm2"] = LABELS["env_mmm"]
MARGIN = 6
for names in GROUPS:
    if not all(os.path.exists(f"{src}/{n}.png") for n in names):
        print("no rendering of", *names)
        continue
    ims = [Image.open(f"{src}/{n}.png").convert("RGBA") for n in names]
    w, h = ims[0].size
    boxes = [im.getchannel("A").point(lambda a: 255 if a > 8 else 0).getbbox() for im in ims]
    x0, y0 = (max(0, min(b[k] for b in boxes) - MARGIN) for k in (0, 1))
    x1, y1 = min(w, max(b[2] for b in boxes) + MARGIN), min(h, max(b[3] for b in boxes) + MARGIN)
    for n, im in zip(names, ims):
        white = Image.new("RGBA", im.size, "white")
        white.alpha_composite(im)
        white.convert("RGB").crop((x0, y0, x1, y1)).save(f"{out}/{n}.png")
        if n in LABELS:
            with open(f"{src}/{n}.json") as f:
                at = json.load(f)
            with open(f"{out}/{n}_labels.tex", "w") as f:
                for text, key in LABELS[n]:
                    x, y = at[key]
                    f.write(text % ((x * w - x0) / (x1 - x0), 1 - ((1 - y) * h - y0) / (y1 - y0)) + "\n")
    print(*names, (x1 - x0, y1 - y0))
