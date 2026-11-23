# lattice.webp and lattice.mp4: the conventional fcc cell, in which an atom jumps into the vacancy and back, and the three neighbor shells
# of a site, in one turn of 7.5 s. The labels stay; the line that marks an edge of the cell turns with the cell.
# Make them in this folder with, for any scratch folder T:
#   blender -b -P render_lattice_turn.py -- T/pictures
#   python3 animate_lattice.py T/pictures T/frames <folder for lattice.webp> <folder for lattice.mp4>
# With numbers of frames after these arguments only those frames are made, for a look. Needs what turn_film.py needs.
import json
import os
import sys

from PIL import Image

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import turn_film as film

pictures, frames, webp, mp4 = sys.argv[1:5]
os.makedirs(frames, exist_ok=True)
FPS, COUNT = 24, 180                                       # of the MP4; the WebP has every second picture
W, H, CM = 1920, 720, 118                                  # size of a frame, and of one centimeter, in pixels
TOP, BOTTOM = 64, 590                                      # the pictures lie between the letters of the panels and the labels
PANELS = {"cell": (227, 358, 221, r"conventional cell,\\ vacancy and jump"),            # the middle of a picture from the upper left,
          "shell1": (703, 310, 234, r"first neighbors:\\ 12 at $d_1=a/\sqrt2$"),        # half the width that its turn takes, its label
          "shell2": (1192, 310, 234, r"second neighbors:\\ 6 at $d_2=a$"),
          "shell3": (1681, 310, 234, r"third neighbors:\\ 24 at $d_3=\sqrt{3/2}\,a$")}
cm = lambda x, y: f"({x / CM:.4f},{(H - y) / CM:.4f})"     # a place given in pixels from the upper left
labels = "".join(f"\\node[capt] at {cm(x, BOTTOM)} {{{text}}};\n\\node[above lab] at {cm(x - half, TOP)} {{\\panel{{{letter}}}}};\n"
                 for (x, y, half, text), letter in zip(PANELS.values(), "abcd"))


def page(k):
    """the labels, and the line with the letter a at the places that render_lattice_turn.py has written for picture k"""
    x0, y0 = PANELS["cell"][:2]
    width, height = Image.open(os.path.join(pictures, f"cell_{k:03d}.png")).size
    with open(os.path.join(pictures, f"cell_{k:03d}.json")) as file:
        at = {name: cm(x0 + (x - 0.5) * width, y0 + (0.5 - y) * height) for name, (x, y) in json.load(file).items()}
    return (f"\\begin{{tikzpicture}}\n\\path[use as bounding box] (0,0) rectangle {cm(W, 0)};\n{labels}"
            f"\\draw[dbl, black!75] {at['from']} -- {at['to']};\n\\node[black!75, inner sep=0pt] at {at['a']} {{$a$}};\n\\end{{tikzpicture}}\n")


wanted = [int(a) for a in sys.argv[5:]] or range(COUNT)
film.pages("".join(page(k) for k in wanted), frames, CM)
for number, k in enumerate(wanted, 1):                     # the line and its label lie above the cell and its atoms in every picture
    turning = [(os.path.join(pictures, f"{name}_{k:03d}.png"), x, y) for name, (x, y, *_) in PANELS.items()]
    film.put(Image.open(os.path.join(frames, f"labels_{number:03d}.png")), turning).save(os.path.join(frames, f"{k:03d}.png"))
if not sys.argv[5:]:
    film.webp(frames, range(0, COUNT, 2), FPS, (760, 285), 66, os.path.join(webp, "lattice.webp"))
    film.mp4(frames, COUNT, FPS, os.path.join(mp4, "lattice.mp4"))
