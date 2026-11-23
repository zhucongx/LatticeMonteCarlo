# cluster_shapes.webp and cluster_shapes.mp4: the eight cluster shapes of the energy model, each in a conventional cell, in one turn of
# 7.5 s, with the number of such clusters per site of an infinite fcc crystal. The names, the numbers and the legend stay.
# Make them in this folder with, for any scratch folder T:
#   blender -b -P render_cluster_shapes_turn.py -- T/pictures
#   python3 animate_cluster_shapes.py T/pictures T/frames <folder for cluster_shapes.webp> <folder for cluster_shapes.mp4>
# With numbers of frames after these arguments only those frames are made, for a look. Needs what turn_film.py needs.
import os
import sys

from PIL import Image

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import turn_film as film

pictures, frames, webp, mp4 = sys.argv[1:5]
os.makedirs(frames, exist_ok=True)
FPS, COUNT = 24, 180                                       # of the MP4; the WebP has every third picture
W, H, CM = 1654, 1080, 105.5                               # size of a frame, and of one centimeter, in pixels
X, Y = (202, 619, 1036, 1453), (166, 668)                  # the middle of the cells of a column and of a row, from the upper left
SHAPES = [("site", "single site", 1), ("pair1", "pair 1", 6), ("pair2", "pair 2", 3), ("pair3", "pair 3", 12),
          ("triplet111", "triplet 1-1-1", 8), ("triplet112", "triplet 1-1-2", 12), ("triplet113", "triplet 1-1-3", 24),
          ("triplet123", "triplet 1-2-3", 24)]             # file, name, number per site
PLACE = {file: (X[i % 4], Y[i // 4]) for i, (file, _, _) in enumerate(SHAPES)}
cm = lambda x, y: f"({x / CM:.4f},{(H - y) / CM:.4f})"     # a place given in pixels from the upper left
labels = "".join(f"\\node[capt] at {cm(PLACE[file][0], PLACE[file][1] + 204)} {{{name}\\\\[-1pt] {number} per site}};\n"
                 for file, name, number in SHAPES)
legend = "".join(f"\\draw[{color}, line width=3.6pt] {cm(x, 1046)} -- ++(0.6,0);\n\\node[leg] at {cm(x + 0.75 * CM, 1046)} {{{text}}};\n"
                 for x, color, text in ((138, "cOne", "1: first neighbors"), (618, "cTwo", "2: second neighbors"),
                                        (1145, "cThree", "3: third neighbors")))
film.pages(f"\\begin{{tikzpicture}}\n\\path[use as bounding box] (0,0) rectangle {cm(W, 0)};\n{labels}{legend}\\end{{tikzpicture}}\n", frames, CM)
page = Image.open(os.path.join(frames, "labels_001.png"))
for k in [int(a) for a in sys.argv[5:]] or range(COUNT):
    film.put(page, [(os.path.join(pictures, f"{file}_{k:03d}.png"), *PLACE[file]) for file in PLACE]).save(os.path.join(frames, f"{k:03d}.png"))
if not sys.argv[5:]:
    film.webp(frames, range(0, COUNT, 3), FPS, (760, 496), 80, os.path.join(webp, "cluster_shapes.webp"))
    film.mp4(frames, COUNT, FPS, os.path.join(mp4, "cluster_shapes.mp4"))
