# Pictures of cluster_shapes.png: the eight cluster shapes of the energy model, drawn in a conventional cell.
#   blender -b -P render_clusters.py -- <output folder>
import itertools
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import scene as s

out = s.output_folder()[0]
os.makedirs(out, exist_ok=True)
A = 4.046
COLOR = {2: "red", 4: "blue", 6: "green"}
TYPES = {
    "site": [(0, 0, 0)],
    "pair1": [(0, 0, 0), (1, 0, 1)],
    "pair2": [(0, 0, 0), (2, 0, 0)],
    "pair3": [(0, 0, 0), (2, 1, 1)],
    "triplet111": [(0, 0, 0), (1, 0, 1), (1, 1, 0)],
    "triplet112": [(0, 0, 0), (2, 0, 0), (1, 0, 1)],
    "triplet113": [(0, 0, 0), (2, 1, 1), (1, 0, 1)],
    "triplet123": [(0, 0, 0), (2, 1, 1), (2, 0, 0)],
}
sites = [v for v in itertools.product(range(3), repeat=3) if sum(v) % 2 == 0]
for name, members in TYPES.items():
    s.reset(); s.new_scene(640, 600)
    cam = s.camera((A / 2, A / 2, A / 2), 18, 14, 7.5)
    s.box((0, 0, 0), (A, A, A), 0.02)
    for v in sites:
        p = tuple(c * A / 2 for c in v)
        s.sphere(p, 0.50 if v in members else 0.22, "7B828C" if v in members else "D3D7DC")
    for p, q in itertools.combinations(members, 2):
        d2 = sum((x - y) ** 2 for x, y in zip(p, q))
        s.rod([c * A / 2 for c in p], [c * A / 2 for c in q], 0.14, COLOR[d2])
    s.lights(cam)
    s.render(os.path.join(out, f"cluster_{name}.png"))
