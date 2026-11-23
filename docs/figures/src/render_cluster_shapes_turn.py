# Pictures of a full turn for the animation cluster_shapes: the eight cluster shapes of the energy model, each in a conventional cell, with
# the bonds in the color of the neighbor shell of their length. The pictures have the size that they take in a frame of
# animate_cluster_shapes.py. Pictures that exist already are kept, so that a run can be continued.
#   blender -b -P render_cluster_shapes_turn.py -- <output folder> [first picture] [last picture] [step]
import itertools
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import scene as s

COUNT = 180                                                # pictures of a turn
out, *numbers = s.output_folder()
first, last, step = [*map(int, numbers), *(0, COUNT - 1, 1)[len(numbers):]]
os.makedirs(out, exist_ok=True)
A = 4.046                                                  # lattice parameter in Angstrom
AZ, EL = 18, 14                                            # the view of picture 0
COLOR = {2: "red", 4: "blue", 6: "green"}                  # of a bond, by its squared length in (a/2)^2
SHAPES = {                                                 # the sites of a cluster, in units of a/2
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
place = lambda v: tuple(c * A / 2 for c in v)


def picture(path, azimuth, members):
    s.reset(); s.new_scene(430, 400, 64)
    cam = s.camera((A / 2, A / 2, A / 2), azimuth, EL, 7.5)
    s.box((0, 0, 0), (A, A, A), 0.02)
    for v in sites:
        s.sphere(place(v), 0.50 if v in members else 0.22, "7B828C" if v in members else "D3D7DC")
    for p, q in itertools.combinations(members, 2):
        s.rod(place(p), place(q), 0.14, COLOR[sum((x - y) ** 2 for x, y in zip(p, q))])
    s.lights(cam)
    s.render(path[:-4] + ".part.png")                      # an interrupted picture does not count
    os.replace(path[:-4] + ".part.png", path)


for k in range(first, last + 1, step):
    for name, members in SHAPES.items():
        path = os.path.join(out, f"{name}_{k:03d}.png")
        if not os.path.exists(path):
            picture(path, AZ + 360 * k / COUNT, members)
