# Pictures of a full turn for the animation lattice: the conventional cubic cell, in which an atom jumps into the vacancy and back, and
# the three neighbor shells of a site in a box of eight cells. The places of the line that marks an edge of the cell, and of its label,
# are written next to the pictures of the cell. The pictures have the size that they take in a frame of animate_lattice.py. Pictures that
# exist already are kept, so that a run can be continued.
#   blender -b -P render_lattice_turn.py -- <output folder> [first picture] [last picture] [step]
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
RADIUS = 0.62                                              # drawn radius of an atom of the cell
JUMP = 24                                                  # pictures that a jump takes
VAC, ATOM = (1, 0, 1), (1, 1, 2)                           # in units of a/2: centers of the front face and of the top face
EDGE = {"from": (0, 0, A + 1.7), "to": (A, 0, A + 1.7),    # the line that marks an edge of the cell, above the upper edge in front
        "a": (A / 2, -0.45, A + 2.18)}                     # in picture 0 and clear of the atoms in all pictures, and the place of its label
site = lambda v: tuple(c * A / 2 for c in v)


def way(k):
    """part of the way from its site to the vacancy that the atom has made in picture k: it jumps when the two sites are seen from the
    side, around the pictures COUNT / 5 (there) and 7 COUNT / 10 (back)"""
    there, back = ((k - c * COUNT) / JUMP + 0.5 for c in (0.2, 0.7))
    w = min(max(there, 0.0), 1.0) - min(max(back, 0.0), 1.0)
    return w * w * (3 - 2 * w)


def cell(path, azimuth, f):
    s.reset(); s.new_scene(498, 448, 64)
    cam = s.camera((A / 2, A / 2, A / 2), azimuth, EL, 8.0)
    s.box((0, 0, 0), (A, A, A), 0.022)
    for v in itertools.product(range(3), repeat=3):
        if sum(v) % 2 == 0 and v not in (VAC, ATOM):
            s.sphere(site(v), RADIUS, "Al")
    start, end = site(ATOM), site(VAC)
    for place, part in ((start, f), (end, 1 - f)):         # the share of the vacancy on a site, drawn once the atom has left the site
        r = min(RADIUS * part ** (1 / 3), part * A / 2 ** 0.5 - RADIUS)
        if r > 0.02:
            s.vacancy(place, r)                            # the ring faces the camera
    s.sphere(tuple(a + (b - a) * f for a, b in zip(start, end)), RADIUS, "Al")
    s.lights(cam)
    s.render(path[:-4] + ".part.png", EDGE)
    os.replace(path[:-4] + ".part.json", path[:-4] + ".json")
    os.replace(path[:-4] + ".part.png", path)              # an interrupted picture does not count


def shell(path, azimuth, color, d2):
    s.reset(); s.new_scene(484, 472, 64)
    cam = s.camera((0, 0, 0), azimuth, EL, 12.2)
    pts = [site(v) for v in itertools.product(range(-2, 3), repeat=3) if sum(v) % 2 == 0 and sum(c * c for c in v) == d2]
    s.box((-A, -A, -A), (A, A, A), 0.028)                  # the eight cells around the site
    for p, q in s.hull(pts, color, 0.13):
        s.rod(p, q, 0.075, s.mix(color, "FFFFFF", 0.15))
    for p in pts:
        s.sphere(p, 0.52, s.mix(color, "FFFFFF", 0.12))
    s.sphere((0, 0, 0), 0.56, "555555")
    s.rod((0, 0, 0), max(pts, key=lambda p: p[0] * 0.8 - p[1] * 0.9 + p[2] * 0.55), 0.09, "555555")
    s.lights(cam)
    s.render(path[:-4] + ".part.png")
    os.replace(path[:-4] + ".part.png", path)


SHELLS = {1: ("red", 2), 2: ("blue", 4), 3: ("green", 6)}  # color, squared distance in (a/2)^2
for k in range(first, last + 1, step):
    azimuth = AZ + 360 * k / COUNT
    path = os.path.join(out, f"cell_{k:03d}.png")
    if not os.path.exists(path):
        cell(path, azimuth, way(k))
    for number, look in SHELLS.items():
        path = os.path.join(out, f"shell{number}_{k:03d}.png")
        if not os.path.exists(path):
            shell(path, azimuth, *look)
