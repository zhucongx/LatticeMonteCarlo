# Pictures of one vacancy jump for the animation jump_path: the atom on its way from its site (picture 000) to the vacancy, and the 18
# first neighbors of the two sites. The four atoms of the window it passes move aside a little; this is a sketch of the relaxation.
# Pictures that exist already are kept, so that a run can be continued.
#   blender -b -P render_jump_path.py -- <output folder> [first picture] [last picture]
import itertools
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import jump_path_motion as motion
import scene as s

out, *numbers = s.output_folder()
first, last = (int(numbers[0]), int(numbers[1])) if numbers else (0, motion.STEPS)
os.makedirs(out, exist_ok=True)
A = 4.046                                                  # lattice parameter in Angstrom
H, R2 = A / 2, 2 ** 0.5
RADIUS, OPENING = 0.64, 0.20                               # drawn radius of Al; shift of the window atoms, in Angstrom
add = lambda p, q: tuple(a + b for a, b in zip(p, q))
d2 = lambda p, q: sum((a - b) ** 2 for a, b in zip(p, q))
shell = [v for v in itertools.product((-1, 0, 1), repeat=3) if sum(c * c for c in v) == 2]
VAC, ATOM = (0, 0, 0), (1, 1, 0)                           # sites in units of a/2
near = sorted({v for v in shell if v != ATOM} | {add(ATOM, v) for v in shell if add(ATOM, v) != VAC})
assert len(near) == 18


def pos(v):
    """x along the jump, from the atom to the vacancy; z along <001>; origin at the midpoint"""
    return (-(v[0] + v[1] - 1) / R2 * H, (v[1] - v[0]) / R2 * H, v[2] * H)


def picture(path, f):
    s.reset(); s.new_scene(1080, 810, 192)
    cam = s.camera((0, 0, 0), 12, 23, 10.1)
    push = OPENING * (1 - (2 * f - 1) ** 2) ** 2
    where = {}
    for v in near:
        x, y, z = pos(v)
        if abs(x) < 0.1:                                   # window
            n = (y * y + z * z) ** 0.5
            y, z = y * (1 + push / n), z * (1 + push / n)
        where[v] = (x, y, z)
    for p, q in itertools.combinations(near, 2):
        if d2(p, q) == 2:
            s.rod(where[p], where[q], 0.08, "A9B4BF")
    for v in near:
        s.sphere(where[v], RADIUS, "Al")
    start, end = pos(ATOM), pos(VAC)
    atom_radius = RADIUS * s.RADIUS["Mg"]
    for site, part in ((start, f), (end, 1 - f)):          # the share of the vacancy on a site, drawn once the atom has left the site
        r = min(RADIUS * part ** (1 / 3), part * (end[0] - start[0]) - atom_radius)
        if r > 0.02:
            s.vacancy(site, r)
    s.sphere(tuple(a + (b - a) * f for a, b in zip(start, end)), atom_radius, "Mg")
    s.lights(cam)
    s.render(path + ".part.png")                           # an interrupted picture does not count
    os.replace(path + ".part.png", path)


for k in range(first, last + 1):
    path = os.path.join(out, f"jump_{k:03d}.png")
    if not os.path.exists(path):
        picture(path, motion.fraction(k))
