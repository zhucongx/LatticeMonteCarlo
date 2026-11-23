# Pictures of a full turn for the animation environment, two for every step of the turn: the 58 sites within the third neighbor shell of
# one of the two sites of a jump, the vacancy (white with a red ring) and the atom that jumps (red). The 18 first neighbors of one of the
# two are colored by their class, under the point group mmm in turn_mmm_*.png and under mm2 in turn_mm2_*.png; the other 40 sites are
# small and pale. Bonds join first neighbors within the 18 sites and within the 40 sites. Pictures that exist already are kept, so that
# a run can be continued.
#   blender -b -P render_environment_turn.py -- <output folder> [first picture] [last picture] [step]
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
H, R2 = A / 2, 2 ** 0.5
add = lambda p, q: tuple(a + b for a, b in zip(p, q))
d2 = lambda p, q: sum((a - b) ** 2 for a, b in zip(p, q))
shells = [v for v in itertools.product(range(-2, 3), repeat=3) if sum(c * c for c in v) in (2, 4, 6) and sum(v) % 2 == 0]
VAC, ATOM = (0, 0, 0), (1, 1, 0)                           # sites in units of a/2
sites = {v for v in shells if v != ATOM} | {add(ATOM, v) for v in shells if add(ATOM, v) != VAC}
near = sorted(v for v in sites if 2 in (d2(v, VAC), d2(v, ATOM)))
far = sorted(sites - set(near))
assert len(near) == 18 and len(far) == 40
# The code tells the sites apart by their distance from the midpoint of the jump and by their coordinate x along the jump, which it
# counts from the vacancy to the atom: under mmm by |x|, under mm2 by x (PositionCompareMMM and PositionCompareMM2 of the program).
# Among the 18 first neighbors |x| alone, in units of half the jump distance, tells the classes apart.
CLASS = {3: "66717E", 2: "4A80BD", 1: "4FA06B", 0: "DFA83A"}


def color(v, group):
    """mmm: one color for every |x|; mm2: the same colors, pale on the side of the vacancy"""
    x = v[0] + v[1] - 1
    return s.mix(CLASS[-x], "FFFFFF", 0.46) if group == "mm2" and x < 0 else CLASS[abs(x)]


def pos(v):
    """x along the jump <110>, from the vacancy to the atom, y along <-110>, z along <001>, origin at the midpoint"""
    return ((v[0] + v[1] - 1) / R2 * H, (v[1] - v[0]) / R2 * H, v[2] * H)


def picture(path, azimuth, group):
    s.reset(); s.new_scene(920, 828, 64)
    cam = s.camera((0, 0, 0), azimuth, 14, 13.0)
    for sites_of_a_size in (near, far):
        for p, q in itertools.combinations(sites_of_a_size, 2):
            if d2(p, q) == 2 and sites_of_a_size is near:
                s.rod(pos(p), pos(q), 0.07, "AEB5BD")
            elif d2(p, q) == 2:                            # pale and without a shadow
                r = s.rod(pos(p), pos(q), 0.025, "CDD2D8", 1.0, True)
                r.visible_shadow = r.visible_diffuse = r.visible_glossy = False
    for v in near:
        s.sphere(pos(v), 0.60, color(v, group))
    for v in far:
        s.sphere(pos(v), 0.27, "C6CBD2")
    s.vacancy(pos(VAC), 0.66)                              # the ring faces the camera
    s.sphere(pos(ATOM), 0.66, "red")
    s.lights(cam)
    s.render(path + ".part.png")                           # an interrupted picture does not count
    os.replace(path + ".part.png", path)


for k in range(first, last + 1, step):
    for group in ("mmm", "mm2"):
        path = os.path.join(out, f"turn_{group}_{k:03d}.png")
        if not os.path.exists(path):
            picture(path, -16 + 360 * k / COUNT, group)
