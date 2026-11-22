# Pictures of environment.png: the sites around a vacancy jump in the local frame of the jump.
#   blender -b -P render_environment.py -- <output folder>
import itertools
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import scene as s

out = s.output_folder()[0]
os.makedirs(out, exist_ok=True)
A = 4.046
H = A / 2
R2 = 2 ** 0.5
vecs = [v for v in itertools.product(range(-2, 3), repeat=3) if sum(v) % 2 == 0 and sum(c * c for c in v) in (2, 4, 6)]
add = lambda p, q: tuple(a + b for a, b in zip(p, q))
d2 = lambda p, q: sum((a - b) ** 2 for a, b in zip(p, q))
VAC, ATOM = (0, 0, 0), (1, 1, 0)
env = sorted({v for v in vecs if v != ATOM} | {add(ATOM, v) for v in vecs if add(ATOM, v) != VAC})
assert len(env) == 58
near = [v for v in env if d2(v, VAC) == 2 or d2(v, ATOM) == 2]
assert len(near) == 18


def pos(v):
    """local frame: x along the jump <110>, y along <-110>, z along <001>, origin at the midpoint"""
    return ((v[0] + v[1] - 1) / R2 * H, (v[1] - v[0]) / R2 * H, v[2] * H)


def xyz(v):
    return (v[0] + v[1] - 1, v[1] - v[0], v[2])         # in units of d/2, d/2, a/2


# numbers of the 18 first neighbors: by increasing x, for the same x in the order (y+,z+), (y-,z+), (y-,z-), (y+,z-)
def number_key(v):
    x, y, z = xyz(v)
    order = {(1, 1): 0, (-1, 1): 1, (-1, -1): 2, (1, -1): 3}
    if y != 0 and z != 0:
        k = order[(1 if y > 0 else -1, 1 if z > 0 else -1)]
    else:
        k = 0 if y >= 0 else 1
    return (x, k)


numbered = {v: i + 1 for i, v in enumerate(sorted(near, key=number_key))}
AZ, EL = -16, 12                                        # y points to the upper left
edges = [(p, q) for p, q in itertools.combinations(near, 2) if d2(p, q) == 2]


def kind(v):
    x, y, z = xyz(v)
    return {3: "end", 2: "ring", 1: "plane", 0: "window"}[abs(x)], x


MMM = {"end": "66717E", "ring": "4A80BD", "plane": "4FA06B", "window": "DFA83A"}
light = lambda c: s.mix(c, "FFFFFF", 0.46)
MM2 = {(k, x): (light(c) if x < 0 else c) for k, c in MMM.items() for x in (-3, -2, -1, 0, 1, 2, 3)}


def panel(name, what):
    s.reset(); s.new_scene(1100, 860, 160)
    cam = s.camera((0, 0, 0), AZ, EL, 14.0 if what == "all" else 10.4)
    if what == "all":
        for p, q in itertools.combinations(env, 2):
            if d2(p, q) == 2:
                s.rod(pos(p), pos(q), 0.05, "C3C7CC")
        for v in env:
            s.sphere(pos(v), 0.46, "7E858F" if v in near else "C4C9CF")
    else:
        for p, q in edges:
            kp, kq = kind(p), kind(q)
            same = kp[0] == kq[0] and (what == "mmm" or kp == kq)
            s.rod(pos(p), pos(q), 0.12, (MMM[kp[0]] if what == "mmm" else MM2[kp]) if same else "B4B9BF")
        for v in near:
            k = kind(v)
            s.sphere(pos(v), 0.68, MMM[k[0]] if what == "mmm" else MM2[k])
    s.vacancy(pos(VAC), 0.68)
    s.sphere(pos(ATOM), 0.68, "C8423B")
    s.lights(cam)
    labels = {str(numbered[v]): pos(v) for v in near}
    labels.update(o=(0, 0, 0), x=(1.5, 0, 0), y=(0, 1.5, 0), z=(0, 0, 1.5))         # the axes of the local frame, 1.5 Angstrom long
    s.render(os.path.join(out, name + ".png"), labels)


panel("env_all", "all")
panel("env_mmm", "mmm")
panel("env_mm2", "mm2")
