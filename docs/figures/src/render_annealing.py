# Pictures of a simulated annealing run for annealing.webp: the box with the solute atoms and the vacancy, without
# the solvent. All pictures are shifted in the periodic box by the same vector, such that the solutes of the reference configuration
# are in the middle. A picture is named after its configuration; pictures that exist already are kept, so that a run can be continued.
#   blender -b -P render_annealing.py -- <output folder> <size in pixels> <samples> <reference configuration> <configuration> ...
import math
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import scene as s


def periodic_center(points, L):
    """mean of positions in a periodic box, component by component on a circle"""
    out = []
    for c in zip(*points):
        angles = [2 * math.pi * v / L for v in c]
        out.append(math.atan2(sum(map(math.sin, angles)), sum(map(math.cos, angles))) * L / (2 * math.pi) % L)
    return out


out, size, samples, reference = s.output_folder()[:4]
pos, species, L = s.read_cfg(reference)
center = periodic_center([p for p, e in zip(pos, species) if e != "Al"], L)
planes = sorted({round(p[0], 3) for p in pos})
d = planes[1] - planes[0]                                  # distance of two lattice planes, half a lattice parameter
shift = [(math.floor((L / 2 - c) / d) + 0.5) * d for c in center]         # no atom lies on a face of the box

for cfg in s.output_folder()[4:]:
    path = os.path.join(out, os.path.basename(cfg).split(".")[0] + ".png")
    if os.path.exists(path):
        continue
    pos, species, L = s.read_cfg(cfg)
    s.reset(); s.new_scene(int(size), int(size), int(samples))
    cam = s.camera((L / 2, L / 2, L / 2), 24, 16, 61)
    s.box((0, 0, 0), (L, L, L), 0.11)
    for p, e in zip(pos, species):
        p = tuple((c + v) % L for c, v in zip(p, shift))
        if e in ("Mg", "Zn"):
            s.sphere(p, 1.45 * s.RADIUS[e], e, 1.0, 0.16)
        elif e == "X":
            s.vacancy(p, 1.45)
    s.lights(cam)
    s.render(path + ".part.png")                           # an interrupted picture does not count
    os.replace(path + ".part.png", path)
