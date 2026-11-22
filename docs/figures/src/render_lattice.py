# Pictures of lattice.png: the conventional cubic cell with a vacancy, and the three neighbor shells.
#   blender -b -P render_lattice.py -- <output folder>
import itertools
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import scene as s

out = s.output_folder()[0]
os.makedirs(out, exist_ok=True)
A = 4.046
AZ, EL = 18, 14

# (a) conventional cell
s.reset(); s.new_scene(900, 860)
cam = s.camera((A / 2, A / 2, A / 2), AZ, EL, 8.0)
s.box((0, 0, 0), (A, A, A), 0.022)
sites = [v for v in itertools.product(range(3), repeat=3) if sum(v) % 2 == 0]
vac = (1, 0, 1)                                          # center of the front face; the jump is drawn in fig_lattice.tex
for v in sites:
    p = tuple(c * A / 2 for c in v)
    if v == vac:
        s.vacancy(p, 0.62)
    else:
        s.sphere(p, 0.62, "Al")
s.lights(cam)
s.render(os.path.join(out, "lattice_cell.png"))

# (b)-(d) neighbor shells
shells = {1: ("red", 2), 2: ("blue", 4), 3: ("green", 6)}     # color, squared length in (a/2)^2
vectors = [v for v in itertools.product(range(-2, 3), repeat=3) if sum(v) % 2 == 0]
for number, (color, d2) in shells.items():
    s.reset(); s.new_scene(900, 860)
    pts = [tuple(c * A / 2 for c in v) for v in vectors if sum(c * c for c in v) == d2]
    s.box((-A, -A, -A), (A, A, A), 0.028)                     # the eight cells around the site
    for p, q in s.hull(pts, color, 0.13):
        s.rod(p, q, 0.075, s.mix(color, "FFFFFF", 0.15))
    for p in pts:
        s.sphere(p, 0.52, s.mix(color, "FFFFFF", 0.12))
    s.sphere((0, 0, 0), 0.56, "555555")
    example = max(pts, key=lambda p: p[0] * 0.8 - p[1] * 0.9 + p[2] * 0.55)
    s.rod((0, 0, 0), example, 0.09, "555555")
    cam = s.camera((0, 0, 0), AZ, EL, 12.2)
    s.lights(cam)
    s.render(os.path.join(out, f"shell{number}.png"))
