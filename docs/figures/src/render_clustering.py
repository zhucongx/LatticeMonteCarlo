# Pictures of a KMC run for clustering.webp and clustering.mp4: the solute atoms and the vacancy, without the solvent. Solute atoms in
# clusters have their full size, all others are small and pale. The clusters are those of Ansys with smallest_cluster_criteria MIN: groups
# of at least MIN solute atoms that are connected as nearest neighbors, with the vacancy counted. The vacancy is a white ball; its ring is
# added by animate_clustering.py, so that it stays visible when atoms are in front of it. The camera turns from one azimuth to another.
# The folders hold the configurations of one trajectory, see animate_clustering.py. Pictures that exist already are kept.
#   blender -b -P render_clustering.py -- <output folder> <size in pixels> <samples> <first azimuth> <last azimuth>
#                                         <first picture> <last picture> <folder> ...
import glob
import itertools
import json
import math
import os
import sys

import bpy
from mathutils import Vector

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import scene as s

PICTURES, EVEN, LAST = 100, 300_000, 130_000_000           # picture i shows the step next to EVEN (exp(c i) - 1), the last one LAST
MIN = 10
VACANCY = 1.75
A = 4.046
NEIGHBORS = [v for v in itertools.product((-1, 0, 1), repeat=3) if sum(c * c for c in v) == 2]

out = s.output_folder()[0]
size, samples, azimuth, last_azimuth, first, final = (int(v) for v in s.output_folder()[1:7])
files = {}
for folder in s.output_folder()[7:]:
    for f in glob.glob(os.path.join(folder, "*.cfg.gz")):
        if os.path.basename(f)[:-7].isdigit():
            files.setdefault(int(os.path.basename(f)[:-7]), f)
c = math.log(1 + LAST / EVEN) / (PICTURES - 1)
steps = [min(files, key=lambda n: abs(n - EVEN * (math.exp(c * i) - 1))) for i in range(PICTURES)]


def read(step):
    """{site: element} of the solute atoms, the site of the vacancy, the length of the box, and the number of its lattice planes; sites
    are in units of half a lattice parameter"""
    pos, species, L = s.read_cfg(files[step])
    planes = round(2 * L / A)
    sites = [tuple(round(2 * x / A) % planes for x in p) for p in pos]
    return {p: e for p, e in zip(sites, species) if e not in ("Al", "X")}, sites[species.index("X")], L, planes


def clusters(solutes, vacancy, planes):
    found, left = [], {*solutes, vacancy}
    while left:
        group, new = set(), {left.pop()}
        while new:
            group |= new
            new = {tuple((a + b) % planes for a, b in zip(p, v)) for p in new for v in NEIGHBORS} & left
            left -= new
        if len(group) >= MIN:
            found.append(group)
    return found


def best_faces():
    """the lattice plane at which the drawn box begins, for x, y and z: the one at which a face cuts the fewest contacts in clusters"""
    cut = {}
    for step in set(steps):
        solutes, vacancy, _, planes = read(step)
        for group in clusters(solutes, vacancy, planes):
            for p, v, k in itertools.product(group, NEIGHBORS, range(3)):
                if v[k] == 1 and tuple((a + b) % planes for a, b in zip(p, v)) in group:
                    face = k, (p[k] + 1) % planes
                    cut[face] = cut.get(face, 0) + 1
    return [min(range(planes), key=lambda n: cut.get((k, n), 0)) for k in range(3)]


os.makedirs(out, exist_ok=True)
note = os.path.join(out, "faces.json")
if not os.path.exists(note):
    with open(note, "w") as f:
        json.dump(best_faces(), f)
with open(note) as f:
    faces = json.load(f)

for i in range(first, final + 1):
    path = os.path.join(out, f"{i:03d}")
    if os.path.exists(path + ".png"):
        continue
    solutes, vacancy, L, planes = read(steps[i])
    in_clusters = set().union(*clusters(solutes, vacancy, planes))
    place = lambda site: Vector([((n - face) % planes + 0.5) * A / 2 for n, face in zip(site, faces)])
    s.reset(); s.new_scene(size, size, samples)
    u = i / (PICTURES - 1)
    cam = s.camera((L / 2, L / 2, L / 2 - 1.5), azimuth + (last_azimuth - azimuth) * u * u * (3 - 2 * u), 16, 88)
    s.box((0, 0, 0), (L, L, L), 0.15)
    large, small = bpy.data.collections.new("large"), bpy.data.collections.new("small")
    large.objects.link(s.sphere(place(vacancy), 0.93 * VACANCY, "white"))
    for site, e in solutes.items():
        if site in in_clusters:
            large.objects.link(s.sphere(place(site), 1.50 * s.RADIUS[e], e, 1.0, 0.22))
        else:
            ball = s.sphere(place(site), 0.70 * s.RADIUS[e], s.mix(e, "white", 0.6))
            ball.visible_shadow = False
            small.objects.link(ball)
    for balls in (large, small):                           # a light for each kind of balls: shadows of large balls on large balls only
        before = set(bpy.data.objects)
        s.lights(cam)
        light, = set(bpy.data.objects) - before
        light.light_linking.receiver_collection = light.light_linking.blocker_collection = balls
    rim = place(vacancy) + cam.rotation_quaternion @ Vector((0.96 * VACANCY, 0, 0))
    s.render(path + ".part.png", {"vacancy": place(vacancy), "rim": rim})
    with open(path + ".part.json") as f:
        labels = json.load(f)
    with open(path + ".json", "w") as f:
        json.dump({"step": steps[i], **labels}, f)
    os.remove(path + ".part.json")
    os.replace(path + ".part.png", path + ".png")          # an interrupted picture does not count
