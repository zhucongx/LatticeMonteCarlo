# Pictures of the animation trapped_vacancy: the solute atoms and the vacancy of a KMC run, without the solvent, in the look of
# solute_box.py. Solute atoms in clusters have their full size, all others are small and pale; the clusters are those of Ansys in
# the configuration of the picture. The red line is the path of the vacancy during the last TRAIL pictures, paler with age. The vacancy
# is a white ball; animate_trapped_vacancy.py draws it again on top, with its ring, so that neither atoms nor the path hide it, and adds
# the numbers of the clusters. Pictures that exist already are kept.
#   blender -b -P render_trapped_vacancy.py -- <data file of trapped_vacancy_data.py> <output folder> <size in pixels> <samples>
#                                              [first picture] [last picture]
import itertools
import json
import os
import sys

import bpy
import numpy as np
from mathutils import Vector

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import scene as s
import trapped_vacancy_data as data
from trapped_vacancy_data import A, NAMED, NEIGHBORS

TRAIL = 5                                                  # pictures after which the path has faded
VACANCY = 1.75
RED = "C8423B"

source, out, size, samples, *some = s.output_folder()
d = data.load(source)
species, planes, shown = d["species"], d["planes"], d["shown"]
first, final = (int(v) for v in some) if some else (0, len(shown) - 1)
L = planes * A / 2


def best_faces():
    """the lattice plane at which the drawn box begins, for x, y and z: the one at which a face cuts the fewest contacts in clusters"""
    cut = np.zeros((3, planes), int)
    for step in sorted(set(shown.tolist())):
        for group, _ in data.clusters(data.sites(d, step), species, planes):
            for p, v, k in itertools.product(group, NEIGHBORS, range(3)):
                if v[k] == 1 and tuple((p + v) % planes) in group:
                    cut[k, (p[k] + 1) % planes] += 1
    return [int(n) for n in cut.argmin(1)]


os.makedirs(out, exist_ok=True)
note = os.path.join(out, "faces.json")
if not os.path.exists(note):
    with open(note, "w") as f:
        json.dump(best_faces(), f)
with open(note) as f:
    faces = np.array(json.load(f))
cell = lambda site: (np.asarray(site) - faces) % planes   # where a site is drawn, in units of A / 2
place = lambda site: (cell(site) + 0.5) * A / 2


def trail(i):
    """{(site, site): age} of the jumps of the vacancy during the last TRAIL pictures; age 0 is the time since the picture before"""
    age = {}
    for back in range(min(TRAIL, i)):
        path = d["path"][shown[i - back - 1]:shown[i - back] + 1] % planes
        for jump in {tuple(sorted((tuple(p), tuple(q)))) for p, q in zip(path[:-1], path[1:])}:
            age.setdefault(jump, back)
    return age


for i in range(first, final + 1):
    name = os.path.join(out, f"{i:03d}")
    if os.path.exists(name + ".png"):
        continue
    site = data.sites(d, shown[i])
    large = set().union(*(group for group, _ in data.clusters(site, species, planes)))
    s.reset(); s.new_scene(int(size), int(size), int(samples))
    cam = s.camera((L / 2, L / 2, L / 2 - 1.5), 24, 16, 88)
    s.box((0, 0, 0), (L, L, L), 0.15)
    balls = {"large": bpy.data.collections.new("large"), "small": bpy.data.collections.new("small")}
    vacancy = Vector(place(site[d["vacancy"]]))
    balls["large"].objects.link(s.sphere(vacancy, 0.93 * VACANCY, "white"))
    for p, e in zip(site, species):
        if e in ("Al", "X"):
            continue
        if tuple(p) in large:
            balls["large"].objects.link(s.sphere(place(p), 1.50 * s.RADIUS[e], e, 1.0, 0.22))
        else:
            ball = s.sphere(place(p), 0.70 * s.RADIUS[e], s.mix(e, "white", 0.6))
            ball.visible_shadow = False
            balls["small"].objects.link(ball)
    for kind in balls.values():                            # a light for each kind of balls: shadows of large balls on large balls only
        before = set(bpy.data.objects)
        s.lights(cam)
        light, = set(bpy.data.objects) - before
        light.light_linking.receiver_collection = light.light_linking.blocker_collection = kind
    for (p, q), age in trail(i).items():
        jump = (np.array(q) - p + planes // 2) % planes - planes // 2
        # a jump through a face of the box is drawn in two halves, each up to its face
        inside = (cell(q) - cell(p) == jump).all()
        for a, b in [(place(p), place(q))] if inside else [(place(p), place(p) + jump * A / 4), (place(q) - jump * A / 4, place(q))]:
            line = s.rod(a, b, 0.20, RED, round(1 - age / TRAIL, 2), True)
            line.visible_shadow = line.visible_diffuse = line.visible_glossy = False
    # for animate_trapped_vacancy.py: where the vacancy, the middle of every numbered cluster and the large balls are in the picture
    labels = {"vacancy": vacancy, "rim": vacancy + cam.rotation_quaternion @ Vector((0.96 * VACANCY, 0, 0))}
    for k in range(1, NAMED + 1):
        labels[f"cluster {k}"] = place(site[d["names"][i] == k]).mean(0)
    for n, p in enumerate(sorted(large)):
        labels[f"ball {n}"] = place(p)
    s.render(name + ".part.png", labels)
    with open(name + ".part.json") as f:
        labels = json.load(f)
    with open(name + ".json", "w") as f:
        json.dump({"step": int(shown[i]), **{k: v for k, v in labels.items() if k[:4] != "ball"},
                   "large": [v for k, v in labels.items() if k[:4] == "ball"]}, f)
    os.remove(name + ".part.json")
    os.replace(name + ".part.png", name + ".png")          # an interrupted picture does not count
