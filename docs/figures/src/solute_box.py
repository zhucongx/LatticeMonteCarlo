# The box of a KMC run for render_clustering.py and render_quench.py: the solute atoms and the vacancy, without the solvent. Solute atoms
# in clusters have their full size, all others are small and pale. The clusters are those of Ansys with smallest_cluster_criteria MIN:
# groups of at least MIN solute atoms that are connected as nearest neighbors, with the vacancy counted. The vacancy is a white ball; its
# ring is added by the animate_*.py script, so that it stays visible when atoms are in front of it.
import itertools
import json
import os

import bpy
from mathutils import Vector

import scene as s

MIN = 10
VACANCY = 1.75
A = 4.046
NEIGHBORS = [v for v in itertools.product((-1, 0, 1), repeat=3) if sum(c * c for c in v) == 2]


def read(file):
    """{site: element} of the solute atoms, the site of the vacancy, the length of the box, and the number of its lattice planes; sites
    are in units of half a lattice parameter"""
    pos, species, L = s.read_cfg(file)
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


def best_faces(files):
    """the lattice plane at which the drawn box begins, for x, y and z: the one at which a face cuts the fewest contacts in clusters"""
    cut = {}
    for file in files:
        solutes, vacancy, _, planes = read(file)
        for group in clusters(solutes, vacancy, planes):
            for p, v, k in itertools.product(group, NEIGHBORS, range(3)):
                if v[k] == 1 and tuple((a + b) % planes for a, b in zip(p, v)) in group:
                    face = k, (p[k] + 1) % planes
                    cut[face] = cut.get(face, 0) + 1
    return [min(range(planes), key=lambda n: cut.get((k, n), 0)) for k in range(3)]


def faces(out, files):
    """best_faces for the configurations of all pictures; kept in the output folder, so that every call draws the same box"""
    os.makedirs(out, exist_ok=True)
    note = os.path.join(out, "faces.json")
    if not os.path.exists(note):
        with open(note, "w") as f:
            json.dump(best_faces(files), f)
    with open(note) as f:
        return json.load(f)


def picture(path, file, step, faces, azimuth, size, samples):
    """path.png, and in path.json the step and where the vacancy and the right end of its ring are in the picture"""
    solutes, vacancy, L, planes = read(file)
    in_clusters = set().union(*clusters(solutes, vacancy, planes))
    place = lambda site: Vector([((n - face) % planes + 0.5) * A / 2 for n, face in zip(site, faces)])
    s.reset(); s.new_scene(size, size, samples)
    cam = s.camera((L / 2, L / 2, L / 2 - 1.5), azimuth, 16, 88)
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
        json.dump({"step": step, **labels}, f)
    os.remove(path + ".part.json")
    os.replace(path + ".part.png", path + ".png")          # an interrupted picture does not count
