# Pictures of the animation vacancy_jumps: the atoms around the vacancy during some single steps of a KMC run.
# Pictures that exist already are kept, and the script stops after the given number of seconds, so that it can be called again.
#   blender -b -P render_vacancy_jumps.py -- <output folder> <run folder> [seconds] [picture numbers, e.g. 0,45]
import os
import sys
import time

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import bpy
import numpy as np
import scene as s
import vacancy_jumps_data as data
from vacancy_jumps_data import H, STEPS, SUB

out, run, *more = s.output_folder()
limit = float(more[0]) if more else 1e9
numbers = [int(n) for n in more[1].split(",")] if len(more) > 1 else range(STEPS * SUB + 1)
RADIUS = 8.85                                              # of the region, in Angstrom; it holds all first neighbors of the path
BALL, SOLVENT = 0.95, 0.36                                 # drawn radii
TRAIL = 7                                                  # jumps after which the trail has faded
RED = "C8423B"
NEIGHBORS = [(i, j, k) for i in (-1, 0, 1) for j in (-1, 0, 1) for k in (-1, 0, 1) if i * i + j * j + k * k == 2]
TURN = np.array([[1, 0, 0], [0, 0, 1], [0, -1, 0]])        # turns the crystal so that the path is not hidden

d = data.load(run)
species, path, selected = d["species"], d["path"], d["selected"]
center = path.mean(0)                                      # ends up as the center of the smallest sphere around the path
for i in range(4000):
    farthest = path[((path - center) ** 2).sum(1).argmax()]
    center = center + (farthest - center) / (i + 2)
near = np.rint(center).astype(int)
site = (d["site"] - near + d["period"] // 2) % d["period"] - d["period"] // 2 + near        # the image next to the center
inside = np.where(((site - center) ** 2).sum(1) ** 0.5 * H < RADIUS)[0]
jumping = sorted(set(selected[:STEPS]))
assert {tuple(p + v) for p in path for v in NEIGHBORS} <= {tuple(p) for p in site[inside]}
xyz = lambda p: TURN @ ((p - center) * H)

s.reset(); scene = s.new_scene(1080, 1080, 48)
scene.cycles.seed = 1                                      # the same noise in every picture
scene.render.use_persistent_data = True
cam = s.camera((0, 0, 0), 21, 10, 19.2)
s.lights(cam)
view = np.array(cam.location.normalized())


def atom(p, element, highlight=0.0):
    """solvent atoms are small, and paler towards the back and the rim; highlight marks the one that jumps"""
    if element != "Al":
        return s.sphere(p, BALL * s.RADIUS[element], element)
    depth, rim = np.clip((6 - p @ view) / 14, 0, 1), np.clip((np.linalg.norm(p) - 4.5) / 5, 0, 1)
    pale = round(20 * float(0.3 + 0.6 * depth + 0.1 * rim)) / 20
    return s.sphere(p, SOLVENT * (1 + 0.65 * highlight), s.mix(s.mix("Al", "FFFFFF", pale), "858D96", round(highlight, 2)))


def line(p, q, radius, alpha):
    s.rod(p, q, radius, RED, round(alpha, 2), True)


def plain(objects):
    """no shadows on other atoms: a dark ball would look like a marked one"""
    for ob in objects:
        ob.visible_shadow = ob.visible_diffuse = False


for i in inside:
    if i not in jumping and species[i] != "X":
        atom(xyz(site[i]), species[i])
plain(o for o in bpy.data.objects if o.type == "MESH")


def picture(n, k, name):
    """picture k of step n of the window; the last one, k = 0 of step STEPS, is the state after all jumps"""
    fixed = set(bpy.data.objects)
    f = data.progress(k)
    here, v = {i: site[i] for i in jumping}, path[0]
    for a in selected[:n]:
        here[a], v = v, here[a]
    assert (v == path[n]).all()
    moving, last = selected[n] if n < STEPS else None, selected[n - 1] if n else None
    marked = {last: max(0.0, 0.5 - k / 4)}
    if moving is not None:
        marked[moving] = min(1.0, 0.5 + min(k, SUB - k) / 4)
    for i in jumping:
        p = xyz(here[i])
        if i == moving:
            p = now = p + (xyz(v) - p) * f
        atom(p, species[i], marked.get(i, 0.0))
    old, new = xyz(v), xyz(here[moving]) if moving is not None else None
    # the vacancy shrinks at its old site and grows at the new one, next to the atom; it is not drawn where it would be a mere speck
    size = 0 if moving is None else BALL * s.RADIUS[species[moving]] if species[moving] != "Al" else SOLVENT * (1 + 0.65 * marked[moving])
    for p, part in ((old, 1 - f), (new, f)):
        r = min(BALL * part ** (1 / 3), part * 2 ** 0.5 * H - size)
        if part > 0 and r > 0.3:
            s.vacancy(p, r)
    # thin lines to the twelve neighbors of the vacancy
    for p, strength in ((v, 1 - 2 * f), (here.get(moving), 2 * f - 1)):
        if strength > 0:
            for w in NEIGHBORS:
                line(xyz(p), xyz(p + w), 0.05, 0.6 * strength)
    # the trail: the last jumps of the vacancy
    pair = lambda m: tuple(sorted((tuple(path[m]), tuple(path[m + 1]))))
    trail = {}
    for age in range(min(TRAIL, n)):
        trail[pair(n - age - 1)] = max(trail.get(pair(n - age - 1), 0), 1 - (age + f) / TRAIL)
    # the jump that is being made: its own trail behind the atom and an older one ahead of it, both only as far as a vacancy is drawn
    if f > 0:
        older = trail.pop(pair(n), 0)
        if older > 0 and f < 0.6:
            line(now, old, 0.11, older)
        if f > 0.4:
            line(new, now, 0.11, 1.0)
    for (p, q), alpha in trail.items():
        line(xyz(np.array(p)), xyz(np.array(q)), 0.11, alpha)
    created = [o for o in bpy.data.objects if o not in fixed]
    plain(created)
    s.render(name + ".part.png")
    os.replace(name + ".part.png", name)                   # an interrupted picture does not count
    for o in created:
        bpy.data.objects.remove(o)


begin = time.time()
for number in numbers:
    name = os.path.join(out, f"{number:03d}.png")
    if not os.path.exists(name) and time.time() - begin < limit:
        picture(*divmod(number, SUB), name)
