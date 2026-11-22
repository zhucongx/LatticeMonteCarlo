# Pictures of frame_fields.png, one frame of Ansys: the solute atoms, the clusters, and the largest cluster colored by the field
# binding_energy.   blender -b -P render_frame.py -- <output folder> <step>.xyz.gz
import gzip
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import scene as s

out, frame = s.output_folder()[:2]
os.makedirs(out, exist_ok=True)
AZ, EL = 114, 16                                       # a side of the box from which the clusters of the frame do not hide each other
CLUSTER = ["8A63A9", "3E9B5C", "D487AC", "E6C343", "2F9AA0", "8B6A4A", "565C66"]     # one color per cluster_id, none is that of a species
BIND, VMAX = ["3F78B5", "F4F4F4", "C8423B"], 0.15      # BIND of style.py, and the range of the color bar in eV


def read_xyz(path):
    """atoms of an extended XYZ file as dictionaries of the fields, and the length of the cubic box"""
    with gzip.open(path, "rt") as f:
        n, comment = int(f.readline()), f.readline()
        fields = comment.split("Properties=")[1].split()[0].split(":")
        length = float(comment.split('Lattice="')[1].split()[0])
        atoms = []
        for _ in range(n):
            words, atom = f.readline().split(), {}
            for name, kind, size in zip(fields[0::3], fields[1::3], map(int, fields[2::3])):
                values = [{"S": str, "I": int, "R": float}[kind](w) for w in words[:size]]
                atom[name], words = (values[0] if size == 1 else values), words[size:]
            atoms.append(atom)
    return atoms, length


atoms, L = read_xyz(frame)
for a in atoms:                                        # a cluster that crosses the boundary of the box is drawn as a whole
    a["xyz"] = a["pos"]
    if a["cluster_id"]:
        a["xyz"] = [p + L * round((c - p) / L) for p, c in zip(a["pos"], a["geometry_center"])]
members = [a for a in atoms if a["cluster_id"] == 1]
vacancy = next(a for a in atoms if a["species"] == "X")


def solute(a):
    e = a["species"]
    if e == "X":
        s.vacancy(a["xyz"], 1.75)
    elif e != "Al" and a["cluster_id"]:
        s.sphere(a["xyz"], 1.50 * s.RADIUS[e], e, 1.0, 0.22)
    elif e != "Al":
        s.sphere(a["xyz"], 0.80 * s.RADIUS[e], s.mix(e, "FFFFFF", 0.55))


def cluster(a):
    if a["cluster_id"]:
        e = a["species"]
        s.sphere(a["xyz"], 0.95 if e == "Al" else 1.50 * s.RADIUS[e], CLUSTER[a["cluster_id"] - 1], 1.0, 0.16 if e == "Al" else 0.22)


def bound(energy):
    t = max(-1.0, min(1.0, energy / VMAX))
    return s.mix(BIND[1], BIND[0 if t < 0 else 2], abs(t))


for name, draw in (("frame_solutes", solute), ("frame_clusters", cluster)):
    s.reset(); s.new_scene(1000, 1000, 160)
    cam = s.camera((L / 2, L / 2, L / 2), AZ, EL, 88)
    s.box((0, 0, 0), (L, L, L), 0.15)
    for a in atoms:
        draw(a)
    s.lights(cam)
    s.render(os.path.join(out, name + ".png"), {"vacancy": vacancy["xyz"]})

s.reset(); s.new_scene(1000, 1000, 160)
cam = s.camera([sum(a["xyz"][k] for a in members) / len(members) for k in range(3)], AZ, 60, 22)
for i, a in enumerate(members):
    s.sphere(a["xyz"], 0.88 * s.RADIUS[a["species"]], bound(a["binding_energy"]), 1.0, 0.06)
    for b in members[:i]:
        if sum((p - q) ** 2 for p, q in zip(a["xyz"], b["xyz"])) < 3.5 ** 2:
            s.rod(a["xyz"], b["xyz"], 0.14, "B4B9BF")
s.lights(cam)
s.render(os.path.join(out, "frame_binding.png"))
