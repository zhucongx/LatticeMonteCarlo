# Pictures of a full turn for frame_fields.webp and frame_fields.mp4: the clusters of one frame of Ansys in the periodic box, colored by
# the field cluster_id (id_###.png) and by the field binding_energy (binding_###.png). The frame is step 10 000 000 of the KMC run of
# docs/README.md, see frame_fields_data.py. Solvent atoms of a cluster are smaller than its solute atoms, and atoms that are in no
# cluster are not drawn. The vacancy is a white ball; its ring is added by animate_frame_fields.py, so that it stays visible when atoms
# are in front of it. The box is shifted so that its faces cut no cluster. Pictures that exist already are kept, so that a run can be
# continued.
#   blender -b -P render_frame_turn.py -- <output folder> <step>.xyz.gz [first picture] [last picture] [step]
import os
import sys

import bpy
from mathutils import Vector

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import frame_fields_data as data
import scene as s

COUNT = 180                                                # pictures of a turn
WIDTH, HEIGHT = 960, 900
AZIMUTH, ELEVATION, SCALE = 156, 16, 88                    # in the first picture no cluster hides another; the rest as in solute_box.py
VACANCY = 1.75
out, frame, *numbers = s.output_folder()
first, last, step = [*map(int, numbers), *(0, COUNT - 1, 1)[len(numbers):]]
os.makedirs(out, exist_ok=True)
members, L, cut = data.load(frame)
print("contacts in clusters that the faces of the box cut:", cut)
vacancy = Vector(next(a["xyz"] for a in members if a["species"] == "X"))
COLOR = {"id": lambda a: data.CLUSTER[a["cluster_id"] - 1], "binding": lambda a: data.bound(a["binding_energy"])}

for name, color in COLOR.items():
    todo = [k for k in range(first, last + 1, step) if not os.path.exists(os.path.join(out, f"{name}_{k:03d}.png"))]
    if not todo:
        continue
    s.reset(); s.new_scene(WIDTH, HEIGHT, 64)
    s.box((0, 0, 0), (L, L, L), 0.15)
    for a in members:
        e = a["species"]
        if e == "X":
            s.sphere(a["xyz"], 0.93 * VACANCY, "white")
        else:                                              # the solvent atoms of a cluster are smaller
            s.sphere(a["xyz"], 0.95 if e == "Al" else 1.50 * s.RADIUS[e], color(a), 1.0, 0.16 if e == "Al" else 0.22)
    fixed = set(bpy.data.objects)
    for k in todo:
        path = os.path.join(out, f"{name}_{k:03d}")
        # the camera looks at a point below the center of the box, so that the turning box is in the middle of the picture
        cam = s.camera((L / 2, L / 2, L / 2 - 3.1), AZIMUTH + 360 * k / COUNT, ELEVATION, SCALE)
        s.lights(cam)                                      # the light turns with the camera
        rim = vacancy + cam.rotation_quaternion @ Vector((0.96 * VACANCY, 0, 0))
        s.render(path + ".part.png", {"vacancy": vacancy, "rim": rim})
        os.replace(path + ".part.json", path + ".json")
        os.replace(path + ".part.png", path + ".png")      # an interrupted picture does not count
        for ob in set(bpy.data.objects) - fixed:
            bpy.data.objects.remove(ob)
