# Pictures of a KMC run for clustering.webp and clustering.mp4: the solute atoms and the vacancy, without the solvent, drawn by
# solute_box.py. The camera turns from one azimuth to another. The folders hold the configurations of one trajectory, see
# animate_clustering.py. Pictures that exist already are kept.
#   blender -b -P render_clustering.py -- <output folder> <size in pixels> <samples> <first azimuth> <last azimuth>
#                                         <first picture> <last picture> <folder> ...
import glob
import math
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import scene as s
import solute_box

PICTURES, EVEN, LAST = 100, 300_000, 130_000_000           # picture i shows the step next to EVEN (exp(c i) - 1), the last one LAST

out = s.output_folder()[0]
size, samples, azimuth, last_azimuth, first, final = (int(v) for v in s.output_folder()[1:7])
files = {}
for folder in s.output_folder()[7:]:
    for f in glob.glob(os.path.join(folder, "*.cfg.gz")):
        if os.path.basename(f)[:-7].isdigit():
            files.setdefault(int(os.path.basename(f)[:-7]), f)
c = math.log(1 + LAST / EVEN) / (PICTURES - 1)
steps = [min(files, key=lambda n: abs(n - EVEN * (math.exp(c * i) - 1))) for i in range(PICTURES)]
faces = solute_box.faces(out, [files[n] for n in set(steps)])

for i in range(first, final + 1):
    path = os.path.join(out, f"{i:03d}")
    if os.path.exists(path + ".png"):
        continue
    u = i / (PICTURES - 1)
    solute_box.picture(path, files[steps[i]], steps[i], faces, azimuth + (last_azimuth - azimuth) * u * u * (3 - 2 * u), size, samples)
