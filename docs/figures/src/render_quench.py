# Pictures of a KMC run that follows a time-temperature table, for quench.webp and quench.mp4: the solute atoms and the vacancy as in
# clustering.webp, drawn by solute_box.py, with the camera at rest. The pictures are equally spaced in the time of the run, and each
# shows the configuration that the run wrote nearest to its time. The run is described in animate_quench.py. Pictures that exist
# already are kept.
#   blender -b -P render_quench.py -- <output folder> <size in pixels> <samples> <first picture> <last picture> <run folder>
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import scene as s
import solute_box

PICTURES, AZIMUTH = 120, 24

out, run = s.output_folder()[0], s.output_folder()[5]
size, samples, first, final = (int(v) for v in s.output_folder()[1:5])
with open(os.path.join(run, "kmc_log.txt")) as f:
    rows = [line.split()[:2] for line in f][1:]            # steps and time
file = lambda step: os.path.join(run, f"{step}.cfg.gz")
time = {int(step): float(t) for step, t in rows if os.path.exists(file(step))}
end = float(rows[-1][1])
steps = [min(time, key=lambda n: abs(time[n] - end * i / (PICTURES - 1))) for i in range(PICTURES)]
faces = solute_box.faces(out, [file(n) for n in set(steps)])

for i in range(first, final + 1):
    path = os.path.join(out, f"{i:03d}")
    if not os.path.exists(path + ".png"):
        solute_box.picture(path, file(steps[i]), steps[i], faces, AZIMUTH, size, samples)
