# environment.webp and environment.mp4: the sites around a jump, one turn in 7.5 s.
# Make them in this folder with, for any scratch folder T:
#   blender -b -P render_environment_turn.py -- T/pictures
#   python3 animate_environment.py T/pictures T/frames <output folder>
import os
import subprocess
import sys

from PIL import Image

pictures, frames, out = sys.argv[1:4]
os.makedirs(frames, exist_ok=True)
FPS, COUNT = 24, 180                                       # of the MP4; the WebP has every second picture
args = []
for k in range(COUNT):
    path = os.path.join(frames, f"{k:03d}.png")
    im = Image.open(os.path.join(pictures, f"turn_{k:03d}.png")).convert("RGBA")
    white = Image.new("RGBA", im.size, "white")
    white.alpha_composite(im)
    white.convert("RGB").save(path)
    if k % 2 == 0:
        small = os.path.join(frames, f"small_{k:03d}.png")
        white.convert("RGB").resize((600, 540), Image.LANCZOS).save(small)
        args += ["-d", str(round(1000 * (k + 2) / FPS) - round(1000 * k / FPS)), small]

# Everything moves from one frame to the next, so every frame is stored in full (-kmax 1).
subprocess.run(["img2webp", "-loop", "0", "-kmin", "0", "-kmax", "1", "-sharp_yuv", "-lossy", "-q", "56", "-m", "6", *args,
                "-o", os.path.join(out, "environment.webp")], check=True)
# the colors of the frames, in the form that players expect at this size
colors = "scale=out_color_matrix=bt709:out_range=tv,setparams=color_primaries=bt709:color_trc=iec61966-2-1"
subprocess.run(["ffmpeg", "-loglevel", "error", "-y", "-framerate", str(FPS), "-i", os.path.join(frames, "%03d.png"), "-vf", colors,
                "-c:v", "libx264", "-pix_fmt", "yuv420p", "-crf", "18", "-movflags", "+faststart", os.path.join(out, "environment.mp4")],
               check=True)
print([f"{os.path.getsize(os.path.join(out, 'environment.' + e)) // 1024} kB {e}" for e in ("webp", "mp4")])
