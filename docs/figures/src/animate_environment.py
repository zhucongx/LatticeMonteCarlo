# environment.webp and environment.mp4: the sites around a jump with the classes of the 18 first neighbors, under mmm on the left and
# under mm2 on the right, one turn in 7.5 s. Make them in this folder with, for any scratch folder T:
#   blender -b -P render_environment_turn.py -- T/pictures
#   pdflatex -output-directory=T fig_environment_turn.tex
#   python3 animate_environment.py T/pictures T/fig_environment_turn.pdf T/frames <output folder>
# Needs gs, cwebp, webpmux, ffmpeg and Pillow.
import os
import subprocess
import sys
import tempfile

from PIL import Image, ImageChops

pictures, pdf, frames, out = sys.argv[1:5]
os.makedirs(frames, exist_ok=True)
FPS, COUNT = 24, 180                                       # of the MP4; the WebP has every second picture
W, H, CM = 1920, 1080, 120                                 # size of a frame, and of one centimeter of fig_environment_turn.tex, in pixels
LEFT, TOP = {"mmm": 20, "mm2": 980}, 10                    # where the pictures of the two panels go
QUALITY = ["-q", "44", "-m", "6", "-sharp_yuv"]            # of the WebP

with tempfile.NamedTemporaryFile(suffix=".png") as f:      # the labels and the legend, the same in every frame
    subprocess.run(["gs", "-q", "-dNOPAUSE", "-dBATCH", "-sDEVICE=pngalpha", f"-r{2.54 * CM * 2}", "-dTextAlphaBits=4",
                    "-dGraphicsAlphaBits=4", f"-sOutputFile={f.name}", pdf], check=True)
    labels = Image.open(f.name).convert("RGBA").resize((W, H), Image.LANCZOS)

# The WebP holds its first frame in full and of every other frame the rectangle in which something turns: the labels are stored once
# and stand still, and a frame takes two thirds of the bytes.
small = lambda k: os.path.join(frames, f"small_{k:03d}.png")
box = None
for k in range(COUNT):
    im = Image.new("RGBA", (W, H), "white")
    for group, left in LEFT.items():
        im.alpha_composite(Image.open(os.path.join(pictures, f"turn_{group}_{k:03d}.png")).convert("RGBA"), (left, TOP))
    im.alpha_composite(labels)
    im.convert("RGB").save(os.path.join(frames, f"{k:03d}.png"))
    if k % 2 == 0:
        im = im.convert("RGB").resize((760, 428), Image.LANCZOS)
        im.save(small(k))
        if k == 0:
            first = im
        else:
            b = ImageChops.difference(first, im).getbbox()
            box = b if box is None else (min(box[0], b[0]), min(box[1], b[1]), max(box[2], b[2]), max(box[3], b[3]))
x, y = box[0] // 2 * 2, box[1] // 2 * 2                    # a rectangle starts at even coordinates
crop = ["-crop", str(x), str(y), str(box[2] - x), str(box[3] - y)]
args = []
for k in range(0, COUNT, 2):
    webp = small(k)[:-3] + "webp"
    subprocess.run(["cwebp", "-quiet", *QUALITY, *(crop if k else []), small(k), "-o", webp], check=True)
    duration = round(1000 * (k + 2) / FPS) - round(1000 * k / FPS)
    args += ["-frame", webp, f"+{duration}+{x if k else 0}+{y if k else 0}+0-b"]                         # replaces what is there
subprocess.run(["webpmux", *args, "-loop", "0", "-bgcolor", "255,255,255,255", "-o", os.path.join(out, "environment.webp")], check=True)
# the colors of the frames, in the form that players expect at this size
colors = "scale=out_color_matrix=bt709:out_range=tv,setparams=color_primaries=bt709:color_trc=iec61966-2-1"
subprocess.run(["ffmpeg", "-loglevel", "error", "-y", "-framerate", str(FPS), "-i", os.path.join(frames, "%03d.png"), "-vf", colors,
                "-c:v", "libx264", "-pix_fmt", "yuv420p", "-crf", "18", "-threads", "4", "-movflags", "+faststart",
                os.path.join(out, "environment.mp4")], check=True)
print([f"{os.path.getsize(os.path.join(out, 'environment.' + e)) // 1024} kB {e}" for e in ("webp", "mp4")])
