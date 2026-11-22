# jump_path.webp and jump_path.mp4: a jump of an atom into the vacancy and back, next to the energy along the path, in 6 s.
# Make them in this folder with, for any scratch folder T:
#   blender -b -P render_jump_path.py -- T/pictures
#   pdflatex -output-directory=T fig_jump_path.tex
#   python3 animate_jump_path.py T/pictures T/fig_jump_path.pdf T/frames <output folder>
import os
import subprocess
import sys
import tempfile

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import jump_path_motion as motion
from PIL import Image

pictures, pdf, frames, out = sys.argv[1:5]
os.makedirs(frames, exist_ok=True)
FPS, COUNT = 30, 180                                       # of the MP4; the WebP has every second frame
THERE = 10 + motion.STEPS                                  # frame in which the atom arrives at F; it leaves I in frame 10
LEAVE = THERE + 42                                         # frame in which it leaves F
BACK = LEAVE + motion.STEPS                                # frame in which it arrives at I
FADE = 14                                                  # frames in which the mark moves from one barrier to the other
W, H, CM = 2160, 1080, 135                                 # size of a frame, and of one centimeter of fig_jump_path.tex, in pixels
T, UNIT = (11.60, 6.22), (3.35, 2.45)                      # transition state, units of x and E in fig_jump_path.tex, in centimeters
FINE = 8                                                   # the ball is placed to this fraction of a pixel


def picture_number(n):
    return min(max(n - THERE + motion.STEPS, 0), motion.STEPS) if n < LEAVE else max(BACK - n, 0)


def fade(n, start):
    """rises from 0 to 1 in FADE frames"""
    w = min(max((n - start) / FADE, 0), 1)
    return w * w * (3 - 2 * w)


def page(number, scale):
    with tempfile.NamedTemporaryFile(suffix=".png") as f:
        subprocess.run(["gs", "-q", "-dNOPAUSE", "-dBATCH", "-sDEVICE=pngalpha", f"-r{2.54 * CM * scale}", "-dTextAlphaBits=4",
                        "-dGraphicsAlphaBits=4", f"-dFirstPage={number}", f"-dLastPage={number}", f"-sOutputFile={f.name}", pdf], check=True)
        return Image.open(f.name).convert("RGBA")


jump, jump_back = (page(number, 2).resize((W, H), Image.LANCZOS) for number in (1, 2))
ball = page(3, FINE)
size = ball.width // FINE + 2


def put_ball(im, x, y):
    """the ball with its center at the pixel coordinates x, y"""
    left, top = round(FINE * x - ball.width / 2), round(FINE * y - ball.height / 2)
    patch = Image.new("RGBA", (size * FINE, size * FINE), (255, 255, 255, 0))
    patch.paste(ball, (left % FINE, top % FINE))
    im.alpha_composite(patch.resize((size, size), Image.BOX), (left // FINE, top // FINE))


for n in range(COUNT):
    k = picture_number(n)
    x = motion.place(k)
    im = Image.new("RGBA", (W, H), "white")
    im.alpha_composite(Image.open(os.path.join(pictures, f"jump_{k:03d}.png")).convert("RGBA"), (0, 108))
    im.alpha_composite(Image.blend(jump, jump_back, fade(n, THERE + FADE / 2) - fade(n, BACK + 1)))
    put_ball(im, CM * (T[0] + UNIT[0] * x), H - CM * (T[1] + UNIT[1] * motion.energy(x)))
    im.convert("RGB").save(os.path.join(frames, f"{n:03d}.png"))

args = []
for n in range(0, COUNT, 2):
    small = os.path.join(frames, f"small_{n:03d}.png")
    Image.open(os.path.join(frames, f"{n:03d}.png")).resize((760, 380), Image.LANCZOS).save(small)
    args += ["-d", str(round(1000 * (n + 2) / FPS) - round(1000 * n / FPS)), small]
subprocess.run(["img2webp", "-loop", "0", "-sharp_yuv", "-lossy", "-q", "98", "-m", "4", *args, "-o", os.path.join(out, "jump_path.webp")],
               check=True)
# the colors of the frames, in the form that players expect at this size
colors = "scale=out_color_matrix=bt709:out_range=tv,setparams=color_primaries=bt709:color_trc=iec61966-2-1"
subprocess.run(["ffmpeg", "-loglevel", "error", "-y", "-framerate", str(FPS), "-i", os.path.join(frames, "%03d.png"), "-vf", colors,
                "-c:v", "libx264", "-pix_fmt", "yuv420p", "-crf", "18", "-movflags", "+faststart", os.path.join(out, "jump_path.mp4")],
               check=True)
print([f"{os.path.getsize(os.path.join(out, 'jump_path.' + e)) // 1024} kB {e}" for e in ("webp", "mp4")])
