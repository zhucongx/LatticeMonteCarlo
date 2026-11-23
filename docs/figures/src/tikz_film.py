# Shared by animate_kmc_step.py, animate_kmc_second_order.py and animate_metropolis.py, the animations that are line drawings:
# draws their frames with TikZ (one LaTeX run for all pages that are new, Ghostscript to PNG) and writes the WebP and the MP4.
# Needs pdflatex, gs, img2webp, ffmpeg and Pillow.
import hashlib
import os
import subprocess

from PIL import Image

FPS = 30                                                   # frames per second of the MP4; the WebP takes every second frame
FIRST, LAST = 0.5, 2.5                                     # seconds on the first and on the last frame
SRC = os.path.dirname(os.path.abspath(__file__))
HEAD = "\\documentclass[11pt,tikz,border=0pt]{standalone}\n\\input{figstyle}\n"


def ease(x):
    """from 0 to 1, slow at both ends"""
    x = min(max(x, 0.0), 1.0)
    return x * x * (3 - 2 * x)


def mix(a, b, w):
    """the color between the RGB colors a (w = 0) and b (w = 1), as TikZ wants it"""
    return "{rgb,255:red,%d;green,%d;blue,%d}" % tuple(round(x + (y - x) * w) for x, y in zip(a, b))


def draw(pages, preamble, folder, width):
    """pages: the TikZ code of every frame, each a picture `width` cm wide and 9/16 of that high. Returns its pictures, 1920 x 1080
    pixels (drawn three times as large and scaled down); frames with the same code share one, and one that exists is not drawn again."""
    os.makedirs(folder, exist_ok=True)
    preamble = HEAD + preamble + "\\begin{document}\n"
    files = [os.path.join(folder, hashlib.md5((preamble + p).encode()).hexdigest()[:12] + ".png") for p in pages]
    todo = {png: p for png, p in zip(files, pages) if not os.path.exists(png)}
    print(f"{len(pages)} frames, {len(set(files))} different pictures, {len(todo)} to draw", flush=True)
    if not todo:
        return files
    tex = os.path.join(folder, "frames.tex")
    with open(tex, "w") as file:
        file.write(preamble + "".join(todo.values()) + "\\end{document}\n")
    done = subprocess.run(["pdflatex", "-interaction=batchmode", "-halt-on-error", f"-output-directory={folder}", tex],
                          stdout=subprocess.DEVNULL, env=dict(os.environ, TEXINPUTS=f"{SRC}:{folder}:"))
    if done.returncode:
        raise SystemExit(f"pdflatex failed, see {tex[:-3]}log")

    for number, png in enumerate(todo, 1):                 # one page after the other: one Ghostscript at a time
        subprocess.run(["gs", "-q", "-dNOPAUSE", "-dBATCH", "-sDEVICE=png16m", f"-r{3 * 1920 * 2.54 / width:.4f}", "-g5760x3240",
                        "-dDownScaleFactor=3", f"-dFirstPage={number}", f"-dLastPage={number}", f"-sOutputFile={png}.part",
                        tex[:-3] + "pdf"], check=True)
        os.replace(png + ".part", png)
    return files


def write(files, out, name):
    """files: the picture of every frame of the film without the two stops. Writes <out>/<name>.webp (760 pixels wide, half the frames,
    a picture that stays gets its time instead of being repeated; lossless, so that nothing that fades leaves a trace on the white)
    and <out>/<name>.mp4."""
    last = len(files) - 1
    half = files[:-1:2] + files[-1:]                       # the last frame is always among them
    args, ms = [], [round(2000 * n / FPS) for n in range(len(half) + 1)]
    for n, png in enumerate(half):
        d = ms[n + 1] - ms[n] + (round(1000 * FIRST) if n == 0 else 0) + (round(1000 * LAST) if n == len(half) - 1 else 0)
        small = png[:-4] + "_small.png"
        if args and args[-1] == small:
            args[-2] = str(int(args[-2]) + d)
            continue
        if not os.path.exists(small):
            Image.open(png).convert("RGB").resize((760, 428), Image.LANCZOS).save(small)
        args += ["-d", str(d), small]
    webp, mp4 = (os.path.join(out, f"{name}.{e}") for e in ("webp", "mp4"))
    subprocess.run(["img2webp", "-loop", "0", "-lossless", "-exact", "-m", "6", *args, "-o", webp], check=True, stdout=subprocess.DEVNULL)
    video = subprocess.Popen(["ffmpeg", "-loglevel", "error", "-y", "-f", "rawvideo", "-pix_fmt", "rgb24", "-s", "1920x1080", "-r", str(FPS),
                              "-i", "-", "-vf", "scale=out_color_matrix=bt709:out_range=tv,setparams=color_primaries=bt709:color_trc=iec61966-2-1",
                              "-c:v", "libx264", "-threads", "4", "-pix_fmt", "yuv420p", "-crf", "18", "-movflags", "+faststart", mp4],
                             stdin=subprocess.PIPE)
    for png in files[:1] * round(FPS * FIRST) + files + files[-1:] * round(FPS * LAST):
        video.stdin.write(Image.open(png).convert("RGB").tobytes())
    video.stdin.close()
    video.wait()
    print(f"{len(args) // 3} pictures, {last / FPS + FIRST + LAST:.1f} s,",
          ", ".join(f"{os.path.getsize(p) // 1024} kB {p[-4:]}" for p in (webp, mp4)))
