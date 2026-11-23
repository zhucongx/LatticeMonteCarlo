# Shared by animate_cluster_definition.py, animate_barrier_groups.py, animate_energy_change.py and animate_restart.py, the animations
# that replace line drawings: draws their frames with TikZ (one LaTeX run for all pages that are new, Ghostscript to PNG) and writes
# <name>.webp and <name>.mp4. It is not started on its own. Needs pdflatex, gs, img2webp, ffmpeg and Pillow.
import hashlib
import os
import subprocess

from PIL import Image

FPS = 30                                                   # frames per second of the film and of the MP4; the WebP takes every second frame
WIDTH = 16.256                                             # width of a picture in cm: 1920 pixels at 300 dpi, and 760 pixels in the WebP
CM = 1920 / WIDTH                                          # pixels of a centimeter
WEBP = ["-near_lossless", "60", "-q", "100", "-m", "6"]     # options of img2webp: no pixel is off by more than 2 of 255
BEGIN = "\\documentclass[11pt,tikz,border=0pt]{standalone}\n\\input{figstyle}\n"


def frames(seconds):
    return round(FPS * seconds)


def ramp(f, start, length=1):
    """0 up to the frame start, 1 from the frame start + length on, linear in between"""
    return min(max((f - start) / length, 0.0), 1.0)


def ease(x):
    return x * x * (3 - 2 * x)


def opacity(x):
    """option of a node or a scope that is drawn with the opacity x"""
    return f"opacity={x:.3f}"


def make(name, preamble, pictures, height, folder, out, only=()):
    """pictures: the TikZ code of every frame of the film, height: of the pictures in pixels (1920 wide). With numbers of frames in
    only, these alone are drawn and their files are printed."""
    here = os.path.dirname(os.path.abspath(__file__))
    os.makedirs(folder, exist_ok=True)
    start = BEGIN + preamble + "\\begin{document}\n"
    # The white mark in the corner cannot be seen. With it every page has something transparent, and Ghostscript draws a page in the
    # same way whatever the other pages of the LaTeX run hold; without it the edges of a picture can differ by single pixels.
    box = (f"\\path[use as bounding box] (0,0) rectangle ({WIDTH},{height / CM:.5f});\n"
           "\\fill[white, opacity=0.5] (0,0) rectangle (0.01,0.01);\n")
    wanted = list(only) or range(len(pictures))
    pages = [f"\\begin{{tikzpicture}}\n{box}{pictures[f]}\\end{{tikzpicture}}\n" for f in wanted]
    files = [os.path.join(folder, hashlib.md5((start + p).encode()).hexdigest()[:12] + ".png") for p in pages]
    todo = {png: p for png, p in zip(files, pages) if not os.path.exists(png)}
    print(f"{name}: {len(pictures)} frames, {len(set(files))} different pictures, {len(todo)} to draw", flush=True)
    if todo:
        tex = os.path.join(folder, "frames.tex")
        with open(tex, "w") as file:
            file.write(start + "".join(todo.values()) + "\\end{document}\n")
        subprocess.run(["pdflatex", "-interaction=batchmode", "-halt-on-error", f"-output-directory={folder}", tex], check=True,
                       stdout=subprocess.DEVNULL, env=dict(os.environ, TEXINPUTS=here + ":"))
        for number, png in enumerate(todo, 1):             # one page after the other, three times as many pixels inside
            subprocess.run(["gs", "-q", "-dNOPAUSE", "-dBATCH", "-sDEVICE=png16m", "-r900", "-dDownScaleFactor=3", f"-dFirstPage={number}",
                            f"-dLastPage={number}", f"-sOutputFile={png}.part", tex[:-3] + "pdf"], check=True)
            os.replace(png + ".part", png)
    if only:
        print("\n".join(f"{f} {png}" for f, png in zip(wanted, files)))
        return
    # the WebP: every second frame, a picture that stays for the sum of its times
    args, ms = [], [round(1000 * f / FPS) for f in range(0, len(files) + 2, 2)]
    small_height = round(height * 760 / 1920)
    for n, png in enumerate(files[::2]):
        d, small = ms[n + 1] - ms[n], png[:-4] + "_small.png"
        if args and args[-1] == small:
            args[-2] = str(int(args[-2]) + d)
            continue
        if not os.path.exists(small):
            Image.open(png).convert("RGB").resize((760, small_height), Image.LANCZOS).save(small)
        args += ["-d", str(d), small]
    webp, mp4 = (os.path.join(out, f"{name}.{e}") for e in ("webp", "mp4"))
    subprocess.run(["img2webp", "-loop", "0", *WEBP, *args, "-o", webp], check=True)
    # the MP4: 1920 x 1080, a lower picture in the middle on white
    video = subprocess.Popen(["ffmpeg", "-loglevel", "error", "-y", "-f", "rawvideo", "-pix_fmt", "rgb24", "-s", "1920x1080", "-r", str(FPS),
                              "-i", "-", "-vf", "scale=out_color_matrix=bt709:out_range=tv,setparams=color_primaries=bt709:color_trc=iec61966-2-1",
                              "-c:v", "libx264", "-threads", "4", "-pix_fmt", "yuv420p", "-crf", "18", "-movflags", "+faststart", mp4],
                             stdin=subprocess.PIPE)
    sheet, last = Image.new("RGB", (1920, 1080), "white"), None
    for png in files:
        if png != last:
            sheet.paste(Image.open(png).convert("RGB"), (0, (1080 - height) // 2))
            data, last = sheet.tobytes(), png
        video.stdin.write(data)
    video.stdin.close()
    video.wait()
    print(f"{len(args) // 3} pictures in the WebP, {len(files) / FPS:.1f} s,",
          ", ".join(f"{os.path.getsize(p) // 1024} kB {p[-4:]}" for p in (webp, mp4)))
