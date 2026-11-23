# What the animations lattice and cluster_shapes share: renderings that turn are put on a page of labels that stay, and the frames are
# written as an animated WebP and as an MP4. Used by animate_lattice.py and animate_cluster_shapes.py.
# Needs pdflatex, gs, cwebp, webpmux, ffmpeg, numpy and Pillow.
import os
import subprocess

import numpy as np
from PIL import Image

HERE = os.path.dirname(os.path.abspath(__file__))
HEAD = "\\documentclass[11pt,tikz,border=0pt]{standalone}\n\\input{figstyle}\n\\begin{document}\n"


def pages(body, folder, cm):
    """the pages of a LaTeX file with this body as pictures on white, labels_001.png and so on in the folder; cm: pixels per centimeter"""
    tex = os.path.join(folder, "labels.tex")
    with open(tex, "w") as file:
        file.write(HEAD + body + "\\end{document}\n")
    subprocess.run(["pdflatex", "-interaction=batchmode", "-halt-on-error", f"-output-directory={folder}", tex], check=True,
                   stdout=subprocess.DEVNULL, env=dict(os.environ, TEXINPUTS=HERE + ":"))
    subprocess.run(["gs", "-q", "-dNOPAUSE", "-dBATCH", "-sDEVICE=png16m", f"-r{2 * 2.54 * cm}", "-dDownScaleFactor=2", "-dTextAlphaBits=4",
                    "-dGraphicsAlphaBits=4", "-sOutputFile=" + os.path.join(folder, "labels_%03d.png"), tex[:-3] + "pdf"], check=True)


def put(page, pictures):
    """the page with the pictures [(file, x, y of its center in pixels)] on it; what lies outside the page is cut off"""
    im = page.convert("RGBA")
    for path, x, y in pictures:
        picture = Image.open(path).convert("RGBA")
        left, top = x - picture.width // 2, y - picture.height // 2
        inside = (max(-left, 0), max(-top, 0), min(picture.width, im.width - left), min(picture.height, im.height - top))
        im.alpha_composite(picture.crop(inside), (left + inside[0], top + inside[1]))
    return im.convert("RGB")


def webp(frames, numbers, fps, size, quality, path):
    """animated WebP of the frames with these numbers (every n-th one of a film with fps frames per second) at this size. The first
    frame is stored in full, with the quality 90 for the labels. The others hold the blocks of 8 x 8 pixels in which something moves in
    the course of the turn, and the blocks around them unless a label reaches into them, so that their edges lie on white; the rest of
    them is transparent. In this way the labels are drawn once and cannot flicker."""
    files = [os.path.join(frames, f"small_{n:03d}") for n in numbers]
    for n, file in zip(numbers, files):
        Image.open(os.path.join(frames, f"{n:03d}.png")).resize(size, Image.LANCZOS).save(file + ".png")
    first = np.asarray(Image.open(files[0] + ".png"))
    moves = np.zeros(first.shape[:2], bool)
    for file in files[1:]:
        moves |= (np.asarray(Image.open(file + ".png")) != first).any(axis=2)
    h, w = moves.shape
    rows, columns = -(-h // 8), -(-w // 8)
    blocks = lambda pixels: np.pad(pixels, ((0, -h % 8), (0, -w % 8))).reshape(rows, 8, columns, 8).any(axis=(1, 3))
    around = np.pad(blocks(moves), 1)
    around = np.any([around[i:i + rows, j:j + columns] for i in range(3) for j in range(3)], axis=0)
    labels = blocks((first < 255).any(axis=2) & ~moves)    # what is drawn and never moves
    drawn = np.kron(blocks(moves) | around & ~labels, np.ones((8, 8), bool))[:h, :w]
    ys, xs = np.flatnonzero(drawn.any(axis=1)), np.flatnonzero(drawn.any(axis=0))
    top, left = ys[0], xs[0]                               # multiples of 8: the place of a frame has to be even
    part = np.s_[top:ys[-1] + 1, left:xs[-1] + 1]
    args = []
    for i, (n, file) in enumerate(zip(numbers, files)):
        if i:
            s = np.asarray(Image.open(file + ".png"))
            Image.fromarray(np.dstack([np.where(drawn[..., None], s, 255), 255 * drawn]).astype(np.uint8)[part]).save(file + ".png")
        subprocess.run(["cwebp", "-quiet", "-q", str(quality if i else 90), "-m", "6", "-sharp_yuv", "-exact", file + ".png", "-o", file + ".webp"],
                       check=True)
        ms = round(1000 * (n + numbers.step) / fps) - round(1000 * n / fps)
        args += ["-frame", file + ".webp", f"+{ms}+{left if i else 0}+{top if i else 0}+0+b"]
    subprocess.run(["webpmux", *args, "-loop", "0", "-o", path], check=True, capture_output=True)
    print(f"{path}: {len(files)} pictures, {os.path.getsize(path) / 1e6:.3f} MB")


def mp4(frames, count, fps, path):
    """MP4 of all frames, 1920 x 1080 with the frames in the middle"""
    # the colors of the frames, in the form that players expect at this size
    colors = "scale=out_color_matrix=bt709:out_range=tv,setparams=color_primaries=bt709:color_trc=iec61966-2-1"
    video = subprocess.Popen(["ffmpeg", "-loglevel", "error", "-y", "-f", "rawvideo", "-pix_fmt", "rgb24", "-s", "1920x1080", "-r", str(fps),
                              "-i", "-", "-vf", colors, "-c:v", "libx264", "-pix_fmt", "yuv420p", "-crf", "18", "-threads", "4",
                              "-movflags", "+faststart", path], stdin=subprocess.PIPE)
    for n in range(count):
        frame = Image.open(os.path.join(frames, f"{n:03d}.png"))
        screen = Image.new("RGB", (1920, 1080), "white")
        screen.paste(frame, ((1920 - frame.width) // 2, (1080 - frame.height) // 2))
        video.stdin.write(screen.tobytes())
    video.stdin.close()
    video.wait()
    print(f"{path}: {count} pictures, {count / fps:.1f} s, {os.path.getsize(path) / 1e6:.3f} MB")
