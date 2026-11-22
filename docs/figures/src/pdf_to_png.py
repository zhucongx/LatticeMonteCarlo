"""A still from its PDF, on white and with a palette of 256 colors.

    python3 pdf_to_png.py <in.pdf> <out.png> [width in pixels, default 1600]

In the choice of the palette a color counts with the number of its pixels to the power 0.4. A plain median cut spends the
palette on what covers most of the picture, and a single orange ball among many gray ones turns brown.
Needs Ghostscript, numpy and Pillow.
"""
import os
import subprocess
import sys
import tempfile

import numpy as np
from PIL import Image

src, dst = sys.argv[1:3]
width = int(sys.argv[3]) if len(sys.argv) > 3 else 1600
with tempfile.TemporaryDirectory() as folder:
    page = os.path.join(folder, "page.png")
    subprocess.run(["gs", "-q", "-dNOPAUSE", "-dBATCH", "-sDEVICE=png16m", "-r900", "-dTextAlphaBits=4",
                    "-dGraphicsAlphaBits=4", "-sOutputFile=" + page, src], check=True)
    im = Image.open(page).convert("RGB")
    im = im.resize((width, round(im.height * width / im.width)), Image.LANCZOS)

colors, inverse, counts = np.unique(np.asarray(im).reshape(-1, 3), axis=0, return_inverse=True, return_counts=True)
weights = counts ** 0.4
sample = np.repeat(colors, np.ceil(weights).astype(int), axis=0)
start = Image.fromarray(sample[None]).quantize(256, Image.Quantize.MEDIANCUT).getpalette()
palette = np.array(start, dtype=float).reshape(-1, 3)


def nearest():
    return np.concatenate([((c[:, None].astype(float) - palette) ** 2).sum(2).argmin(1) for c in np.array_split(colors, 20)])


for _ in range(10):
    index = nearest()
    for k in np.unique(index):
        palette[k] = np.average(colors[index == k], axis=0, weights=weights[index == k])
palette = palette.round()
palette[np.abs(palette - 255).sum(1).argmin()] = 255                # the background is exactly white
out = Image.fromarray(nearest()[inverse].reshape(im.height, im.width).astype(np.uint8))
out.putpalette(palette.astype(np.uint8).ravel().tolist())
out.save(dst, optimize=True)
print(dst, im.size, f"{os.path.getsize(dst) / 1000:.0f} kB")
