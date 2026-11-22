# PDF of one page -> PNG on white with a palette of 256 colors, for figures with renderings and color bars. The palette is fitted to
# the picture, so that the shades of a color that covers a small area are kept as well; smooth shades between two colors of the palette
# are dithered, edges and text are not. Needs Ghostscript, numpy, Pillow and scikit-learn.
#   python3 pdf_to_png_dithered.py <in.pdf> <out.png> [width in pixels, default 1600]
import os
import subprocess
import sys
import tempfile

import numpy as np
from PIL import Image, ImageFilter
from sklearn.cluster import KMeans

src, dst = sys.argv[1:3]
width = int(sys.argv[3]) if len(sys.argv) > 3 else 1600
with tempfile.TemporaryDirectory() as folder:
    page = os.path.join(folder, "page.png")
    subprocess.run(["gs", "-q", "-dNOPAUSE", "-dBATCH", "-sDEVICE=png16m", "-r900", "-dTextAlphaBits=4", "-dGraphicsAlphaBits=4",
                    f"-sOutputFile={page}", src], check=True)
    im = Image.open(page).convert("RGB")
im = im.resize((width, round(im.height * width / im.width)), Image.LANCZOS)


def lab(rgb):
    """CIELAB, in which colors at the same distance look equally different"""
    c = np.asarray(rgb) / 255
    c = np.where(c > 0.04045, ((c + 0.055) / 1.055) ** 2.4, c / 12.92)
    c = c @ np.array([[0.4124, 0.2126, 0.0193], [0.3576, 0.7152, 0.1192], [0.1805, 0.0722, 0.9505]]) / [0.9505, 1.0, 1.089]
    c = np.where(c > 0.008856, np.cbrt(c), 7.787 * c + 16 / 116)
    return np.stack([116 * c[:, 1] - 16, 500 * (c[:, 0] - c[:, 1]), 200 * (c[:, 1] - c[:, 2])], axis=1)


pixels = np.asarray(im)
colors, number, count = np.unique(pixels.reshape(-1, 3), axis=0, return_inverse=True, return_counts=True)
number = number.reshape(pixels.shape[:2])
weight = np.sqrt(count)
group = KMeans(255, n_init=1, max_iter=40, random_state=0).fit(lab(colors), sample_weight=weight).labels_
mean = np.stack([np.bincount(group, weight * colors[:, k]) for k in range(3)], axis=1) / np.bincount(group, weight)[:, None]
palette = Image.new("P", (1, 1))
palette.putpalette(np.vstack([[255, 255, 255], np.round(mean)]).astype(np.uint8).tobytes())
index = 1 + group[number]                                                         # the nearest color
off = np.linalg.norm(lab(colors) - lab(mean)[group], axis=1)[number] > 1
smooth = (np.asarray(im.filter(ImageFilter.MaxFilter(3))).astype(int) - np.asarray(im.filter(ImageFilter.MinFilter(3)))).max(axis=2) < 24
index = np.where(off & smooth, np.asarray(im.quantize(palette=palette, dither=Image.Dither.FLOYDSTEINBERG)), index)
index[(pixels == 255).all(axis=2)] = 0                                            # the background stays white
picture = Image.fromarray(index.astype(np.uint8))
picture.putpalette(palette.getpalette())
picture.save(dst, optimize=True)
print(dst, im.size, f"{os.path.getsize(dst) / 1000:.0f} kB")
