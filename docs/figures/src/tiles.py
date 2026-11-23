# The six small animations in the table of README.md: the part with the atoms of six animations of the documentation.
#   python3 tiles.py <folder docs/figures>
# Reads <name>.webp there and writes tile_<name>.webp. Needs Pillow and webpmux.
import os
import re
import subprocess
import sys

from PIL import Image

WIDTH, HEIGHT = 380, 364                                   # of the white page that the part is put on, in pixels of the animation
SHOWN = 320                                                # width of a tile in pixels
TILES = {                                                  # name: the box that is cut out, and every which picture is kept
    "quench": ((12, 30, 392, 394), 2),
    "vacancy_jumps": ((4, 27, 388, 395), 2),
    "trapped_vacancy": ((12, 30, 392, 394), 2),
    "frame_fields": ((20, 6, 372, 344), 2),
    "annealing": ((4, 26, 392, 404), 2),
    "jump_path": ((10, 30, 390, 342), 1),
}
WHITE = {"frame_fields": (0, 0, 44, 30)}                   # a panel label that is painted over


def durations(path):
    """ the time in ms for which every picture of an animated WebP file is shown """
    info = subprocess.run(["webpmux", "-info", path], check=True, capture_output=True, text=True).stdout
    return [int(line.split()[6]) for line in info.split("\n") if re.match(r"\s*\d+:", line)]


def tile(folder, name, box, step):
    source = os.path.join(folder, name + ".webp")
    shown, im, pictures, times = durations(source), Image.open(source), [], []
    assert len(shown) == im.n_frames
    for first in range(0, im.n_frames, step):
        im.seek(first)
        picture = im.convert("RGB")
        if name in WHITE:
            picture.paste("white", WHITE[name])
        cut = picture.crop(box)
        cut.thumbnail((WIDTH, HEIGHT), Image.LANCZOS)      # a larger part is scaled down to fit
        page = Image.new("RGB", (WIDTH, HEIGHT), "white")
        page.paste(cut, ((WIDTH - cut.width) // 2, (HEIGHT - cut.height) // 2))
        pictures.append(page.resize((SHOWN, round(SHOWN * HEIGHT / WIDTH)), Image.LANCZOS))
        times.append(sum(shown[first:first + step]))
    output = os.path.join(folder, f"tile_{name}.webp")
    pictures[0].save(output, save_all=True, append_images=pictures[1:], duration=times, loop=0, quality=72, method=6, minimize_size=True)
    print(f"tile_{name}.webp: {len(pictures)} pictures, {sum(times) / 1000:.1f} s, {os.path.getsize(output) / 1000:.0f} kB")


if __name__ == "__main__":
    for name, (box, step) in TILES.items():
        tile(sys.argv[1], name, box, step)
