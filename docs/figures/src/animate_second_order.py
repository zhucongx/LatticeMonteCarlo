# second_order.webp and second_order.mp4: the same flickering vacancy in first-order and in second-order KMC, side by side.
# The two runs are those of second_order_model.py with fixed seeds; every frame is drawn with TikZ.
#   python3 animate_second_order.py <folder for the frames> <output folder> [numbers of single frames: draw only these, for a look]
# Needs pdflatex, gs, img2webp, ffmpeg and Pillow. A frame that is already in its folder is not drawn again.
import hashlib
import math
import os
import random
import subprocess
import sys
from collections import namedtuple
from concurrent.futures import ThreadPoolExecutor

from PIL import Image

from second_order_model import I, J, K, L, M, MG, exchange, first_order, neighbors, rates, second_order, select

folder, out = sys.argv[1:3]
SEED_FIRST, SEED_SECOND = 22820, 41                        # chosen so that both examples take the way K, I, L, M; they end there
FPS, HOLD = 30, 2.0                                        # frames per second (the WebP takes every second frame), end in seconds

# ------------------------------------------------------------------ the two runs
Step = namedtuple("Step", "time start end mg probabilities")               # mg: the Mg atoms before the step
first, second = runs = ([], [])
rng, k, mg, t = random.Random(SEED_FIRST), K, MG, 0.0
while k != M:
    r = rates(k, mg)
    dt, i = first_order(k, mg, rng)
    t += dt
    first.append(Step(t, k, i, mg, [x / sum(r) for x in r]))
    k, mg = i, exchange(k, i, mg)
rng, j, k, mg, t = random.Random(SEED_SECOND), J, K, MG, 0.0
while k != M:
    probabilities, t_2 = second_order(k, j, mg)
    t += t_2
    i = select(neighbors(k), probabilities, rng)
    second.append(Step(t, k, i, mg, probabilities))
    j, k, mg = k, i, exchange(k, i, mg)
flicker = next(n for n, s in enumerate(first) if s.end not in (K, I))      # number of jumps between K and I

# ------------------------------------------------------------------ the film: both panels run through the same time
Scene = namedtuple("Scene", "panel n a b c t0 t1 arrows")  # step n: from frame a to b the time runs up to it, from b to c the atom jumps


def length(panel, n):
    """frames of the two parts of a scene, and whether the arrows are drawn"""
    if panel == 1:
        return 12, 12, True
    if n >= flicker:
        return 14 if n == flicker else 10, 10, True
    return (10, 6, True) if n < 2 else (4, 4, False) if n < 6 else (2, 4, False)


film, f, t = [], 20, 0.0
for time, panel, n in sorted((s.time, panel, n) for panel in (0, 1) for n, s in enumerate(runs[panel])):
    run_up, jump, arrows = length(panel, n)
    film.append(Scene(panel, n, f, f + run_up, f + run_up + jump, t, time, arrows))
    f, t = f + run_up + jump, time
last = f


def ease(x):
    x = min(max(x, 0.0), 1.0)
    return x * x * (3 - 2 * x)


def state(panel, f):
    """at frame f: number of steps done, part of the jump under way or None, time, frame at which the last jump ended"""
    done, since, time = 0, -99, 0.0
    for s in film:
        if f < s.a:
            break
        if f < s.c and s.panel == panel:
            return done, ease((f - s.b) / (s.c - s.b)) if f > s.b else None, s.t0 + (s.t1 - s.t0) * min(1, (f - s.a) / (s.b - s.a)), since
        if f < s.c:
            time = s.t0 + (s.t1 - s.t0) * (f - s.a) / (s.c - s.a)
            break
        time = s.t1
        if s.panel == panel:
            done, since = s.n + 1, s.c
    return done, None, min(time, runs[panel][-1].time), since


# ------------------------------------------------------------------ the picture, lengths in cm
A = 1.2                                                    # distance of two sites
SITES = [(a, b) for b in (-1, 0, 1, 2) for a in range(-2 - b // 2, 4 - b // 2)]
NAME = {J: "j", K: "k", I: "i", L: "l", M: "m"}
LEVEL = {K: 0.28, I: 0.80, L: 1.32, M: 1.84}               # heights of the sites in the lower panels
TIME = 12.8                                                # length of a microsecond there
PREAMBLE = r"""\documentclass[11pt,tikz,border=0pt]{standalone}
\input{figstyle}
\tikzset{al/.style={atom, minimum size=20pt}, mg/.style={atom=cMg, minimum size=22pt},
  va/.style={vac, minimum size=18.4pt, line width=1.6pt}, txt/.style={inner sep=0pt, font=\normalsize}}
\newcommand{\axes}{%
  \foreach \y/\t in {LEVELS} {
    \draw[black!14, line width=0.5pt] (0,\y) -- (6.4,\y);
    \draw[line width=0.6pt] (0,\y) -- (-0.08,\y) node[anchor=mid east, inner sep=3pt] {$\t$};}
  \draw[line width=0.6pt] (0,2.15) -- (0,0) -- (6.4,0);
  \foreach \t in {0, 0.1, 0.2, 0.3, 0.4, 0.5} {
    \draw[line width=0.6pt] ({SCALE*\t},0) -- ++(0,-0.08) node[font=\footnotesize, anchor=north, inner sep=2.5pt] {\t};}
  \node[font=\footnotesize, anchor=north, inner sep=0pt] at (3.2,-0.5) {time (\textmu s)};
  \node[font=\footnotesize, rotate=90, inner sep=0pt] at (-0.62,1.06) {site};}
\newcommand{\legend}{%
  \node[al, minimum size=10pt] at (12.05,8.65) {};
  \node[mg, minimum size=11pt] at (13.20,8.65) {};
  \node[va, minimum size=9.5pt, line width=1.1pt] at (14.45,8.65) {};
  \foreach \x/\t in {12.05/Al, 13.20/Mg, 14.45/vacancy} \node[leg, anchor=mid west, yshift=3.2pt] at (\x+0.2,8.65) {\t};}
\begin{document}
""".replace("LEVELS", ", ".join(f"{y:.2f}/{NAME[s]}" for s, y in LEVEL.items())).replace("SCALE", str(TIME))


def place(site):
    return A * (site[0] + site[1] / 2), A * site[1] * math.sqrt(3) / 2


def node(style, x, y, opacity=1.0):
    return f"\\node[{style}{f', opacity={opacity:.2f}' if opacity < 0.995 else ''}] at ({x:.3f},{y:.3f}) {{}};"


def panel_code(panel, f):
    run = runs[panel]
    done, jump, now, since = state(panel, f)
    k, mg, to = (run[done - 1].end, exchange(run[done - 1].start, run[done - 1].end, run[done - 1].mg), None) if done else (K, MG, None)
    if jump is not None:
        to = run[done].end
    code = [f"\\begin{{scope}}[xshift={8.4 * panel}cm]",
            f"\\node[txt, anchor=base west] at (0.1,8.55) {{\\panel{{{'ab'[panel]}}} {('first', 'second')[panel]} order}};",
            "\\begin{scope}[shift={(2.9,5.75)}]"]
    # the atoms and the vacancy; the atom that jumps covers the old place of the vacancy and uncovers the new one
    code += [node("mg" if s in mg else "al", *place(s)) for s in SITES if s not in (k, to)]
    (x0, y0), (x1, y1) = place(k), place(to or k)
    if to:
        code += [node("va", x0, y0, 1 - jump), node("va", x1, y1, jump),
                 node("mg" if to in mg else "al", x1 + (x0 - x1) * jump, y1 + (y0 - y1) * jump)]
    else:
        code.append(node("va", x0, y0))
    # arrows: the probabilities of the six sites; the one that was drawn turns red
    scene = next((s for s in film if (s.panel, s.n) == (panel, done)), None)
    if scene and scene.arrows and not to and f > since:
        for i, p in zip(neighbors(k), run[done].probabilities):
            x, y = place(i)
            a, b = (f"({x0 + (x - x0) * d / A:.3f},{y0 + (y - y0) * d / A:.3f})" for d in (0.42, 0.82 if i in mg else 0.85))
            code.append(f"\\draw[-{{Stealth[length={4 + 4.5 * p:.2f}pt, width={4 + 10 * p:.2f}pt]}}, line width={0.5 + 6.5 * p:.2f}pt, line cap=butt, "
                        f"{'cOne' if i == run[done].end and f > scene.b - 7 else 'black!70'}, opacity={ease((f - since) / 5):.2f}] {a} -- {b};")
    # the names of the sites: below those of the lower row and above those of the upper one, where no atom passes over them
    code += [f"\\node[txt, font=\\small, anchor=base] at ({place(s)[0]:.3f},{place(s)[1] + (0.44 if s[1] else -0.69):.3f}) {{${NAME[s]}$}};" for s in NAME]
    code.append("\\end{scope}")
    # counter and clock
    clock = run[done - 1].time if done else 0.0
    code += [f"\\node[txt, anchor=base west] at (0.85,3.62) {{steps: \\textbf{{{done}}}}};",
             f"\\node[txt, anchor=base east] at (7.25,3.62) {{{('', 'mean ')[panel]}time: \\textbf{{{clock:.3f}}} \\textmu s}};"]
    # the site of the vacancy against time
    way = [(0.0, LEVEL[K])] + [(s.time, LEVEL[site]) for s in run[:done] for site in (s.start, s.end)] + [(now, LEVEL[k])]
    if to:
        way.append((now, LEVEL[k] + (LEVEL[to] - LEVEL[k]) * jump))
    way = [f"({TIME * x:.3f},{y:.3f})" for x, y in way]
    code += ["\\begin{scope}[shift={(0.85,0.95)}]\\axes",
             f"\\draw[cOne, line width=1.0pt, line cap=butt, line join=miter] {' -- '.join(way)};",
             f"\\node[vac, minimum size=5pt, line width=1.0pt] at {way[-1]} {{}};"]
    if done > (flicker, 0)[panel]:                         # the steps that the flicker takes
        n, word = ((flicker, "steps"), (1, "step"))[panel]
        code.append(f"\\draw[decorate, decoration={{brace, amplitude=3.5pt}}, line width=0.6pt] (0,0.90) -- ({TIME * run[n - 1].time - 0.04:.3f},0.90) "
                    f"node[midway, above=3.5pt, font=\\footnotesize, fill=white, inner sep=1.2pt] {{{n} {word}}};")
    return "\n".join(code + ["\\end{scope}", "\\end{scope}", ""])


def page(f):
    return (f"\\begin{{tikzpicture}}\n\\path[use as bounding box] (0,0) rectangle (16,9);\n{panel_code(0, f)}{panel_code(1, f)}"
            "\\legend\n\\end{tikzpicture}\n")


# ------------------------------------------------------------------ the frames
os.makedirs(folder, exist_ok=True)
wanted = [int(a) for a in sys.argv[3:]] or range(last + 1)
pages = [page(f) for f in wanted]
files = [os.path.join(folder, hashlib.md5((PREAMBLE + p).encode()).hexdigest()[:12] + ".png") for p in pages]
todo = {png: p for png, p in zip(files, pages) if not os.path.exists(png)}
print(f"{flicker} jumps between k and i, {len(first)} and {len(second)} steps; frames 0 to {last}, {len(todo)} to draw", flush=True)
if todo:
    tex = os.path.join(folder, "frames.tex")
    with open(tex, "w") as file:
        file.write(PREAMBLE + "".join(todo.values()) + "\\end{document}\n")
    subprocess.run(["pdflatex", "-interaction=batchmode", "-halt-on-error", f"-output-directory={folder}", tex], check=True,
                   stdout=subprocess.DEVNULL, env=dict(os.environ, TEXINPUTS=os.path.dirname(os.path.abspath(__file__)) + ":"))

    def draw(job):
        number, png = job                                  # 1920 x 1080 pixels, three times as many inside
        subprocess.run(["gs", "-q", "-dNOPAUSE", "-dBATCH", "-sDEVICE=png16m", "-r914.4", "-dDownScaleFactor=3", f"-dFirstPage={number}",
                        f"-dLastPage={number}", f"-sOutputFile={png}.part", tex[:-3] + "pdf"], check=True)
        os.replace(png + ".part", png)

    with ThreadPoolExecutor(4) as pool:
        list(pool.map(draw, enumerate(todo, 1)))
if sys.argv[3:]:
    print("\n".join(f"{f} {png}" for f, png in zip(wanted, files)))
    sys.exit()

# ------------------------------------------------------------------ the two files
args, ms = [], [round(1000 * f / FPS) for f in range(0, last + 3, 2)]
for n, png in enumerate(files[::2]):
    d = ms[n + 1] - ms[n] + (round(1000 * HOLD) if 2 * n == last else 0)
    small = png[:-4] + "_small.png"
    if args and args[-1] == small:
        args[-2] = str(int(args[-2]) + d)                  # the same picture again
        continue
    if not os.path.exists(small):
        Image.open(png).convert("RGB").resize((760, 428), Image.LANCZOS).save(small)
    args += ["-d", str(d), small]
webp, mp4 = (os.path.join(out, "second_order." + e) for e in ("webp", "mp4"))
subprocess.run(["img2webp", "-loop", "0", "-sharp_yuv", "-lossy", "-q", "85", "-m", "6", *args, "-o", webp], check=True)
video = subprocess.Popen(["ffmpeg", "-loglevel", "error", "-y", "-f", "rawvideo", "-pix_fmt", "rgb24", "-s", "1920x1080", "-r", str(FPS), "-i", "-",
                          "-vf", "scale=out_color_matrix=bt709:out_range=tv,setparams=color_primaries=bt709:color_trc=iec61966-2-1",
                          "-c:v", "libx264", "-pix_fmt", "yuv420p", "-crf", "18", "-movflags", "+faststart", mp4], stdin=subprocess.PIPE)
for png in files + files[-1:] * round(FPS * HOLD):
    video.stdin.write(Image.open(png).convert("RGB").tobytes())
video.stdin.close()
video.wait()
print(f"{len(args) // 3} pictures, {last / FPS + HOLD:.1f} s,", ", ".join(f"{os.path.getsize(p) // 1024} kB {p[-4:]}" for p in (webp, mp4)))
