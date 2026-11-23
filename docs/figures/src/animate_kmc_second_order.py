# kmc_second_order.webp and kmc_second_order.mp4: one step of second-order KMC (KineticMcSecondOmpi). A small vacancy runs through the
# chart of the paths along two of them, one without and one with two returns, while panel (a) shows what these paths do to the vacancy;
# then the sum over all paths lights up, which is what the program computes (KineticMcSecondOmpi::CalculateTime). Drawn with TikZ.
#   python3 animate_kmc_second_order.py <folder for the frames> <output folder> [numbers of single frames: draw only these, for a look]
import sys
from collections import namedtuple

from tikz_film import draw, ease, write

folder, out = sys.argv[1:3]
# ------------------------------------------------------------------ the two paths
# A path: its name, and the arrows of the chart it takes, each with the site to which the vacancy jumps meanwhile (None: it stays).
# K is the site of the vacancy at the start, J the one it came from, I the neighbor site of the drawing.
PATHS = (("path 1: no return", (("AB", "I"), ("BC", None), ("CS", None))),
         ("path 2: two returns", (("AB", "I"), ("BD", "K"), ("DE", "J"), ("ED", "K"), ("DE", "I"), ("EF", None), ("FS", None))))
ARROWS = {"AB": ("(A.east)", "(B.west)"), "BD": ("(B.east)", "(D.west)"), "BC": ("(B.east)", "(C.west)"), "DE": ("(D.east)", "(E.west)"),
          "EF": ("(E.south)", "(F.north)"), "ED": ("(E.north)", "($(E.north)+(0,0.5)$)", "($(D.north)+(0,0.5)$)", "(D.north)"),
          "CS": ("(C.south)", "(C.south |- S.north)"), "FS": ("(F.south)", "(F.south |- S.north)")}
LENGTHS = {"ED": (0.5, 4.3, 0.5)}                          # of the parts of an arrow with corners, in cm


def frames(arrow, site):
    """frames in which the vacancy runs along an arrow of the chart, and for which it stays at its end"""
    return (9, 14) if arrow in "CS FS" else (15, 14) if site else (11, 20)


# ------------------------------------------------------------------ the film
Scene = namedtuple("Scene", "path n start run end")        # arrow n of a path: from `start` the vacancy runs, from `run` to `end` it stays
BEGIN, RESET = (6, 12), 8                                  # frames at the first box of each path, and between two paths
film, f = [], 0
for path, (_, arrows) in enumerate(PATHS):
    f += BEGIN[path]
    for n, (arrow, site) in enumerate(arrows):
        run, stay = frames(arrow, site)
        film.append(Scene(path, n, f, f + run, f + run + stay))
        f += run + stay
    f += RESET
last = f


def state(f):
    """at frame f: the path (None after the last one), the number of its arrows that are done, the part of the next one that is done
    (None if the vacancy stays), and how far the picture has gone back to the start after the path (0 to 1)"""
    for path in range(len(PATHS)):
        scenes = [s for s in film if s.path == path]
        if f < scenes[-1].end + RESET:
            under_way = next((s for s in scenes if s.start < f < s.run), None)
            return (path, sum(f >= s.run for s in scenes), (f - under_way.start) / (under_way.run - under_way.start) if under_way else None,
                    max(0, f - scenes[-1].end) / RESET)
    return None, 0, None, 0.0


# ------------------------------------------------------------------ the picture, lengths in cm
WIDTH = 453.81 / 28.4528                                   # of the picture with its border
X0, Y1 = -38.289 / 28.4528, 39.711 / 28.4528
BOX = f"({X0:.4f},{Y1 - WIDTH * 9 / 16:.4f}) rectangle ({X0 + WIDTH:.4f},{Y1:.4f})"
SITE = {"J": 0.0, "K": 2.6, "I": 7.5}                      # x of the sites in (a)
BOXES = {"A": ("box, text width=2.1cm", "0.1,0", r"vacancy\\ came from $j$"),
         "B": ("box, text width=2.55cm", "4.25,0", r"first attempt:\\ to a site $i\neq j$"),
         "D": ("backbox, text width=2.8cm", "8.7,1.0", r"return:\\ the step goes on"),
         "C": ("endbox, text width=2.8cm", "8.7,-1.0", r"no return:\\ the step ends at $i$"),
         "E": ("box, text width=2.8cm", "13.0,1.0", r"next attempt:\\ to any site $i$"),
         "F": ("endbox, text width=2.8cm", "13.0,-1.0", r"no return:\\ the step ends at $i$"),
         "S": ("sumbox, text width=15.35cm", "6.725,-2.65",
               r"sum over all paths: the probability that the step ends at $i$, and the mean time of the step")}
PREAMBLE = r"""\tikzset{
  st/.style={atom=#1, minimum size=15pt}, st/.default=cAl,
  box/.style={draw=black!70, rounded corners=3pt, align=center, inner sep=4pt, fill=black!5, minimum height=1.15cm, line width=0.6pt,
              outer sep=0.3pt},
  endbox/.style={box, fill=cThree!22},
  backbox/.style={box, fill=cOne!18},
  sumbox/.style={box, fill=cFour!26, minimum height=0.85cm},
  site/.style={font=\footnotesize, align=center, anchor=north, inner sep=0pt},
  pr/.style={font=\small, inner sep=3pt},
  here/.style={line width=1.6pt, draw=cDark},
  thick arr/.style={-{Stealth[length=7.5pt]}, line width=1.9pt},
  walker/.style={vac, minimum size=9.5pt, line width=1.1pt},
  txt/.style={anchor=base west, inner sep=0pt, font=\normalsize}}
% (a): the sites, whose places the arrows keep while the vacancy is away, and all that stays (also what stands at the arrows)
\newcommand{\sites}{%
  \node[circle, inner sep=0pt, minimum size=15pt] (j) at (0,0) {};
  \node[circle, inner sep=0pt, minimum size=14.5pt, line width=1.1pt] (k) at (2.6,0) {};
  \node[circle, inner sep=0pt, minimum size=15pt] (i) at (7.5,0) {};
  \node[st] (l1) at (9.4,0.95) {};
  \node[st] (l2) at (9.8,0) {};
  \node[st] (l3) at (9.4,-0.95) {};
  \node[text=cThree!60!black, anchor=west, align=left] at (10.3,0) {$1-q_i$: on to one of\\ 11 other sites};
  \path[line width=0.8pt] (k) to[bend left=22] node[above] {$p_i$: jump to $i$} (i);
  \path[line width=0.8pt] (i) to[bend left=22] node[below, text=cOne!85!black] {$q_i$: immediate return} (k);
  \node[site] at (0,-0.40) {\strut previous\\ \strut site $j$};
  \node[site] at (7.5,-0.40) {\strut neighbor\\ \strut site $i$};
  \node[lab] at (-1.75,1.35) {\panel{a}};}
% (b): the arrows of the chart, and what stands at them (drawn last, over the small vacancy)
\newcommand{\chart}{%
  \draw[arr] (A) -- (B);
  \draw[arr] (B.east) -- (D.west);
  \draw[arr] (B.east) -- (C.west);
  \draw[arr] (D) -- (E);
  \draw[arr] (E) -- (F);
  \draw[arr, cOne] (E.north) -- ++(0,0.5) -| (D.north);
  \draw[arr] (C) -- (C.south |- S.north);
  \draw[arr] (F) -- (F.south |- S.north);
  \node[lab] (b) at (-1.3,2.5) {\panel{b}};}
\newcommand{\chartlabels}{%
  \begin{scope}[line width=0.8pt]                          % that of the arrows, which sets the distance of a label
    \path (A) -- node[pr, above] {$\dfrac{p_i}{1-p_j}$} (B);
    \path (B.east) -- node[pr, above, sloped] {$q_i$} (D.west);
    \path (B.east) -- node[pr, below, sloped] {$1-q_i$} (C.west);
    \path (D) -- node[pr, above] {$p_i$} (E);
    \path (E) -- node[pr, right] {$1-q_i$} (F);
    \path (E.north) -- ++(0,0.5) -| node[pos=0.25, above, text=cOne!85!black, inner sep=2.5pt] {$q_i$: return again} (D.north);
  \end{scope}}
"""


def arrow(style, width, path):
    """an arrow of (a) that is thick while the vacancy takes it"""
    return f"\\draw[{'thick arr' if width else 'arr'}, {style}] {path};"


def along(arrow, x):
    """the point at the part x of an arrow of the chart"""
    points, lengths = ARROWS[arrow], LENGTHS.get(arrow, (1,))
    x *= sum(lengths)
    for a, b, length in zip(points, points[1:], lengths):
        if x <= length:
            break
        x -= length
    return f"($ {a}!{x / length:.3f}!{b} $)"


def page(f):
    path, done, part, back = state(f)
    name, arrows = PATHS[path] if path is not None else ("", ())

    def site(n):
        """where the vacancy is after n arrows of the path"""
        return ([s for _, s in arrows[:n] if s] or ["K"])[-1]

    jump = arrows[done] if part is not None else None      # the arrow of the chart that is under way
    now = jump or (arrows[done - 1] if done else None)     # and the one that the picture shows
    before = site(done if jump else done - 1)
    lit = now is not None and back == 0
    code = [f"\\begin{{tikzpicture}}\n\\path[use as bounding box] {BOX};"]
    # ---------------- (a) sites and probabilities; an arrow is thick while the path is on it
    code += ["\\begin{scope}[shift={(0.45,0)}]\\sites",
             arrow("dashed, black!60", path is not None and now is None, "(j) -- (k)"),
             arrow("black", lit and now[1] == "I", "(k) to[bend left=22] (i)"),
             arrow("cOne", lit and now[1] == "K" and before == "I", "(i) to[bend left=22] (k)")]
    code += [arrow("cThree!80!black", lit and now[1] is None, f"(i) -- (l{n})") for n in (1, 2, 3)]
    # the vacancy and the atom with which it changes places; after a path all go back to the start
    v = SITE[site(done)]
    if jump and jump[1]:
        v += (SITE[jump[1]] - v) * ease(part)
    al = SITE["J"] + SITE["K"] - v if v < SITE["K"] else SITE["J"]
    mg = SITE["I"] + SITE["K"] - v if v > SITE["K"] else SITE["I"]
    word = max(0.0, 1 - abs(v - SITE["K"])) * (1 - back) + back            # the word below the vacancy belongs to its site
    code.append(f"\\node[site, opacity={word:.2f}] at (2.6,-0.40) {{\\strut vacancy}};")
    for style, x, home in (("vac, minimum size=14.5pt", v, "K"), ("st", al, "J"), ("st=cMg", mg, "I")):
        for x, opacity in ((x, 1 - back), (SITE[home], back)) if back and x != SITE[home] else ((x, 1),):
            code.append(f"\\node[{style}, opacity={opacity:.2f}] at ({x:.3f},0) {{}};")
    code.append("\\end{scope}")
    # ---------------- (b) the chart, the boxes and arrows of the path, and the small vacancy on its way
    whole = path is None                                   # the end: the sum over all paths
    visited = {"A"} | {a[1] for a, _ in arrows[:done]} if lit or now is None and not whole else {"S"} if whole else set()
    code.append("\\begin{scope}[shift={(0,-4.3)}]")
    code += [f"\\node[{style}{', here' if b in visited else ''}{', fill=cFour!55' if whole and b == 'S' else ''}] ({b}) at ({where}) {{{text}}};"
             for b, (style, where, text) in BOXES.items()]
    code.append("\\chart")
    for a in ("CS", "FS") if whole else dict.fromkeys(a for a, _ in arrows[:done]) if back == 0 else ():
        code.append(f"\\draw[thick arr{', cOne' if a == 'ED' else ''}] {' -- '.join(ARROWS[a])};")
    if path is not None and back == 0:
        code.append(f"\\node[walker] at {along(jump[0], ease(part)) if jump else along(now[0], 1) if now else '(A.east)'} {{}};")
    code.append("\\chartlabels")
    title = None if whole else name if back == 0 else PATHS[path + 1][0] if path + 1 < len(PATHS) else None
    if title:
        code.append(f"\\node[txt] at ($(b.base east)+(0.22,0)$) {{{title}}};")
    else:
        code += ["\\node[txt] at ($(b.base east)+(0.22,0)$) {the program does not walk these paths:};",
                 "\\node[txt] at ($(b.base east)+(0.22,-0.5)$) {it sums over all of them};"]
    code += ["\\end{scope}", "\\end{tikzpicture}\n"]
    return "\n".join(code)


wanted = [int(a) for a in sys.argv[3:]] or range(last + 1)
files = draw([page(f) for f in wanted], PREAMBLE, folder, WIDTH)
if sys.argv[3:]:
    print("\n".join(f"{f} {png}" for f, png in zip(wanted, files)))
else:
    write(files, out, "kmc_second_order")
