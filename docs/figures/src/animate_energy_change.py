# energy_change.webp and energy_change.mp4: the clusters that change when two sites exchange their atoms, first for two distant sites,
# then for a jump of the vacancy, drawn on a triangular lattice. No run is needed. It follows
# EnergyChangePredictorPairSite::GetDeFromLatticeIdPair (lmc/pred/src/EnergyChangePredictorPairSite.cpp). Every frame is drawn with
# TikZ. Make both files in this folder with, for any scratch folder T:
#   python3 animate_energy_change.py T/frames <output folder> [numbers of single frames: draw only these, for a look]
# Needs what tikz_frames.py needs.
import math
import sys

from tikz_frames import ease, frames, make, opacity, ramp

folder, out = sys.argv[1:3]
HEIGHT = 624                                               # of a frame in pixels
H = math.sqrt(3) / 2                                       # distance of two rows, in units of the distance of two sites
SHIFT = 16.6                                               # from (a) to (b)
REACH = 2.32                                               # a circle of this radius holds the first three neighbor shells (1, 1.73 and 2)
I, J = ((3.5, 3), (3.5, 3)), ((10.5, 3), (4.5, 3))         # the two sites in (a) and in (b): x and row
# what i and j hold at the start; the exchange leaves i with Mg and j with Zn in (a), i with the vacancy and j with Zn in (b)
START = (("zn", "mg"), ("zn", "va"))
# examples of clusters, as x and row of their sites: with exactly one of the two sites, and with none or both of them
CHANGED = (([(3.5, 3), (2.5, 3), (3, 4)], [(10.5, 3), (11, 4)]), ([(3.5, 3), (2.5, 3), (3, 4)], [(4.5, 3), (5, 4)]))
SAME = (([(6.5, 5), (7.5, 5), (7, 4)], [(11.5, 1), (12, 2)]), ([(3.5, 3), (4.5, 3), (4, 2)],))
for p in (0, 1):
    assert all((I[p] in c) + (J[p] in c) == 1 for c in CHANGED[p]) and all((I[p] in c) + (J[p] in c) != 1 for c in SAME[p])

# ------------------------------------------------------------------ the film: frames at which its parts begin
FADE = frames(0.4)
starts, t = [], frames(0.5)
for move in (frames(0.9), frames(0.5)):                    # the exchange over a distance takes longer than the jump
    circles, t = t, t + frames(1.0)
    arrow, t = t, t + frames(0.4)
    leave, t = t, t + move + frames(0.3)
    changed, t = t, t + frames(1.2)
    same, t = t, t + frames(1.2)
    starts.append((circles, arrow, leave, move, changed, same))
LAST = same + FADE + frames(2.5)                           # the end is held


# ------------------------------------------------------------------ the picture: the lattices in units of the distance of two sites
def at(x, row, panel=0):
    return f"({x + SHIFT * panel:.3f},{H * row:.3f})"


def shape(style, sites, panel):
    return f"\\draw[{style}] " + " -- ".join(at(x, row, panel) for x, row in sites) + (" -- cycle;" if len(sites) > 2 else ";")


def arc(p, q, part, side):
    """the point at a part of the way from p to q on a curve that leaves p at 26 degrees to one side: that of the arrow, or its mirror image"""
    d, c, s = 0.3915 * (q[0] - p[0]), math.cos(math.radians(26)), math.sin(math.radians(26)) * side
    curve = [p, (p[0] + d * c, p[1] + abs(d) * s / H), (q[0] - d * c, q[1] + abs(d) * s / H), q]
    weights = [(1 - part) ** 3, 3 * part * (1 - part) ** 2, 3 * part ** 2 * (1 - part), part ** 3]
    return [sum(w * point[k] for w, point in zip(weights, curve)) for k in (0, 1)]


def dots(last_even, last_odd, panel):
    sites = [(x + 0.5 * (row % 2), row) for row in range(7) for x in range((last_odd if row % 2 else last_even) + 1)]
    return "".join(f"\\node[dot] at {at(x, row, panel)} {{}};" for x, row in sites if (x, row) not in (I[panel], J[panel]))


PREAMBLE = r"""
\tikzset{dot/.style={atom=cAl!80, minimum size=7.5pt}, big/.style={minimum size=14pt, font=\small},
  mg/.style={atom=cMg, big}, zn/.style={atom=cZn, big}, va/.style={vac, big, minimum size=13pt},
  on mg/.style={font=\small, inner sep=0pt}, on va/.style={on mg}, on zn/.style={on mg, text=white},
  reach/.style={fill=cFour!28, draw=cFour!80!black, line width=0.5pt},
  changed/.style={draw=cThree, fill=cThree!38, line width=2.2pt},
  same/.style={draw=cDark!75, fill=cDark!22, line width=2.2pt}}
\pgfmathsetmacro{\h}{sqrt(3)/2}
\newcommand{\lattices}{%
DOTS}
\newcommand{\fixed}{%
  \node[above lab] at (-0.35,6*\h+0.5) {\panel{a}}; \node[above lab] at (SHIFT-0.35,6*\h+0.5) {\panel{b}};
  \begin{scope}[shift={(-0.35,-1.25)}, x=1cm, y=1cm]
    \specieslegend{0.25}{0}
    \vacancylegend{3.85}{0}
    \draw[reach] (6.15,0) circle (0.2); \node[leg] at (6.45,0) {up to third neighbors};
    \draw[changed] (10.45,-0.14) -- ++(0.36,0) -- ++(-0.18,0.31) -- cycle; \node[leg] at (10.95,0) {changed};
    \draw[same] (13.05,-0.14) -- ++(0.36,0) -- ++(-0.18,0.31) -- cycle; \node[leg] at (13.55,0) {unchanged};
  \end{scope}
  % the places of what the two sites hold at the end, for the ends of the arrows
  \node[big, circle, inner sep=0pt] (ia) at (3.5,3*\h) {}; \node[big, circle, inner sep=0pt] (ja) at (10.5,3*\h) {};
  \node[va, draw opacity=0, fill opacity=0] (ib) at (SHIFT+3.5,3*\h) {}; \node[big, circle, inner sep=0pt] (jb) at (SHIFT+4.5,3*\h) {};}
""".replace("DOTS", dots(14, 13, 0) + "%\n" + dots(8, 7, 1)).replace("SHIFT", str(SHIFT))


def picture(f):
    under, atoms, over = [], [], []
    for panel, (circles, arrow, leave, move, changed, same) in enumerate(starts):
        i, j = I[panel], J[panel]
        if f > circles:                                    # the sites up to the third neighbors of each of the two
            both = f"{at(*i, panel)} circle ({REACH}) {at(*j, panel)} circle ({REACH});"
            under.append(f"\\begin{{scope}}[transparency group, {opacity(ramp(f, circles, FADE))}] \\fill[reach] {both}"
                         f" \\draw[reach, fill=none] {both} \\end{{scope}}")
        for style, begin, shapes in (("changed", changed, CHANGED[panel]), ("same", same, SAME[panel])):
            if f > begin:
                under.append(f"\\begin{{scope}}[{opacity(ramp(f, begin, FADE))}] " + " ".join(shape(style, s, panel) for s in shapes)
                             + " \\end{scope}")
        part = ease(ramp(f, leave, move))                  # the two sites exchange what they hold
        first, second = START[panel]
        if panel == 0:                                     # two atoms on two curves: that of the arrow and its mirror image
            atoms += [f"\\node[{first}] at {at(*arc(i, j, part, -1), panel)} {{}};", f"\\node[{second}] at {at(*arc(j, i, part, 1), panel)} {{}};"]
        else:                                              # the atom covers the old place of the vacancy and uncovers the new one
            atoms += [f"\\node[va, {opacity(1 - part)}] at {at(*j, panel)} {{}};", f"\\node[va, {opacity(part)}] at {at(*i, panel)} {{}};",
                      f"\\node[{first}] at {at(i[0] + (j[0] - i[0]) * part, i[1], panel)} {{}};"]
        # the names of the sites stay where they are; a name is white on a Zn atom
        on_i, on_j = (first, second) if part < 0.12 else (second, first) if part > 0.88 else ("va", "va")
        atoms += [f"\\node[on {on_i}] at {at(*i, panel)} {{$i$}};", f"\\node[on {on_j}] at {at(*j, panel)} {{$j$}};"]
        if f > arrow:
            shown = opacity(ramp(f, arrow, frames(0.3)))
            over.append(f"\\draw[dbl, {shown}] (ia) to[bend right=26] (ja);"
                        f" \\node[fill=white, inner xsep=1.5pt, inner ysep=2.5pt, {shown}] at (7,\\h) {{exchange}};" if panel == 0 else
                        f"\\draw[{{Stealth[length=3.4pt]}}-{{Stealth[length=3.4pt]}}, line width=0.8pt, {shown}]"
                        " (ib.76) to[out=58, in=122, looseness=1.0] (jb.104);")
    return "\n".join(["\\begin{scope}[shift={(0.45,1.18)}, x=0.62cm, y=0.62cm]", "\\fixed"] + under + ["\\lattices"] + atoms + over
                     + ["\\end{scope}", ""])


make("energy_change", PREAMBLE, [picture(f) for f in range(LAST + 1)], HEIGHT, folder, out, [int(a) for a in sys.argv[3:]])
