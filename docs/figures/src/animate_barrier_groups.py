# barrier_groups.webp and barrier_groups.mp4: the four groups of jumps around a cluster, shown by a vacancy that makes two jumps of
# each group, drawn on a triangular lattice. No run is needed: the cluster and the jumps come from analysis_schematics.py, and their
# groups are checked by the rules of ExitTime::GetAverageBarriers (lmc/ansys/src/ExitTime.cpp). As that function does, the film puts
# the vacancy on the first site of a jump in place of the atom there, and puts the atoms back afterwards. Every frame is drawn with
# TikZ. Make both files in this folder with, for any scratch folder T:
#   python3 animate_barrier_groups.py T/frames <output folder> [numbers of single frames: draw only these, for a look]
# Needs what tikz_frames.py needs.
import sys

from analysis_schematics import JUMPS, barrier_groups, element, find_clusters, group_of_jump, neighbors
from tikz_frames import ease, frames, make, opacity, ramp

folder, out = sys.argv[1:3]
HEIGHT = 768                                               # of a frame in pixels
RULES = {"in": ("cFour", "between two sites of the cluster"),
         "to": ("cOne", "from a site of the cluster to a nearest neighbor site outside"),
         "on": ("cThree", "between two sites outside that are both nearest neighbors of the cluster"),
         "off": ("cTwo", "from a nearest neighbor site to one that is not a nearest neighbor of the cluster")}
BOTH = ("in", "on")                                        # drawn with two heads: the jump back belongs to the same group
cluster = find_clusters()[2][0]
PATCH = sorted((r, c) for r in range(-1, 7) for c in range(-(r % 2), 7))                               # the sites of the drawing
STYLE = {s: {"M": "mg", "Z": "zn"}.get(element(s), "al") if s in cluster else "al" for s in PATCH}     # its vacancy is a solvent atom here

# ------------------------------------------------------------------ the film: frames at which its parts begin
PLACE, JUMP, STAY, BACK = 5, 10, 4, 4                      # frames of the four parts of a jump
STEP = PLACE + JUMP + STAY + BACK
TEXT = frames(0.2)                                         # a rule comes and goes in this time
turns, jumps, t = {}, [], frames(0.5)                      # turns: first and last frame of a group; jumps: first frame, sites, group, heads
for name in RULES:
    begin, t = t, t + frames(0.4)
    for first, second in JUMPS[name]:
        for a, b in ((first, second), (second, first))[:1 + (name in BOTH)]:
            assert b in neighbors(a) and group_of_jump(a, b, cluster) == name, (name, a, b)
            jumps.append((t, a, b, name, 1 + (a != first)))
            t += STEP
    turns[name] = (begin, t)
    t += TEXT
LAST = t + frames(2.5)                                     # the end is held

# ------------------------------------------------------------------ the picture: the lattice in units of the distance of two sites
PREAMBLE = barrier_groups() + r"""
\tikzset{al/.style={atom=cAl, minimum size=13.4pt}, mg/.style={atom=cMg, minimum size=15.4pt}, zn/.style={atom=cZn, minimum size=13.4pt},
  va/.style={vac, minimum size=13.2pt, line width=1.3pt},
  rim/.style={cTwo!75, line width=27.5pt}, region/.style={cTwo!20, line width=25.9pt},
  near rim/.style={cGray!75, line width=27.5pt}, near region/.style={cGray!18, line width=25.9pt},
  one/.style={-{Stealth[length=8pt, width=7.6pt]}, line width=2.3pt, shorten <=3pt, shorten >=1pt},
  half/.style={one, shorten <=1pt},
  both/.style={{Stealth[length=8pt, width=7.6pt]}-{Stealth[length=8pt, width=7.6pt]}, line width=2.3pt, shorten <=1pt, shorten >=1pt},
  swatch/.style={line width=0.8pt, rounded corners=3pt},
  rule/.style={anchor=west, inner sep=0pt, align=flush left, text width=6.3cm, yshift=-0.8pt}}
\newcommand{\lines}[2]{\foreach \x/\y/\u/\v in #1 \draw[#2] (\x,\y) -- (\u,\v);}
\newcommand{\regions}{\lines{\bondsNear}{near rim} \lines{\bondsNear}{near region} \lines{\bondsCluster}{rim} \lines{\bondsCluster}{region}}
\newcommand{\legend}{%
  \begin{scope}[shift={(7.4,0.55)}, x=1cm, y=1cm]
    \foreach \t/\s/\c [count=\i from 0] in {in/both/cFour, to/one/cOne, on/both/cThree, off/one/cTwo} {
      \draw[\s, \c] (0,-0.7*\i) -- ++(0.92,0); \node[leg] at (1.15,-0.7*\i) {\texttt{\t}};}
    \foreach \t/\c/\f [count=\i from 0] in {cluster/cTwo/20, nearest neighbors/cGray/18} {
      \draw[swatch, \c!75, fill=\c!\f] (0.1,-3.30-0.75*\i) rectangle ++(0.72,0.46); \node[leg] at (1.15,-3.07-0.75*\i) {\t};}
    \foreach \t [count=\i from 0] in {Al, Mg, Zn} {
      \node[atom=c\t, minimum size=12pt] at (0.22+1.18*\i,-4.75) {}; \node[leg] at (0.46+1.18*\i,-4.75) {\t};}
    \node[va, minimum size=11.4pt, line width=1.1pt] at (3.76,-4.75) {}; \node[leg] at (4.0,-4.75) {vacancy};
  \end{scope}}
"""


def place(site):
    return site[1] + 0.5 * (site[0] % 2), -site[0] * 3 ** 0.5 / 2


def node(style, site, shown=1.0, to=None, part=0.0):
    """an atom or the vacancy on a site, or at a part of its way to another site"""
    (x, y), (u, v) = place(site), place(to or site)
    return f"\\node[{style}, {opacity(shown)}] at ({x + (u - x) * part:.3f},{y + (v - y) * part:.3f}) {{}};" if shown > 0 else ""


def arrow(first, second, style, color, shown):
    line = "(%.3f,%.3f) -- (%.3f,%.3f);" % (*place(first), *place(second))
    return (f"\\begin{{scope}}[transparency group, {opacity(shown)}] \\draw[{style}, white, line width=4.2pt] {line}"
            f" \\draw[{style}, {color}] {line} \\end{{scope}}")


def picture(f):
    now = next((j for j in jumps if j[0] <= f < j[0] + STEP), None)
    code = ["\\begin{scope}[shift={(0.80,5.25)}, x=0.92cm, y=0.92cm]", "\\regions"]
    code += [node(STYLE[s], s) for s in PATCH if not now or s not in now[1:3]]
    if now:
        start, a, b = now[:3]
        put, jump, back = ramp(f, start, PLACE), ease(ramp(f, start + PLACE, JUMP)), ramp(f, start + PLACE + JUMP + STAY, BACK)
        # the vacancy takes the place of the atom of the first site, the atom of the second site jumps into it, and both atoms are put back
        code += [node(STYLE[a], a, 1 - put if jump == 0 else back), node("va", a, put * (1 - jump)), node("va", b, jump * (1 - back)),
                 node(STYLE[b], b, back), node(STYLE[b], b, 1 - back, a, jump)]
    for start, a, b, name, heads in jumps:                 # the arrows of the jumps that were made or are under way
        shown = ramp(f, start + PLACE, JUMP)
        if shown > 0 and name not in BOTH:
            code.append(arrow(a, b, "one", RULES[name][0], shown))
        elif shown > 0 and heads == 2:                     # the jump back: the arrow gets its second head
            code.append(arrow(b, a, "both", RULES[name][0], shown))
        elif shown > 0 and f < start + STEP + PLACE + JUMP:
            code.append(arrow(a, b, "half", RULES[name][0], shown))
    code.append("\\legend")
    for i, (name, (color, rule)) in enumerate(RULES.items()):                                          # the rule of the group whose turn it is
        shown = min(ramp(f, turns[name][0], TEXT), 1 - ramp(f, turns[name][1], TEXT))
        if shown > 0:
            code.append(f"\\node[rule, {opacity(shown)}] at ($(7.4,0.55)+(2.2cm,{-0.7 * i:.2f}cm)$) {{\\strut {rule}\\strut}};")
    return "\n".join([c for c in code if c] + ["\\end{scope}", ""])


make("barrier_groups", PREAMBLE, [picture(f) for f in range(LAST + 1)], HEIGHT, folder, out, [int(a) for a in sys.argv[3:]])
