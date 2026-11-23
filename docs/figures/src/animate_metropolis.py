# metropolis.webp and metropolis.mp4: five trial exchanges of CMC, or of SA at its current temperature, in the order of the code
# (GenerateLatticeIdJumpPair and SelectEvent in CanonicalMcAbstract.cpp): two sites that hold different elements, the energy change,
# for an energy change that is not negative a random number against exp(-dE/kT), and the exchange if it is accepted.
# The sites, energy changes and random numbers are an example; who is accepted follows from them here. Drawn with TikZ over the plot
# of the acceptance probability, which is made first.
#   python3 animate_metropolis.py <folder for the frames> <output folder> [numbers of single frames: draw only these, for a look]
import math
import os
import sys
from collections import namedtuple

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import style
from style import *
from tikz_film import draw, ease, write

folder, out = sys.argv[1:3]
os.makedirs(folder, exist_ok=True)
# ------------------------------------------------------------------ the example
# The box at the start, and the trials: two sites, dE/kT, the random number (None where the code draws none), the bend of the arrow.
START = {(0, 1): "Mg", (6, 0): "Mg", (7, 2): "Mg", (1, 4): "Zn", (8, 5): "Zn", (5, 3): "Zn", (3, 5): "X"}
END = {(1, 4): "Mg", (5, 5): "Mg", (6, 0): "Mg", (7, 2): "Zn", (2, 1): "Zn", (8, 5): "Zn", (4, 2): "X"}        # the box at the end
Trial = namedtuple("Trial", "a b x u bend")
TRIALS = [Trial((0, 1), (5, 5), -2.4, None, "bend left=25"),
          Trial((8, 5), (6, 2), 2.6, 0.38, "bend left=30"),
          Trial((3, 5), (4, 2), 0.4, 0.52, "bend right=30"),
          Trial((5, 3), (2, 1), -0.5, None, "bend right=30"),
          Trial((1, 4), (7, 2), 0.9, 0.18, "bend left=42, looseness=1.25")]
NAME = {"Al": "Al", "Mg": "Mg", "Zn": "Zn", "X": "vacancy"}
boxes, accepted = [dict(START)], []                        # the box before every trial
for t in TRIALS:
    box = dict(boxes[-1])
    one, two = box.get(t.a, "Al"), box.get(t.b, "Al")
    assert one != two and (t.u is None) == (t.x < 0)       # as in the code: different elements, a random number unless dE < 0
    accepted.append(t.x < 0 or t.u < math.exp(-t.x))       # SelectEvent
    if accepted[-1]:
        box[t.a], box[t.b] = two, one
    boxes.append({site: e for site, e in box.items() if e != "Al"})
assert boxes[-1] == END

# ------------------------------------------------------------------ panel (b), and where its axes are
style.OUT = folder
XLIM, YLIM, AXES, SCALE = (-3, 5.2), (0, 1.15), (4.8, 2.7), 0.6            # the axes in inches, and the scale of \plot in figstyle.tex
x = np.linspace(*XLIM, 821)
fig, ax = plt.subplots()
ax.plot(x, np.minimum(1, np.exp(-x)), color=RED, lw=2.2)
ax.axvline(0, color="k", lw=1, ls="--")
ax.text(-1.5, 0.915, "always", ha="center", va="baseline", fontsize=14)
ax.text(-1.5, 0.83, "accepted", ha="center", va="baseline", fontsize=14)
ax.text(1.15, 0.5, r"$\exp(-\Delta E/k_{\text{B}}T)$", fontsize=14)
ax.set_xlim(*XLIM)
ax.set_ylim(*YLIM)
ax.set_xlabel(r"$\Delta E/k_{\text{B}}T$")
ax.set_ylabel("Acceptance probability")
grid(ax)
set_size(*AXES, ax)
edge = fig.get_tightbbox(fig._get_renderer()).padded(mpl.rcParams["savefig.pad_inches"])       # what save() cuts out, in inches
corner = ax.get_position()
ORIGIN = [2.54 * SCALE * (c * s - e) for c, s, e in zip((corner.x0, corner.y0), fig.get_size_inches(), (edge.x0, edge.y0))]
UNIT = [2.54 * SCALE * a / (hi - lo) for a, (lo, hi) in zip(AXES, (XLIM, YLIM))]          # cm of the figure for 1 on an axis
save(fig, "plot_metropolis.pdf")

# ------------------------------------------------------------------ the film: the stages of a trial and their frames
STAGES = ("pair", "energy", "random", "result")


def frames(n):
    return {"pair": 14, "energy": 18, "random": 20 if TRIALS[n].u is not None else 0, "result": 26}


film = [(n, name, frames(n)[name]) for n in range(len(TRIALS)) for name in STAGES]
last = sum(count for _, _, count in film)


def state(f):
    """at frame f: the trial, and for each of its stages the number of frames that have passed (0: not begun)"""
    trial, passed, start = 0, dict.fromkeys(STAGES, 0), 0
    for n, name, count in film:
        if f <= start:
            break
        if n != trial:
            trial, passed = n, dict.fromkeys(STAGES, 0)
        passed[name], start = min(f - start, count), start + count
    return trial, passed


# ------------------------------------------------------------------ the picture, lengths in cm
WIDTH = 425.142 / 28.4528                                  # of the picture with its border
X0, Y1 = -2.723 / 28.4528, 163.964 / 28.4528
BOX = f"({X0:.4f},{Y1 - WIDTH * 9 / 16:.4f}) rectangle ({X0 + WIDTH:.4f},{Y1:.4f})"
COLUMNS = ((0.55, "trial"), (2.6, "sites"), (5.5, r"$\Delta E/k_{\mathrm{B}}T$"), (8.35, "acceptance probability"), (11.4, "random number"),
           (13.75, "exchange"))
ROW = 0.39                                                 # distance of two rows of the table
PREAMBLE = r"""\tikzset{dot/.style={atom=cAl!80, minimum size=6.4pt}, am/.style={atom=cMg, minimum size=11pt},
  az/.style={atom=cZn, minimum size=11pt}, ax/.style={vac, minimum size=10.5pt},
  site/.style={circle, inner sep=0pt, minimum size=11pt},
  cell/.style={font=\footnotesize, inner sep=0pt, anchor=base},
  point/.style={circle, inner sep=0pt, minimum size=5.5pt, fill=cFour, draw=black, line width=0.5pt},
  number/.style={diamond, inner sep=0pt, minimum size=7pt, fill=cTwo, draw=black, line width=0.5pt}}
\newcommand{\fixed}{%
  \plot{b}{6.3}{0}{plot_metropolis.pdf}
  \path (b.north) ++(0,-4.200) coordinate (base);
  \specieslegend{0.2}{0.3}
  \vacancylegend{3.8}{0.3}
  \plotlabel{b}{b}
  \node[above lab] at ($(-0.05,0 |- b.north)+(0,0.02)$) {\panel{a}};
  \foreach \x/\t in {COLUMNS} \node[cell] at (\x,-0.36) {\t};
  \draw[black!45, line width=0.5pt] (0.05,-0.51) -- (14.75,-0.51);}
""".replace("COLUMNS", ", ".join(f"{x}/{{{t}}}" for x, t in COLUMNS))
STYLE = {"Al": "dot", "Mg": "am", "Zn": "az", "X": "ax"}


def fade(w):
    return f", opacity={w:.2f}" if w < 0.995 else ""


def row(n, passed):
    """the cells of the table for trial n, as far as the trial has come"""
    t, box = TRIALS[n], boxes[n]
    cells = [str(n + 1), f"{NAME[box.get(t.a, 'Al')]} and {NAME[box.get(t.b, 'Al')]}"]
    if passed["energy"]:
        cells += [f"${t.x:.1f}$", "1" if t.x < 0 else f"{math.exp(-t.x):.2f}"]
    if passed["random"] or passed["result"]:
        cells.append("none" if t.u is None else f"{t.u:.2f}")
    if passed["result"]:
        cells.append("\\textcolor{cThree!70!black}{accepted}" if accepted[n] else "\\textcolor{cOne!85!black}{rejected}")
    return [f"\\node[cell] at ({x},{-0.87 - ROW * n:.3f}) {{{c}}};" for (x, _), c in zip(COLUMNS, cells)]


def page(f):
    n, passed = state(f)
    t, box = TRIALS[n], boxes[n]
    done = dict.fromkeys(STAGES, 99)
    swap = ease(passed["result"] / 14) if accepted[n] else 0.0
    mark = min(passed["pair"] / 4, 1) * (1 if accepted[n] else 1 - ease((passed["result"] - 4) / 12))    # a rejected pair fades
    code = [f"\\begin{{tikzpicture}}\n\\path[use as bounding box] {BOX};\n\\fixed"]
    # ---------------- (a) the box, the two sites of the trial and their exchange
    code += ["\\begin{scope}[shift={(0,0 |- base)}, x=0.5879cm, y=0.5879cm]", "\\draw[black!70, line width=0.6pt] (0,0) rectangle (9,7);",
             "\\begin{scope}[shift={(0.5,0.5)}]"]
    if mark:
        code += [f"\\fill[cOne!30{fade(mark * (1 - swap))}] ({x},{y}) circle[radius=7.5pt];" for x, y in (t.a, t.b)]
    code += [f"\\node[{STYLE[box.get((x, y), 'Al')]}] at ({x},{y}) {{}};" for x in range(9) for y in range(7) if (x, y) not in (t.a, t.b)]
    code += [f"\\node[site] (a) at ({t.a[0]},{t.a[1]}) {{}}; \\node[site] (b) at ({t.b[0]},{t.b[1]}) {{}};"]
    if mark:
        code.append(f"\\draw[dbl, line width=1.0pt, cOne{fade(mark)}] (a) to[{t.bend}] (b);")
    code += [f"\\path (a.center) to[{t.bend}] node[pos={w:.3f}, {STYLE[box.get(site, 'Al')]}] {{}} (b.center);"
             for site, w in sorted(((t.a, swap), (t.b, 1 - swap)), key=lambda s: box.get(s[0], "Al") != "Al")]      # an Al atom passes below
    code += ["\\end{scope}", "\\end{scope}"]
    # ---------------- (b) the trial on the curve, and the random number against the probability
    p = min(1.0, math.exp(-t.x))
    code.append(f"\\begin{{scope}}[shift={{({6.3 + ORIGIN[0] - XLIM[0] * UNIT[0]:.4f},{ORIGIN[1]:.4f})}}, x={UNIT[0]:.4f}cm, y={UNIT[1]:.4f}cm]")
    code += ["\\node[point] at (2.35,1.055) {}; \\node[cell, anchor=base west] at (2.55,1.025) {trial};",
             "\\node[number] at (2.35,0.925) {}; \\node[cell, anchor=base west] at (2.55,0.895) {random number};"]
    if passed["random"]:
        w = fade(min(passed["random"] / 4, 1))
        code += [f"\\draw[black!70, line width=0.6pt, densely dashed{w}] ({t.x},{p:.4f}) -- ({t.x},{t.u});", f"\\node[number{w}] at ({t.x},{t.u}) {{}};"]
    if passed["energy"]:
        code.append(f"\\node[point{fade(min(passed['energy'] / 4, 1))}] at ({t.x},{p:.4f}) {{}};")
    code.append("\\end{scope}")
    # ---------------- the table of the trials
    for m in range(n + 1 if passed["pair"] else n):
        code += row(m, passed if m == n else done)
    code.append("\\end{tikzpicture}\n")
    return "\n".join(code)


wanted = [int(a) for a in sys.argv[3:]] or range(last + 1)
files = draw([page(f) for f in wanted], PREAMBLE, folder, WIDTH)
if sys.argv[3:]:
    print("\n".join(f"{f} {png}" for f, png in zip(wanted, files)))
else:
    print("accepted:", accepted)
    write(files, out, "metropolis")
