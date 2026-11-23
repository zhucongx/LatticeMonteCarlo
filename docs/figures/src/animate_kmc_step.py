# kmc_step.webp and kmc_step.mp4: three steps of first-order KMC (KineticMcFirstSerial), each in the order of the code: the twelve rates,
# the cumulative probabilities, the first random number for the time step, the second one for the jump, the time, the exchange.
# The rates and the random numbers are an example; everything else is computed from them here. Every frame is drawn with TikZ.
#   python3 animate_kmc_step.py <folder for the frames> <output folder> [numbers of single frames: draw only these, for a look]
import math
import sys
from collections import namedtuple

from tikz_film import draw, ease, mix, write

folder, out = sys.argv[1:3]
# ------------------------------------------------------------------ the example
# A step: the atoms on the twelve neighbor sites, the rates of their jumps in 1e6/s, the random numbers for the time and for the jump.
# The ring is a sketch, but it keeps one thing of the lattice: the atom that has jumped is found opposite to the direction of its jump.
Step = namedtuple("Step", "atoms rates u x")
STEPS = [Step("Mg Al Al Al Al Al Al Al Al Al Al Al".split(), [2.6, 0.7, 0.9, 0.5, 0.6, 0.5, 0.8, 0.6, 1.0, 0.5, 0.7, 0.6], 0.67, 0.39),
         Step("Al Al Al Al Al Al Mg Al Al Zn Al Al".split(), [0.64, 0.48, 0.80, 0.96, 0.48, 0.80, 2.24, 0.64, 0.48, 7.04, 0.80, 0.64], 0.85, 0.68),
         Step("Al Al Al Zn Al Mg Al Mg Al Zn Al Al".split(), [0.4, 1.0, 0.4, 6.0, 0.4, 1.6, 0.4, 4.4, 0.6, 3.0, 0.8, 1.0], 0.55, 0.62)]
total = [sum(s.rates) for s in STEPS]                      # 1e6/s
prob = [[r / t for r in s.rates] for s, t in zip(STEPS, total)]
stop = [[sum(p[:i + 1]) for i in range(12)] for p in prob]                 # cumulative probabilities
chosen = [next(i for i in range(12) if c[i] >= s.x) for s, c in zip(STEPS, stop)]          # as SelectEvent: the first that reaches x
dt = [-math.log(s.u) / t * 1000 for s, t in zip(STEPS, total)]             # as CalculateTime, in ns
clock = [sum(dt[:n]) for n in range(len(STEPS) + 1)]
for n in range(len(STEPS) - 1):
    assert STEPS[n + 1].atoms[(chosen[n] + 6) % 12] == STEPS[n].atoms[chosen[n]]
    assert abs(dt[n] - round(dt[n])) < 0.25 and abs(clock[n + 1] - round(clock[n + 1])) < 0.25           # shown without decimals

# ------------------------------------------------------------------ the film: the stages of a step and their frames
STAGES = (("rates", 22), ("bar", 22), ("time", 22), ("draw", 24), ("clock", 14), ("jump", 12), ("move", 10))
film = [(n, name, count) for n in range(len(STEPS)) for name, count in STAGES if n < len(STEPS) - 1 or name not in ("jump", "move")]
last = sum(count for _, _, count in film)


def state(f):
    """at frame f: the step, and for each of its stages the number of frames that have passed (0: not begun)"""
    step, passed, start = 0, dict.fromkeys(dict(STAGES), 0), 0
    for n, name, count in film:
        if f <= start:
            break
        if n != step:
            step, passed = n, dict.fromkeys(dict(STAGES), 0)
        passed[name], start = min(f - start, count), start + count
    return step, passed


# ------------------------------------------------------------------ the picture, lengths in cm
WIDTH = 474.015 / 28.4528                                  # of the picture with its border
X0, Y1 = -82.668 / 28.4528, 85.2355 / 28.4528
BOX = f"({X0:.4f},{Y1 - WIDTH * 9 / 16:.4f}) rectangle ({X0 + WIDTH:.4f},{Y1:.4f})"
COLOR = {"Al": "cAl", "Mg": "cMg", "Zn": "cZn"}
GRAY, RED = (224, 224, 224), (241, 208, 206)               # black!12 and cOne!25: a time step of the clock, and the newest one
BAR, NS = 9.6, 0.096                                       # length of the bar of (b), and of a nanosecond in (c)
PREAMBLE = r"""\tikzset{txt/.style={anchor=base west, inner sep=0pt}}
\newcommand{\fixed}{%
  \foreach \num in {1,...,12} {
    \pgfmathsetmacro{\ang}{90 - (\num-1)*30}
    \pgfmathsetmacro{\opp}{\ang + 180}
    \node[font=\footnotesize, anchor=\opp, inner sep=2pt] at (\ang:2.38cm) {\num};}
  \node at (0,-3.2) {arrow width: rate $r_i$};
  \begin{scope}[shift={(3.95,-0.45)}]
    \foreach \t in {0, 0.25, 0.5, 0.75, 1} {
      \draw[line width=0.6pt] (\t*BAR,0.7) -- ++(0,0.13) node[font=\footnotesize, anchor=south, inner sep=2pt] {\t};}
    \node[anchor=south west, inner xsep=0pt] at (0,1.25) {cumulative probability};
  \end{scope}
  \specieslegend{4.15}{-3.2}
  \vacancylegend{7.75}{-3.2}
  \node[lab] at (-2.85,2.95) {\panel{a}};
  \node[lab] at (3.40,2.95) {\panel{b}};
  \node[lab] at (-2.85,-3.92) {\panel{c}};
  \begin{scope}[shift={(3.95,-5.42)}]
    \foreach \t in {0, 20, ..., 100} {
      \draw[line width=0.6pt] (\t*NS,0) -- ++(0,-0.13) node[font=\footnotesize, anchor=north, inner sep=2pt] {\t};}
    \node[anchor=north, inner sep=0pt] at (0.5*BAR,-0.50) {time (ns)};
  \end{scope}}
""".replace("BAR", str(BAR)).replace("NS", str(NS))


def place(i, radius=2.1):
    a = math.radians(90 - 30 * i)
    return radius * math.cos(a), radius * math.sin(a)


def fade(w):
    return f", opacity={w:.2f}" if w < 0.995 else ""


def page(f):
    n, passed = state(f)
    s, k = STEPS[n], chosen[n]
    jump, move = ease(passed["jump"] / 11), ease(passed["move"] / 10)
    code = [f"\\begin{{tikzpicture}}\n\\path[use as bounding box] {BOX};\n\\fixed"]
    # ---------------- (a) the vacancy, its neighbors and the arrows of the jumps
    kx, ky = place(k)
    sx, sy = -kx * move, -ky * move                        # the picture follows the vacancy to its new site
    if move > 0.65:                                        # the old neighbors go out of sight, then those of the new site come
        code += [f"\\node[atom={COLOR[a]}, minimum size=14pt{fade((move - 0.65) / 0.35)}] at "
                 f"({place(i)[0] + kx + sx:.3f},{place(i)[1] + ky + sy:.3f}) {{}};" for i, a in enumerate(STEPS[n + 1].atoms) if i != (k + 6) % 12]
    elif move < 0.35:
        code += [f"\\node[atom={COLOR[a]}, minimum size=14pt{fade(1 - move / 0.35)}] (n{i}) at ({place(i)[0] + sx:.3f},{place(i)[1] + sy:.3f}) {{}};"
                 for i, a in enumerate(s.atoms) if i != k]
    code += [f"\\node[vac, minimum size=14pt] (v) at ({kx * jump + sx:.3f},{ky * jump + sy:.3f}) {{}};",
             f"\\node[atom={COLOR[s.atoms[k]]}, minimum size=14pt] (n{k}) at ({kx * (1 - jump) + sx:.3f},{ky * (1 - jump) + sy:.3f}) {{}};"]
    marked = passed["draw"] > 6                            # the random number has arrived on the bar
    if not passed["jump"]:
        for i in range(min(passed["rates"], 12)):          # the events are built one after the other
            r = s.rates[i] / 20                            # one scale for all steps: 17 pt for the sum of the rates of the last step
            code.append(f"\\draw[-{{Stealth[length={4 + 18 * r:.2f}pt, width={0.85 * (4 + 18 * r):.2f}pt]}}, line width={17 * r:.3f}pt, "
                        f"{'cOne' if marked and i == k else 'black!70'}, line cap=butt] (v) -- (n{i});")
    # ---------------- (b) the probabilities, laid end to end, and the second random number
    code.append("\\begin{scope}[shift={(3.95,-0.45)}" + (f", opacity={1 - move:.2f}" if move else "") + "]")
    for i in range(min(passed["bar"], 12)):
        a, b = BAR * (stop[n][i] - prob[n][i]), BAR * stop[n][i]
        code += [f"\\fill[{COLOR[s.atoms[i]]}!85] ({a:.4f},0) rectangle ({b:.4f},0.7);",
                 f"\\draw[white, line width=0.8pt] ({b:.4f},0) -- ({b:.4f},0.7);",
                 f"\\node[font=\\footnotesize, anchor=north, inner sep=3pt] at ({(a + b) / 2:.4f},0) {{{i + 1}}};"]
    code.append(f"\\draw[black!65, line width=0.6pt, opacity=1] (0,0) rectangle ({BAR},0.7);")
    if passed["draw"]:
        x, a, b = BAR * s.x, BAR * (stop[n][k] - prob[n][k]), BAR * stop[n][k]
        if marked:
            code += [f"\\draw[cOne, line width=1.4pt] ({a:.4f},0) rectangle ({b:.4f},0.7);",
                     f"\\draw[decorate, decoration={{brace, mirror, amplitude=4pt}}, line width=0.5pt] ({a:.4f},-0.50) -- ({b:.4f},-0.50) "
                     f"node[midway, below=5pt] {{$p_{{{k + 1}}} = r_{{{k + 1}}}/\\sum_i r_i$}};"]
        code += [f"\\draw[jump] ({x:.4f},1.85) -- ({x:.4f},{1.85 - 1.11 * min(passed['draw'], 6) / 6:.4f});",
                 f"\\node[anchor=south] at ({x:.4f},1.85) {{second random number {s.x:.2f}}};"]
    code.append("\\end{scope}")
    # ---------------- (c) the time step and the clock
    if passed["rates"] >= 12:
        code.append(f"\\node[txt{fade(1 - move)}] at (-2.15,-4.62) {{$\\sum_i r_i = {total[n] / 10:.1f}\\times10^{{7}}\\ \\mathrm{{s}}^{{-1}}$}};")
    if passed["time"]:
        w = fade(min(passed["time"] / 4, 1 - move))
        code += [f"\\node[txt{w}] at (-2.15,-5.17) {{first random number $u = {s.u:.2f}$}};",
                 f"\\node[txt{w}] at (-2.15,-5.72) {{time step $-\\ln(u)/\\sum_i r_i = {dt[n]:.0f}$ ns}};"]
    grown = min(passed["clock"] / 10, 1)
    now = clock[n] + dt[n] * grown
    code.append("\\begin{scope}[shift={(3.95,-5.42)}]")
    for m in range(n + 1 if grown else n):
        a, b = NS * clock[m], NS * (clock[m] + dt[m] * (grown if m == n else 1))
        newest = 1 - move if m == n else 0.0
        code += [f"\\fill[fill={mix(GRAY, RED, newest)}] ({a:.4f},0) rectangle ({b:.4f},0.45);",
                 f"\\draw[white, line width=0.8pt] ({b:.4f},0) -- ({b:.4f},0.45);"]
        if m < n or grown == 1:
            code.append(f"\\node[font=\\footnotesize, inner sep=0pt] at ({(a + b) / 2:.4f},0.225) {{{dt[m]:.0f} ns}};")
    code.append(f"\\draw[black!65, line width=0.6pt] (0,0) rectangle ({BAR},0.45);")
    if grown:
        a, b = NS * clock[n], NS * now
        code.append(f"\\draw[cOne, line width=1.4pt{fade(1 - move)}] ({a:.4f},0) rectangle ({b:.4f},0.45);")
    code += [f"\\node[txt, anchor=base east] at ({BAR},0.68) {{time: \\makebox[1.2em][r]{{\\textbf{{{now:.0f}}}}} ns}};",     # of a fixed width
             "\\end{scope}", "\\end{tikzpicture}\n"]
    return "\n".join(code)


wanted = [int(a) for a in sys.argv[3:]] or range(last + 1)
files = draw([page(f) for f in wanted], PREAMBLE, folder, WIDTH)
if sys.argv[3:]:
    print("\n".join(f"{f} {png}" for f, png in zip(wanted, files)))
else:
    print(f"jumps {[k + 1 for k in chosen]}, time steps {[round(t, 2) for t in dt]} ns")
    write(files, out, "kmc_step")
