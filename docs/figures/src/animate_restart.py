# restart.webp and restart.mp4: a run that stops, what has to be done with its log, and its restart, as a time line. No run is needed.
# It follows Dump and the handling of restart_steps in lmc/mc/src/KineticMcAbstract.cpp and McAbstract.cpp. Every frame is drawn with
# TikZ. Make both files in this folder with, for any scratch folder T:
#   python3 animate_restart.py T/frames <output folder> [numbers of single frames: draw only these, for a look]
# Needs what tikz_frames.py needs.
import sys

from tikz_frames import frames, make, opacity, ramp

folder, out = sys.argv[1:3]
HEIGHT = 624                                               # of a frame in pixels
XO, XN, XS, DC, DL = 1.85, 6.60, 8.25, 2.82, 0.47          # x of step 0, of step N and of the stop; config_dump_steps and log_dump_steps, in cm
FIRST_ROWS = [XO] + [XN + k * DL for k in range(-7, 4)]    # the rows of the first run; those of its first steps are left out
FIRST_SHEETS = {XO: "0.cfg.gz", XN - DC: "", XN: "N.cfg.gz"}                   # its configurations
ROWS = [XN + k * DL for k in range(1, 19)]                 # the rows of the restart: none for step N
SHEETS = [XN + k * DC for k in (1, 2, 3)]                  # its configurations
KEYS = ["\\file{config_filename N.cfg.gz}", "\\file{restart_steps N}", "\\file{restart_energy}, \\file{restart_time}"]

# ------------------------------------------------------------------ the film: frames at which its parts begin
SPEED = (XS - XO - 0.2) / frames(2.0)                      # of the front of both runs, in cm per frame
POP = 3                                                    # frames in which a row or a configuration appears
t = frames(0.5)
BEGIN, t = t, t + frames(0.2)                              # the first run begins: its first row and configuration
RUN, t = t, t + frames(2.0)                                # it advances
STOP, t = t, t + frames(0.9)                               # and stops
MARK, t = t, t + frames(1.1)                               # the rows after step N are marked
REMOVE, t = t, t + frames(0.9)                             # and removed
KEY, t = t, t + frames(2.2)                                # the keys of the restart, one after the other
AGAIN, t = t, t + frames(0.9)                              # the restart begins; the row of step N is not written again
GO, t = t, t + round((3 * DC + 0.1) / SPEED)               # it advances to maximum_steps
END, t = t, t + frames(0.2) + POP                          # and writes end.cfg.gz
LAST = t + frames(2.5)                                     # the end is held


def reached(f, start, x0, x):
    """how far a row or a configuration at x has appeared in frame f, if the front left x0 in the frame start"""
    return ramp(f, start + (x + 0.1 - x0) / SPEED, POP)


# ------------------------------------------------------------------ the picture, lengths in cm
PREAMBLE = r"""
\newcommand{\file}[1]{{\footnotesize\ttfamily\detokenize{#1}}}
\colorlet{gone}{cOne!70}
\tikzset{bar/.style={black!50, line width=0.6pt, rounded corners=2.5pt},
  row/.style={circle, fill=cDark, minimum size=4.8pt, inner sep=0pt},
  norow/.style={circle, draw=cDark, densely dotted, line width=0.7pt, fill=white, minimum size=5pt, inner sep=0pt},
  lbl/.style={inner sep=1.5pt, anchor=base},
  side/.style={anchor=base west, inner sep=0pt, yshift=-3.2pt},
  tick/.style={line width=0.6pt},
  sheet/.style={draw=cFour!70!black, fill=cFour!30, line width=0.6pt}}
% a configuration: a sheet with a folded corner at (#1,#2), drawn with the option #3
\newcommand{\sheet}[3]{%
  \begin{scope}[shift={(#1,#2)}, transparency group, #3]
    \filldraw[sheet] (-0.14,0) -- (-0.14,0.38) -- (0.04,0.38) -- (0.14,0.28) -- (0.14,0) -- cycle;
    \draw[sheet, fill=none] (0.04,0.38) -- (0.04,0.28) -- (0.14,0.28);
  \end{scope}}
\newcommand{\cut}[3]{% a cut through a bar or an axis at x = #1 from y = #2 to #3
  \draw[white, line width=3pt, line cap=butt] (#1-0.08,#2) -- (#1+0.08,#3);
  \draw[black!60, line width=0.6pt] (#1-0.13,#2) -- (#1+0.03,#3) (#1-0.03,#2) -- (#1+0.13,#3);}
\def\xo{XO}\def\xN{XN}\def\xs{XS}\def\dc{DC}\def\dl{DL}
% the axis of the steps, the names of the rows and the legend
\newcommand{\fixed}{%
  \draw[arr] (\xo-0.2,-1.5) -- (\xN+3*\dc+0.75,-1.5);
  \cut{\xo+0.5}{-1.62}{-1.38}
  \foreach \x in {\xo,\xN-\dc,\xN,\xN+\dc,\xN+2*\dc,\xN+3*\dc} \draw[tick] (\x,-1.5) -- ++(0,0.1);
  \node[lbl] at (\xo,-1.9) {0};
  \node[lbl] at (\xN,-1.9) {\file{N}};
  \node[lbl] at (\xN+3*\dc,-1.9) {\file{maximum_steps}};
  \node[side] at (0,-1.5) {step};
  \node[side] at (0,1.7) {first run};
  \node[side] at (0,-0.3) {restart};
  \sheet{11.35}{2.16}{} \node[anchor=base west, inner sep=0pt] at (11.65,2.26) {configuration};
  \node[row] at (11.35,1.82) {}; \node[side] at (11.65,1.82) {row of the log};}
""".replace("XO", str(XO)).replace("XN", str(XN)).replace("XS", str(XS)).replace("DC", str(DC)).replace("DL", str(DL))


def picture(f):
    code = ["\\begin{scope}[shift={(0.03,2.27)}]", "\\fixed"]
    if f > BEGIN:                                          # ---------------- the first run
        begin = opacity(ramp(f, BEGIN, frames(0.2)))
        x = min(XO + 0.2 + SPEED * max(f - RUN, 0), XS)    # its front
        mark, remove = opacity(ramp(f, MARK, frames(0.3))), ramp(f, REMOVE, frames(0.4))
        code.append(f"\\begin{{scope}}[transparency group, {begin}] \\begin{{scope}} \\clip[rounded corners=2.5pt] (\\xo-0.2,1.5) rectangle ({x:.3f},1.9);"
                    f" \\fill[cGray!22] (\\xo-0.2,1.5) rectangle ({x:.3f},1.9);"
                    + (f" \\fill[cOne!16, {mark}] (\\xN+0.5*\\dl,1.5) rectangle (\\xs,1.9);" if f > MARK else "")
                    + f" \\end{{scope}} \\draw[bar] (\\xo-0.2,1.5) rectangle ({x:.3f},1.9);"
                    + (" \\cut{\\xo+0.5}{1.42}{1.98}" if x > XO + 0.7 else "") + " \\end{scope}")
        for r in FIRST_ROWS:                               # a row after step N turns pale when it is removed
            if (shown := reached(f, RUN, XO + 0.2, r) if r > XO else 1) > 0:
                color = f"fill=gone!{100 * remove:.0f}!cDark, " if r > XN + 0.1 else ""
                code.append(f"\\node[row, {color}{opacity(shown) if r > XO else begin}] at ({r:.3f},1.7) {{}};")
        for s, name in FIRST_SHEETS.items():
            if (shown := reached(f, RUN, XO + 0.2, s) if s > XO else 1) > 0:
                shown = opacity(shown) if s > XO else begin
                code.append(f"\\sheet{{{s:.3f}}}{{1.98}}{{{shown}}}" + (f" \\node[lbl, {shown}] at ({s:.3f},2.48) {{\\file{{{name}}}}};" if name else ""))
        if f > STOP:
            shown = opacity(ramp(f, STOP, frames(0.2)))
            code.append(f"\\draw[cOne, line width=1.6pt, {shown}] (\\xs,1.38) -- (\\xs,2.02);"
                        f" \\node[side, text=cOne!85!black, {shown}] at (\\xs+0.2,1.7) {{run stops}};")
        if f > MARK:
            code.append(f"\\begin{{scope}}[{mark}] \\draw[decorate, decoration={{brace, mirror, amplitude=3.5pt}}, line width=0.6pt, cOne!85!black]"
                        " (\\xN+0.6*\\dl,1.4) -- (\\xs-0.08,1.4);"
                        " \\node[anchor=north west, inner sep=0pt, text=cOne!85!black] at (\\xN+0.35,1.13) {remove these rows}; \\end{scope}")
    if f > AGAIN:                                          # ---------------- the restart
        x = min(XN + 0.2 + SPEED * max(f - GO, 0), XN + 3 * DC + 0.2)
        code.append(f"\\draw[bar, fill=cTwo!16, {opacity(ramp(f, AGAIN, frames(0.3)))}] (\\xN-0.2,-0.5) rectangle ({x:.3f},-0.1);")
        code += [f"\\node[row, {opacity(reached(f, GO, XN + 0.2, r))}] at ({r:.3f},-0.3) {{}};" for r in ROWS if reached(f, GO, XN + 0.2, r) > 0]
        if f > END:                                        # at maximum_steps also end.cfg.gz, which lies behind the configuration of the step
            shown = opacity(ramp(f, END + frames(0.2), POP))
            code.append(f"\\sheet{{\\xN+3*\\dc+0.1}}{{0.05}}{{{shown}}} \\node[lbl, {shown}] at (\\xN+3*\\dc,0.55) {{\\file{{end.cfg.gz}}}};")
        code += [f"\\sheet{{{s:.3f}}}{{-0.02}}{{{opacity(reached(f, GO, XN + 0.2, s))}}}" for s in SHEETS if reached(f, GO, XN + 0.2, s) > 0]
    if f > KEY:                                            # ---------------- from the one to the other
        code.append(f"\\draw[arr, cTwo, line width=1pt, {opacity(ramp(f, KEY, frames(0.3)))}] (\\xN,1.62) -- (\\xN,-0.2); \\node[row] at (\\xN,1.7) {{}};")
        gray = [100 * ramp(f, KEY + frames(0.3 + 0.4 * k), frames(0.2)) for k in range(3)]             # a line appears: from white to black
        code.append("\\node[anchor=east, align=right, inner sep=0pt] at (\\xN-0.35,0.42) {"
                    + "\\\\ ".join(f"\\textcolor{{black!{g:.0f}}}{{{key}}}" for g, key in zip(gray, KEYS)) + "};")
    if f > AGAIN:
        shown = opacity(ramp(f, AGAIN, frames(0.3)))
        code.append(f"\\node[norow, {shown}] (again) at (\\xN,-0.3) {{}};"
                    f" \\node[anchor=base east, inner sep=0pt, yshift=-3.2pt, {shown}] (note) at (\\xN-0.35,-0.93) {{not written again}};"
                    f" \\draw[black!70, line width=0.6pt, {shown}] ($(note.east)+(0.06,0.03)$) to[out=20, in=-100] (again);")
        # the two intervals, once the front has passed them
        for x, width, name in ((XN + 3 * DL, DL, "log_dump_steps"), (XN + 2 * DC, DC, "config_dump_steps")):
            if (shown := reached(f, GO, XN + 0.2, x + width)) > 0:
                code.append(f"\\draw[dbl, {opacity(shown)}] ({x:.3f},-0.76) -- ({x + width:.3f},-0.76);"
                            f" \\node[lbl, {opacity(shown)}] at ({x + width / 2:.3f},-1.14) {{\\file{{{name}}}}};")
    return "\n".join(code + ["\\end{scope}", ""])


make("restart", PREAMBLE, [picture(f) for f in range(LAST + 1)], HEIGHT, folder, out, [int(a) for a in sys.argv[3:]])
