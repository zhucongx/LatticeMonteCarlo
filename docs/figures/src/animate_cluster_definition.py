# cluster_definition.webp and cluster_definition.mp4: how Ansys finds clusters, one rule after the other, drawn on a triangular lattice
# with 6 nearest neighbors instead of 12. No run is needed: the sites, and what is connected, dropped and added, come from
# analysis_schematics.py, which follows SoluteCluster::FindAtomListOfClusters (lmc/ansys/src/SoluteCluster.cpp). Every frame is drawn
# with TikZ. Make both files in this folder with, for any scratch folder T:
#   python3 animate_cluster_definition.py T/frames <output folder> [numbers of single frames: draw only these, for a look]
# Needs what tikz_frames.py needs.
import sys

from analysis_schematics import SMALLEST_CLUSTER_CRITERIA, SOLVENT_BOND_CRITERIA, cluster_definition, find_clusters
from tikz_frames import frames, make, opacity, ramp

folder, out = sys.argv[1:3]
HEIGHT = 1080                                              # of a frame in pixels
clusters = find_clusters()[2]                              # the largest one first: the order of cluster_id
# the pointers of the two numbers end above the site (1, 2) and below the site (6, 6)
assert len(clusters) == 2 and len(clusters[0]) > len(clusters[1]) and (1, 2) in clusters[0] and (6, 6) in clusters[1]

# ------------------------------------------------------------------ the film: frames at which its parts begin
t = frames(0.5)
CONNECT, t = t, t + frames(2.4)                            # rule 1, and after a moment the bonds appear
DROP, t = t, t + frames(2.9)                               # rule 2: the small groups fade, the others get their outline
COUNT, t = t, t + frames(2.4)                              # rule 3: the numbers appear,
JOIN, t = t, t + frames(2.3)                               # those above the limit get a ring, and the outline takes them in
NUMBER, t = t, t + frames(0.7)                             # the clusters are numbered
LAST = t + frames(2.5)                                     # the end is held
FADE = frames(0.4)
BONDS, FADING, OUTLINE = CONNECT + frames(0.6), DROP + frames(0.6), DROP + frames(1.5)
COUNTS, GROW, LABELS = COUNT + frames(1.0), JOIN + frames(0.9), NUMBER + frames(0.3)
RULES = [(CONNECT, "1", "Solute atoms and vacancies that are\\\\ nearest neighbors are connected.", ""),
         (DROP, "2", f"A group of fewer than {SMALLEST_CLUSTER_CRITERIA} of them\\\\ is dropped.",
          f"smallest\\char95 cluster\\char95 criteria {SMALLEST_CLUSTER_CRITERIA}"),
         (COUNT, "3", f"A solvent atom joins a group\\\\ if more than {SOLVENT_BOND_CRITERIA} of its nearest\\\\ neighbors are solute atoms\\\\ "
          "or vacancies of the group.", f"solvent\\char95 bond\\char95 criteria {SOLVENT_BOND_CRITERIA}"),
         (NUMBER, "", "The clusters are numbered\\\\ from the largest to the smallest.", ""),
         (LAST + 99,)]

# ------------------------------------------------------------------ the picture: the lattice in units of the distance of two sites
PREAMBLE = cluster_definition() + r"""
\tikzset{al/.style={atom=cAl, minimum size=21.4pt}, mg/.style={atom=cMg, minimum size=24.8pt}, zn/.style={atom=cZn, minimum size=21.4pt},
  va/.style={vac, minimum size=21pt, line width=2pt},
  joined/.style={circle, draw=cThree, line width=3.2pt, minimum size=22.4pt, inner sep=0pt},
  bond/.style={cDark, line width=4.1pt}, rim/.style={cTwo!75, line width=30pt}, region/.style={cTwo!20, line width=27.8pt},
  number/.style={font=\Large, inner sep=0pt}, pointer/.style={black!70, line width=0.6pt},
  txt/.style={inner sep=0pt, anchor=north west, align=flush left, text width=5.45cm}}
\newcommand{\lines}[2]{\foreach \x/\y/\u/\v in #1 \draw[#2] (\x,\y) -- (\u,\v);}
\newcommand{\nodes}[2]{\foreach \x/\y in #1 \node[#2] at (\x,\y) {};}
\newcommand{\outline}[2]{\begin{scope}[transparency group, #2] \lines{#1}{rim} \lines{#1}{region} \end{scope}}
\newcommand{\bonds}[2]{\begin{scope}[transparency group, #2] \lines{#1}{bond} \end{scope}}
\newcommand{\legend}{%
  \begin{scope}[shift={(10.42,4.78)}, al/.append style={minimum size=11.4pt}, mg/.append style={minimum size=13.2pt},
      zn/.append style={minimum size=11.4pt}, va/.append style={minimum size=11.2pt, line width=1.1pt},
      joined/.append style={minimum size=12pt, line width=1.8pt}, bond/.append style={line width=2.2pt}]
    \node[al] at (0,0) {}; \node[mg] at (0,-0.56) {}; \node[zn] at (0,-1.12) {}; \node[va] at (0,-1.68) {};
    \draw[bond] (-0.27,-2.24) -- ++(0.54,0);
    \begin{scope}[transparency group, opacity=0.3] \node[mg] at (0,-2.80) {}; \end{scope}
    \node[al] at (0,-3.36) {}; \node[joined] at (0,-3.36) {};
    \foreach \t [count=\i from 0] in {Al (solvent), Mg, Zn, vacancy, nearest neighbors, dropped, joins the group}
      \node[leg] at (0.42,-0.56*\i) {\t};
  \end{scope}}
"""


def picture(f):
    bond, drop = opacity(ramp(f, BONDS, FADE)), opacity(1 - 0.7 * ramp(f, FADING, FADE))
    outline, grow = ramp(f, OUTLINE, FADE), ramp(f, GROW, FADE)
    code = ["\\begin{scope}[shift={(0.72,7.78)}, x=1.05cm, y=1.05cm]"]
    if outline > 0 and grow < 1:                           # around the groups that are left,
        code.append(f"\\outline{{\\bondsKept}}{{{opacity(outline * (1 - grow))}}}")
    if grow > 0:                                           # and around these with the solvent atoms that joined
        code.append(f"\\outline{{\\bondsClusters}}{{{opacity(grow)}}}")
    if f > BONDS:
        code.append(f"\\bonds{{\\bondsKept}}{{{bond}}}")
    code.append("\\nodes{\\sitesAl}{al} \\nodes{\\keptMg}{mg} \\nodes{\\keptZn}{zn} \\nodes{\\keptX}{va}")
    code.append(f"\\begin{{scope}}[transparency group, {drop}] " + (f"\\bonds{{\\bondsDropped}}{{{bond}}} " if f > BONDS else "")
                + "\\nodes{\\droppedMg}{mg} \\nodes{\\droppedZn}{zn} \\end{scope}")
    if f > JOIN:
        code.append(f"\\nodes{{\\joinedAl}}{{joined, {opacity(ramp(f, JOIN, FADE))}}}")
    if f > COUNTS:                                         # the neighbors of a solvent atom that are solute atoms or vacancies of the group
        code.append(f"\\foreach \\x/\\y/\\n in \\counts \\node[number, {opacity(ramp(f, COUNTS, FADE))}] at (\\x,\\y) {{\\n}};")
    if f > LABELS:
        code.append(f"\\begin{{scope}}[{opacity(ramp(f, LABELS, FADE))}]\n"
                    "\\draw[pointer] (2.5,0.62) -- (2.5,-0.35); \\node[anchor=south, inner sep=2pt] at (2.9,0.62) {\\texttt{cluster\\char95 id} 1};\n"
                    "\\draw[pointer] (6,-6.70) -- (6,-5.71); \\node[anchor=north, inner sep=2pt] at (6,-6.70) {\\texttt{cluster\\char95 id} 2};\n"
                    "\\end{scope}")
    code.append("\\end{scope}")
    # the rule that is applied: it comes in a fifth of a second, and goes in as much before the next one
    for (start, number, text, key), (end, *_) in zip(RULES, RULES[1:]):
        shown = min(ramp(f, start, frames(0.2)), 1 - ramp(f, end - frames(0.2), frames(0.2)))
        if shown > 0:
            key = f"\\\\[3pt] \\texttt{{{key}}}" if key else ""
            code.append(f"\\begin{{scope}}[{opacity(shown)}] \\node[txt, text width=0.4cm] at (10.15,8.14) {{\\textbf{{{number}}}}};"
                        f" \\node[txt] at ({10.6 if number else 10.15},8.14) {{{text}{key}}}; \\end{{scope}}")
    return "\n".join(code + ["\\legend", ""])


make("cluster_definition", PREAMBLE, [picture(f) for f in range(LAST + 1)], HEIGHT, folder, out, [int(a) for a in sys.argv[3:]])
