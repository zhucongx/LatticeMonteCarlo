# banner.svg, the picture at the top of README.md: a close-packed layer of Al, Mg and Zn atoms in which a vacancy makes a round of
# jumps. It is a vector drawing that moves by itself (SMIL), without a background, so that it fits a light and a dark page.
#   python3 banner.py <output folder>
# The vacancy goes around a ring of Al atoms and changes places twice with a Zn atom and once with a Mg atom on the way, there and back.
# After the round every site holds the element it had at the start, so the animation repeats without a break.
import os
import sys

WIDTH, HEIGHT = 760, 128
A, RADIUS = 38.0, 15.6                                     # distance of two neighbors and radius of an atom, in pixels
ROWS, COLUMNS = 3, 20
MOVE, REST = 0.30, 0.45                                    # seconds that a jump takes, and that the vacancy stays afterwards
BALL = {"Al": ("#F3F4F6", "#AEB3BA", "#6C737D"), "Mg": ("#FBE6CF", "#DE8D3A", "#A25B17"), "Zn": ("#D3E3F5", "#3F78B5", "#234A78")}
RED = "#C8423B"
SOLUTE = {(0, 2): "Zn", (1, 2): "Mg", (2, 3): "Zn", (1, 3): "Zn", (2, 2): "Mg",                 # a cluster on the left
          (2, 7): "Zn", (0, 12): "Mg", (2, 11): "Mg", (1, 14): "Zn", (0, 15): "Mg", (1, 15): "Zn", (2, 15): "Mg", (0, 16): "Zn",
          (2, 18): "Zn", (0, 0): "Mg", (0, 5): "Mg"}
START = (1, 5)                                             # the site of the vacancy
ROUND = [(1, 6), (1, 7), ((2, 7), 2), (1, 8), (1, 9), (0, 10), (0, 9), (0, 8), (0, 7), (0, 6), ((0, 5), 1), (1, 5)]


def place(site):
    row, column = site
    return 9.0 + A * column + (A / 2 if row % 2 else 0), 31.0 + A * 0.866 * row


def jumps():
    """ the jumps of the round: (site that the atom leaves, site of the vacancy that it enters) """
    vacancy, found = START, []
    for step in ROUND:
        if isinstance(step[0], tuple):                     # an exchange with a neighbor that is not on the ring, there and back
            neighbor, times = step
            found += [(neighbor, vacancy), (vacancy, neighbor)] * times
        else:
            found.append((step, vacancy))
            vacancy = step
    assert vacancy == START
    return found


def moving(points, total):
    """ the animation of an atom or the vacancy that is at points[k][1] from the time points[k][0] on, after a jump of MOVE seconds """
    x0, y0 = points[0][1]
    values, times = ["0 0"], [0.0]
    for (time, (x, y)), (_, (px, py)) in zip(points[1:], points):
        values += [f"{px - x0:.1f} {py - y0:.1f}", f"{x - x0:.1f} {y - y0:.1f}"]
        times += [time - MOVE, time]
    values.append(values[-1]); times.append(total)
    return ('<animateTransform attributeName="transform" type="translate" dur="%gs" repeatCount="indefinite" values="%s" keyTimes="%s"/>'
            % (total, ";".join(values), ";".join(f"{t / total:.4f}" for t in times)))


def main(folder):
    sites = [(row, column) for row in range(ROWS) for column in range(COLUMNS) if (row, column) != START]
    at = {site: site for site in sites}                    # where every atom is, by the site of its start
    atom_at = {site: site for site in sites}               # which atom every site holds
    track = {site: [(0.0, place(site))] for site in sites}
    vacancy = [(0.0, place(START))]
    found = jumps()
    total = len(found) * (MOVE + REST)
    for k, (source, target) in enumerate(found):
        time = k * (MOVE + REST) + REST + MOVE              # the jump ends then
        atom = atom_at.pop(source)
        atom_at[target], at[atom] = atom, target
        track[atom].append((time, place(target)))
        vacancy.append((time, place(source)))
    kind = lambda site: SOLUTE.get(site, "Al")
    assert all(kind(atom) == kind(site) for site, atom in atom_at.items()), "the round has to end as it begins"
    out = ['<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 %d %d" width="%d" height="%d" role="img">' % (WIDTH, HEIGHT, WIDTH, HEIGHT),
           "<title>A vacancy jumps through a close-packed layer of Al, Mg and Zn atoms</title>", "<defs>"]
    for name, (light, middle, dark) in BALL.items():
        out.append('<radialGradient id="%s" cx="0.36" cy="0.30" r="0.74"><stop offset="0" stop-color="%s"/><stop offset="0.42" stop-color="%s"/>'
                   '<stop offset="1" stop-color="%s"/></radialGradient>' % (name, light, middle, dark))
    out += ['<linearGradient id="ends" x1="0" x2="1" y1="0" y2="0"><stop offset="0" stop-color="#fff" stop-opacity="0"/>'
            '<stop offset="0.16" stop-color="#fff"/><stop offset="0.84" stop-color="#fff"/><stop offset="1" stop-color="#fff" stop-opacity="0"/>'
            "</linearGradient>", '<mask id="fade"><rect width="%d" height="%d" fill="url(#ends)"/></mask>' % (WIDTH, HEIGHT), "</defs>",
            '<g mask="url(#fade)">']
    x, y = place(START)                                    # the vacancy lies under the atoms: an atom that jumps passes over its ring
    out.append('<circle cx="%.1f" cy="%.1f" r="%.1f" fill="#fff" stroke="%s" stroke-width="2.8">%s</circle>'
               % (x, y, RADIUS - 1.9, RED, moving(vacancy, total)))
    for site in sorted(sites, key=lambda s: len(track[s])):    # the atoms that jump are drawn last, on top
        x, y = place(site)
        circle = '<circle cx="%.1f" cy="%.1f" r="%.1f" fill="url(#%s)"' % (x, y, RADIUS, kind(site))
        out.append(circle + "/>" if len(track[site]) == 1 else circle + ">" + moving(track[site], total) + "</circle>")
    out += ["</g>", "</svg>"]
    with open(os.path.join(folder, "banner.svg"), "w", encoding="utf-8") as f:
        f.write("\n".join(out) + "\n")
    print(f"{len(found)} jumps, {total:g} s, {len(sites)} atoms")


if __name__ == "__main__":
    main(sys.argv[1])
