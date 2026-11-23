"""The sites of the animations cluster_definition and barrier_groups, on a triangular lattice with 6 nearest neighbors instead of 12.
Imported by animate_cluster_definition.py and animate_barrier_groups.py; not started on its own.

    cluster_definition()   what is connected, dropped and added, by the rules of SoluteCluster::FindAtomListOfClusters
                           (lmc/ansys/src/SoluteCluster.cpp)
    barrier_groups()       cluster 1 of the first animation, its nearest neighbor sites, and some jumps. Their groups are checked by
                           the rules of ExitTime::GetAverageBarriers (lmc/ansys/src/ExitTime.cpp).
"""

# The atoms: "." solvent (Al), M Mg, Z Zn, X vacancy. Every other row is moved to the right by half a site.
ROWS = """
. . . . . . . . .
 . . M Z . Z M .
. . Z X . . M . .
 . M . M . . . .
. . Z . . . Z . .
 . . . . . M Z .
. . . M . . . M .
 . . . . . . . .
"""
SMALLEST_CLUSTER_CRITERIA = 4
SOLVENT_BOND_CRITERIA = 2
SOLVENT = "."
ATOMS = {(r, c): e for r, row in enumerate(ROWS.strip("\n").split("\n")) for c, e in enumerate(row.split())}

# Jumps of the vacancy, as pairs of sites, for the second drawing
JUMPS = {"in": [((2, 2), (2, 3)), ((3, 2), (4, 2))],
         "to": [((2, 4), (1, 4)), ((3, 3), (4, 3))],
         "on": [((1, 1), (2, 1)), ((5, 1), (5, 2))],
         "off": [((0, 2), (-1, 1)), ((3, 4), (4, 5))]}


def element(site):
    return ATOMS.get(site, SOLVENT)


def neighbors(site):
    r, c = site
    left = c - 1 + r % 2
    return [(r, c - 1), (r, c + 1), (r - 1, left), (r - 1, left + 1), (r + 1, left), (r + 1, left + 1)]


def around(sites):
    """The sites and all their nearest neighbors."""
    return set(sites) | {n for s in sites for n in neighbors(s)}


# ---------------- rules of the program ----------------
def connected_groups(sites):
    """FindAtomListOfClustersBFSHelper: the groups of sites that are connected by nearest neighbors."""
    unvisited, groups = set(sites), []
    while unvisited:
        queue, group = [unvisited.pop()], []
        while queue:
            site = queue.pop(0)
            group.append(site)
            for n in neighbors(site):
                if n in unvisited:
                    unvisited.remove(n)
                    queue.append(n)
        groups.append(group)
    return groups


def bonds(site, group):
    """Nearest neighbors of the site that belong to the group and are not solvent atoms."""
    return sum(n in group and element(n) != SOLVENT for n in neighbors(site))


def find_clusters():
    """FindAtomListOfClusters: the dropped groups, the other groups, and the clusters they become, the largest one first."""
    groups = connected_groups([s for s in ATOMS if element(s) != SOLVENT])
    dropped = [g for g in groups if len(g) < SMALLEST_CLUSTER_CRITERIA]
    kept = [g for g in groups if g not in dropped]
    clusters = [set(g) for g in kept]
    for cluster in clusters:
        if len(cluster) <= SOLVENT_BOND_CRITERIA:
            continue
        for site in around(cluster):
            if bonds(site, cluster) > SOLVENT_BOND_CRITERIA:
                cluster.add(site)
    return dropped, kept, sorted(clusters, key=len, reverse=True)


def group_of_jump(first, second, cluster):
    """GetAverageBarriers: the group of the jump of the vacancy from the first to the second site, or None."""
    near = around(cluster)
    if first not in near:
        return None
    if first in cluster and second in cluster:
        return "in"
    if first in cluster and second not in cluster and second in near:
        return "to"
    if first not in cluster and second not in cluster and second in near:
        return "on"
    if first not in cluster and second not in near:
        return "off"
    return None


# ---------------- output: \newcommand{\name}{x/y, x/y, ...} for \foreach ----------------
def point(site):
    return "%.3f/%.3f" % (site[1] + 0.5 * (site[0] % 2), -site[0] * 3 ** 0.5 / 2)


def command(name, entries):
    return "\\newcommand{\\%s}{%s}\n" % (name, ", ".join(entries))


def species(prefix, sites, names):
    return "".join(command(prefix + names[e], [point(s) for s in sorted(sites) if element(s) == e]) for e in names)


def links(name, sites):
    """All pairs of nearest neighbors among the sites."""
    return command(name, [point(s) + "/" + point(n) for s in sorted(sites) for n in neighbors(s) if n in sites and s < n])


def cluster_definition():
    dropped, kept, clusters = find_clusters()
    counts = {}                          # the solvent atoms next to a group: their nearest neighbors in it
    for group in kept:
        for s in around(group) - set(group):
            assert s not in counts and s in ATOMS, "next to two groups, or outside of the drawing"
            counts[s] = bonds(s, group)
    dropped, kept, clusters = (set().union(*lists) for lists in (dropped, kept, clusters))
    text = species("sites", ATOMS, {SOLVENT: "Al"})
    text += species("kept", kept, {"M": "Mg", "Z": "Zn", "X": "X"})
    text += species("dropped", dropped, {"M": "Mg", "Z": "Zn"})
    text += species("joined", clusters, {SOLVENT: "Al"})
    text += links("bondsKept", kept) + links("bondsDropped", dropped) + links("bondsClusters", clusters)
    return text + command("counts", ["%s/%d" % (point(s), n) for s, n in sorted(counts.items())])


def barrier_groups():
    cluster = find_clusters()[2][0]
    near = around(cluster)
    patch = {(r, c) for r in range(-1, 7) for c in range(-(r % 2), 7)}
    assert around(near) <= patch, "a neighbor outside of the drawing"
    text = ""
    for name, jumps in JUMPS.items():
        for first, second in jumps:
            assert second in neighbors(first) and group_of_jump(first, second, cluster) == name, (name, first, second)
            if name in ("in", "on"):     # drawn with two heads: the jump back belongs to the same group
                assert group_of_jump(second, first, cluster) == name
        text += command("jumps" + name.capitalize(), [point(a) + "/" + point(b) for a, b in jumps])
    text += links("bondsCluster", cluster) + links("bondsNear", near)
    text += species("cluster", cluster, {"M": "Mg", "Z": "Zn"})
    # GetAverageBarriers first puts a solvent atom on the site of the vacancy
    return text + command("sitesAl", [point(s) for s in sorted(patch) if s not in cluster or element(s) in (SOLVENT, "X")])
