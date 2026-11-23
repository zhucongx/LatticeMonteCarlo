# The run of the animation trapped_vacancy, for render_trapped_vacancy.py and animate_trapped_vacancy.py. Reads the log once, finds for
# every step where the vacancy is and the time in the alloy, and writes what the two scripts need to one file:
#   python3 trapped_vacancy_data.py <run folder> <file of Ansys for the start of the KMC example> <data file>
# The run continues the KMC example of docs/README.md: its last configuration, 10000000.cfg.gz, is the start.cfg.gz of a run with the same
# parameter file, but maximum_steps 3000000 and log_dump_steps 1, so that kmc_log.txt has a row for every step. The run folder holds
# start.cfg.gz and kmc_log.txt of this run. The alloy is that of the KMC example: it keeps the vacancies that the solid solution at the
# start of the example has in equilibrium at the temperature of the run, and every time of the run is multiplied by the same factor.
# Needs pandas; takes about two minutes.
import gzip
import itertools
import os
import sys

import numpy as np
from alloy_time import factor

A = 4.046                                                  # lattice parameter in Angstrom; sites are counted in units of A / 2
MIN, BONDS = 10, 3                                         # smallest_cluster_criteria and solvent_bond_criteria of Ansys
NAMED = 4                                                  # clusters that keep their cluster_id of the start; all others are OTHER
MATRIX, OTHER = 0, NAMED + 1
PICTURES = 111                                             # equally spaced in time
NEIGHBORS = np.array([v for v in itertools.product((-1, 0, 1), repeat=3) if sum(c * c for c in v) == 2])


def read_cfg(path):
    """elements, sites and the number of lattice planes along an edge of a cubic configuration file"""
    with gzip.open(path, "rt") as f:
        lines = f.read().split("\n")
    number = int(lines[0].split("=")[1])
    length = float(next(l for l in lines if l.startswith("H0(1,1)")).split("=")[1].split()[0])
    body = [l for l in lines[[l[:11] for l in lines].index("entry_count") + 1:] if not l.startswith("auxiliary")]
    planes = round(2 * length / A)
    site = np.rint([[float(x) * planes for x in body[3 * i + 2].split()[:3]] for i in range(number)]).astype(int) % planes
    return np.array([body[3 * i + 1].strip() for i in range(number)]), site, planes


def sites(d, step):
    """site of every atom at a step: an atom is where the vacancy was when they last changed places"""
    last = np.full(len(d["species"]), -1)
    np.maximum.at(last, d["selected"][:step], np.arange(step))
    site = np.where(last[:, None] < 0, d["start"], d["path"][last] % d["planes"])
    site[d["vacancy"]] = d["path"][step] % d["planes"]
    return site


def clusters(site, species, planes):
    """the clusters of Ansys with Al as the solvent, largest first: of each one the sites of the solute atoms and the vacancy, and the
    sites of the solvent atoms that join"""
    left, found = {tuple(p) for p in site[species != "Al"]}, []
    around = lambda p: [tuple((a + b) % planes for a, b in zip(p, v)) for v in NEIGHBORS]
    while left:
        group, new = set(), {left.pop()}
        while new:
            group |= new
            new = {q for p in new for q in around(p)} & left
            left -= new
        if len(group) >= MIN:
            near = [q for p in group for q in around(p) if q not in group]
            found.append((group, {q for q in near if near.count(q) > BONDS}))
    return sorted(found, key=lambda c: -len(c[0]) - len(c[1]))


def locate(d):
    """Where the vacancy is in every row of the log: MATRIX, or at a cluster, which means on one of its sites or on a nearest neighbor
    site of one. The clusters are those of the configuration of the row. The NAMED largest clusters of the start keep their numbers:
    when solute atoms leave or join, the number stays with the largest part. Every other cluster is OTHER.
    Also returns, for every picture, the number of the cluster of each solute atom, 0 for one in a cluster without number or in none."""
    species, planes, selected, vacancy = d["species"], d["planes"], d["selected"], d["vacancy"]
    number = lambda p: (p[..., 0] % planes * planes + p[..., 1] % planes) * planes + p[..., 2] % planes
    grid = np.stack(np.meshgrid(*[np.arange(planes)] * 3, indexing="ij"), -1).reshape(-1, 3)
    around = number(grid[:, None] + NEIGHBORS).tolist()   # the twelve neighbors of every site
    path = number(d["path"]).tolist()
    atom = np.full(planes ** 3, -1)
    atom[number(d["start"])] = np.arange(len(species))
    atom = atom.tolist()
    solute = (species != "Al") & (species != "X")
    # groups of solute atoms that are connected as nearest neighbors, without the vacancy: the group of every site, the sites of every
    # group, and the groups that carry the number of a cluster
    label, members, title, count = [-1] * planes ** 3, {}, {}, 0

    def connect(some):
        """groups the solute sites `some` anew and returns the new groups"""
        nonlocal count
        before = count
        while some:
            new, group = [some.pop()], set()
            while new:
                p = new.pop()
                group.add(p)
                for q in around[p]:
                    if q in some:
                        some.discard(q)
                        new.append(q)
            count += 1
            members[count] = group
            for p in group:
                label[p] = count
        return range(before + 1, count + 1)

    def name(groups):
        return min((title[g] for g in groups if g in title), default=OTHER)

    everything = connect(set(number(d["start"][solute]).tolist()))
    for k, (cluster, _) in enumerate(clusters(d["start"], species, planes)[:NAMED], 1):
        cluster = set(number(np.array(sorted(cluster))).tolist())
        title[max(everything, key=lambda g: len(members[g] & cluster))] = k
    where, names = np.zeros(len(path), np.int8), np.zeros((len(d["shown"]), len(species)), np.int8)
    for n, x in enumerate(path):
        touched = {label[q] for q in around[x]} - {-1}
        if 1 + sum(len(members[g]) for g in touched) >= MIN:                  # the vacancy belongs to a cluster
            where[n] = name(touched)
        else:                                                                 # or is next to a solvent atom that has joined one
            for q in around[x]:
                bonds = [label[w] for w in around[q] if label[w] >= 0 and len(members[label[w]]) >= MIN] if label[q] < 0 else []
                for g in set(bonds):
                    if bonds.count(g) > BONDS:
                        where[n] = name([g])
        for i in np.flatnonzero(d["shown"] == n):
            for g, k in title.items():
                names[i, [atom[p] for p in members[g]]] = k
        if n + 1 < len(path):                                                 # the jump of this step: the atom at y comes to x
            y, jumping = path[n + 1], selected[n]
            assert atom[x] == vacancy and atom[y] == jumping and y in around[x], f"step {n}: the log does not fit the configuration"
            atom[x], atom[y] = jumping, vacancy
            if solute[jumping]:                                               # its group and those it meets at x are grouped anew
                old = {g: members.pop(g) for g in touched}
                label[y] = -1
                new = connect(set().union(*old.values()) - {y} | {x})
                for k, g in sorted((title.pop(g), g) for g in old if g in title):
                    title.setdefault(max(new, key=lambda h: len(members[h] & old[g])), k)
        if n % 500_000 == 0:
            print(f"row {n}", flush=True)
    return where, names


def prepare(run, frame, output):
    import pandas as pd
    species, start, planes = read_cfg(os.path.join(run, "start.cfg.gz"))
    log = pd.read_csv(os.path.join(run, "kmc_log.txt"), sep=r"\s+",
                      usecols=["steps", "time", "temperature", "selected", "vac1", "vac2", "vac3"])
    steps, time, temperature = (log[column].values for column in ("steps", "time", "temperature"))
    assert (steps == np.arange(len(log))).all(), "the log needs a row for every step"
    assert (temperature == temperature[0]).all(), "the run has to be at one temperature"
    d = dict(species=species, start=start, planes=planes, vacancy=int(np.flatnonzero(species == "X")[0]), time=time,
             temperature=temperature[0], selected=log["selected"].values.astype(np.int32),
             path=np.rint(log[["vac1", "vac2", "vac3"]].values / (A / 2)).astype(np.int32))        # not wrapped back into the box
    # the time in the alloy: the time of the run times the factor of alloy_time.py
    d["factor"] = factor(frame, temperature[0])
    d["alloy"] = time * d["factor"]
    # picture i shows the row that holds at the time i / (PICTURES - 1) of the total time of the run
    d["shown"] = np.searchsorted(time, np.linspace(0, time[-1], PICTURES), side="right") - 1
    d["where"], d["names"] = locate(d)
    np.savez(output, **d)


def load(path):
    d = dict(np.load(path))
    return {**d, "planes": int(d["planes"]), "vacancy": int(d["vacancy"]), "temperature": float(d["temperature"])}


if __name__ == "__main__":
    prepare(*sys.argv[1:4])
