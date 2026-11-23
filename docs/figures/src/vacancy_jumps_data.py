# The steps shown by the animation vacancy_jumps, for render_vacancy_jumps.py and animate_vacancy_jumps.py.
# The run is the KMC example of docs/README.md with log_dump_steps 1 and maximum_steps 200000. Its folder holds kmc_log.txt and the start
# configuration end_of_cmc.cfg.gz, which is end.cfg.gz of the CMC example.
import gzip
import math
import os

import numpy as np

FIRST, STEPS = 1446, 26                                    # the first step and the number of steps
SUB = 10                                                   # pictures per step: one at rest and nine during the jump
H = 4.046 / 2                                              # unit of the integer coordinates: half a lattice parameter, in Angstrom


def progress(k):
    """how far the jump has come in picture k of a step, with a smooth start and end"""
    return 0.5 * (1 - math.cos(math.pi * k / SUB))


def load(run):
    """the configuration at step FIRST and the rows FIRST ... FIRST + STEPS of the log"""
    with gzip.open(os.path.join(run, "end_of_cmc.cfg.gz"), "rt") as f:
        lines = f.read().split("\n")
    number = int(lines[0].split("=")[1])
    length = float(next(l for l in lines if l.startswith("H0(1,1)")).split("=")[1].split()[0])
    body = [l for l in lines[[l[:11] for l in lines].index("entry_count") + 1:] if not l.startswith("auxiliary")]
    species = [body[3 * i + 1].strip() for i in range(number)]
    period = round(length / H)                             # sites along an edge of the box
    site = np.rint([[float(x) * period for x in body[3 * i + 2].split()[:3]] for i in range(number)]).astype(int) % period
    with open(os.path.join(run, "kmc_log.txt")) as f:
        names = f.readline().split()
        log = np.array([[float(x) for x in f.readline().split()] for _ in range(FIRST + STEPS + 1)])
    column = dict(zip(names, log.T))
    assert (column["steps"] == np.arange(len(log))).all()
    selected = column["selected"].astype(int)
    path = np.rint(log[:, -3:] / H).astype(int)            # vac1 vac2 vac3, not wrapped back into the box
    # The configuration follows from the start by exchanging the vacancy and the atom `selected` in every step.
    vacancy = species.index("X")
    for n in range(len(log)):
        if n == FIRST:
            start = site.copy()
        assert ((site[vacancy] - path[n]) % period == 0).all(), f"step {n}: the vacancy is not where the log says"
        jump = (site[selected[n]] - site[vacancy] + period // 2) % period - period // 2
        assert (jump ** 2).sum() == 2, f"step {n}: not a jump to a first neighbor"
        site[[vacancy, selected[n]]] = site[[selected[n], vacancy]]
    w = slice(FIRST, FIRST + STEPS + 1)
    return dict(species=species, period=period, site=start, vacancy=vacancy, path=path[w], selected=selected[w], step=column["steps"][w].astype(int),
                time=column["time"][w], barrier=column["Ea"][w], temperature=column["temperature"][FIRST])
