# The time in an alloy that belongs to the time of a KMC run (page Methods, "Time and vacancy concentration").
# The box of a run holds one vacancy in all its sites. The alloy of the animations holds the vacancies that the solid solution at the
# start of the run has in equilibrium at a temperature T, and keeps their number while clusters form, as the box keeps its one. Clusters
# then trap the same share of the vacancies in both, and a time of the run times a constant is the time in the alloy:
#   python3 alloy_time.py <frame> <T in K>
# prints this constant. The frame is the file that Ansys wrote for the start of the run, with the field profile_energy of every site.
# The animations and plots.py use the functions.
import gzip
import re
import sys

import numpy as np

kB = 8.617333262145e-5                                     # eV/K


def fields(frame, *names):
    """ the fields names of all sites in a file of Ansys, as text """
    with gzip.open(frame, "rt") as f:
        f.readline()
        parts = re.search(r"Properties=(\S+)", f.readline()).group(1).split(":")
        column = {name: sum(int(width) for width in parts[2:3 * i:3]) for i, name in enumerate(parts[::3])}
        rows = [line.split() for line in f]
    return [np.array([row[column[name]] for row in rows]) for name in names]


def boltzmann(energy, temp):
    return np.exp(-energy[:, None] / (kB * np.atleast_1d(temp)))


def c_equilibrium(energy, temp):
    """ all vacancies in equilibrium: the concentration in pure Al times the mean of exp(-E/kT) over the sites """
    return 1.67 * np.exp(0.7 - 0.66 / (kB * temp)) * boltzmann(energy, temp).mean(axis=0)


def factor(frame, temp):
    """ the time in the alloy divided by the time of the run: one vacancy in the box, against those of the frame in equilibrium at temp """
    energy = fields(frame, "profile_energy")[0].astype(float)
    return 1 / (len(energy) * c_equilibrium(energy, float(temp))[0])


if __name__ == "__main__":
    print(f"{factor(*sys.argv[1:3]):.4g}")
