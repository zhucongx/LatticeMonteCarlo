# Model of the animation second_order: one vacancy on a triangular lattice (6 neighbor sites, a sketch of the 12 of fcc) with rates that
# do not change with time: the vacancy exchanges places with a Mg atom on a neighbor site 100 times faster than with an Al atom.
# first_order() does one step of KineticMcFirstSerial, second_order() is the calculation of a second-order step (the GetSecondOrder
# functions of KineticMcSecondAbstract), with the names used there.
# Compare the formulas of the second order with many first-order runs:
#   python3 second_order_model.py [number of runs, default 100000]
import math
import random
import sys

RATE = {"Al": 1.0, "Mg": 100.0}                            # 1/microsecond
NEIGHBORS = ((1, 0), (0, 1), (-1, 1), (-1, 0), (0, -1), (1, -1))
# the example: the vacancy has come from J to K; all sites but those of MG hold Al
J, K, I, L, M = (-1, 0), (0, 0), (1, 0), (1, 1), (2, 1)
MG = frozenset({I, (-2, 2), (3, -1)})


def neighbors(v):
    return [(v[0] + a, v[1] + b) for a, b in NEIGHBORS]


def rates(v, mg):
    return [RATE["Mg" if i in mg else "Al"] for i in neighbors(v)]


def exchange(v, i, mg):
    """the Mg atoms after the vacancy has jumped from v to i"""
    return (mg - {i}) | {v} if i in mg else mg


def select(sites, probabilities, rng):
    """as SelectEvent: the first event whose cumulative probability reaches the random number"""
    x, cumulative = rng.random(), 0.0
    for i, p in zip(sites, probabilities):
        cumulative += p
        if cumulative >= x:
            break
    return i


def first_order(k, mg, rng):
    """time of the step and the site to which the vacancy jumps"""
    r = rates(k, mg)
    total = sum(r)
    dt = -math.log(rng.random()) / total
    return dt, select(neighbors(k), [x / total for x in r], rng)


def second_order(k, j, mg):
    """probabilities that the step ends on the neighbors of k, and the mean time of the step; the vacancy has come from j"""
    total_rate_k = sum(rates(k, mg))
    t_1 = 1 / total_rate_k
    beta_bar_k = beta_k = gamma_bar_k_j = gamma_k_j = beta_k_j = alpha_k_j = ts_numerator = ts_j_numerator = 0.0
    events = []
    for i, rate_k_i in zip(neighbors(k), rates(k, mg)):
        total_rate_i = sum(rates(i, exchange(k, i, mg)))
        probability_k_i = rate_k_i / total_rate_k
        probability_i_k = rate_k_i / total_rate_i           # the same atom jumps back: the same rate in this model
        beta_bar_k_i = probability_k_i * probability_i_k
        beta_k_i = probability_k_i * (1 - probability_i_k)
        is_previous_event = i == j
        t_i = 1 / total_rate_i
        beta_bar_k += beta_bar_k_i
        beta_k += beta_k_i
        gamma_bar_k_j += 0 if is_previous_event else beta_bar_k_i
        gamma_k_j += 0 if is_previous_event else beta_k_i
        beta_k_j += beta_k_i if is_previous_event else 0
        alpha_k_j += probability_k_i if is_previous_event else 0
        ts_numerator += (t_1 + t_i) * beta_bar_k_i
        ts_j_numerator += 0 if is_previous_event else (t_1 + t_i) * beta_bar_k_i
        events.append((is_previous_event, beta_k_i))
    ts = ts_numerator / beta_bar_k
    ts_j = ts_j_numerator / gamma_bar_k_j
    one_over_one_minus_a_j = 1 / (1 - alpha_k_j)
    t_2 = one_over_one_minus_a_j * (gamma_k_j * t_1 + gamma_bar_k_j * (ts_j + t_1 + beta_bar_k / beta_k * ts))
    return [one_over_one_minus_a_j * (gamma_bar_k_j / beta_k) * beta_k_j if is_previous_event else
            one_over_one_minus_a_j * (1 + gamma_bar_k_j / beta_k) * beta_k_i for is_previous_event, beta_k_i in events], t_2


def returns_end(k, j, mg, rng):
    """the same step jump by jump: the site i that the vacancy does not leave to k at once, the time of its arrival there, the jumps"""
    t, jumps = 0.0, 0
    while True:
        dt, i = first_order(k, mg, rng)
        if jumps == 0 and i == j:                          # then the step before would not have ended on k
            continue
        t, jumps = t + dt, jumps + 1
        dt, l = first_order(i, exchange(k, i, mg), rng)
        if l != k:
            return i, t, jumps
        t, jumps = t + dt, jumps + 1


if __name__ == "__main__":
    n = int(sys.argv[1]) if len(sys.argv) > 1 else 100000
    rng = random.Random(1)
    mg = exchange(K, I, MG)
    for k, j, atoms in ((K, J, MG), (I, K, mg), (L, I, mg)):
        probabilities, t_2 = second_order(k, j, atoms)
        ends, times, jumps = zip(*(returns_end(k, j, atoms, rng) for _ in range(n)))
        mean = sum(times) / n
        error = math.sqrt(sum((t - mean) ** 2 for t in times)) / n
        print(f"step from {k}, previous site {j}: {n} runs, {sum(jumps) / n:.2f} jumps on average")
        print(f"  mean time  {mean:.5f} +- {error:.5f}   formula {t_2:.5f}")
        for i, p in zip(neighbors(k), probabilities):
            f = ends.count(i) / n
            print(f"  ends on {str(i):8} {f:.5f} +- {math.sqrt(p * (1 - p) / n):.5f}   formula {p:.5f}")
