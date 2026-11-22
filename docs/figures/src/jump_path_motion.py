# The motion of the animation jump_path, for render_jump_path.py and animate_jump_path.py.
# The path is E(x) = 16 x^4 - 4.8 x^3 - 7.595 x^2, that is D = 1, K_s = 64 and Delta E = -1.2 in units in which K_s D^2 / 64 = 1:
# minima I and F, transition state T at x = 0.
X_I, X_T, X_F = -0.3875, 0.0, 0.6125
STEPS, STEP_T = 56, 26                                     # picture 0 is at I, picture STEP_T at T, picture STEPS at F


def energy(x):
    return 16 * x ** 4 - 4.8 * x ** 3 - 7.595 * x ** 2


def place(k, slope_t=0.011):
    """x of picture k: at rest at I and F, slow at T; two cubic polynomials that meet at T with the same slope"""
    if k <= STEP_T:
        x0, x1, v0, v1, n, w = X_I, X_T, 0.0, slope_t, STEP_T, k / STEP_T
    else:
        x0, x1, v0, v1, n, w = X_T, X_F, slope_t, 0.0, STEPS - STEP_T, (k - STEP_T) / (STEPS - STEP_T)
    return x0 * (1 + 2 * w) * (1 - w) ** 2 + x1 * w * w * (3 - 2 * w) + n * (v0 * w * (1 - w) ** 2 - v1 * w * w * (1 - w))


def fraction(k):
    """place of the atom on the way from its site (0) to the vacancy (1)"""
    return (place(k) - X_I) / (X_F - X_I)
