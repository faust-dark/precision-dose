"""Deterministic numerical primitives.

No third-party dependencies. Everything here is either exact in closed form
(matrix exponential for the small linear systems we use) or a documented
algorithm with a convergence guarantee. The matrix exponential is
scaling-and-squaring with a Taylor series, which is accurate to machine
precision for the small (<= 4x4) matrices of a compartment model; the tests
compare it against the closed-form biexponential solution.
"""

import math

# --------------------------------------------------------------------------
# linear algebra (small dense matrices, lists of lists)
# --------------------------------------------------------------------------


def mat_mul(a, b):
    n, m, p = len(a), len(b), len(b[0])
    out = [[0.0] * p for _ in range(n)]
    for i in range(n):
        ai = a[i]
        for k in range(m):
            aik = ai[k]
            if aik == 0.0:
                continue
            bk = b[k]
            for j in range(p):
                out[i][j] += aik * bk[j]
    return out


def mat_add(a, b):
    return [[a[i][j] + b[i][j] for j in range(len(a[0]))] for i in range(len(a))]


def mat_scale(a, s):
    return [[v * s for v in row] for row in a]


def identity(n):
    return [[1.0 if i == j else 0.0 for j in range(n)] for i in range(n)]


def mat_vec(a, x):
    return [sum(a[i][j] * x[j] for j in range(len(x))) for i in range(len(a))]


def norm_inf(a):
    return max(sum(abs(v) for v in row) for row in a)


def solve(a, b):
    """Gaussian elimination with partial pivoting. Raises on singular input."""
    n = len(a)
    m = [list(row) + [b[i]] for i, row in enumerate(a)]
    for col in range(n):
        pivot = max(range(col, n), key=lambda r: abs(m[r][col]))
        if abs(m[pivot][col]) < 1e-300:
            raise ZeroDivisionError("singular matrix")
        if pivot != col:
            m[col], m[pivot] = m[pivot], m[col]
        pv = m[col][col]
        for r in range(col + 1, n):
            f = m[r][col] / pv
            if f == 0.0:
                continue
            for c in range(col, n + 1):
                m[r][c] -= f * m[col][c]
    x = [0.0] * n
    for i in range(n - 1, -1, -1):
        s = m[i][n] - sum(m[i][j] * x[j] for j in range(i + 1, n))
        x[i] = s / m[i][i]
    return x


def mat_exp(m):
    """Matrix exponential via scaling-and-squaring with a Taylor series."""
    n = len(m)
    nrm = norm_inf(m)
    squarings = 0
    if nrm > 0.5:
        squarings = max(0, int(math.ceil(math.log2(nrm))) + 1)
    scale = 2.0 ** squarings
    a = mat_scale(m, 1.0 / scale)
    result = identity(n)
    term = identity(n)
    for k in range(1, 20):
        term = mat_scale(mat_mul(term, a), 1.0 / k)
        result = mat_add(result, term)
        if norm_inf(term) < 1e-18:
            break
    for _ in range(squarings):
        result = mat_mul(result, result)
    return result


def cholesky(a):
    """Lower-triangular Cholesky factor of a small SPD matrix."""
    n = len(a)
    low = [[0.0] * n for _ in range(n)]
    for i in range(n):
        for j in range(i + 1):
            s = a[i][j] - sum(low[i][k] * low[j][k] for k in range(j))
            if i == j:
                if s <= 0.0:
                    raise ValueError("matrix is not positive definite")
                low[i][j] = math.sqrt(s)
            else:
                low[i][j] = s / low[j][j]
    return low


# --------------------------------------------------------------------------
# deterministic normal sampling
# --------------------------------------------------------------------------

# Wichura AS241 inverse normal CDF. Used instead of random.gauss so that the
# same seed produces the same draws on any Python build (random.gauss has
# changed implementation across versions).
_A = (
    3.3871328727963666080, 1.3314166789178437745e2, 1.9715909503065514427e3,
    1.3731693765509461125e4, 4.5921953931549871457e4, 6.7265770927008700853e4,
    3.3430575583588128105e4, 2.5090809287301226727e3,
)
_B = (
    4.2313330701600911252e1, 6.8718700749205790830e2, 5.3941960214247511077e3,
    2.1213794301586595867e4, 3.9307895800092710610e4, 2.8729085735721942674e4,
    5.2264952788528545610e3,
)
_C = (
    1.42343711074968357734, 4.63033784615654529590, 5.76949722146069140550,
    3.64784832476320460504, 1.27045825245236838258, 2.41780725177450611770e-1,
    2.27238449892691845833e-2, 7.74545014278341407640e-4,
)
_D = (
    2.05319162663775882187, 1.67638483018380384940, 6.89767334985100004550e-1,
    1.48103976427480074590e-1, 1.51986665636164571966e-2, 5.47593808499534494600e-4,
    1.05075007164441684324e-9,
)
_E = (
    6.65790464350110377720, 5.46378491116411436990, 1.78482653991729133580,
    2.96560571828504891230e-1, 2.65321895265761230930e-2, 1.24266094738807843860e-3,
    2.71155556874348757815e-5, 1.00725088747240825720e-7,
)
_F = (
    5.99832206555887937690e-1, 1.36929880922735805310e-1, 1.48753612908506148525e-2,
    7.86869131145613259100e-4, 1.84631831751005468180e-5, 1.42151175831644588870e-7,
    2.04426310338993978564e-15,
)


def _poly(coefs, r):
    """Horner over coefficients given in ASCENDING power order (c0 + c1 r + ...)."""
    acc = 0.0
    for c in reversed(coefs):
        acc = acc * r + c
    return acc


def norm_ppf(p):
    """Inverse standard normal CDF, Wichura AS241 (accuracy ~1e-15)."""
    if not (0.0 < p < 1.0):
        raise ValueError("p must be in (0,1)")
    q = p - 0.5
    if abs(q) <= 0.425:
        r = 0.180625 - q * q
        num = _poly(_A, r) * q
        den = _poly(_B, r) * r + 1.0
        return num / den
    r = p if q < 0 else 1.0 - p
    r = math.sqrt(-math.log(r))
    if r <= 5.0:
        r -= 1.6
        val = _poly(_C, r) / (_poly(_D, r) * r + 1.0)
    else:
        r -= 5.0
        val = _poly(_E, r) / (_poly(_F, r) * r + 1.0)
    return -val if q < 0 else val


def norm_cdf(x):
    return 0.5 * math.erfc(-x / math.sqrt(2.0))


class DeterministicRNG:
    """MT19937 (Python's Random) + AS241 inverse transform.

    Draws depend only on (seed, draw index) in a documented order, so the same
    seed gives the same population on any build of Python 3.
    """

    def __init__(self, seed):
        import random

        self._rng = random.Random(seed)
        self.count = 0

    def uniform(self):
        self.count += 1
        # 53-bit double in (0,1): identical across CPython versions.
        return (self._rng.getrandbits(53) + 0.5) / 9007199254740992.0

    def normal(self):
        return norm_ppf(self.uniform())

    def standard_normals(self, n):
        return [self.normal() for _ in range(n)]


# --------------------------------------------------------------------------
# deterministic optimiser (Nelder-Mead)
# --------------------------------------------------------------------------


def nelder_mead(f, x0, step=0.3, max_iter=2000, x_tol=1e-9, f_tol=1e-12):
    """Deterministic Nelder-Mead simplex minimiser.

    Returns (x_best, f_best, info). info carries iterations, simplex diameter
    and a converged flag, which callers are expected to gate on.
    """
    n = len(x0)
    simplex = [list(x0)]
    for i in range(n):
        p = list(x0)
        p[i] += step
        simplex.append(p)
    values = [f(p) for p in simplex]

    def order():
        pairs = sorted(zip(values, range(len(values))), key=lambda t: (t[0], t[1]))
        new_s, new_v = [], []
        for v, i in pairs:
            new_s.append(simplex[i])
            new_v.append(v)
        simplex[:] = new_s
        values[:] = new_v

    order()
    it = 0
    for it in range(1, max_iter + 1):
        order()
        diameter = max(
            math.sqrt(sum((p[j] - simplex[0][j]) ** 2 for j in range(n))) for p in simplex
        )
        if diameter < x_tol and abs(values[-1] - values[0]) < f_tol:
            return simplex[0], values[0], {
                "iterations": it,
                "simplex_diameter": diameter,
                "converged": True,
            }
        centroid = [sum(p[j] for p in simplex[:-1]) / n for j in range(n)]
        worst = simplex[-1]
        xr = [centroid[j] + 1.0 * (centroid[j] - worst[j]) for j in range(n)]
        fr = f(xr)
        if fr < values[0]:
            xe = [centroid[j] + 2.0 * (centroid[j] - worst[j]) for j in range(n)]
            fe = f(xe)
            if fe < fr:
                simplex[-1], values[-1] = xe, fe
            else:
                simplex[-1], values[-1] = xr, fr
        elif fr < values[-2]:
            simplex[-1], values[-1] = xr, fr
        else:
            xc = [centroid[j] + 0.5 * (worst[j] - centroid[j]) for j in range(n)]
            fc = f(xc)
            if fc < values[-1]:
                simplex[-1], values[-1] = xc, fc
            else:
                for i in range(1, len(simplex)):
                    simplex[i] = [
                        simplex[0][j] + 0.5 * (simplex[i][j] - simplex[0][j])
                        for j in range(n)
                    ]
                    values[i] = f(simplex[i])
    diameter = max(
        math.sqrt(sum((p[j] - simplex[0][j]) ** 2 for j in range(n))) for p in simplex
    )
    return simplex[0], values[0], {
        "iterations": it,
        "simplex_diameter": diameter,
        "converged": False,
    }
