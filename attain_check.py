"""Theorem 6.6 (attainment): for a general-position technology and admissible rational
capacities, the eliminant of the Nash / Leontief-Fisher critical system is irreducible
over Q of degree exactly C(m-1, l-1).

The critical map  c_j = Psi_j(lambda) = sum_T Abar[j][T] W_T / (Abar^T lambda)_T
is homogeneous of degree -1 in lambda, so the ratios c_j/c_1 depend only on the
projective class of lambda.  Dehomogenizing (lambda_1 = 1) turns the system into l-1
equations in l-1 unknowns, which keeps the elimination small.

Instances: Vandermonde aggregated technology on distinct positive nodes, so every l
columns are independent (general position); positive rational weights; and capacities
built from an admissible multiplier vector, so an interior solution serving all m types
exists by construction.
"""
from math import comb
from sympy import (Rational, symbols, Poly, QQ, factor_list, prod, resultant, expand,
                   nsimplify, nroots, Float)

# ---------------------------------------------------------------- instances


def build(m, l, nodes, weights, lam0, shift=None):
    Abar = [[Rational(r) ** j for r in nodes] for j in range(l)]
    W = [Rational(w) for w in weights]
    lam0 = [Rational(x) for x in lam0]
    den = [sum(Abar[j][t] * lam0[j] for j in range(l)) for t in range(m)]
    assert all(d > 0 for d in den), "multiplier not admissible"
    X = [W[t] / den[t] for t in range(m)]
    c = [sum(Abar[j][t] * X[t] for t in range(m)) for j in range(l)]
    if shift is not None:
        # move the capacities off the thin set of vectors with a rational multiplier
        c = [c[j] + Rational(shift[j], 10**6) for j in range(l)]
    return Abar, W, c, lam0


def ratio_equations(Abar, W, c, m, l, mu):
    """l-1 polynomial equations in mu = (mu_2,...,mu_l), from c_j/c_1 = Psi_j/Psi_1."""
    lam = [Rational(1)] + list(mu)
    den = [sum(Abar[j][t] * lam[j] for j in range(l)) for t in range(m)]
    # Psi_j * prod(den) is a polynomial
    P = [sum(Abar[j][t] * W[t] * prod(den[s] for s in range(m) if s != t)
             for t in range(m)) for j in range(l)]
    return [expand(c[0] * P[j] - c[j] * P[0]) for j in range(1, l)]


# ---------------------------------------------------------------- checks


def check(m, l, nodes, weights, lam0, shift, verbose=True):
    beta = comb(m - 1, l - 1)
    Abar, W, c, lam0r = build(m, l, nodes, weights, lam0, shift)
    mu = symbols(f"m2:{l + 1}")
    eqs = ratio_equations(Abar, W, c, m, l, mu)
    target = [float(x / lam0r[0]) for x in lam0r[1:]]   # unperturbed solution, dehomogenized

    if l == 2:
        p = Poly(eqs[0], mu[0], domain=QQ)
    else:
        p = eqs[0]
        for k in range(len(mu) - 1):
            p = resultant(p, eqs[k + 1], mu[k])
        p = Poly(expand(p), mu[-1], domain=QQ)

    facs = [(f, e) for (f, e) in factor_list(p.as_expr())[1]
            if Poly(f, p.gens[0], domain=QQ).degree() > 0]
    # the economic solution is the real root near the unperturbed multiplier ratio
    best, bestgap, bestexp = None, None, None
    for f, e in facs:
        for r in nroots(Poly(f, p.gens[0], domain=QQ), n=25):
            if abs(complex(r).imag) < 1e-18:
                gap = abs(float(complex(r).real) - target[-1])
                if bestgap is None or gap < bestgap:
                    best, bestgap, bestexp = f, gap, e
    deg = Poly(best, p.gens[0], domain=QQ).degree() if best is not None else -1
    served = True
    if l == 2 and best is not None:
        root = min((r for r in nroots(Poly(best, p.gens[0], domain=QQ), n=25)
                    if abs(complex(r).imag) < 1e-18),
                   key=lambda r: abs(float(complex(r).real) - target[-1]))
        rr = float(complex(root).real)
        served = all(float(Abar[0][t] + Abar[1][t] * rr) > 0 for t in range(m))
    ok = (deg == beta) and bestexp == 1 and served and \
        bestgap is not None and bestgap < 1e-3
    if verbose:
        tag = "PASS" if ok else "FAIL"
        print(f"  {tag}  m={m} l={l}: C({m-1},{l-1})={beta}, minimal polynomial of the "
              f"solution has degree {deg}, simple factor={bestexp == 1}, "
              f"all types served={served}, root gap={bestgap:.2e}")
    return ok


CASES = [
    # m, l, Vandermonde nodes, weights, base multiplier, capacity shift (thousandths)
    (3, 2, [1, 2, 3], [1, 1, 1], [1, 1], [1, 0]),   # shifts are millionths
    (4, 2, [1, 2, 3, 5], [1, 1, 1, 1], [1, 1], [1, 0]),
    (4, 2, [1, 3, 7, 11], [2, 3, 5, 7], [3, 1], [0, 1]),
    (5, 2, [1, 2, 3, 4, 6], [1, 1, 1, 1, 1], [1, 1], [1, 0]),
    (6, 2, [1, 2, 3, 4, 5, 7], [1, 1, 1, 1, 1, 1], [1, 1], [2, 1]),
    (7, 2, [1, 2, 3, 4, 5, 6, 8], [1, 1, 1, 1, 1, 1, 1], [1, 1], [1, 0]),
    (4, 3, [1, 2, 3, 5], [1, 1, 1, 1], [1, 1, 1], [1, 0, 0]),
    (5, 3, [1, 2, 3, 4, 6], [1, 1, 1, 1, 1], [1, 1, 1], [1, 0, 0]),
]

if __name__ == "__main__":
    print("Attainment of the bound C(m-1, l-1): random general-position instances")
    passed = failed = 0
    for m, l, nodes, w, lam0, shift in CASES:
        try:
            ok = check(m, l, nodes, w, lam0, shift)
        except Exception as exc:                                    # noqa: BLE001
            print(f"  FAIL  m={m} l={l}: {type(exc).__name__}: {exc}")
            ok = False
        passed, failed = passed + bool(ok), failed + (not ok)
    print(f"{passed} passed, {failed} failed")
