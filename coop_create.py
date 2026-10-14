"""Checks for the two additions that give Sections 3 and 8 their own results.

(1) Section 3, the boundary inside cooperative theory.  A four-player TU game with
    integer coalition values whose CORE is a polytope over Q, but whose Nash bargaining
    solution on that core has every coordinate of degree three over Q.  Three constraints
    bind (efficiency and two coalitions), the binding incidence columns take m = 4
    distinct values, and C(m-1, l-1) = C(3,2) = 3, so the bound of Theorem 6.2 is met.

(2) Section 8, general-l comparative statics.  With q_i = (A^T p)_i, u_i = w_i / q_i and
    G = sum_i (w_i / q_i^2) a_i a_i^T, the claim is
        du_i/dc_k = (w_i / q_i^2) * (G^{-1} a_i)_k,
        sum_i a_{ji} du_i/dc_k = delta_{jk},
    so agent i gains from an expansion of resource k if and only if (G^{-1} a_i)_k > 0.
    Checked against exact finite differences of the solved equilibrium.
"""
import itertools
from sympy import (Rational, symbols, Poly, QQ, real_roots, minimal_polynomial, nsimplify,
                   Matrix, diff, solve, simplify, S)

PASS, FAIL = [], []


def check(name, ok):
    (PASS if ok else FAIL).append(name)
    print(f"  {'PASS' if ok else 'FAIL'}  {name}")


# --------------------------------------------------------------- (1) core Nash
def core_nash():
    n = 4
    V = Rational(12)
    v = {(0,): 0, (1,): 0, (2,): 0, (3,): 0,
         (0, 1): 8, (1, 2): 7, (0, 2): 0, (0, 3): 0, (1, 3): 0, (2, 3): 0,
         (0, 1, 2): 8, (0, 1, 3): 8, (1, 2, 3): 7, (0, 2, 3): 0}
    v = {k: Rational(x) for k, x in v.items()}

    # superadditivity on disjoint pairs
    sup = True
    allS = [S_ for k in range(1, n) for S_ in itertools.combinations(range(n), k)]
    for A in allS:
        for B in allS:
            if set(A) & set(B):
                continue
            U = tuple(sorted(set(A) | set(B)))
            vU = V if len(U) == n else v[U]
            if vU < v[A] + v[B]:
                sup = False
    check("core-Nash game is superadditive", sup)

    # the two binding coalitions are {0,1} and {1,2}; parametrize the face
    s = symbols("s")
    x = [8 - s, s, 7 - s, s - 3]
    obj_num = Poly(4 * s**3 - 54 * s**2 + 202 * s - 168, s, domain=QQ)
    # first-order condition
    foc = simplify(sum(diff(xi, s) / xi for xi in x))
    check("first-order condition matches the cubic",
          simplify(foc - obj_num.as_expr() / (s * (8 - s) * (7 - s) * (s - 3))) == 0)
    check("cubic is irreducible over Q", obj_num.is_irreducible)

    roots = [r for r in real_roots(obj_num) if 3 < r < 7]
    check("exactly one feasible root", len(roots) == 1)
    r = roots[0]
    xs = [8 - r, r, 7 - r, r - 3]
    check("all four coordinates strictly positive", all(xi > 0 for xi in xs))

    # in the core: every coalition gets at least its value
    incore = all(sum(xs[i] for i in A) >= v[A] for A in allS) and simplify(sum(xs) - V) == 0
    check("solution lies in the core", incore)

    # every coordinate has degree three
    t = symbols("t")
    degs = [Poly(minimal_polynomial(xi, t), t).degree() for xi in xs]
    check("every coordinate has degree 3 = C(3,2)", degs == [3, 3, 3, 3])

    # the bound of Theorem 6.2: l = 3 binding rows, m = 4 distinct incidence columns
    cols = {(1, 1, 0), (1, 1, 1), (1, 0, 1), (1, 0, 0)}
    from math import comb
    check("m = 4 distinct columns, l = 3, bound C(3,2) = 3",
          len(cols) == 4 and comb(len(cols) - 1, 3 - 1) == 3)

    # the Nash value beats the equal split and the vertices of the face
    prod = lambda y: y[0] * y[1] * y[2] * y[3]
    check("Nash product exceeds the value at the face endpoints",
          all(prod(xs) > prod([8 - e, e, 7 - e, e - 3])
              for e in [Rational(31, 10), Rational(69, 10)]))


# --------------------------------------------------- (2) general-l comparative statics
def comparative_statics():
    # three resources, five agents, rational technology and capacities
    A = Matrix([[1, 1, 1, 1, 1],
                [1, 2, 3, 4, 6],
                [1, 4, 9, 16, 36]])
    w = Matrix([1, 1, 1, 1, 1])
    lam = Matrix([Rational(3), Rational(1), Rational(1, 2)])       # admissible multipliers
    q = (A.T * lam)
    u = Matrix([w[i] / q[i] for i in range(5)])
    c = A * u

    G = sum(((w[i] / q[i]**2) * A[:, i] * A[:, i].T for i in range(5)),
            Matrix.zeros(3, 3))
    Ginv = G.inv()

    # predicted derivatives
    pred = [[(w[i] / q[i]**2) * (Ginv * A[:, i])[k] for k in range(3)] for i in range(5)]

    # independent numerical check: Newton-solve the KKT system at c and at c + eps e_k
    import mpmath as mp
    mp.mp.dps = 40
    Af = [[mp.mpf(int(A[j, i])) for i in range(5)] for j in range(3)]
    wf = [mp.mpf(int(w[i])) for i in range(5)]
    cf = [mp.mpf(str(c[j])) for j in range(3)]
    lam0 = [mp.mpf(str(lam[j])) for j in range(3)]

    def resid(L, target):
        q_ = [sum(Af[j][i] * L[j] for j in range(3)) for i in range(5)]
        return [sum(Af[j][i] * wf[i] / q_[i] for i in range(5)) - target[j] for j in range(3)]

    def newton(target, start):
        L = list(start)
        for _ in range(60):
            q_ = [sum(Af[j][i] * L[j] for j in range(3)) for i in range(5)]
            F = mp.matrix(resid(L, target))
            J = mp.matrix(3, 3)
            for j in range(3):
                for jj in range(3):
                    J[j, jj] = -sum(Af[j][i] * Af[jj][i] * wf[i] / q_[i]**2 for i in range(5))
            step = mp.lu_solve(J, F)
            L = [L[t] - step[t] for t in range(3)]
            if max(abs(x) for x in step) < mp.mpf(10)**(-30):
                break
        return L

    eps = mp.mpf(10)**(-12)
    ok_all = True
    for k in range(3):
        tgt = [cf[j] + (eps if j == k else 0) for j in range(3)]
        Lp = newton(tgt, lam0)
        qp = [sum(Af[j][i] * Lp[j] for j in range(3)) for i in range(5)]
        for i in range(5):
            du_num = (wf[i] / qp[i] - mp.mpf(str(u[i]))) / eps
            du_pred = mp.mpf(str(pred[i][k]))
            if abs(du_num - du_pred) > mp.mpf(10)**(-6) * max(1, abs(du_pred)):
                ok_all = False
    check("du_i/dc_k = (w_i/q_i^2)(G^{-1}a_i)_k, against finite differences", ok_all)

    # adding-up: sum_i a_{ji} du_i/dc_k = delta_{jk}
    adding = all(simplify(sum(A[j, i] * pred[i][k] for i in range(5))
                          - (1 if j == k else 0)) == 0
                 for j in range(3) for k in range(3))
    check("sum_i a_{ji} du_i/dc_k = delta_{jk}", adding)

    # gainers and losers coexist for every expanded resource
    mixed = all(any(pred[i][k] > 0 for i in range(5)) and any(pred[i][k] < 0 for i in range(5))
                for k in range(3))
    check("every expansion has both gainers and losers", mixed)

    # with l >= 3 the gain/loss split is a halfspace in technology space and need not be
    # monotone in any single index; check on a technology that is not on a moment curve
    B = Matrix([[1, 2, 3, 1, 4],
                [2, 1, 3, 5, 1],
                [5, 4, 1, 2, 1]])
    lamB = Matrix([Rational(1), Rational(1), Rational(1)])
    qB = B.T * lamB
    GB = sum(((1 / qB[i]**2) * B[:, i] * B[:, i].T for i in range(5)), Matrix.zeros(3, 3))
    GBinv = GB.inv()
    predB = [[(1 / qB[i]**2) * (GBinv * B[:, i])[k] for k in range(3)] for i in range(5)]
    rowsigns = [set(1 if GBinv[k, j] > 0 else -1 for j in range(3)) for k in range(3)]
    check("rows of G^{-1} carry both signs, so the cutoff is a halfspace, not a ranking",
          any(len(r) == 2 for r in rowsigns))

    indices = {"a1": lambda i: B[0, i], "a2": lambda i: B[1, i], "a3": lambda i: B[2, i],
               "a2/a1": lambda i: B[1, i] / B[0, i], "a3/a1": lambda i: B[2, i] / B[0, i],
               "a3/a2": lambda i: B[2, i] / B[1, i]}
    nonmono = {}
    for k in range(3):
        gain = [i for i in range(5) if predB[i][k] > 0]
        lose = [i for i in range(5) if predB[i][k] < 0]
        if not gain or not lose:
            continue
        bad = []
        for name, f in indices.items():
            # a scalar index orders the split if all gainers sit on one side of all losers
            if not (max(f(i) for i in gain) < min(f(j) for j in lose)
                    or min(f(i) for i in gain) > max(f(j) for j in lose)):
                bad.append(name)
        nonmono[k] = bad
    failing_somewhere = set().union(*[set(v) for v in nonmono.values()]) if nonmono else set()
    check("every candidate scalar index fails to order the split for some resource",
          failing_somewhere == set(indices))
    for k, bad in nonmono.items():
        print(f"      resource {k + 1}: indices that fail to order the split: {sorted(bad)}")


def two_resource_consistency():
    """With l = 2 the halfspace cutoff is g12/g11, which must equal the C(t) of Section 8."""
    from sympy import symbols as syms
    t = syms("t", positive=True)
    r = [Rational(1), Rational(2), Rational(4), Rational(7, 2)]
    w = [Rational(1), Rational(2), Rational(1), Rational(3)]
    p2 = Rational(1)
    q = [p2 * (t + ri) for ri in r]
    A = Matrix([[1] * 4, r])
    G = sum(((w[i] / q[i]**2) * A[:, i] * A[:, i].T for i in range(4)), Matrix.zeros(2, 2))
    cutoff = simplify(G[0, 1] / G[0, 0])
    C_of_t = simplify(sum(w[i] * r[i] / (t + r[i])**2 for i in range(4))
                      / sum(w[i] / (t + r[i])**2 for i in range(4)))
    check("l = 2 halfspace cutoff g12/g11 equals the C(t) of Section 8",
          simplify(cutoff - C_of_t) == 0)
    # and the sign rule matches: agent i gains iff r_i > cutoff
    tv = Rational(3, 2)
    Gv = G.subs(t, tv)
    Ginv = Gv.inv()
    rule = all(((Ginv * A[:, i].subs(t, tv))[1] > 0) == (r[i] > cutoff.subs(t, tv))
               for i in range(4))
    check("gain rule (G^{-1}a_i)_2 > 0 coincides with r_i > C(t)", rule)


if __name__ == "__main__":
    print("Section 3: a cooperative concept defined by a product creates complexity")
    core_nash()
    print("Section 8: general-l comparative statics of a capacity expansion")
    comparative_statics()
    two_resource_consistency()
    print(f"{len(PASS)} passed, {len(FAIL)} failed")
