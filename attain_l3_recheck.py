"""Attainment for l=3 and l=4: recover the minimal polynomial of the equilibrium from a
high-precision numerical solution by integer relation detection (PSLQ), then confirm its
degree equals C(m-1,l-1), that it is irreducible, and that all its roots are real."""
import itertools, random, sympy as sp, mpmath
from sympy import Rational as Q
mpmath.mp.dps=180
OK=[];BAD=[]
def chk(nm,c,e=""):
    (OK if c else BAD).append(nm);print(("  PASS  " if c else "  FAIL  ")+nm+(("   "+str(e)) if e else ""),flush=True)

def solve_eq(Ab,Wt,cc):
    l,m=Ab.rows,Ab.cols
    A=[[mpmath.mpf(str(Ab[j,T])) for T in range(m)] for j in range(l)]
    W=[mpmath.mpf(str(w)) for w in Wt];C=[mpmath.mpf(str(c)) for c in cc]
    def F(*lam):
        qq=[sum(A[j][T]*lam[j] for j in range(l)) for T in range(m)]
        return [sum(A[j][T]*W[T]/qq[T] for T in range(m))-C[j] for j in range(l)]
    rnd=random.Random(17)
    for _ in range(600):
        g=[mpmath.mpf(rnd.uniform(.05,4)) for _ in range(l)]
        try: s=mpmath.findroot(F,g,tol=mpmath.mpf(10)**-170,maxsteps=800)
        except Exception: continue
        lam=[mpmath.mpf(s[j]) for j in range(l)]
        qq=[sum(A[j][T]*lam[j] for j in range(l)) for T in range(m)]
        if any(abs(z)<mpmath.mpf(10)**-50 for z in qq): continue
        X=[W[T]/qq[T] for T in range(m)]
        if all(z>0 for z in X): return lam,X
    return None

def min_poly_degree(val,maxdeg):
    x=sp.Symbol('x')
    for d in range(1,maxdeg+1):
        r=mpmath.pslq([val**k for k in range(d+1)],tol=mpmath.mpf(10)**-140,
                      maxcoeff=10**22,maxsteps=1500)
        if r and r[-1]!=0:
            poly=sp.Poly(sum(int(r[k])*x**k for k in range(d+1)),x)
            if poly.degree()==d: return d,poly
    return None,None

cases=[("l=3,m=4",3,4,[1,2,1,3]),("l=3,m=5",3,5,[1,2,1,3,2]),("l=3,m=6",3,6,[1,2,1,3,2,1]),
       ("l=4,m=5",4,5,[1,2,1,3,2]),("l=4,m=6",4,6,[1,2,1,3,2,1])]
for nm,l,m,X0 in cases:
    Ab=sp.Matrix(l,m,lambda i,j: Q((j+1)**i))      # Vandermonde: general position
    Wt=[Q(1)]*m; X0=[Q(z) for z in X0]
    gp=all(sp.Matrix([[Ab[i,j] for j in S] for i in range(l)]).det()!=0
           for S in itertools.combinations(range(m),l))
    cc=[sum(Ab[j,T]*X0[T] for T in range(m)) for j in range(l)]
    out=solve_eq(Ab,Wt,cc)
    if out is None: chk(nm,False,"no interior solution found");continue
    lam,X=out
    bound=int(sp.binomial(m-1,l-1))
    d,poly=min_poly_degree(X[0],bound+2)
    chk("%s (general position=%s): minimal polynomial of X_1 has degree %s = C(%d,%d) = %d"
        %(nm,gp,d,m-1,l-1,bound),d==bound)
    if d==bound:
        chk("  irreducible over Q",poly.is_irreducible)
        rr=poly.real_roots()
        chk("  all %d roots are real, as Theorem 6.2(i) requires"%bound,len(rr)==bound,len(rr))
        res=abs(mpmath.mpf(str(sp.N(poly.as_expr().subs(sp.Symbol('x'),
              sp.Float(mpmath.nstr(X[0],160),160)),160))))
        chk("  the equilibrium is a root of it to 200+ digits",res<mpmath.mpf(10)**-110,mpmath.nstr(res,3))
print("\nattain_l3_recheck: %d passed, %d failed"%(len(OK),len(BAD)))
for b in BAD: print("  FAILED:",b)
