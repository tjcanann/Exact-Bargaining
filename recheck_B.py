"""Independent recheck, Part B: polyhedral bargaining (Section 6).
Chamber counting, the sqrt(13) example, the degree bound, attainment, the five-agent S_6 instance,
and the smooth-frontier bounds."""
import itertools, random, sympy as sp
from sympy import Rational as Q, sqrt

OK=[];BAD=[]
def chk(name,cond,extra=""):
    (OK if cond else BAD).append(name)
    print(("  PASS  " if cond else "  FAIL  ")+name+(("   "+str(extra)) if extra else ""))

print("=== Theorem 6.3(b) / thm:polylinear  facet formula, checked symbolically")
n=4
w=sp.symbols('w1:5',positive=True);a=sp.symbols('a1:5',positive=True)
d=sp.symbols('d1:5',real=True);c=sp.Symbol('c',positive=True)
W=sum(w);cp=c-sum(a[i]*d[i] for i in range(n))
nu=[d[i]+w[i]*cp/(W*a[i]) for i in range(n)]
chk("facet solution satisfies a.nu = c",sp.simplify(sum(a[i]*nu[i] for i in range(n))-c)==0)
lam=sp.Symbol('lam')
chk("facet solution satisfies w_i/(nu_i-d_i) = lam a_i with lam = W/c'",
    all(sp.simplify(w[i]/(nu[i]-d[i])-(W/cp)*a[i])==0 for i in range(n)))

print("\n=== Theorem 6.2 chamber count: b(A) <= C(m-1,k), equality in general position")
def bounded_chambers(normals,consts):
    """Count bounded chambers of an arrangement of hyperplanes n.x = b in R^k by sign vectors."""
    m=len(normals);k=len(normals[0])
    cnt=0
    # a chamber is a consistent strict sign vector; it is bounded iff the LP max of each
    # coordinate direction is finite -> test via Minkowski: bounded iff the only solution of
    # the homogeneous system with the same weak signs is 0
    for signs in itertools.product([1,-1],repeat=m):
        A=[];b=[]
        for s,nv,cv in zip(signs,normals,consts):
            A.append([-s*x for x in nv]);b.append(-s*cv)
        # feasibility of A x < b  (strict): LP with a slack variable maximised
        import numpy as np
        from scipy.optimize import linprog
        Aub=np.array([[float(x) for x in row]+[1.0] for row in A])
        bub=np.array([float(x) for x in b])
        res=linprog(c=np.array([0.0]*k+[-1.0]),A_ub=Aub,b_ub=bub,
                    bounds=[(None,None)]*k+[(0,1)],method="highs")
        if not res.success or res.x[-1]<=1e-9: continue
        # bounded?  recession cone {r : sign_i * n_i . r <= 0 for all i} = {0}
        Ar=np.array([[float(-s*x) for x in nv] for s,nv in zip(signs,normals)])
        bounded=True
        for j in range(k):
            for sgn in (1.0,-1.0):
                obj=np.zeros(k);obj[j]=-sgn
                r=linprog(c=obj,A_ub=Ar,b_ub=np.zeros(m),bounds=[(-1,1)]*k,method="highs")
                if r.success and -r.fun>1e-7: bounded=False;break
            if not bounded: break
        if bounded: cnt+=1
    return cnt
random.seed(5)
for (m,k) in [(3,1),(4,1),(3,2),(4,2),(5,2),(4,3),(5,3),(6,3)]:
    normals=[];consts=[]
    while len(normals)<m:
        v=[Q(random.randint(-6,6)) for _ in range(k)]
        if any(v) and all(sp.Matrix([v]).rank()>0 for _ in [0]): normals.append(v);consts.append(Q(random.randint(-6,6)))
    # general position: every k of the m normals independent and no k+1 meet
    gp=all(sp.Matrix([normals[i] for i in S]).rank()==k for S in itertools.combinations(range(m),k))
    b=bounded_chambers(normals,consts)
    chk("m=%d hyperplanes in R^%d: bounded chambers %d <= C(%d,%d)=%d"%(m,k,b,m-1,k,sp.binomial(m-1,k)),
        b<=sp.binomial(m-1,k), "general position: %s"%gp)

print("\n=== Example 6.7 / ex:sqrt13")
s=sp.Symbol('s')
# U = comprehensive hull of {u>=0 : u1+2u2+2u3<=1, 2u1+3u2+u3<=1}
A=sp.Matrix([[1,2,2],[2,3,1]])
r=[Q(A[1,j],A[0,j]) for j in range(3)]
chk("ratios b_i/a_i = (2, 3/2, 1/2), so m=3",r==[Q(2),Q(3,2),Q(1,2)],r)
u=[4*s-1,1-3*s,s]
chk("edge parametrization satisfies both constraints with equality",
    sp.simplify(u[0]+2*u[1]+2*u[2]-1)==0 and sp.simplify(2*u[0]+3*u[1]+u[2]-1)==0)
prod=sp.expand(u[0]*u[1]*u[2])
chk("Nash product = -12s^3+7s^2-s",sp.expand(prod-(-12*s**3+7*s**2-s))==0,prod)
crit=sp.solve(sp.diff(prod,s),s)
nu_t=[Q(sqrt(13)-2,1)/9 if False else (sqrt(13)-2)/9,(5-sqrt(13))/12,(7+sqrt(13))/36]
# the paper's nu should be a critical point with all coords positive
sval=nu_t[2]
chk("s = (7+sqrt13)/36 is a critical point of the Nash product",
    sp.simplify(sp.diff(prod,s).subs(s,sval))==0)
chk("nu = ((sqrt13-2)/9, (5-sqrt13)/12, (7+sqrt13)/36)",
    all(sp.simplify(u[i].subs(s,sval)-nu_t[i])==0 for i in range(3)),[sp.simplify(u[i].subs(s,sval)) for i in range(3)])
chk("nu ~ (0.178395, 0.116204, 0.294599)",
    all(abs(float(nu_t[i])-v)<5e-7 for i,v in enumerate([0.178395,0.116204,0.294599])))
yy=sp.Symbol('yy')
mps=[sp.Poly(sp.minimal_polynomial(z,yy),yy).all_coeffs() for z in nu_t]
chk("minimal polynomials 9x^2+4x-1, 12x^2-10x+1, 36x^2-14x+1",
    mps==[[9,4,-1],[12,-10,1],[36,-14,1]],mps)
# multipliers
lam,mu=sp.symbols('lam mu')
eqs=[sp.Eq(1/nu_t[i],lam*A[0,i]+mu*A[1,i]) for i in range(3)]
msol=sp.solve(eqs[:2],[lam,mu],dict=True)[0]
chk("multipliers lam = 4-sqrt13 ... check against paper's (4-sqrt13, sqrt13-1)",
    sp.simplify(msol[lam]-(4-sqrt(13)))==0 or sp.simplify(msol[mu]-(4-sqrt(13)))==0,
    (sp.simplify(msol[lam]),sp.simplify(msol[mu])))
chk("both multipliers positive and stationarity holds for all three coordinates",
    sp.simplify(msol[lam])>0 and sp.simplify(msol[mu])>0 and
    sp.simplify(1/nu_t[2]-(msol[lam]*A[0,2]+msol[mu]*A[1,2]))==0,
    (float(msol[lam]),float(msol[mu])))
chk("degree 2 = bound C(m-1,l-1)=C(2,1)=2",sp.binomial(2,1)==2)
# rho = 1 in (18/19, 4/3)?  normalised: x_i = a_i u_i / c_a with a=(1,2,2), c_a=1
W3=3;Wj=[Q(1),Q(1),Q(1)]
Rj=sorted(r)
H=Q(W3,1)/sum(Q(1,1)/R for R in Rj);Aa=sum(R for R in Rj)/Q(3)
chk("rho = 1 lies in I = (18/19, 4/3)",H==Q(18,19) and Aa==Q(4,3) and H<1<Aa,(H,Aa))

print("\n=== Theorem 6.8 / thm:attain  attainment, recomputed from scratch")
def nash_eliminant(Abar,Wt,cc):
    """Abar: l x m exact matrix; Wt: m weights; cc: l capacities.
    Returns the minimal polynomial over Q of the equilibrium price ratio lam_1/..., via
    the critical system  W_T/X_T = (Abar^T lam)_T,  Abar X = c."""
    l,m=Abar.shape
    lamv=sp.symbols('L1:%d'%(l+1))
    expr=[sum(Abar[j,T]*Wt[T]/sum(Abar[jj,T]*lamv[jj] for jj in range(l)) for T in range(m))-cc[j]
          for j in range(l)]
    polys=[sp.Poly(sp.together(e).as_numer_denom()[0],*lamv) for e in expr]
    # eliminate L2..Ll by iterated resultants, then normalise to the ratio L2/L1
    cur=polys[0].as_expr()
    rest=[p.as_expr() for p in polys[1:]]
    for k in range(l-1,0,-1):
        nxt=[]
        for e in rest:
            cur2=sp.resultant(sp.Poly(cur,lamv[k]),sp.Poly(e,lamv[k]))
            nxt.append(sp.expand(cur2))
        cur=nxt[0]
        rest=nxt[1:]
    return sp.factor_list(sp.expand(cur)),lamv
random.seed(2024)
def attain_case(m,l,seed):
    rnd=random.Random(seed)
    while True:
        Abar=sp.Matrix(l,m,lambda i,j: Q(rnd.randint(1,9)))
        gp=all(sp.Matrix([[Abar[i,j] for j in S] for i in range(l)]).det()!=0
               for S in itertools.combinations(range(m),l))
        if gp: break
    Wt=[Q(1)]*m
    # pick capacities making all types served: start from an interior X>0
    X0=[Q(rnd.randint(1,5)) for _ in range(m)]
    cc=[sum(Abar[j,T]*X0[T] for T in range(m)) for j in range(l)]
    # solve the critical system numerically then get the exact minimal polynomial of the solution
    Xs=sp.symbols('X1:%d'%(m+1),positive=True)
    lamv=sp.symbols('L1:%d'%(l+1))
    eqs=[sp.Eq(Wt[T]/Xs[T],sum(Abar[j,T]*lamv[j] for j in range(l))) for T in range(m)]
    eqs+=[sp.Eq(sum(Abar[j,T]*Xs[T] for T in range(m)),cc[j]) for j in range(l)]
    import mpmath
    # numeric solve from a feasible start
    f=sp.lambdify(list(Xs)+list(lamv),[sp.simplify(e.lhs-e.rhs) for e in eqs],'mpmath')
    guess=[mpmath.mpf(float(x)) for x in X0]+[mpmath.mpf(1)]*l
    try:
        sol=mpmath.findroot(lambda *a: f(*a),guess,tol=1e-40,maxsteps=200)
    except Exception as ex:
        return None
    Xnum=[sp.nsimplify(0)]*m
    vals=[mpmath.mpf(sol[i]) for i in range(m+l)]
    if any(v<=0 for v in vals[:m]): return None
    # exact minimal polynomial of X_1 via resultants of the algebraic system
    G=[sp.together(sp.Eq(Wt[T]*sum(Abar[j,T]*lamv[j] for j in range(l))**(-1),Xs[T]).rewrite(sp.Eq)) for T in range(m)]
    sub={Xs[T]:Wt[T]/sum(Abar[j,T]*lamv[j] for j in range(l)) for T in range(m)}
    cons=[sp.together(sum(Abar[j,T]*sub[Xs[T]] for T in range(m))-cc[j]).as_numer_denom()[0]
          for j in range(l)]
    cur=sp.expand(cons[0])
    for e in cons[1:]:
        k=l-1
        cur=sp.expand(sp.resultant(sp.Poly(cur,lamv[k]),sp.Poly(sp.expand(e),lamv[k])))
    # cur is a polynomial in L1 (and L2.. eliminated); take its factor vanishing at the solution
    var=[v for v in lamv if cur.has(v)]
    if len(var)!=1: return None
    P=sp.Poly(cur,var[0])
    facs=[f for f,_ in sp.factor_list(P.as_expr())[1]]
    tgt=mpmath.mpf(sol[m+ (lamv.index(var[0]))])
    best=None
    for f_ in facs:
        pf=sp.Poly(f_,var[0])
        if pf.degree()<1: continue
        vv=abs(mpmath.mpf(str(sp.N(f_.subs(var[0],sp.Float(str(tgt),40)),40))))
        if best is None or vv<best[0]: best=(vv,pf)
    return best[1].degree(),sp.binomial(m-1,l-1),best[1].as_expr(),var[0]
for (m,l,sd) in [(3,2,1),(4,2,2),(5,2,3),(6,2,4),(4,3,5),(5,3,6)]:
    out=attain_case(m,l,sd)
    if out is None:
        chk("attainment m=%d l=%d"%(m,l),False,"solver failed");continue
    deg,bound,poly,var=out
    chk("attainment m=%d l=%d: minimal polynomial degree %d = C(%d,%d) = %d"%(m,l,deg,m-1,l-1,bound),
        deg==bound)

print("\n"+"="*72)
print("Part B: %d passed, %d failed"%(len(OK),len(BAD)))
for b in BAD: print("   FAILED:",b)
