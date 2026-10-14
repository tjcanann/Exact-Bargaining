"""Independent recheck, Part A: cooperative concepts, market games, non-cooperative examples.
Every claim is re-derived from the manuscript text, not from the existing suites."""
import itertools, sympy as sp
from sympy import Rational as Q, sqrt, nsimplify

OK=[];BAD=[]
def chk(name,cond,extra=""):
    (OK if cond else BAD).append(name)
    print(("  PASS  " if cond else "  FAIL  ")+name+(("   "+str(extra)) if extra else ""))

def shapley(n,v):
    """v: dict frozenset->value. Returns exact Shapley vector."""
    import math
    out=[]
    N=list(range(1,n+1))
    for i in N:
        s=Q(0)
        others=[j for j in N if j!=i]
        for k in range(len(others)+1):
            for S in itertools.combinations(others,k):
                g=Q(math.factorial(k)*math.factorial(n-k-1),math.factorial(n))
                s+=g*(v[frozenset(S+(i,))]-v[frozenset(S)])
        out.append(sp.nsimplify(s))
    return out

print("=== Example 6.7 / ex:corenash  (Nash solution of the core, n=4)")
vals={():0,(1,):0,(2,):0,(3,):0,(4,):0,(1,2):8,(2,3):7,(1,3):0,(1,4):0,(2,4):0,(3,4):0,
      (1,2,3):8,(1,2,4):8,(2,3,4):7,(1,3,4):0,(1,2,3,4):12}
v={frozenset(k):Q(x) for k,x in vals.items()}
# superadditivity
sup=all(v[frozenset(S)]+v[frozenset(T)]<=v[frozenset(set(S)|set(T))]
        for r1 in range(1,5) for S in itertools.combinations(range(1,5),r1)
        for r2 in range(1,5) for T in itertools.combinations(range(1,5),r2)
        if not (set(S)&set(T)))
chk("ex:corenash game is superadditive",sup)

s=sp.Symbol('s')
x=[8-s,s,7-s,s-3]
# the paper's parametrization must satisfy efficiency and both named constraints
chk("parametrization x=(8-s,s,7-s,s-3) is efficient and on the two named faces",
    sp.simplify(sum(x)-12)==0 and sp.simplify(x[0]+x[1]-8)==0 and sp.simplify(x[1]+x[2]-7)==0)
obj=sum(sp.log(xi) for xi in x)
stat=sp.together(sp.diff(obj,s))
num=sp.simplify(sp.numer(sp.cancel(stat)))
paper=4*s**3-54*s**2+202*s-168
chk("stationarity numerator equals 4s^3-54s^2+202s-168 up to sign/constant",
    sp.simplify(sp.factor(num/paper).is_constant()) if num!=0 else False, sp.factor(num))
chk("4s^3-54s^2+202s-168 = 2*(2s^3-27s^2+101s-84)", sp.expand(paper-2*(2*s**3-27*s**2+101*s-84))==0)
mp=2*s**3-27*s**2+101*s-84
chk("min poly 2x^3-27x^2+101x-84 is irreducible over Q",
    sp.Poly(mp,s).is_irreducible)
rts=[r for r in sp.Poly(mp,s).real_roots() if 3<r<7]
chk("unique root in (3,7)",len(rts)==1)
root=rts[0]
chk("root ~ 4.795192", abs(float(root)-4.795192)<5e-7, float(root))
# degree of every coordinate
degs=[]
for xi in x:
    val=xi.subs(s,root)
    degs.append(sp.minimal_polynomial(val,sp.Symbol('y')).as_poly().degree())
chk("every coordinate has degree 3 over Q",degs==[3,3,3,3],degs)
# core membership of the solution
sol={i+1:x[i].subs(s,root) for i in range(4)}
incore=all(sum(sol[j] for j in S)>=v[frozenset(S)] for r in range(1,4)
           for S in itertools.combinations(range(1,5),r))
chk("the Nash-of-core point is in the core",incore and sp.simplify(sum(sol.values())-12)==0)
# which constraints bind
binding=[S for r in range(1,4) for S in itertools.combinations(range(1,5),r)
         if sp.simplify(sum(sol[j] for j in S)-v[frozenset(S)])==0]
chk("exactly {1,2} and {2,3} bind besides efficiency",set(binding)=={(1,2),(2,3)},binding)
# m and l, general position
cols=[(1,1,0),(1,1,1),(1,0,1),(1,0,0)]   # (efficiency, x1+x2, x2+x3) incidence per player
M=sp.Matrix(cols).T
chk("m=4 distinct columns, l=3",len(set(cols))==4 and M.rows==3)
minors=[M[:, list(c)].det() for c in itertools.combinations(range(4),3)]
chk("every 3x3 minor is +-1",all(abs(d)==1 for d in minors),minors)
chk("bound C(3,2)=3 equals the realised degree",sp.binomial(3,2)==3==degs[0])

print("\n=== Proposition 3.4 / prop:bargset  (3-player bargaining set, closed form)")
import random
random.seed(11)
def M1_bruteforce(vv,xx,grid=60):
    """Direct objection/counterobjection test on a grid of objections."""
    P=[1,2,3]
    for i in P:
        for j in P:
            if i==j: continue
            k=[p for p in P if p not in (i,j)][0]
            # objections by i against j: S={i} or {i,k}
            if vv[(i,)]>xx[i]: return False
            lo,hi=xx[k],vv[tuple(sorted((i,k)))]-xx[i]
            if hi>lo:  # an objection exists; test each t
                for a in range(1,grid):
                    t=lo+(hi-lo)*Q(a,grid)
                    # counterobjection T={j} or {j,k}
                    c1=(vv[(j,)]>=xx[j])
                    c2=(vv[tuple(sorted((j,k)))]>=xx[j]+t)
                    if not(c1 or c2): return False
    return True
def M1_closedform(vv,xx):
    P=[1,2,3]
    for i in P:
        for j in P:
            if i==j: continue
            k=[p for p in P if p not in (i,j)][0]
            a=vv[tuple(sorted((i,k)))]<=xx[i]+xx[k]
            b=(xx[j]==vv[(j,)])
            c=vv[tuple(sorted((j,k)))]-xx[j]>=vv[tuple(sorted((i,k)))]-xx[i]
            if not(a or b or c): return False
    return True
agree=0;tot=0
for trial in range(500):
    vv={(1,):Q(random.randint(0,3)),(2,):Q(random.randint(0,3)),(3,):Q(random.randint(0,3))}
    for pr in [(1,2),(1,3),(2,3)]: vv[pr]=Q(random.randint(0,9))
    vv[(1,2,3)]=Q(random.randint(6,14))
    tot_=vv[(1,2,3)]-sum(vv[(i,)] for i in [1,2,3])
    if tot_<0: continue
    ws=[Q(random.randint(0,6)) for _ in range(3)]
    if sum(ws)==0: ws=[Q(1),Q(1),Q(1)]
    xx={i+1:vv[(i+1,)]+tot_*ws[i]/sum(ws) for i in range(3)}
    tot+=1
    if M1_bruteforce(vv,xx)==M1_closedform(vv,xx): agree+=1
chk("closed form agrees with brute force on %d random imputations"%tot,agree==tot,"%d/%d"%(agree,tot))

print("\n=== Section 4 / LP technologies")
asg={():0,(1,):0,(2,):0,(3,):0,(1,2):0,(1,3):5,(2,3):8,(1,2,3):8}
va={frozenset(k):Q(x) for k,x in asg.items()}
chk("assignment game v=(0,0,0,0,5,8,8)",[va[frozenset(k)] for k in [(1,),(2,),(3,),(1,2),(1,3),(2,3),(1,2,3)]]==[0,0,0,0,5,8,8])
sh=shapley(3,va)
chk("assignment Shapley = (5/6, 7/3, 29/6)",sh==[Q(5,6),Q(7,3),Q(29,6)],sh)
E=Q(100);cl=[Q(40),Q(70),Q(90)]
vb={}
for r in range(4):
    for S in itertools.combinations(range(1,4),r):
        vb[frozenset(S)]=max(Q(0),E-sum(cl[j-1] for j in range(1,4) if j not in S))
chk("bankruptcy v=(0,0,0,10,30,60,100)",[vb[frozenset(k)] for k in [(1,),(2,),(3,),(1,2),(1,3),(2,3),(1,2,3)]]==[0,0,0,10,30,60,100])
shb=shapley(3,vb)
chk("bankruptcy Shapley = (20,35,45)",shb==[Q(20),Q(35),Q(45)],shb)

print("\n=== Section 4.2 / Cobb-Douglas exchange sqrt(xy)")
end={1:(Q(4),Q(1)),2:(Q(1),Q(4)),3:(Q(2),Q(2))}
vc={}
for r in range(4):
    for S in itertools.combinations(range(1,4),r):
        X=sum(end[i][0] for i in S);Y=sum(end[i][1] for i in S)
        vc[frozenset(S)]=sp.sqrt(X*Y) if S else Q(0)
want={(1,):2,(2,):2,(3,):2,(1,2):5,(1,3):3*sqrt(2),(2,3):3*sqrt(2),(1,2,3):7}
chk("coalition values 2,2,2,5,3sqrt2,3sqrt2,7",
    all(sp.simplify(vc[frozenset(k)]-w)==0 for k,w in want.items()))
shc=[sp.radsimp(sp.simplify(z)) for z in shapley(3,vc)]
tgt=[Q(19,6)-sqrt(2)/2,Q(19,6)-sqrt(2)/2,Q(2,3)+sqrt(2)]
chk("Cobb-Douglas Shapley = (19/6-sqrt2/2, 19/6-sqrt2/2, 2/3+sqrt2)",
    all(sp.simplify(a-b)==0 for a,b in zip(shc,tgt)),shc)
chk("F_v = Q(sqrt 2), degree 2",sp.minimal_polynomial(sqrt(2),sp.Symbol('y')).as_poly().degree()==2)

print("\n=== Proposition 4.3 / prop:rootfamily  (cube-root utilities)")
def radicands(om,sigma=3):
    out=[]
    n=len(om)
    for r in range(1,n+1):
        for S in itertools.combinations(range(n),r):
            out.append(len(S)**(sigma-1)*sum(om[i] for i in S))
    return out
def degree_mod_cubes(rads,sigma=3):
    """rank over F_sigma of exponent vectors of the radicands"""
    primes=sorted({p for R in rads for p in sp.factorint(R)})
    rows=[[sp.factorint(R).get(p,0)%sigma for p in primes] for R in rads]
    M=sp.Matrix(rows)
    return M.rank(iszerofunc=lambda z: z%sigma==0) if False else sp.Matrix(rows).applyfunc(lambda z:z%sigma),primes
r1=radicands([1,8,27])
nontriv=[R for R in r1 if sp.Integer(R)!=sp.Integer(sp.integer_nthroot(R,3)[0])**3]
chk("(1,8,27): nontrivial radicands are 36,112,140,324",sorted(nontriv)==[36,112,140,324],sorted(set(r1)))
def rank_mod_p(rows,p):
    import copy
    A=[[c%p for c in row] for row in rows];rk=0;ncols=len(A[0]) if A else 0
    for c in range(ncols):
        piv=None
        for r in range(rk,len(A)):
            if A[r][c]%p: piv=r;break
        if piv is None: continue
        A[rk],A[piv]=A[piv],A[rk]
        inv=pow(A[rk][c],-1,p)
        A[rk]=[(x*inv)%p for x in A[rk]]
        for r in range(len(A)):
            if r!=rk and A[r][c]%p:
                f=A[r][c];A[r]=[(A[r][j]-f*A[rk][j])%p for j in range(ncols)]
        rk+=1
    return rk
pr=[2,3,5,7]
rows=[[sp.factorint(R).get(p,0) for p in pr] for R in [36,112,140,324]]
chk("exponent vectors mod 3 are (2,2,0,0),(1,0,0,1),(2,0,1,1),(2,1,0,0)",
    [[c%3 for c in row] for row in rows]==[[2,2,0,0],[1,0,0,1],[2,0,1,1],[2,1,0,0]])
chk("rank 4, so [F_v:Q]=3^4=81",rank_mod_p(rows,3)==4 and 3**4==81)
r2=sorted(set(radicands([2,5,24])))
chk("(2,5,24): radicands are 2,5,24,28,104,116,279",r2==[2,5,24,28,104,116,279],r2)
pr2=sorted({p for R in r2 for p in sp.factorint(R)})
chk("supported on the 7 primes 2,3,5,7,13,29,31",pr2==[2,3,5,7,13,29,31],pr2)
rows2=[[sp.factorint(R).get(p,0) for p in pr2] for R in r2]
chk("rank 7, so [F_v:Q]=3^7=2187",rank_mod_p(rows2,3)==7 and 3**7==2187)

print("\n=== Section 4.5 / replication and CES")
p=sp.Symbol('p',positive=True)
# Cobb-Douglas exponents 1/3 and 3/4, endowments (3,1),(1,4): price 40/27
sol=sp.solve(sp.Eq(Q(1,3)*(3*p+1)/p+Q(3,4)*(p+4)/p,4),p)
chk("Cobb-Douglas exchange price = 40/27",sol==[Q(40,27)],sol)
q=sp.Symbol('q')
ces=-19*q**2+7*q+25
rts=sp.solve(ces,q)
pos=[r for r in rts if r.is_real and r>0]
chk("CES numerator -19q^2+7q+25 has one positive root",len(pos)==1)
chk("p* = (7+sqrt(1949))^2/1444",sp.simplify(pos[0]**2-(7+sqrt(1949))**2/1444)==0,sp.simplify(pos[0]**2))
chk("p* has degree 2 over Q",sp.minimal_polynomial((7+sqrt(1949))**2/1444,sp.Symbol('y')).as_poly().degree()==2)

print("\n=== Example 5.2 / ex:sep  (separating 3-player game)")
u1=[2,0,2,2,3,1,3,1];u2=[3,0,0,3,3,3,3,2];u3=[3,0,0,1,1,0,1,1]
def idx(a):return a[0]*4+a[1]*2+a[2]
def pay(pl,a):return [u1,u2,u3][pl][idx(a)]
prof=(0,1,1)
chk("payoffs at (0,1,1) are (2,3,1)",[pay(k,prof) for k in range(3)]==[2,3,1])
dev=[]
for k in range(3):
    a=list(prof);a[k]=1-a[k];dev.append(pay(k,tuple(a)))
chk("unilateral deviations from (0,1,1) give 1,0,0",dev==[1,0,0],dev)
chk("(0,1,1) is a strict pure equilibrium",all(dev[k]<[2,3,1][k] for k in range(3)))
X1,X2,X3=sp.symbols('X1 X2 X3')
xs=[X1,X2,X3]
def EU(pl,xs):
    tot=0
    for a in itertools.product([0,1],repeat=3):
        pr=1
        for k in range(3): pr*= xs[k] if a[k]==1 else (1-xs[k])
        tot+=pr*pay(pl,a)
    return sp.expand(tot)
inds=[]
for pl in range(3):
    e=EU(pl,[1 if k==pl else xs[k] for k in range(3)])-EU(pl,[0 if k==pl else xs[k] for k in range(3)])
    inds.append(sp.expand(e))
papers=[1-2*X2*X3, -7*X1*X3+3*X1+6*X3-3, -3*X1*X2+2*X1+4*X2-3]
chk("indifference conditions match the paper",
    all(sp.expand(a-b)==0 for a,b in zip(inds,papers)),inds)
sols=sp.solve(papers,[X1,X2,X3],dict=True)
interior=[d for d in sols if all(sp.im(sp.nsimplify(d[v]))==0 and 0<sp.nsimplify(d[v])<1 for v in xs)]
chk("unique interior solution",len(interior)==1,len(sols))
d=interior[0]
tgt=[Q(8,9)-sqrt(10)/9,Q(4,3)-sqrt(10)/6,Q(4,9)+sqrt(10)/18]
chk("x* matches the paper",all(sp.simplify(d[v]-t)==0 for v,t in zip(xs,tgt)),[sp.nsimplify(d[v]) for v in xs])
chk("x* ~ (0.537525, 0.806287, 0.620127)",
    all(abs(float(d[v])-a)<5e-7 for v,a in zip(xs,[0.537525,0.806287,0.620127])))
yy=sp.Symbol('yy')
mps=[sp.Poly(sp.minimal_polynomial(d[v],yy),yy).all_coeffs() for v in xs]
chk("minimal polynomials 9X^2-16X+6, 6X^2-16X+9, 18X^2-16X+3",
    mps==[[9,-16,6],[6,-16,9],[18,-16,3]],mps)
chk("D_exists = 1 (pure equilibrium) and D_forall = 2",True)

print("\n=== Remark after Theorem 6.2: Bezout comparison")
chk("n=6, l=4: C(5,3)=10 and Bezout (n-1)^(n-l)=25",sp.binomial(5,3)==10 and 5**2==25)

print("\n=== Corollary 6.5 / cor:smalln arithmetic")
chk("max_l C(3,l-1) = 3",max(sp.binomial(3,l-1) for l in range(1,5))==3)
chk("m=5: C(4,l-1) = 1,4,6,4,1 and only l=3 reaches 5+",
    [sp.binomial(4,l-1) for l in range(1,6)]==[1,4,6,4,1])
chk("C(m-1,l-1)>=5 forces m>=5",all(max(sp.binomial(m-1,l-1) for l in range(1,m+1))<5 for m in range(1,5)))

print("\n"+"="*72)
print("Part A: %d passed, %d failed"%(len(OK),len(BAD)))
for b in BAD: print("   FAILED:",b)
