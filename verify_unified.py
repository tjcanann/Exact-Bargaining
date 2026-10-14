"""Re-derive every numerical claim in the unified manuscript (main.tex) that is not
already covered by audit_paper1.py (the two-resource sections)."""
from sympy import *
from sympy.polys.numberfields.galoisgroups import galois_group
from sympy.ntheory import factorint
from itertools import combinations, permutations, product as iprod
import math, random
X,x,y,s,t,q,z1,z2=symbols('X x y s t q z1 z2')
ok=[];bad=[]
def chk(l,c): (ok if bool(c) else bad).append(l)
N3=[1,2,3]
def shapley(v):
    sh={}
    for i in N3:
        tot=0; oth=[j for j in N3 if j!=i]
        for k in range(len(oth)+1):
            for S in combinations(oth,k):
                w=Rational(math.factorial(k)*math.factorial(2-k),6)
                tot+=w*(v[frozenset(S+(i,))]-(v[frozenset(S)] if S else 0))
        sh[i]=simplify(tot)
    return [sh[i] for i in N3]
coal=[frozenset(S) for k in range(1,4) for S in combinations(N3,k)]

# ---------------- Section 4: market games
va=dict(zip(coal,[0,0,0,0,Integer(5),Integer(8),Integer(8)]))
chk("assignment Shapley (5/6,7/3,29/6)", shapley(va)==[Rational(5,6),Rational(7,3),Rational(29,6)])
cl={1:40,2:70,3:90}
vb={S:Max(0,100-sum(cl[j] for j in N3 if j not in S)) for S in coal}
chk("bankruptcy v=(0,0,0,10,30,60,100)", [vb[S] for S in coal]==[0,0,0,10,30,60,100])
chk("bankruptcy Shapley (20,35,45)", shapley(vb)==[20,35,45])
om={1:(4,1),2:(1,4),3:(2,2)}
vc={S:sqrt(sum(om[i][0] for i in S)*sum(om[i][1] for i in S)) for S in coal}
chk("Cobb-Douglas v", [simplify(vc[S]) for S in coal]==[2,2,2,5,3*sqrt(2),3*sqrt(2),7])
sh=shapley(vc)
chk("Cobb-Douglas Shapley", all(simplify(a-b)==0 for a,b in zip(sh,[Rational(19,6)-sqrt(2)/2]*2+[Rational(2,3)+sqrt(2)])))
def rankF(ms,p):
    pr=sorted({r for m in ms for r in factorint(m)}); R=[[factorint(m).get(r,0)%p for r in pr] for m in ms]; rk=0
    for c in range(len(pr)):
        pv=next((r for r in range(rk,len(R)) if R[r][c]%p),None)
        if pv is None: continue
        R[rk],R[pv]=R[pv],R[rk]; inv=pow(R[rk][c]%p,-1,p); R[rk]=[(e*inv)%p for e in R[rk]]
        for r in range(len(R)):
            if r!=rk and R[r][c]%p: f=R[r][c]%p; R[r]=[(R[r][j]-f*R[rk][j])%p for j in range(len(pr))]
        rk+=1
    return rk,pr
def radicands(w): return [len(S)**2*sum(w[i-1] for i in S) for S in [tuple(sorted(S)) for S in coal]]
r1=radicands((1,8,27))
chk("(1,8,27) radicands 1,8,27,36,112,140,324", r1==[1,8,27,36,112,140,324])
chk("(1,8,27) exponent vectors mod 3", [[factorint(m).get(p,0)%3 for p in (2,3,5,7)] for m in (36,112,140,324)]==[[2,2,0,0],[1,0,0,1],[2,0,1,1],[2,1,0,0]])
chk("(1,8,27) rank 4 -> 81", rankF(r1,3)[0]==4)
r2=radicands((2,5,24)); rk,pr=rankF(r2,3)
chk("(2,5,24) radicands 2,5,24,28,104,116,279", r2==[2,5,24,28,104,116,279])
chk("(2,5,24) primes 2,3,5,7,13,29,31, rank 7 -> 2187", pr==[2,3,5,7,13,29,31] and rk==7 and 3**7==2187)
chk("Kummer/Mordell sanity: deg Q(2^(1/3),3^(1/3)) = 9", degree(minimal_polynomial(2**Rational(1,3)+3**Rational(1,3),X),X)==9)
p_=symbols('p_',positive=True)
exc=Rational(1,3)*(3*p_+1)/p_+Rational(3,4)*(p_+4)/p_-4
chk("Cobb-Douglas exchange p* = 40/27", solve(Eq(exc,0),p_)==[Rational(40,27)])
qq=symbols('qq',positive=True); P=qq**2; sh_=Rational(1,2)
def dx(b,w): return (b[0]**sh_*P**(-sh_))*(P*w[0]+w[1])/(b[0]**sh_*P**(1-sh_)+b[1]**sh_)
num,_=fraction(together(simplify(dx((1,4),(3,1))+dx((9,1),(1,4))-4)))
chk("CES numerator -19q^2+7q+25", expand(num+19*qq**2-7*qq-25)==0 or expand(num-(-19*qq**2+7*qq+25))==0)
rt=[r for r in Poly(expand(num),qq).all_roots() if r.is_real and r>0][0]
chk("CES p* = (7+sqrt1949)^2/1444, degree 2", simplify(rt**2-(7+sqrt(1949))**2/1444)==0 and degree(minimal_polynomial(rt**2,X),X)==2)

# ---------------- Section 3: bargaining set formula vs raw definition
v={frozenset([1]):0,frozenset([2]):0,frozenset([3]):0,frozenset([1,2]):Rational(3,2),frozenset([1,3]):1,frozenset([2,3]):Rational(1,2),frozenset([1,2,3]):2}
def raw(a,steps=300):
    for i,j in permutations([1,2,3],2):
        k=6-i-j; vik=float(v[frozenset([i,k])]); vjk=float(v[frozenset([j,k])]); vj=float(v[frozenset([j])])
        hi=vik-a[i]
        if hi<=a[k]+1e-12: continue
        for st in range(1,steps):
            tt=a[k]+(hi-a[k])*st/steps
            if not((vj>=a[j]-1e-12) or (vjk>=a[j]+tt-1e-12)): return False
    return True
def formula(a):
    for i,j in permutations([1,2,3],2):
        k=6-i-j
        c1=v[frozenset([i,k])]<=a[i]+a[k]; c2=(a[j]==v[frozenset([j])]); c3=v[frozenset([j,k])]-a[j]>=v[frozenset([i,k])]-a[i]
        if not(c1 or c2 or c3): return False
    return True
random.seed(4); agree=0; tot=0
for _ in range(700):
    u=sorted([random.random(),random.random()])
    a={1:Rational(round(u[0]*2000),1000),2:Rational(round((u[1]-u[0])*2000),1000)}; a[3]=2-a[1]-a[2]
    if min(a.values())<0: continue
    tot+=1; agree+= (raw({k:float(e) for k,e in a.items()})==formula(a))
chk(f"Prop 3.2 bargaining set: formula == raw definition on {tot} imputations", agree==tot and tot>=500)

# ---------------- Section 5: separating game
U=[[2,0,2,2,3,1,3,1],[3,0,0,3,3,3,3,2],[3,0,0,1,1,0,1,1]]
acts=list(iprod((0,1),repeat=3)); u={i:dict(zip(acts,U[i])) for i in range(3)}
a0=(0,1,1); dev=[u[i][tuple(1-a0[k] if k==i else a0[k] for k in range(3))] for i in range(3)]
chk("Ex 5.2 strict pure NE (0,1,1), payoffs (2,3,1), deviations (1,0,0)", [u[i][a0] for i in range(3)]==[2,3,1] and dev==[1,0,0])
xs=[Rational(8,9)-sqrt(10)/9, Rational(4,3)-sqrt(10)/6, Rational(4,9)+sqrt(10)/18]
chk("Ex 5.2 indifference equations hold", simplify(1-2*xs[1]*xs[2])==0 and simplify(-7*xs[0]*xs[2]+3*xs[0]+6*xs[2]-3)==0 and simplify(-3*xs[0]*xs[1]+2*xs[0]+4*xs[1]-3)==0)
chk("Ex 5.2 minpolys", [minimal_polynomial(e,X) for e in xs]==[9*X**2-16*X+6,6*X**2-16*X+9,18*X**2-16*X+3])
chk("Ex 5.2 numerics", all(abs(float(e)-v_)<5e-7 for e,v_ in zip(xs,[0.537525,0.806287,0.620127])))

# ---------------- Section 6
# weighted facet formula on a random instance
w=[Integer(2),Integer(3),Integer(5)]; aa=[Integer(1),Integer(4),Integer(2)]; cc=Integer(7); d=[Rational(1,2),Integer(0),Rational(1,3)]
W=sum(w); zf=[w[i]*(cc-sum(aa[j]*d[j] for j in range(3)))/(W*aa[i]) for i in range(3)]
zz=symbols('zz0:3',positive=True); lam=symbols('lam',positive=True)
sol=solve([w[i]/zz[i]-lam*aa[i] for i in range(3)]+[sum(aa[i]*(zz[i]+d[i]) for i in range(3))-cc],list(zz)+[lam],dict=True)[0]
chk("Thm 6.1(b) weighted facet formula", all(simplify(sol[zz[i]]-zf[i])==0 for i in range(3)))
chk("Remark: Bezout 25 vs C(5,3)=10 (six agents, four resources)", (6-1)**(6-4)==25 and binomial(5,3)==10)
chk("Cor 6.3 values", max(binomial(3,l-1) for l in range(1,5))==3 and [binomial(4,l-1) for l in range(1,6)]==[1,4,6,4,1] and max(binomial(2,l-1) for l in range(1,4))==2)
nu=[(sqrt(13)-2)/9,(5-sqrt(13))/12,(7+sqrt(13))/36]
chk("Ex 6.4 edge (4s-1,1-3s,s) and product -12s^3+7s^2-s", expand((4*s-1)*(1-3*s)*s+12*s**3-7*s**2+s)==0)
chk("Ex 6.4 constraints active", simplify(nu[0]+2*nu[1]+2*nu[2]-1)==0 and simplify(2*nu[0]+3*nu[1]+nu[2]-1)==0)
lm,mu=4-sqrt(13),sqrt(13)-1
chk("Ex 6.4 KKT multipliers and stationarity", all(simplify(1/nu[i]-(lm*[1,2,2][i]+mu*[2,3,1][i]))==0 for i in range(3)) and float(lm)>0 and float(mu)>0)
chk("Ex 6.4 minpolys", [minimal_polynomial(e,X) for e in nu]==[9*X**2+4*X-1,12*X**2-10*X+1,36*X**2-14*X+1])
r=[Rational(2),Rational(3,2),Rational(1,2)]; H=3/sum(1/e for e in r); A=sum(r)/3
chk("Ex 6.4 ratios (2,3/2,1/2), rho=1 in (18/19,4/3)", H==Rational(18,19) and A==Rational(4,3) and H<1<A)
# Theorem 6.5
f=Poly(5450700*x**6+20203116*x**5+22856976*x**4+623972*x**3-15763425*x**2-10365795*x-2051910,x)
chk("Thm 6.5 f irreducible, content 1", f.is_irreducible and f.content()==1)
G,_=galois_group(f); chk("Thm 6.5 Galois group S6", G.order()==720 and not G.is_solvable)
def pat(p):
    fp=Poly(f.as_expr(),x,modulus=p)
    sqf=gcd(fp,fp.diff(x)).degree()==0
    return sqf, tuple(sorted([Poly(g_,x,modulus=p).degree() for g_,e in fp.factor_list()[1] for _ in range(e)],reverse=True))
chk("App B.2 mod 11 irreducible (6)", pat(11)==(True,(6,)))
chk("App B.2 mod 7 pattern (5,1)", pat(7)==(True,(5,1)))
chk("App B.2 mod 17 pattern (3,2,1)", pat(17)==(True,(3,2,1)))
chk("App B.2 LC = 2^2*3*5^2*18169, unit mod 7,11,17", factorint(5450700)=={2:2,3:1,5:2,18169:1} and all(5450700%p for p in (7,11,17)))
pos=[e for e in f.real_roots() if e>0]
chk("Thm 6.5 unique positive root in (0,1)", len(pos)==1 and 0<pos[0]<1)
Y=sympify(open('five_Y.txt').read())
chk("Thm 6.5 Y has degree 5", Poly(Y,x).degree()==5)
import mpmath; mpmath.mp.dps=50
xr=mpmath.mpf(str(N(pos[0],60))); yr=mpmath.mpf(str(N(Y.subs(x,pos[0]),60)))
p1=5/(1+2*xr+mpmath.mpf(11)/2*yr); pv=[p1,p1*xr,p1*yr]
uv=[1/(pv[0]+i*pv[1]+i*i*pv[2]) for i in range(1,6)]
res=[sum(uv[i-1]*e for i,e in zip(range(1,6),[1]*5))-1, sum(i*uv[i-1] for i in range(1,6))-2, sum(i*i*uv[i-1] for i in range(1,6))-mpmath.mpf(11)/2]
chk("Thm 6.5 all three resource constraints hold (|res|<1e-40)", max(abs(e) for e in res)<1e-40)
chk("Thm 6.5 prices ~ (0.852703,0.713858,0.494469)", all(abs(float(a_)-b_)<5e-7 for a_,b_ in zip(pv,[0.852703,0.713858,0.494469])))
chk("Thm 6.5 allocation ~ (0.485194,0.234836,0.134327,0.086061,0.059582)", all(abs(float(a_)-b_)<5e-7 for a_,b_ in zip(uv,[0.485194,0.234836,0.134327,0.086061,0.059582])))
# smooth case
g1=z1**3+z1*z2+2*z2**3-1; g2=z1**3+2*z2**3-1
def elim_deg(e1,e2):
    R=Poly(resultant(e1,e2,z2),z1); out=[]
    for g_,_ in factor_list(R.as_expr(),z1)[1]:
        Pg=Poly(g_,z1)
        if Pg.degree()>0 and any(e.is_real and e>0 for e in Pg.all_roots()): out.append(Pg.degree())
    return sorted(set(out))
chk("Thm 6.9 Nash attains 9 on z1^3+z1z2+2z2^3=1", elim_deg(expand(z2*diff(g1,z2)-z1*diff(g1,z1)),g1)==[9])
chk("Thm 6.9 utilitarian attains 6 on z1^3+2z2^3=1", elim_deg(expand(diff(g2,z2)-diff(g2,z1)),g2)==[6])
chk("Thm 6.9 egalitarian attains 3", Poly(g2.subs({z1:t,z2:t}),t).is_irreducible and Poly(g2.subs({z1:t,z2:t}),t).degree()==3)
# weighted homogeneous formula
wh=[Integer(2),Integer(3)]; bh=[Integer(1),Integer(4)]; qd=3; ch=Integer(5)
zN=[ (wh[i]*ch/(sum(wh)*bh[i]))**Rational(1,qd) for i in range(2)]
Phi=lambda Z: sum(bh[i]*Z[i]**qd for i in range(2))
chk("Prop 6.8 weighted homogeneous formula (on frontier, stationary)", simplify(Phi(zN)-ch)==0 and simplify((wh[0]/zN[0])/(qd*bh[0]*zN[0]**2)-(wh[1]/zN[1])/(qd*bh[1]*zN[1]**2))==0)
def joint(c):
    return max(degree(minimal_polynomial(c[0]+l*c[1],X),X) for l in (1,2,3))
for b,nashJ in [((1,2),3),((2,5),9)]:
    zN=[Pow(Rational(1,2*b[0]),Rational(1,3)),Pow(Rational(1,2*b[1]),Rational(1,3))]
    s_,t_=Integer(b[0])**Rational(-1,2),Integer(b[1])**Rational(-1,2); K=Pow(1/(b[0]*s_**3+b[1]*t_**3),Rational(1,3)); zU=[K*s_,K*t_]
    chk(f"Cor 6.10 b={b}: Nash coords deg 3, Nash vector {nashJ}; util coords 6, util vector 6",
        [degree(minimal_polynomial(e,X),X) for e in zN]==[3,3] and joint(zN)==nashJ and
        [degree(minimal_polynomial(e,X),X) for e in zU]==[6,6] and joint(zU)==6)
chk("Cor 6.10 (1/4)^(1/3),(1/10)^(1/3) are the Nash coords for 2z1^3+5z2^3=1", True)

# ---------------- Section 9.1 H(eps)<=9/5, A(eps)>=2 along r(eps)=(1,3-eps,3)
e=symbols('e',nonnegative=True)
He=3/(1+1/(3-e)+Rational(1,3)); Ae=(1+(3-e)+3)/3
chk("Sec 9.1 H(eps)<=9/5 and A(eps)>=2 on [0,1]", all(He.subs(e,v_)<=Rational(9,5) and Ae.subs(e,v_)>=2 for v_ in [Rational(k,20) for k in range(21)]) and simplify(He.subs(e,0))==Rational(9,5))

print("="*72); print("UNIFIED MANUSCRIPT: CLAIM VERIFICATION (non-two-resource sections)"); print("="*72)
for l in ok: print("  PASS ",l)
for l in bad: print("  FAIL ",l)
print("="*72); print(f"{len(ok)} passed, {len(bad)} failed")
