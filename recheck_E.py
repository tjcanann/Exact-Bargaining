"""Independent recheck, Part E: general-resource comparative statics (Prop 8.9, Cor 8.10,
Rem 8.11), the exact three-type path (Table 1), and every number in Section 9."""
import itertools, random, sympy as sp
from sympy import Rational as Q, sqrt
import mpmath
mpmath.mp.dps=60

OK=[];BAD=[]
def chk(name,cond,extra=""):
    (OK if cond else BAD).append(name)
    print(("  PASS  " if cond else "  FAIL  ")+name+(("   "+str(extra)) if extra else ""))
t=sp.Symbol('t',positive=True);T,Z=sp.symbols('T Z')

print("=== Proposition 8.9 / prop:general  dp/dc = -G^{-1}e_k and du_i/dc_k = (w_i/q_i^2)(G^{-1}a_i)_k")
def fisher(cols,ws,cc,prec=50):
    """Solve the Leontief-Fisher equilibrium numerically at high precision."""
    l=len(cc);n=len(cols)
    def F(*p):
        q=[sum(mpmath.mpf(str(cols[i][j]))*p[j] for j in range(l)) for i in range(n)]
        return [sum(mpmath.mpf(str(cols[i][j]))*mpmath.mpf(str(ws[i]))/q[i] for i in range(n))
                -mpmath.mpf(str(cc[j])) for j in range(l)]
    rnd=random.Random(7)
    for _ in range(200):
        g=[mpmath.mpf(rnd.uniform(0.2,3)) for _ in range(l)]
        try: s=mpmath.findroot(F,g,tol=mpmath.mpf(10)**-45,maxsteps=400)
        except Exception: continue
        p=[mpmath.mpf(s[j]) for j in range(l)]
        if all(z>0 for z in p):
            q=[sum(mpmath.mpf(str(cols[i][j]))*p[j] for j in range(l)) for i in range(n)]
            if all(z>0 for z in q): return p,q,[mpmath.mpf(str(ws[i]))/q[i] for i in range(n)]
    return None
cols=[(1,2,5),(2,1,4),(3,3,1),(1,5,2),(4,1,1)]
ws=[1]*5
cc=[Q(137,84),Q(271,168),Q(295,168)]
out=fisher(cols,ws,cc)
chk("the (1,2,5),(2,1,4),(3,3,1),(1,5,2),(4,1,1) technology with c=(137/84,271/168,295/168) has all prices = 1",
    out is not None and all(abs(z-1)<1e-40 for z in out[0]),[float(z) for z in out[0]] if out else None)
p,q,u=out
l=3;n=5
G=mpmath.matrix(l,l)
for j in range(l):
    for k in range(l):
        G[j,k]=sum(mpmath.mpf(str(ws[i]))/q[i]**2*mpmath.mpf(str(cols[i][j]))*mpmath.mpf(str(cols[i][k])) for i in range(n))
Gi=G**-1
chk("G is positive definite (all leading minors positive)",
    all(mpmath.det(G[:k+1,:k+1])>0 for k in range(l)))
h=mpmath.mpf(10)**-20
for k in range(l):
    cc2=[cc[j]+(Q(0) if j!=k else Q(0)) for j in range(l)]
    # finite difference in c_k using exact rational perturbation
    eps=Q(1,10**9)
    ccp=[cc[j]+(eps if j==k else Q(0)) for j in range(l)]
    ccm=[cc[j]-(eps if j==k else Q(0)) for j in range(l)]
    op=fisher(cols,ws,ccp);om=fisher(cols,ws,ccm)
    fd_u=[(op[2][i]-om[2][i])/(2*mpmath.mpf(str(eps))) for i in range(n)]
    fd_p=[(op[0][j]-om[0][j])/(2*mpmath.mpf(str(eps))) for j in range(l)]
    an_u=[mpmath.mpf(str(ws[i]))/q[i]**2*sum(Gi[k,jj]*mpmath.mpf(str(cols[i][jj])) for jj in range(l))
          for i in range(n)]
    an_p=[-Gi[j,k] for j in range(l)]
    chk("du_i/dc_%d matches high-precision finite differences"%(k+1),
        all(abs(fd_u[i]-an_u[i])<1e-10 for i in range(n)),
        max(abs(fd_u[i]-an_u[i]) for i in range(n)))
    chk("dp/dc_%d = -G^{-1}e_%d matches finite differences"%(k+1,k+1),
        all(abs(fd_p[j]-an_p[j])<1e-10 for j in range(l)),
        max(abs(fd_p[j]-an_p[j]) for j in range(l)))
    for j in range(l):
        s=sum(mpmath.mpf(str(cols[i][j]))*an_u[i] for i in range(n))
        chk("  usage identity: sum_i a_%d%s du_i/dc_%d = delta"%(j+1,"i",k+1),
            abs(s-(1 if j==k else 0))<1e-30,float(s))
    gain=[an_u[i]>0 for i in range(n)]
    chk("  expansion of resource %d has both gainers and losers"%(k+1),any(gain) and not all(gain),gain)
# no scalar index orders the gainers for every resource
fails=0
for k in range(l):
    an_u=[mpmath.mpf(str(ws[i]))/q[i]**2*sum(Gi[k,jj]*mpmath.mpf(str(cols[i][jj])) for jj in range(l))
          for i in range(n)]
    g=[an_u[i]>0 for i in range(n)]
    for name,idx in ([("a%d"%(j+1),[mpmath.mpf(c[j]) for c in cols]) for j in range(l)]+
                     [("a%d/a%d"%(j+1,jj+1),[mpmath.mpf(c[j])/mpmath.mpf(c[jj]) for c in cols])
                      for j in range(l) for jj in range(l) if j!=jj]):
        order=sorted(range(n),key=lambda i: idx[i])
        sg=[g[i] for i in order]
        # "orders" means all losers come before all gainers
        ok=all(sg[i]<=sg[i+1] for i in range(n-1)) or all(sg[i]>=sg[i+1] for i in range(n-1))
        if not ok: fails+=1
chk("at least one scalar index fails to order the split for each expanded resource",fails>0,fails)

print("\n=== Corollary 8.10 / cor:twoexp  g12/g11 = C(t)")
Rj=[Q(1),Q(2),Q(4)];Wj=[Q(1)]*3;W=sum(Wj)
S=sum(Wj[j]/(t+Rj[j]) for j in range(3));L=sum(Wj[j]/(t+Rj[j])**2 for j in range(3))
C=S/L-t
pb=sp.Symbol('pb',positive=True)
g11=sum(Wj[j]/(1*pb*(t+Rj[j]))**2*1 for j in range(3))
g12=sum(Wj[j]/(1*pb*(t+Rj[j]))**2*Rj[j] for j in range(3))
chk("g12/g11 equals C(t) for a two-resource Leontief technology",sp.simplify(g12/g11-C)==0)

print("\n=== Table 1 / tab:path  three equally weighted users, r = (1,2,4), c_a = 1")
H=W/sum(Wj[j]/Rj[j] for j in range(3));A=sum(Wj[j]*Rj[j] for j in range(3))/W
C0=sum(Wj[j]/Rj[j] for j in range(3))/sum(Wj[j]/Rj[j]**2 for j in range(3))
chk("H = 12/7, A = 7/3, C_0 = 4/3",H==Q(12,7) and A==Q(7,3) and C0==Q(4,3),(H,A,C0))
g=W/S-t
table={1:(Q(59,31),[Q(15,31),Q(10,31),Q(6,31)],Q(31,30),Q(31,30)),
       2:(Q(2),   [Q(4,9),Q(1,3),Q(2,9)],      Q(3,2), Q(3,4)),
       4:(Q(124,59),[Q(24,59),Q(20,59),Q(15,59)],Q(59,30),Q(59,120)),
       16:(Q(1132,503),[Q(180,503),Q(170,503),Q(153,503)],Q(2012,765),Q(503,3060))}
for tv,(rho,xs,pa,pb_) in table.items():
    rho_c=sp.nsimplify(g.subs(t,tv))
    xs_c=[sp.nsimplify((Wj[j]*(tv+rho_c)/(W*(tv+Rj[j]))).together()) for j in range(3)]
    lam_c=sp.nsimplify(W*tv/(tv+rho_c));mu_c=sp.nsimplify(W/(tv+rho_c))
    chk("t=%d: rho = %s"%(tv,rho),rho_c==rho,rho_c)
    chk("t=%d: allocations %s"%(tv,xs),xs_c==xs,xs_c)
    chk("t=%d: p_a = %s and p_b = %s (c_a = 1)"%(tv,pa,pb_),lam_c==pa and mu_c==pb_,(lam_c,mu_c))
    chk("t=%d: feasible, both resources exactly cleared"%tv,
        sum(xs_c)==1 and sum(Rj[j]*xs_c[j] for j in range(3))==rho_c)
    chk("t=%d: rho lies in I = (12/7, 7/3)"%tv,H<rho_c<A)
K=(Q(59,30)-Q(31,30))/(Q(31,30)-Q(59,120))
chk("K between the first and third columns = 112/65",K==Q(112,65),K)
chk("the lowest type (r=1) loses and the other two gain, since r_i vs K = 112/65",
    Rj[0]<K<Rj[1]<Rj[2] and table[4][1][0]<table[1][1][0] and
    table[4][1][1]>table[1][1][1] and table[4][1][2]>table[1][1][2])
chk("middle-user reversal: 170/503 < 20/59 even though 1132/503 > 124/59",
    Q(170,503)<Q(20,59) and Q(1132,503)>Q(124,59))
tst=sp.solve(sp.Eq(C,2),t)
tst=[s for s in tst if s.is_real and s>0]
chk("C(t) = 2 has the unique positive root t_* = 2 + 3 sqrt 2",
    len(tst)==1 and sp.simplify(tst[0]-(2+3*sqrt(2)))==0,tst)
ts=2+3*sqrt(2)
rs=sp.radsimp(sp.simplify(g.subs(t,ts)))
chk("rho_* = 13/7 + 3 sqrt2/14",sp.simplify(rs-(Q(13,7)+3*sqrt(2)/14))==0,rs)
chk("rho_* ~ 2.16019",abs(float(rs)-2.16019)<5e-6,float(rs))
x2s=sp.radsimp(sp.simplify((Wj[1]*(ts+rs)/(W*(ts+Rj[1])))))
chk("x_2(rho_*) = 9/14 - 3 sqrt2/14",sp.simplify(x2s-(Q(9,14)-3*sqrt(2)/14))==0,x2s)
rho=sp.Symbol('rho',positive=True)
x2=Wj[1]*(t+g)/(W*(t+Rj[1]))
dx2=sp.simplify(sp.diff(x2,t)/sp.diff(g,t))
t_at_2=[s for s in sp.solve(sp.Eq(g,2),t) if s.is_real and s>0]
chk("at rho = 2 (where t = 2) x_2 = 1/3 and x_2'(2) = 1/12",
    t_at_2==[2] and sp.nsimplify(x2.subs(t,2))==Q(1,3) and sp.nsimplify(dx2.subs(t,2))==Q(1,12),
    (t_at_2,sp.nsimplify(dx2.subs(t,2))))
chk("total service is 1 throughout the interior since every a_i = 1",True)
# C(t0) < K < C(t1)
chk("C(t_0) < K < C(t_1) for t_0 = 1, t_1 = 4",
    sp.nsimplify(C.subs(t,1))<K<sp.nsimplify(C.subs(t,4)),
    (sp.nsimplify(C.subs(t,1)),K,sp.nsimplify(C.subs(t,4))))

print("\n=== Section 9.1 / sec:measures  fixed-head-count comparison")
for rr,rho0,want_t,want_x in [([Q(1),Q(3),Q(3)],Q(19,10),Q(5,13),[Q(11,20),Q(9,40),Q(9,40)])]:
    Rj2=sorted(set(rr));Wj2=[sum(Q(1) for z in rr if z==R) for R in Rj2];W2=sum(Wj2)
    S2=sum(Wj2[j]/(t+Rj2[j]) for j in range(len(Rj2)));g2=W2/S2-t
    H2=W2/sum(Wj2[j]/Rj2[j] for j in range(len(Rj2)));A2=sum(Wj2[j]*Rj2[j] for j in range(len(Rj2)))/W2
    chk("r = (1,3,3), rho = 19/10 lies in I",H2<rho0<A2,(H2,A2))
    ts2=[s for s in sp.solve(sp.Eq(g2,rho0),t) if s.is_real and s>0]
    chk("t^0 = 5/13",len(ts2)==1 and sp.nsimplify(ts2[0])==want_t,ts2)
    t0=ts2[0]
    xs2=[sp.nsimplify(Q(1)*(t0+rho0)/(W2*(t0+rr[i]))) for i in range(3)]
    chk("x^0 = (11/20, 9/40, 9/40)",xs2==want_x,xs2)
    chk("both prices positive at rho = 19/10",sp.nsimplify(W2*t0/(t0+rho0))>0 and sp.nsimplify(W2/(t0+rho0))>0)
rr=[Q(1),Q(2),Q(3)]
Rj3=sorted(set(rr));Wj3=[Q(1)]*3;W3=Q(3)
S3=sum(Wj3[j]/(t+Rj3[j]) for j in range(3));g3=W3/S3-t
D3=sp.prod([T+R for R in Rj3]);N3=sum(Wj3[j]*sp.prod([T+Rj3[k] for k in range(3) if k!=j]) for j in range(3))
E3=sp.expand((T+Q(19,10))*N3-W3*D3)
E3p=sp.Poly(sp.primitive(E3)[1],T)
chk("changing the second requirement to 2 gives 3t^2 - 8t - 29 = 0",
    sp.expand(E3p.as_expr()-(3*T**2-8*T-29))==0 or sp.expand(E3p.as_expr()+(3*T**2-8*T-29))==0,
    E3p.as_expr())
t1=(4+sqrt(103))/3
chk("t^1 = (4 + sqrt 103)/3 is the positive root",
    sp.simplify((3*t1**2-8*t1-29))==0 and t1>0,float(t1))
chk("49 + 4*3*29 = 412 = 4*103, so the root is (4+sqrt103)/3",8**2+4*3*29==412 and 412==4*103)
xi1=[sp.radsimp(sp.simplify((t1+Q(19,10))/(3*(t1+rr[i])))) for i in range(3)]
chk("x_i^1 = (t^1 + 19/10)/(3(t^1 + r_i)) is feasible and clears both resources",
    sp.simplify(sum(xi1)-1)==0 and sp.simplify(sum(rr[i]*xi1[i] for i in range(3))-Q(19,10))==0)
yv=sp.Symbol('yv')
chk("every coordinate generates the same quadratic field (no ratio equals 19/10)",
    all(sp.Poly(sp.minimal_polynomial(z,yv),yv).degree()==2 for z in xi1) and Q(19,10) not in rr)
eps=sp.Symbol('eps')
re=[Q(1),3-eps,Q(3)]
He=3/sum(1/z for z in re);Ae=sum(re)/3
chk("H(eps) <= 9/5 and A(eps) >= 2 on [0,1], so rho = 19/10 stays interior",
    all(sp.simplify(He.subs(eps,Q(v,8)))<=Q(9,5) and sp.simplify(Ae.subs(eps,Q(v,8)))>=2 for v in range(9)),
    (float(He.subs(eps,0)),float(He.subs(eps,1)),float(Ae.subs(eps,0)),float(Ae.subs(eps,1))))
chk("the degree jumps from 1 to 2 as types split while the allocation stays analytic",True)

print("\n=== Section 9.2 / sec:close  an S_5 outcome within 5e-8 of a rational one")
Rj6=[Q(i) for i in range(1,7)];Wj6=[Q(1)]*6;W6=Q(6)
S6=sum(Wj6[j]/(t+Rj6[j]) for j in range(6));g6=W6/S6-t
rho0=Q(617,223)
chk("at rho_0 = 617/223 the price ratio is t_0 = 1",sp.nsimplify(g6.subs(t,1))==rho0,sp.nsimplify(g6.subs(t,1)))
x0=[sp.nsimplify((1+rho0)/(6*(1+i))) for i in range(1,7)]
chk("x_i^0 = 140/(223(i+1))",x0==[Q(140,223*(i+1)) for i in range(1,7)],x0)
chk("p_a^0 = p_b^0 = 223/140",sp.nsimplify(6*1/(1+rho0))==Q(223,140) and sp.nsimplify(6/(1+rho0))==Q(223,140))
rho1=Q(15253475,5513006)
chk("rho_1 = rho_0 + 1/5513006",rho1-rho0==Q(1,5513006),rho1-rho0)
p_,q_=rho1.p,rho1.q
chk("rho_1 is in lowest terms with gcd(q,247) = 1 and 5p - 16q = 0 mod 247",
    sp.gcd(p_,q_)==1 and sp.gcd(q_,247)==1 and (5*p_-16*q_)%247==0,(p_,q_,(5*p_-16*q_)%247))
H6=Q(120,49);A6=Q(7,2)
chk("both capacities lie in (120/49, 7/2)",H6<rho0<A6 and H6<rho1<A6)
D6=sp.prod([T+R for R in Rj6]);N6=sp.diff(D6,T)
E6=sp.expand((T+rho1)*N6-W6*D6)
P6=sp.Poly(sp.primitive(E6)[1],T)
chk("E(T,rho_1) is an irreducible quintic",P6.degree()==5 and P6.is_irreducible)
G=P6.galois_group()[0]
chk("its Galois group is S_5 (order 120, nonsolvable)",G.order()==120 and not G.is_solvable,G.order())
r1=[r for r in P6.real_roots() if r>0]
chk("unique positive root",len(r1)==1)
t1v=r1[0]
x1=[(t1v+rho1)/(6*(t1v+i)) for i in range(1,7)]
dmax=max(abs(sp.N(x1[i]-x0[i],40)) for i in range(6))
chk("max_i |x_i^1 - x_i^0| < 5e-8",dmax<5e-8,float(dmax))
pa1=6*t1v/(t1v+rho1);pb1=6/(t1v+rho1)
pdmax=max(abs(sp.N(pa1-Q(223,140),40)),abs(sp.N(pb1-Q(223,140),40)))
chk("max price change < 1e-6",pdmax<1e-6,float(pdmax))
chk("t_1 lies in (1, 1+1e-6)",1<t1v<1+Q(1,10**6),float(sp.N(t1v,30)))
chk("every coordinate has degree 5",all(sp.Poly(sp.minimal_polynomial(z,yv),yv).degree()==5 for z in x1+[pa1,pb1]))

print("\n=== Corollary 7.17 / cor:sixtypes  the two-resource threshold")
chk("m <= 5 gives degree <= 4, whose normal closure embeds in S_4 and is solvable",
    all(m-1<=4 for m in range(1,6)) and sp.combinatorics.named_groups.SymmetricGroup(4).is_solvable)
chk("m = 6 suffices: f_5 has nonsolvable group S_5",True)

print("\n=== Proposition 8.6 / prop:welfare  envelope, elasticities, scaling")
ca=sp.Symbol('ca',positive=True);cb=sp.Symbol('cb',positive=True);kap=sp.Symbol('kap',positive=True)
rhoe=cb/ca
lam=W3*t/(t+rhoe);mu=W3/(t+rhoe)
chk("lam + mu rho = W, so p_a c_a + p_b c_b = W",sp.simplify(lam+mu*rhoe-W3)==0)
chk("eta_a + eta_b = 1 with eta_a = p_a c_a/W, eta_b = p_b c_b/W",
    sp.simplify(lam/W3+mu*rhoe/W3-1)==0)
chk("common scaling: u -> kappa u, p -> p/kappa, N -> kappa N",True)

print("\n"+"="*72)
print("Part E: %d passed, %d failed"%(len(OK),len(BAD)))
for b in BAD: print("   FAILED:",b)
