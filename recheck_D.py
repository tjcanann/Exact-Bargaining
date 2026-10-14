"""Independent recheck, Part D: the two-resource characterisation (Section 7),
capacity expansion (Section 8), and every number in Sections 8.6 and 9."""
import itertools, random, sympy as sp
from sympy import Rational as Q, sqrt

OK=[];BAD=[]
def chk(name,cond,extra=""):
    (OK if cond else BAD).append(name)
    print(("  PASS  " if cond else "  FAIL  ")+name+(("   "+str(extra)) if extra else ""))

T,Z,t=sp.symbols('T Z t')

def setup(Rj,Wj):
    m=len(Rj);W=sum(Wj)
    D=sp.prod([T+R for R in Rj])
    N=sum(Wj[j]*sp.prod([T+Rj[k] for k in range(m) if k!=j]) for j in range(m))
    E=sp.expand((T+Z)*N-W*D)
    H=W/sum(Q(Wj[j],1)/Rj[j] for j in range(m))
    A=sum(Wj[j]*Rj[j] for j in range(m))/W
    S=sum(Q(Wj[j],1)/(t+Rj[j]) for j in range(m))
    g=W/S-t
    L=sum(Q(Wj[j],1)/(t+Rj[j])**2 for j in range(m))
    return dict(m=m,W=W,D=D,N=N,E=E,H=H,A=A,S=S,g=g,L=L,Rj=Rj,Wj=Wj)

print("=== Section 7.1--7.2  scalar reduction, identities checked symbolically")
for Rj,Wj in [([Q(1),Q(2),Q(4)],[Q(1)]*3),([Q(1),Q(3),Q(7),Q(11)],[Q(2),Q(1),Q(3),Q(5)]),
              ([Q(i) for i in range(1,7)],[Q(1)]*6)]:
    P=setup(Rj,Wj);m=P['m'];W=P['W']
    chk("m=%d: R_1 < H < A < R_m"%m,Rj[0]<P['H']<P['A']<Rj[-1],(P['H'],P['A']))
    # g maps (0,inf) onto (H,A) and is increasing
    gp=sp.simplify(sp.diff(P['g'],t))
    chk("  g'(t) = W L/S^2 - 1 and is positive on (0,inf)",
        sp.simplify(gp-(W*P['L']/P['S']**2-1))==0 and all(sp.simplify(gp.subs(t,Q(v,7)))>0 for v in [1,3,7,20,100]))
    chk("  g(0) = H and g(t) -> A",sp.simplify(sp.limit(P['g'],t,0,'+')-P['H'])==0 and
        sp.simplify(sp.limit(P['g'],t,sp.oo)-P['A'])==0)
    # E(t,g(t)) = 0 identically: the scalar reduction is the eliminant
    chk("  E(t, g(t)) = 0 identically",sp.simplify(P['E'].subs({T:t,Z:P['g']}))==0)
    # leading coefficient W(Z-A), degree m-1
    EP=sp.Poly(P['E'],T)
    chk("  deg_T E = m-1 = %d with leading coefficient W(Z-A)"%(m-1),
        EP.degree()==m-1 and sp.simplify(EP.LC()-W*(Z-P['A']))==0,(EP.degree(),sp.factor(EP.LC())))
    chk("  E is irreducible in Q(Z)[T]",sp.Poly(P['E'],T,domain=sp.QQ.frac_field(Z)).is_irreducible)
    chk("  N(-R_j) != 0 for every j, so gcd(N,D)=1",
        all(sp.simplify(P['N'].subs(T,-R))!=0 for R in Rj) and sp.gcd(sp.Poly(P['N'],T),sp.Poly(P['D'],T)).degree()==0)
    # eq:param  and the KKT conditions
    rho=sp.Symbol('rho',positive=True)
    xs=[Wj[j]*(t+rho)/(W*(t+Rj[j])) for j in range(m)]
    mu=W/(t+rho);lam=W*t/(t+rho)
    chk("  sum x_i = 1 on the branch rho = g(t)",
        sp.simplify((sum(xs)-1).subs(rho,P['g']))==0)
    chk("  sum R_j x_j = rho on the branch",
        sp.simplify((sum(Rj[j]*xs[j] for j in range(m))-rho).subs(rho,P['g']))==0)
    chk("  stationarity W_j/x_j = lam + mu R_j",
        all(sp.simplify(Wj[j]/xs[j]-(lam+mu*Rj[j]))==0 for j in range(m)))
    chk("  lam + mu rho = W and t = lam/mu",
        sp.simplify(lam+mu*rho-W)==0 and sp.simplify(lam/mu-t)==0)
    chk("  x_i = w_i/(S(t)(t+r_i)) (eq:scalarprices), lam = tS, mu = S",
        all(sp.simplify((xs[j]-Wj[j]/(P['S']*(t+Rj[j]))).subs(rho,P['g']))==0 for j in range(m)) and
        sp.simplify((mu-P['S']).subs(rho,P['g']))==0 and sp.simplify((lam-t*P['S']).subs(rho,P['g']))==0)
    # Proposition 8.1: the cutoff identity and signs
    C=sum(Wj[j]*Rj[j]/(t+Rj[j])**2 for j in range(m))/P['L']
    chk("  C(t) = S/L - t",sp.simplify(C-(P['S']/P['L']-t))==0)
    chk("  R_1 < C(t) < g(t) < A at sample t",
        all(Rj[0]<sp.simplify(C.subs(t,Q(v,5)))<sp.simplify(P['g'].subs(t,Q(v,5)))<P['A'] for v in [1,4,13,60]))
    dxdrho=[sp.simplify(sp.diff(xs[j].subs(rho,P['g']),t)/sp.diff(P['g'],t)) for j in range(m)]
    paper=[sp.simplify(xs[j].subs(rho,P['g'])/sp.diff(P['g'],t)*(P['L']/P['S']-1/(t+Rj[j]))) for j in range(m)]
    chk("  dx_i/drho matches eq:dxdrho",all(sp.simplify(a-b)==0 for a,b in zip(dxdrho,paper)))
    chk("  sign(dx_i/drho) = sign(r_i - C(t)) at sample t",
        all(sp.sign(sp.simplify(dxdrho[j].subs(t,Q(v,5))))==sp.sign(sp.simplify(Rj[j]-C.subs(t,Q(v,5))))
            for j in range(m) for v in [1,4,13]))
    lamb_=lam.subs(rho,P['g']);mu_=mu.subs(rho,P['g'])
    chk("  lambda rises and mu falls with rho",
        all(sp.simplify(sp.diff(lamb_,t).subs(t,Q(v,5)))>0 and sp.simplify(sp.diff(mu_,t).subs(t,Q(v,5)))<0
            for v in [1,4,13,60]))
    chk("  mu_t = -L and lam_t = S - tL = L C",
        sp.simplify(sp.diff(mu_,t)+P['L'])==0 and sp.simplify(sp.diff(lamb_,t)-(P['S']-t*P['L']))==0
        and sp.simplify(P['S']-t*P['L']-P['L']*C)==0)
    # Proposition 8.3: C' > 0, C(0) = C_0, C(inf) = A
    M3=sum(Wj[j]/(t+Rj[j])**3 for j in range(m))
    chk("  C'(t) = 2(S M_3 - L^2)/L^2 > 0",
        sp.simplify(sp.diff(C,t)-2*(P['S']*M3-P['L']**2)/P['L']**2)==0 and
        all(sp.simplify(sp.diff(C,t).subs(t,Q(v,5)))>0 for v in [1,4,13,60]))
    C0=sum(Q(Wj[j],1)/Rj[j] for j in range(m))/sum(Q(Wj[j],1)/Rj[j]**2 for j in range(m))
    chk("  C(0) = C_0 and C(t) -> A",
        sp.simplify(sp.limit(C,t,0,'+')-C0)==0 and sp.simplify(sp.limit(C,t,sp.oo)-P['A'])==0)
    # matching capacity derivative
    for j in range(m):
        if P['H']<Rj[j]<P['A']:
            tm=sp.solve(sp.Eq(P['g'],Rj[j]),t)
            tm=[s for s in tm if s.is_real and s>0]
            if tm:
                val=sp.simplify(dxdrho[j].subs(t,tm[0]))
                chk("  at a matching capacity rho=R_%d, dx/drho = w/(W(t+r))"%(j+1),
                    sp.simplify(val-Wj[j]/(W*(tm[0]+Rj[j])))==0 and val>0,val)
                chk("  and x_j = w_j/W there",sp.simplify(xs[j].subs({rho:Rj[j],t:tm[0]})-Wj[j]/W)==0)
                break
    print("")

print("=== Proposition 7.6 / prop:twotype")
R1,R2,rh=sp.symbols('R1 R2 rh')
X1=(R2-rh)/(R2-R1);X2=(rh-R1)/(R2-R1)
chk("X1+X2 = 1 and R1 X1 + R2 X2 = rho",sp.simplify(X1+X2-1)==0 and sp.simplify(R1*X1+R2*X2-rh)==0)

print("\n=== Section 7.4 / the six-type technology and f_5")
P6=setup([Q(i) for i in range(1,7)],[Q(1)]*6)
chk("H = 120/49 and A = 7/2",P6['H']==Q(120,49) and P6['A']==Q(7,2),(P6['H'],P6['A']))
chk("N(T) = D'(T) for equal weights",sp.expand(P6['N']-sp.diff(P6['D'],T))==0)
f5=9*T**5+70*T**4-175*T**3-2800*T**2-7868*T-6624
chk("f_5 = -5 E(T,16/5)",sp.expand(f5+5*P6['E'].subs(Z,Q(16,5)))==0)
chk("f_5(6) = -31728 and f_5(7) = 60408",f5.subs(T,6)==-31728 and f5.subs(T,7)==60408,
    (f5.subs(T,6),f5.subs(T,7)))
rr=sp.Poly(f5,T).real_roots()
chk("f_5 has five real roots, exactly one positive",len(rr)==5 and sum(1 for r in rr if r>0)==1)
ivs=[(-6,-5),(-5,-4),(-3,-2),(-2,-1),(6,7)]
chk("the five roots lie one apiece in (-6,-5),(-5,-4),(-3,-2),(-2,-1),(6,7)",
    all(any(a<r<b for r in rr) for a,b in ivs) and
    all(sum(1 for r in rr if a<r<b)==1 for a,b in ivs),[float(r) for r in rr])
chk("16/5 lies in I = (120/49, 7/2)",P6['H']<Q(16,5)<P6['A'])
tr=[r for r in rr if r>0][0]
u6=[(tr+Q(16,5))/(6*(tr+i)) for i in range(1,7)]
chk("allocation formula u_i = (t+16/5)/(6(t+i)) is feasible and clears both resources",
    sp.simplify(sum(u6)-1)==0 and sp.simplify(sum(i*u6[i-1] for i in range(1,7))-Q(16,5))==0)
pa=6*tr/(tr+Q(16,5));pb=6/(tr+Q(16,5))
chk("prices p_a = 6t/(t+16/5), p_b = 6/(t+16/5) are positive with p_a + (16/5)p_b = 6",
    pa>0 and pb>0 and sp.simplify(pa+Q(16,5)*pb-6)==0)
yv=sp.Symbol('yv')
chk("every allocation coordinate and both prices have degree 5",
    all(sp.Poly(sp.minimal_polynomial(z,yv),yv).degree()==5 for z in u6+[pa,pb]))
chk("t recovered from the prices: t = 6/p_b - 16/5 = (16/5)p_a/(6-p_a)",
    sp.simplify(6/pb-Q(16,5)-tr)==0 and sp.simplify(Q(16,5)*pa/(6-pa)-tr)==0)

print("\n=== Theorem 7.13 / thm:family  congruence family")
chk("247 = 13 * 19",sp.factorint(247)=={13:1,19:1})
def Fpq(p_,q_): return sp.expand((q_*T+p_)*P6['N']-6*q_*P6['D'])
tested=0;good=0
for (p_,q_) in [(16,5),(5223,5),(16+247*3,5),(3*247+16,5),(1731,538),(4193,1292)]:
    g_=sp.gcd(p_,q_)
    p0,q0=p_//g_,q_//g_
    if sp.gcd(q0,247)!=1 or (5*p0-16*q0)%247!=0: continue
    rho=Q(p0,q0)
    if not (P6['H']<rho<P6['A']): continue
    tested+=1
    poly=sp.Poly(sp.primitive(Fpq(p0,q0))[1],T)
    irr=poly.is_irreducible and poly.degree()==5
    G=poly.galois_group()[0] if irr else None
    ok=irr and G.order()==120 and not G.is_solvable
    good+= 1 if ok else 0
    chk("rho = %s satisfies the congruence: irreducible quintic with S_5"%rho,ok,
        "degree %d, |G| = %s"%(poly.degree(),G.order() if G else "-"))
chk("every congruence capacity tested gave S_5",tested>0 and good==tested,"%d/%d"%(good,tested))
# density construction
for z in [Q(5,2),Q(13,4),Q(26,10)]:
    worst=0
    for h in range(1,60):
        s_=247*h+1
        k_=sp.floor(s_*(z-Q(16,5))/247+Q(1,2))
        val=Q(16,5)+247*Q(k_,s_)
        worst=max(worst,abs(float(val-z)))
        if h>40: break
    chk("density: |16/5 + 247 k/s - %s| -> 0 along s = 247h+1"%z,abs(float(val-z))<=247/(2*(247*40+1)))

print("\n=== Appendix A / prop:matching and lem:realroots")
for Rj,Wj in [([Q(1),Q(2),Q(4)],[Q(1)]*3),([Q(1),Q(3),Q(7),Q(11)],[Q(2),Q(1),Q(3),Q(5)])]:
    P=setup(Rj,Wj);m=P['m']
    for k in range(m):
        if not (P['H']<Rj[k]<P['A']): continue
        Ek=sp.Poly(P['E'].subs(Z,Rj[k]),T)
        gg=sp.gcd(Ek,sp.Poly(P['D'],T))
        chk("m=%d: at a matching rho = R_%d, gcd(E,D) = T + R_%d"%(m,k+1,k+1),
            sp.simplify(gg.as_expr()-(T+Rj[k]))==0,gg.as_expr())
    # non-matching: gcd = 1, and the sign pattern of lem:realroots
    rho=Q(P['H']+P['A'],2) if False else sp.nsimplify((P['H']+P['A'])/2)
    Ek=sp.Poly(sp.expand(P['E'].subs(Z,rho)),T)
    chk("m=%d: at a non-matching rational rho, gcd(E,D) = 1"%m,sp.gcd(Ek,sp.Poly(P['D'],T)).degree()==0)
    signs=[sp.sign(sp.simplify(P['E'].subs({T:-R,Z:rho}))) for R in Rj]
    want=[sp.sign(rho-Rj[j])*(-1)**j for j in range(m)]
    chk("m=%d: sgn E(-R_j,rho) = sgn(rho-R_j)(-1)^(j-1)"%m,signs==want,(signs,want))
    rts=Ek.real_roots()
    chk("m=%d: E(T,rho) has m-1 = %d distinct real roots, one positive"%(m,m-1),
        len(rts)==m-1 and len(set(rts))==m-1 and sum(1 for r in rts if r>0)==1,[float(r) for r in rts])

print("\n"+"="*72)
print("Part D: %d passed, %d failed"%(len(OK),len(BAD)))
for b in BAD: print("   FAILED:",b)
