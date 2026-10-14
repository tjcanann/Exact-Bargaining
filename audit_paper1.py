"""Audit of every numerical claim in 'Rationality and Heterogeneity in Two-Resource
Nash Bargaining' (the attached manuscript)."""
from sympy import *
from sympy.polys.numberfields.galoisgroups import galois_group
T,Z,X,t=symbols('T Z X t')
ok=[];bad=[]
def chk(l,c): (ok if bool(c) else bad).append(l)

def theory(R,Wt):
    W=sum(Wt)
    S=sum(Wt[j]/(t+R[j]) for j in range(len(R)))
    L=sum(Wt[j]/(t+R[j])**2 for j in range(len(R)))
    g=W/S-t
    H=W/sum(Wt[j]/R[j] for j in range(len(R))); A=sum(Wt[j]*R[j] for j in range(len(R)))/W
    D=prod([T+r for r in R]); N=sum(Wt[j]*prod([T+R[l] for l in range(len(R)) if l!=j]) for j in range(len(R)))
    E=expand((T+Z)*N-W*D)
    return W,S,L,g,H,A,D,N,E

# ---- six-type family (Sec 4)
R6=[Integer(i) for i in range(1,7)]; W6=[Integer(1)]*6
W,S,L,g,H,A,D,N,E=theory(R6,W6)
chk("H=120/49, A=7/2", H==Rational(120,49) and A==Rational(7,2))
chk("N(T)=D'(T) for equal weights", expand(N-diff(D,T))==0)
f=Poly(9*T**5+70*T**4-175*T**3-2800*T**2-7868*T-6624,T)
chk("f(T) = -5 E(T,16/5)", expand(f.as_expr()+5*E.subs(Z,Rational(16,5)))==0)
chk("coefficient of T^m in E cancels (degree m-1=5)", Poly(E,T).degree()==5)
chk("f(6)=-31728, f(7)=60408", f.eval(6)==-31728 and f.eval(7)==60408)
chk("f irreducible over Q", f.is_irreducible)
h19=Poly(f.as_expr(),T,modulus=19).monic()
chk("mod 19: f/LC = T^5-7T^4+8T^3-5T^2+4T+5", h19==Poly(T**5-7*T**4+8*T**3-5*T**2+4*T+5,T,modulus=19))
chk("mod 19 irreducible", h19.is_irreducible)
f13=Poly(f.as_expr(),T,modulus=13)
chk("mod 13: 9(T^2+5T+5)(T^3-3T^2+5T+1)",
    f13==Poly(9*(T**2+5*T+5)*(T**3-3*T**2+5*T+1),T,modulus=13))
chk("mod 13 factors irreducible & squarefree",
    Poly(T**2+5*T+5,T,modulus=13).is_irreducible and Poly(T**3-3*T**2+5*T+1,T,modulus=13).is_irreducible
    and gcd(f13,f13.diff(T)).degree()==0)
G,_=galois_group(f); chk("Galois group S5 (order 120)", G.order()==120)
rts=f.real_roots()
chk("all five roots real", len(rts)==5)
iv=[(-6,-5),(-5,-4),(-3,-2),(-2,-1),(6,7)]
chk("one root in each of (-6,-5),(-5,-4),(-3,-2),(-2,-1),(6,7)",
    all(sum(1 for r in rts if a<r<b)==1 for a,b in iv))
tt=[r for r in rts if r>0][0]; rho=Rational(16,5)
u=[(tt+rho)/(6*(tt+i)) for i in range(1,7)]
chk("(4.3) allocation satisfies both constraints at t",
    abs(float(sum(u))-1)<1e-12 and abs(float(sum(i*u[i-1] for i in range(1,7)))-16/5)<1e-12)
pa=6*tt/(tt+rho); pb=6/(tt+rho)
chk("stationarity w_i/u_i = p_a + i p_b", max(abs(float(1/u[i-1]-(pa+i*pb))) for i in range(1,7))<1e-10)
# dense family Theorem 4.2 & Sec 6.3 capacity
rho1=Rational(15253475,5513006); p_,q_=rho1.p,rho1.q
chk("rho1 = rho0 + 1/5513006 with rho0=617/223", rho1==Rational(617,223)+Rational(1,5513006))
chk("gcd(q,247)=1 and 5p-16q = 0 mod 247", gcd(q_,247)==1 and (5*p_-16*q_)%247==0)
F1=Poly(E.subs(Z,rho1),T); F1=F1.primitive()[1]
chk("E(T,rho1) irreducible quintic", F1.degree()==5 and F1.is_irreducible)
G1,_=galois_group(F1); chk("E(T,rho1) has Galois group S5", G1.order()==120)
chk("rho0=617/223 at t0=1", g.subs(t,1)==Rational(617,223))
x0=[(1+Rational(617,223))/(6*(1+i)) for i in range(1,7)]
chk("x0_i = 140/(223(i+1))", all(x0[i-1]==Rational(140,223*(i+1)) for i in range(1,7)))
chk("p0 = 223/140 both", 6*1/(1+Rational(617,223))==Rational(223,140) and 6/(1+Rational(617,223))==Rational(223,140))
# nearby S5 outcome within 5e-8
import mpmath; mpmath.mp.dps=50
t1=[r for r in Poly(E.subs(Z,rho1),T).nroots(n=40) if im(r)==0 and re(r)>0][0]
x1=[(t1+rho1)/(6*(t1+i)) for i in range(1,7)]
chk("max|x1-x0| < 5e-8", max(abs(float(x1[i]-x0[i])) for i in range(6))<5e-8)
pa1=6*t1/(t1+rho1); pb1=6/(t1+rho1)
chk("max price diff < 1e-6", max(abs(float(pa1-Rational(223,140))),abs(float(pb1-Rational(223,140))))<1e-6)
chk("rho0, rho1 both in (120/49, 7/2)", all(Rational(120,49)<r<Rational(7,2) for r in (Rational(617,223),rho1)))

# ---- three-type experiment (Sec 5.4)
R3=[Integer(1),Integer(2),Integer(4)]; W3=[Integer(1)]*3
W,S,L,g,H,A,D,N,E3=theory(R3,W3)
chk("H=12/7, A=7/3", H==Rational(12,7) and A==Rational(7,3))
C=S/L-t
chk("C(0)=C0=4/3", simplify(C.subs(t,0))==Rational(4,3))
def row(tv):
    rho=simplify(g.subs(t,tv)); x=[simplify((tv+rho)/(3*(tv+r))) for r in R3]
    return rho,x,simplify(3*tv/(tv+rho)),simplify(3/(tv+rho))
table={1:(Rational(59,31),[Rational(15,31),Rational(10,31),Rational(6,31)],Rational(31,30),Rational(31,30)),
       2:(Integer(2),[Rational(4,9),Rational(1,3),Rational(2,9)],Rational(3,2),Rational(3,4)),
       4:(Rational(124,59),[Rational(24,59),Rational(20,59),Rational(15,59)],Rational(59,30),Rational(59,120)),
       16:(Rational(1132,503),[Rational(180,503),Rational(170,503),Rational(153,503)],Rational(2012,765),Rational(503,3060))}
for tv,exp in table.items():
    got=row(tv)
    chk(f"Table 1 column t={tv}", got[0]==exp[0] and got[1]==exp[1] and got[2]==exp[2] and got[3]==exp[3])
K=(Rational(59,30)-Rational(31,30))/(Rational(31,30)-Rational(59,120))
chk("K = 112/65", K==Rational(112,65))
chk("170/503 < 20/59 and 1132/503 > 124/59", Rational(170,503)<Rational(20,59) and Rational(1132,503)>Rational(124,59))
ts=2+3*sqrt(2)
chk("C(t*)=2 at t*=2+3sqrt2", simplify(C.subs(t,ts)-2)==0)
chk("rho* = 13/7 + 3sqrt2/14", simplify(g.subs(t,ts)-(Rational(13,7)+3*sqrt(2)/14))==0)
x2s=simplify((ts+g.subs(t,ts))/(3*(ts+2)))
chk("x2(rho*) = 9/14 - 3sqrt2/14", simplify(x2s-(Rational(9,14)-3*sqrt(2)/14))==0)
# x2'(rho) at rho=2 (t=2): dx/drho = (dx/dt)/(drho/dt)
x2t=(t+g)/(3*(t+2))
chk("x2'(2) = 1/12", simplify((diff(x2t,t)/diff(g,t)).subs(t,2))==Rational(1,12))
chk("Figure 2: p_a + rho p_b = 3 (=W) on interior", simplify(3*t/(t+g)+g*3/(t+g))==3)

# ---- Sec 6.1
rho=Rational(19,10)
W,S,L,g,H,A,D,N,E=theory([Integer(1),Integer(3)],[Integer(1),Integer(2)])   # r=(1,3,3) aggregated
t0=solve(Eq(g,rho),t); t0=[v for v in t0 if v>0][0]
chk("r=(1,3,3), rho=19/10: t0=5/13", t0==Rational(5,13))
x=[(t0+rho)/(3*(t0+r)) for r in (1,3,3)]
chk("x0 = (11/20, 9/40, 9/40)", x==[Rational(11,20),Rational(9,40),Rational(9,40)])
W,S,L,g,H,A,D,N,E=theory([Integer(1),Integer(2),Integer(3)],[Integer(1)]*3)
pol=Poly(E.subs(Z,rho),T).primitive()[1]
chk("r=(1,2,3): 3t^2-8t-29=0", expand(pol.as_expr()-(3*T**2-8*T-29))==0 or expand(pol.as_expr()+(3*T**2-8*T-29))==0)
chk("t1=(4+sqrt103)/3", simplify((3*T**2-8*T-29).subs(T,(4+sqrt(103))/3))==0)

# ---- Appendix C table
for m,rho,poly,(a_,b_) in [(3,Rational(9,5),3*T**2+2*T-9,(1,2)),(4,Rational(9,4),2*T**3+5*T**2-15*T-33,(2,3)),
                           (5,Rational(8,3),5*T**4+30*T**3-15*T**2-312*T-392,(3,4))]:
    W,S,L,g,H,A,D,N,E=theory([Integer(i) for i in range(1,m+1)],[Integer(1)]*m)
    P=Poly(E.subs(Z,rho),T).primitive()[1]
    if P.LC()<0: P=-P
    chk(f"App C m={m}, rho={rho}: {poly}", expand(P.as_expr()-poly)==0 and P.is_irreducible)
    pos=[r for r in Poly(poly,T).real_roots() if r>0]
    chk(f"App C m={m}: unique positive root in ({a_},{b_})", len(pos)==1 and a_<pos[0]<b_)

print("="*70); print("AUDIT OF THE ATTACHED TWO-RESOURCE PAPER"); print("="*70)
for x_ in ok: print("  PASS ",x_)
for x_ in bad: print("  FAIL ",x_)
print("="*70); print(f"{len(ok)} passed, {len(bad)} failed")
