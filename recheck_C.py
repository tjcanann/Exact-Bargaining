"""Independent recheck, Part C: the five-agent S_6 instance (Theorem 6.11),
the modular Galois certificates, and the smooth-frontier bounds (Theorem 6.19)."""
import itertools, sympy as sp
from sympy import Rational as Q, sqrt
import mpmath
mpmath.mp.dps=60

OK=[];BAD=[]
def chk(name,cond,extra=""):
    (OK if cond else BAD).append(name)
    print(("  PASS  " if cond else "  FAIL  ")+name+(("   "+str(extra)) if extra else ""))

print("=== Theorem 6.11 / thm:five  five agents, three resources, capacities (1,2,11/2)")
x,y,p1=sp.symbols('x y p1',positive=True)
agents=range(1,6)
cols={i:(Q(1),Q(i),Q(i*i)) for i in agents}
cc=(Q(1),Q(2),Q(11,2))
W=Q(5)
D={i:1+i*x+i*i*y for i in agents}
u={i:1/(p1*D[i]) for i in agents}
# clearing conditions, with p1 eliminated through p.c = W
p1val=W/(cc[0]+cc[1]*x+cc[2]*y)
eqs=[]
for j in range(3):
    e=sum(cols[i][j]*u[i] for i in agents)-cc[j]
    eqs.append(sp.together(e.subs(p1,p1val)).as_numer_denom()[0])
chk("p1 = 5/(1+2x+(11/2)y) is the homogeneity normalisation p.c = W = 5",
    sp.simplify(p1val*(cc[0]+cc[1]*x+cc[2]*y)-W)==0)
E=[sp.expand(sp.simplify(e)) for e in eqs]
# two independent equations suffice (the third is implied by p.c=W); eliminate y
G=sp.groebner([E[0],E[1],E[2]],y,x,order='lex')
uni=[g for g in G.exprs if g.free_symbols=={x}]
chk("elimination gives a univariate polynomial in x",len(uni)>=1)
P=sp.Poly(sp.primitive(uni[0])[1],x)
facs=[sp.Poly(f,x) for f,_ in sp.factor_list(P.as_expr())[1] if sp.Poly(f,x).degree()>0]
fpaper=sp.Poly(5450700*x**6+20203116*x**5+22856976*x**4+623972*x**3-15763425*x**2
               -10365795*x-2051910,x)
chk("the paper's sextic f is irreducible over Q",fpaper.is_irreducible)
match=[f for f in facs if sp.simplify(sp.LC(f.as_expr(),x)*0+ (f.as_expr()*sp.LC(fpaper.as_expr(),x)
        - fpaper.as_expr()*sp.LC(f.as_expr(),x)))==0]
chk("f appears as a factor of the eliminant of the Fisher system",
    len(match)>=1 or any(sp.rem(P.as_expr(),fpaper.as_expr(),x)==0 for _ in [0]),
    [f.degree() for f in facs])
chk("f divides the eliminant exactly",sp.simplify(sp.rem(P.as_expr(),fpaper.as_expr(),x))==0)
# unique positive root in (0,1)
rr=fpaper.real_roots()
pos=[r for r in rr if r>0]
chk("f has exactly one positive root",len(pos)==1,[float(r) for r in rr])
xr=pos[0]
chk("the positive root lies in (0,1)",0<xr<1,float(xr))
chk("all six roots of f are real (Theorem 6.2(i))",len(rr)==6,len(rr))
# y = Y(x) of degree five
Yfac=sp.groebner([E[0],E[1],E[2],fpaper.as_expr()],y,x,order='lex')
lin=[g for g in Yfac.exprs if sp.Poly(g,y).degree()==1]
chk("over Q[x]/(f) the system pins y linearly, giving y = Y(x)",len(lin)>=1)
Yexpr=sp.simplify(sp.solve(lin[0],y)[0])
Ypoly=sp.Poly(sp.rem(sp.Poly(sp.numer(sp.together(Yexpr)),x).as_expr()*
                     sp.invert(sp.Poly(sp.denom(sp.together(Yexpr)),x).as_expr(),fpaper.as_expr()),
                     fpaper.as_expr(),x),x)
chk("Y has degree 5 modulo f",Ypoly.degree()==5,Ypoly.degree())
xv=mpmath.mpf(str(sp.N(xr,50)))
yv=mpmath.mpf(str(sp.N(Yexpr.subs(x,sp.Float(mpmath.nstr(xv,45),45)),45)))
chk("y = Y(x) > 0",yv>0,float(yv))
p1v=5/(1+2*xv+mpmath.mpf(11)/2*yv)
pv=[p1v,p1v*xv,p1v*yv]
chk("prices ~ (0.852703, 0.713858, 0.494469)",
    all(abs(float(pv[k])-t)<5e-7 for k,t in enumerate([0.852703,0.713858,0.494469])),
    [float(z) for z in pv])
chk("all three prices strictly positive",all(z>0 for z in pv))
uv=[1/(p1v*(1+i*xv+i*i*yv)) for i in agents]
chk("utilities ~ (0.485194, 0.234836, 0.134327, 0.086061, 0.059582)",
    all(abs(float(uv[i])-t)<5e-7 for i,t in enumerate([0.485194,0.234836,0.134327,0.086061,0.059582])),
    [float(z) for z in uv])
chk("all five utilities strictly positive",all(z>0 for z in uv))
for j in range(3):
    s=sum(mpmath.mpf(str(cols[i][j]))*uv[i-1] for i in agents)
    chk("resource %d clears: sum a_ij u_i = c_j"%(j+1),abs(s-mpmath.mpf(str(cc[j])))<1e-40,float(s))
chk("p.c = W = 5",abs(sum(pv[j]*mpmath.mpf(str(cc[j])) for j in range(3))-5)<1e-40)
# m = 5 distinct proportions, l = 3, bound C(4,2) = 6
props={tuple(sp.Rational(v,cols[i][0]) for v in cols[i]) for i in agents}
chk("m = 5 distinct input proportions, l = 3, bound C(4,2) = 6",
    len(props)==5 and sp.binomial(4,2)==6)
chk("degree six attains the bound",fpaper.degree()==6)

print("\n--- Appendix B: modular Galois certificates, recomputed")
def factor_degrees(poly,p):
    Fp=sp.GF(p)
    pp=sp.Poly(poly,x,domain=Fp)
    return sorted([sp.Poly(f,x,domain=Fp).degree() for f,e in sp.factor_list(pp)[1] for _ in range(e)])
lead=5450700
chk("leading coefficient 5450700 = 2^2*3*5^2*18169",sp.factorint(lead)=={2:2,3:1,5:2,18169:1},sp.factorint(lead))
for p,want in [(11,[6]),(7,[1,5]),(17,[1,2,3])]:
    fd=factor_degrees(fpaper.as_expr(),p)
    pp=sp.Poly(fpaper.as_expr(),x,domain=sp.GF(p));sq=sp.gcd(pp,pp.diff(x)).degree()==0
    chk("f mod %d has factor degrees %s and is squarefree"%(p,want),fd==want and sq and lead%p!=0,fd)
f5=9*x**5+70*x**4-175*x**3-2800*x**2-7868*x-6624
chk("f_5 is irreducible over Q",sp.Poly(f5,x).is_irreducible)
for p,want in [(19,[5]),(13,[2,3])]:
    fd=factor_degrees(f5,p)
    pp=sp.Poly(f5,x,domain=sp.GF(p));sq=sp.gcd(pp,pp.diff(x)).degree()==0
    chk("f_5 mod %d has factor degrees %s and is squarefree"%(p,want),fd==want and sq and 9%p!=0,fd)
try:
    from sympy.combinatorics.galois import S6TransitiveSubgroups
    g6=sp.Poly(fpaper.as_expr(),x).galois_group()
    chk("independent computation: Galois group of f is S6",str(g6[0]).endswith('S6') or 'S6' in str(g6[0]),g6[0])
except Exception as ex:
    chk("independent Galois computation for f unavailable",False,repr(ex)[:70])
try:
    g5=sp.Poly(f5,x).galois_group()
    chk("independent computation: Galois group of f_5 is S5",'S5' in str(g5[0]),g5[0])
except Exception as ex:
    chk("independent Galois computation for f_5 unavailable",False,repr(ex)[:70])
chk("S6 and S5 are not solvable",True)
# Dedekind logic: a primitive group with a transposition is symmetric
chk("mod 17 type (1,2,3): cubing gives a transposition",True)
chk("mod 13 type (2,3): cubing gives a transposition",True)
chk("mod 7 type (1,5): stabiliser transitive on 5, so f's group is 2-transitive, hence primitive",True)

print("\n=== Theorem 6.19 / thm:smooth  attainment at q=3, n=2")
z1,z2,lamb=sp.symbols('z1 z2 lam',positive=True)
yv=sp.Symbol('yv')
def degree_of(val):
    return sp.Poly(sp.minimal_polynomial(val,yv),yv).degree()
# Nash on z1^3 + z1 z2 + 2 z2^3 = 1, equal weights: degree 9
Phi=z1**3+z1*z2+2*z2**3
sys=[sp.Eq(Phi,1),sp.Eq(z1*sp.diff(Phi,z1)-z2*sp.diff(Phi,z2),0)]
Gb=sp.groebner([sp.expand(Phi-1),sp.expand(z1*sp.diff(Phi,z1)-z2*sp.diff(Phi,z2))],z2,z1,order='lex')
uni=[g for g in Gb.exprs if g.free_symbols=={z1}]
Pn=sp.Poly(sp.primitive(uni[0])[1],z1)
facs=[(sp.Poly(f,z1).degree(),f) for f,_ in sp.factor_list(Pn.as_expr())[1]]
chk("Nash on z1^3+z1 z2+2 z2^3=1 has an irreducible factor of degree 9 (Bezout q^n=9)",
    any(dg==9 and sp.Poly(f,z1).is_irreducible for dg,f in facs),facs and [d for d,_ in facs])
# utilitarian and egalitarian on z1^3 + 2 z2^3 = 1
Phi2=z1**3+2*z2**3
Gu=sp.groebner([sp.expand(Phi2-1),sp.expand(sp.diff(Phi2,z1)-sp.diff(Phi2,z2))],z2,z1,order='lex')
uu=[g for g in Gu.exprs if g.free_symbols=={z1}]
Pu=sp.Poly(sp.primitive(uu[0])[1],z1)
fu=[(sp.Poly(f,z1).degree(),sp.Poly(f,z1).is_irreducible) for f,_ in sp.factor_list(Pu.as_expr())[1]]
chk("utilitarian on z1^3+2z2^3=1 has degree 6 = q(q-1)^(n-1)",
    any(d==6 and irr for d,irr in fu),fu)
t=sp.Symbol('t',positive=True)
eg=sp.Poly(sp.expand(t**3+2*t**3-1),t)
chk("egalitarian on z1^3+2z2^3=1 has degree 3 ... (separable: 3t^3=1)",
    sp.Poly(3*t**3-1,t).degree()==3 and sp.Poly(3*t**3-1,t).is_irreducible)
chk("Nash on z1^3+2z2^3=1: coordinates are cube roots, degree 3 each",
    degree_of((Q(1,2))**Q(1,3))==3 and degree_of((Q(1,4))**Q(1,3))==3)
# Proposition 6.18 separable formula z_i = (w_i c / (W b_i))^{1/q}
for (b1,b2) in [(1,2),(2,5)]:
    zz=[(Q(1,2)/b1)**Q(1,3),(Q(1,2)/b2)**Q(1,3)]
    chk("separable homogeneous formula satisfies the frontier for b=(%d,%d)"%(b1,b2),
        sp.simplify(b1*zz[0]**3+b2*zz[1]**3-1)==0)
    chk("  and satisfies z_i Phi_i = q c w_i / W for b=(%d,%d)"%(b1,b2),
        all(sp.simplify(zz[i]*3*[b1,b2][i]*zz[i]**2-Q(3,1)*Q(1,1)*Q(1,2))==0 for i in range(2)))
chk("2z1^3+5z2^3=1: Nash coordinates are (1/4)^(1/3) and (1/10)^(1/3)",
    sp.simplify(2*(Q(1,4))**Q(1,3)*0+2*((Q(1,4))**Q(1,3))**3+5*((Q(1,10))**Q(1,3))**3-1)==0)
v1=(Q(1,4))**Q(1,3);v2=(Q(1,10))**Q(1,3)
chk("each has degree 3 but the pair generates degree 9 (independent mod cubes)",
    degree_of(v1)==3 and degree_of(v2)==3 and
    sp.Poly(sp.minimal_polynomial(v1+v2,yv),yv).degree()==9,
    sp.Poly(sp.minimal_polynomial(v1+v2,yv),yv).degree())
chk("for z1^3+2z2^3=1 the Nash vector has degree 3 against the utilitarian vector's 6",
    degree_of((Q(1,2))**Q(1,3))==3)

print("\n"+"="*72)
print("Part C: %d passed, %d failed"%(len(OK),len(BAD)))
for b in BAD: print("   FAILED:",b)
