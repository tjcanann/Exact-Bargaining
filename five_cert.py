from sympy import *
from sympy.polys.numberfields.galoisgroups import galois_group
import mpmath
mpmath.mp.dps=60
x,y,X,T=symbols('x y X T')
f=Poly(5450700*x**6 + 20203116*x**5 + 22856976*x**4 + 623972*x**3 - 15763425*x**2 - 10365795*x - 2051910, x)
print("f(x) =",f.as_expr())
print("content:",f.content(),"  irreducible over Q:",f.is_irreducible)

# --- modular certificate: find primes giving cycle types (6), (5,1), (2,1,1,1,1) or (2,3,1)
def pattern(p):
    if f.LC()%p==0: return None
    fp=Poly(f.as_expr(),x,modulus=p)
    if gcd(fp,fp.diff(x)).degree()>0: return None     # need squarefree mod p
    fl=fp.factor_list()[1]
    return tuple(sorted([Poly(g,x,modulus=p).degree() for g,e in fl for _ in range(e)],reverse=True))
need={(6,):None,(5,1):None,(2,1,1,1,1):None,(3,2,1):None}
for p in primerange(5,400):
    pt=pattern(p)
    if pt in need and need[pt] is None: need[pt]=p
print("\nfactorization patterns mod p (first prime found):")
for k,v in need.items(): print(f"   cycle type {k}: p = {v}")
print("""   (6): f irreducible mod p  -> irreducible over Q, Galois group transitive
   (5,1): a 5-cycle fixing a point; with transitivity the group is 2-transitive, hence primitive
   (2,1,1,1,1) is a transposition; (3,2,1) cubed is a transposition
   primitive + transposition  =>  full symmetric group S6 (Jordan)""")
G,_=galois_group(f)
print(f"sympy galois_group: order {G.order()}, solvable {G.is_solvable}")

# --- exact-feasibility check at high precision
a=[(1,i,i*i) for i in range(1,6)]; c=(1,2,Rational(11,2)); W=5
root=[r for r in Poly(f.as_expr(),x).nroots(n=50) if abs(im(r))<1e-40 and re(r)>0]
print(f"\npositive real roots of f: {[N(r,20) for r in root]}")
xs=re(root[0])
D=lambda xx,yy:[1+ai[1]*xx+ai[2]*yy for ai in a]
# y from the first ratio equation at x = xs (unique positive root)
S=lambda xx,yy,j: sum(Integer(a[i][j])/D(xx,yy)[i] for i in range(5))
eq=lambda yy: c[0]*S(xs,yy,1)-c[1]*S(xs,yy,0)
ys=mpmath.findroot(lambda yy: eq(yy), 0.58)
p1=W/(c[0]+c[1]*xs+c[2]*ys)
p=[p1,p1*xs,p1*ys]
u=[1/(p[0]*ai[0]+p[1]*ai[1]+p[2]*ai[2]) for ai in a]
res=[sum(a[i][j]*u[i] for i in range(5))-c[j] for j in range(3)]
print(f"prices   p = {[mpmath.nstr(v,15) for v in p]}")
print(f"allocation u = {[mpmath.nstr(v,15) for v in u]}")
print(f"resource residuals A u - c = {[mpmath.nstr(v,3) for v in res]}")
print(f"all prices > 0: {all(v>0 for v in p)};  all u > 0: {all(v>0 for v in u)}")
print("=> KKT with strictly positive prices; concave program => unique global Nash solution.")
