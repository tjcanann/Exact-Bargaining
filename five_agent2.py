"""Five agents, three resources, technology a_i = (1, i, i^2), equal weights.
Dual (price) formulation of weighted Nash / Eisenberg-Gale:
   min_{p>=0}  p.c - sum_i w_i log(p.a_i),   u_i = w_i/(p.a_i).
All three prices positive  <=>  solution in the relative interior of the 2-face."""
import numpy as np
from scipy.optimize import minimize
from sympy import *
from sympy.polys.numberfields.galoisgroups import galois_group
import itertools
x,y,X=symbols('x y X')
n=5
a=[(1,i,i*i) for i in range(1,n+1)]
w=[1]*n

def dual_solve(c):
    c=np.array(c,float)
    def f(p): 
        d=np.array([p@np.array(ai,float) for ai in a])
        return p@c-np.sum(np.log(d))
    def g(p):
        d=np.array([p@np.array(ai,float) for ai in a])
        return c-sum(np.array(ai,float)/di for ai,di in zip(a,d))
    best=None
    for p0 in ([1,0.1,0.01],[1,1,1],[0.5,0.2,0.05]):
        r=minimize(f,np.array(p0,float),jac=g,bounds=[(0,None)]*3,method='L-BFGS-B',
                   options={'ftol':1e-15,'gtol':1e-13,'maxiter':20000})
        if best is None or r.fun<best.fun: best=r
    return best.x

cands=[]
for c2,c3 in itertools.product([Rational(k,4) for k in range(6,17)],[Rational(k,2) for k in range(8,40)]):
    p=dual_solve([1,float(c2),float(c3)])
    if min(p)>1e-4: cands.append(((Integer(1),c2,c3),p))
print(f"capacity vectors (1,c2,c3) on the grid with all three prices positive: {len(cands)}")
if not cands: raise SystemExit
# pick the one with smallest denominators / most interior
cands.sort(key=lambda t:(-min(t[1]/max(t[1]))))
c,p=cands[0]
print(f"chosen capacities c = {c},  numeric prices p ~ {p}")
xs,ys=p[1]/p[0],p[2]/p[0]

# exact system in x=p2/p1, y=p3/p1 after eliminating p1 by ratios of constraints
D=[1+ai[1]*x+ai[2]*y for ai in a]
S=[sum(Integer(w[i])*a[i][j]/D[i] for i in range(n)) for j in range(3)]
e1=together(c[0]*S[1]-c[1]*S[0]); e2=together(c[0]*S[2]-c[2]*S[0])
N1=expand(fraction(e1)[0]); N2=expand(fraction(e2)[0])
R=Poly(resultant(N1,N2,y),x)
print(f"eliminant degree in x: {R.degree()}")
target=None
for fct,_ in factor_list(R.as_expr(),x)[1]:
    Pf=Poly(fct,x)
    if Pf.degree()<1: continue
    val=abs(float(Pf.eval(xs)))/max(1.0,float(max(abs(cc) for cc in Pf.all_coeffs())))
    print(f"   factor degree {Pf.degree():2d}: |f(x*)| scaled = {val:.3e}")
    if val<1e-8: target=Pf
fz=target.primitive()[1]
print(f"\nMINIMAL POLYNOMIAL of the price ratio x = p2/p1 (degree {fz.degree()}):")
print(f"   {fz.as_expr()}")
print(f"   irreducible over Q: {fz.is_irreducible}")
G,alt=galois_group(fz)
print(f"   Galois group: order {G.order()}, solvable: {G.is_solvable}, transitive degree {G.degree}")
open('five_agent_poly.txt','w').write(str(fz.as_expr())+"\n"+str(c)+"\n")
