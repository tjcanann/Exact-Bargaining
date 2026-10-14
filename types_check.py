"""Six agents, three resources, but agent 6's input column is proportional to agent 1's
(same 'type'), so m = 5 distinct types. Prediction: degree C(m-1,K-1)=C(4,2)=6, not C(5,2)=10."""
import numpy as np
from scipy.optimize import minimize
from sympy import *
x,y=symbols('x y')
cols=[(1,i,i*i) for i in range(1,6)]+[(2,2,2)]      # agent 6 = 2 x agent 1
w=[1]*6; n=6
def dual(c):
    c=np.array(c,float)
    f=lambda p: p@c-sum(np.log(p@np.array(a,float)) for a in cols)
    gr=lambda p: c-sum(np.array(a,float)/(p@np.array(a,float)) for a in cols)
    best=None
    for p0 in ([1,.5,.2],[1,1,1],[.6,.3,.1]):
        r=minimize(f,np.array(p0,float),jac=gr,bounds=[(0,None)]*3,method='L-BFGS-B',
                   options={'ftol':1e-15,'gtol':1e-13,'maxiter':20000})
        if best is None or r.fun<best.fun: best=r
    return best.x
import itertools
hit=None
for c2,c3 in itertools.product([Rational(k,4) for k in range(6,20)],[Rational(k,2) for k in range(6,40)]):
    p=dual([2,float(c2),float(c3)])
    if min(p)>1e-4: hit=((Integer(2),c2,c3),p); break
c,p=hit; xs=p[1]/p[0]
print("capacities",c," prices",p)
D=[a[0]+a[1]*x+a[2]*y for a in cols]
S=[sum(Integer(a[j])/Di for a,Di in zip(cols,D)) for j in range(3)]
N1=expand(fraction(together(c[0]*S[1]-c[1]*S[0]))[0])
N2=expand(fraction(together(c[0]*S[2]-c[2]*S[0]))[0])
R=Poly(resultant(N1,N2,y),x)
for fct,_ in factor_list(R.as_expr(),x)[1]:
    Pf=Poly(fct,x)
    if Pf.degree()<1: continue
    v=abs(float(Pf.eval(xs)))/max(1.0,float(max(abs(cc) for cc in Pf.all_coeffs())))
    if v<1e-7: print(f"factor containing the Nash solution: degree {Pf.degree()}  (C(m-1,2) with m=5 is 6; C(n-1,2) with n=6 is 10)")
