"""Test: Nash solution on a k-face of a polyhedral problem with n agents has
degree <= C(n-1,k) = number of bounded chambers of the coordinate arrangement
restricted to the face (Varchenko 1995; Orlik-Terao 1995).
Built from genuine Leontief data: K resources, n agents, positive coefficients,
face V = {A u = c} of dimension k = n - K."""
from sympy import *
from math import comb
import random, cmath
t1,t2=symbols('t1 t2')

def build(n,K,seed):
    random.seed(seed)
    A=Matrix(K,n,lambda i,j: Integer(random.randint(1,7)))
    u0=[Rational(random.randint(1,5),random.randint(1,3)) for _ in range(n)]
    c=A*Matrix(u0)                      # positive chamber nonempty (contains u0)
    w=[Integer(1)]*n                    # equal weights
    # parametrize V by the last k=n-K coordinates as free variables
    k=n-K; free=symbols(f's0:{k}')
    Ab=A[:, :K]; Af=A[:, K:]
    if Ab.det()==0: return None
    ub=Ab.inv()*(c-Af*Matrix(free))
    L=[expand(ub[i]) for i in range(K)]+list(free)
    return A,c,w,L,free

def crit_count_k2(n,K,seed):
    b=build(n,K,seed)
    if b is None: return None
    A,c,w,L,free=b
    L=[e.subs({free[0]:t1,free[1]:t2}) for e in L]
    E=[]
    for var in (t1,t2):
        e=0
        for i in range(n):
            pr=1
            for l in range(n):
                if l!=i: pr*=L[l]
            e+=w[i]*diff(L[i],var)*pr
        E.append(expand(e))
    R=Poly(resultant(E[0],E[1],t2),t1)
    genuine=[]
    for f,_ in factor_list(R.as_expr(),t1)[1]:
        Pf=Poly(f,t1)
        if Pf.degree()<1: continue
        hits=0; realhits=0
        for r in Pf.nroots(n=30,maxsteps=200):
            p2=Poly(E[0].subs(t1,r),t2)
            if p2.degree()<1: continue
            for s in p2.nroots(n=30,maxsteps=200):
                if abs(complex(E[1].subs({t1:r,t2:s})))>1e-8: continue
                if min(abs(complex(l.subs({t1:r,t2:s}))) for l in L)<1e-8: continue
                hits+=1
                if abs(complex(r).imag)<1e-9 and abs(complex(s).imag)<1e-9: realhits+=1
                break
        if hits: genuine.append((Pf.degree(),hits,realhits))
    # count critical points in the positive chamber (the actual Nash candidate)
    return genuine

print("="*76)
print(f"{'n':>2} {'K':>2} {'k':>2} {'C(n-1,k)':>9} {'Bezout (n-1)^k':>15}   genuine factors (deg, #pts, #real)")
print("="*76)
for n,K in [(4,2),(5,3),(6,4)]:
    k=n-K
    for seed in (1,2,3):
        g=crit_count_k2(n,K,seed)
        if g is None: continue
        tot=sum(h for _,h,_ in g); real=sum(r for *_,r in g)
        print(f"{n:>2} {K:>2} {k:>2} {comb(n-1,k):>9} {(n-1)**k:>15}   {g}  total={tot}, real={real}")
