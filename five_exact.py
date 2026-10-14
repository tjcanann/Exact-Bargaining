from sympy import *
import mpmath
from mpmath import iv
x,y=symbols('x y')
n=5; a=[(1,i,i*i) for i in range(1,6)]; c=(Integer(1),Integer(2),Rational(11,2)); W=5
f=Poly(5450700*x**6+20203116*x**5+22856976*x**4+623972*x**3-15763425*x**2-10365795*x-2051910,x)
D=[1+ai[1]*x+ai[2]*y for ai in a]
S=[sum(Integer(ai[j])/Di for ai,Di in zip(a,D)) for j in range(3)]
N1=expand(fraction(together(c[0]*S[1]-c[1]*S[0]))[0])
N2=expand(fraction(together(c[0]*S[2]-c[2]*S[0]))[0])
# work in Q[x]/(f): find y as polynomial in x from the subresultant chain in y
K=QQ.algebraic_field  # not used; do it by hand with polynomial remainders mod f
def red(p): return Poly(p,x).rem(f).as_expr()
P1=Poly(N1,y); P2=Poly(N2,y)
def reduce_poly(P): return Poly([red(cf) for cf in P.all_coeffs()],y)
# Euclid in (Q[x]/f)[y]
def inv_mod_f(e):
    s,t,g=gcdex(Poly(e,x),f)          # s*e + t*f = g
    return (s.as_expr()/g.as_expr())
A_,B_=reduce_poly(P1),reduce_poly(P2)
while B_.degree()>0 and any(simplify(cf)!=0 for cf in B_.all_coeffs()):
    # remainder of A_ by B_ over Q[x]/f
    Q,Rm=A_,Poly(0,y)
    r=A_
    lb=B_.LC(); ilb=inv_mod_f(lb)
    while r.degree()>=B_.degree() and not r.is_zero:
        coeff=red(r.LC()*ilb)
        r=reduce_poly(r-Poly(coeff*y**(r.degree()-B_.degree()),y)*B_)
        if all(simplify(cf)==0 for cf in r.all_coeffs()): r=Poly(0,y); break
    A_,B_=B_,r
lin=A_
print("degree in y of gcd over Q[x]/(f):",lin.degree())
cA,cB=lin.all_coeffs()
Y=red(-cB*inv_mod_f(cA))
Y=Poly(Y,x)
print("y = Y(x) mod f, with Y of degree",Y.degree())
# exact check: N1(x,Y(x)) and N2(x,Y(x)) are 0 mod f
r1=Poly(expand(N1.subs(y,Y.as_expr())),x).rem(f); r2=Poly(expand(N2.subs(y,Y.as_expr())),x).rem(f)
print("N1(x,Y(x)) mod f == 0:",r1.is_zero,"   N2(x,Y(x)) mod f == 0:",r2.is_zero)
# isolate the positive root of f exactly
ints=[I for I in f.intervals() if I[0][1]>0 or I[0][0]>0]
print("positive root isolating intervals:",ints)
lo,hi=ints[0][0]
for _ in range(60):   # refine by bisection with exact sign tests
    mid=(lo+hi)/2
    if f.eval(mid)==0: lo=hi=mid; break
    if sign(f.eval(lo))*sign(f.eval(mid))<0: hi=mid
    else: lo=mid
print(f"refined isolating interval width ~ {float(hi-lo):.2e}")
iv.dps=40
X=iv.mpf([mpmath.mpf(lo.p)/lo.q, mpmath.mpf(hi.p)/hi.q])
def ev(poly,Xi):
    r=iv.mpf(0)
    for cf in Poly(poly,x).all_coeffs(): r=r*Xi+iv.mpf(Rational(cf).p)/Rational(cf).q
    return r
Yi=ev(Y.as_expr(),X)
Di=[1+ai[1]*X+ai[2]*Yi for ai in a]
den=1+2*X+iv.mpf(11)/2*Yi
p1=iv.mpf(W)/den
print("interval x  =",X)
print("interval y  =",Yi)
print("interval p1 =",p1)
print("interval D_i=",Di)
pos=all(v.a>0 for v in [X,Yi,p1]+Di)
print("CERTIFIED: x>0, y>0, p1>0, all D_i>0 ->  all three prices positive and u>0:",pos)
open('five_Y.txt','w').write(str(Y.as_expr()))
