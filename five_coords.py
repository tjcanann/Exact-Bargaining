from sympy import *
x,U=symbols('x U')
f=Poly(5450700*x**6+20203116*x**5+22856976*x**4+623972*x**3-15763425*x**2-10365795*x-2051910,x)
Y=sympify(open('five_Y.txt').read())
den1=1+2*x+Rational(11,2)*Y                      # p1 = 5/den1
coords={}
for i in range(1,6):
    Di=1+i*x+i*i*Y
    coords[f"u{i}"]=(den1, 5*Di)                  # u_i = den1/(5 D_i)
coords["p1"]=(Integer(5),den1)
coords["p2"]=(5*x,den1)
coords["p3"]=(5*Y,den1)
for name,(num,den) in coords.items():
    num=Poly(num,x).rem(f).as_expr(); den=Poly(den,x).rem(f).as_expr()
    R=Poly(resultant(f.as_expr(),expand(U*den-num),x),U)
    fl=[(Poly(g,U).degree(),e) for g,e in factor_list(R.as_expr(),U)[1]]
    print(f"{name}: minimal polynomial degree(s) {fl}")
