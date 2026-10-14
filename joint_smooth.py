from sympy import *
X=symbols('X')
def joint_deg(coords):
    # degree of Q(coords) via a random-ish primitive element
    best=0
    for lam in (1,2,3,5):
        e=coords[0]+lam*coords[1]
        best=max(best,degree(minimal_polynomial(e,X),X))
    return best
print("Separable homogeneous cubic frontier b1 z1^3 + b2 z2^3 = 1, n=2, q=3\n")
for b in [(1,2),(2,5),(3,7)]:
    b1,b2=map(Integer,b)
    # Nash: z_i^3 = c/(n b_i)
    zN=[Pow(Rational(1,2*b1),Rational(1,3)),Pow(Rational(1,2*b2),Rational(1,3))]
    # utilitarian: z_i = K * b_i^(-1/2),  K^3 * sum b_j^(-1/2)^3 ... solve directly
    s=b1**Rational(-1,2); t=b2**Rational(-1,2)
    K=Pow(1/(b1*s**3+b2*t**3),Rational(1,3))
    zU=[K*s,K*t]
    chkU=simplify(b1*zU[0]**3+b2*zU[1]**3-1)
    chkUfoc=simplify(3*b1*zU[0]**2-3*b2*zU[1]**2)
    cN=[degree(minimal_polynomial(z,X),X) for z in zN]
    cU=[degree(minimal_polynomial(z,X),X) for z in zU]
    print(f"b={b}:  NASH coord degrees {cN}, JOINT {joint_deg(zN)}   |  "
          f"UTIL coord degrees {cU}, JOINT {joint_deg(zU)}   (checks {chkU},{chkUfoc})")
print("""
Reading: per-coordinate, Nash <= q = 3 < utilitarian 6. Jointly, the Nash vector can reach
q^n = 9 (when the classes of c/(n b_i) are independent mod cubes) and then EXCEEDS the
utilitarian vector's 6. So the draft's 'reversal' holds for individual payoffs only.""")
