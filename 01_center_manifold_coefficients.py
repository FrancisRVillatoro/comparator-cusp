"""Symbolic verification of the center-manifold coefficients of the symmetric
cross-repressive comparator.

Model (symmetric structural limit):
    X' = a1 f(Y) - X + g,   Y' = a2 f(X) - Y + g,   f(z) = 1/(1+z^n).

At the symmetry-breaking (pitchfork) point:  a_c f'(q_c) = -1, i.e.
    n (q_c - g) q_c^{n-1} = 1 + q_c^n,   a_c = (1+q_c^n)^2 / (n q_c^{n-1}).

Claimed reduced dynamics on the center manifold:
    d' = mu d - b d^3 + h Delta a,
    mu = A_mu (abar - a_c),  A_mu = (n-1) / (2 q_c (1+q_c^n)),
    b   = (n^2 - 1) / (12 q_c^2),
    h   = 1 / (2 (1+q_c^n)).
"""
import sympy as sp

q, n, g = sp.symbols('q n g', positive=True)
z = sp.Symbol('z', positive=True)

f = 1/(1+z**n)
f1 = sp.diff(f, z)
f2 = sp.diff(f, z, 2)
f3 = sp.diff(f, z, 3)

sub = {z: q}
f0q = sp.simplify(f.subs(sub))
f1q = sp.simplify(f1.subs(sub))
f2q = sp.simplify(f2.subs(sub))
f3q = sp.simplify(f3.subs(sub))

# critical value of the mean strength, from a_c f'(q_c) = -1
a_c = sp.simplify(-1/f1q)
print("a_c  =", sp.simplify(a_c))
print("a_c should equal (1+q^n)^2/(n q^(n-1)):",
      sp.simplify(a_c - (1+q**n)**2/(n*q**(n-1))) == 0)

# ---- center manifold: u = u_rho * rho + u_2 * d^2 ----
u_rho = f0q/2                     # from  0 = rho f(q) + (a_c f'(q) - 1) u
u_2 = a_c*f2q/4                   # from  0 = -2u + (a_c f''/2) d^2
print("u_rho - 1/(2(1+q^n)) = ",
      sp.simplify(u_rho - 1/(2*(1+q**n))))
print("u_2 - ((n+1)q^n-(n-1))/(4q(1+q^n)) = ",
      sp.simplify(sp.together(u_2 - ((n+1)*q**n-(n-1))/(4*q*(1+q**n)))))

# ---- linear coefficient A_mu ----
A_mu = sp.simplify(1/a_c - a_c*f2q*u_rho)   # rho-coefficient of d
print("A_mu simplified:", sp.simplify(sp.factor(sp.simplify(A_mu))))
print("A_mu - (n-1)/(2q(1+q^n)) = ",
      sp.simplify(A_mu - (n-1)/(2*q*(1+q**n))))

# ---- cubic coefficient b ----
b_expr = sp.simplify(a_c*f3q/6 + a_c*f2q*u_2)
print("b - (n^2-1)/(12 q^2) = ",
      sp.simplify(sp.together(b_expr - (n**2-1)/(12*q**2))))

# ---- forcing coefficient h ----
h_expr = f0q/2
print("h - 1/(2(1+q^n)) = ", sp.simplify(h_expr - 1/(2*(1+q**n))))

# ---- cusp wedge constant ----
C_expr = sp.simplify(2*((n-1)/(2*q*(1+q**n)))**sp.Rational(3,2)
                     / (3*sp.sqrt(3*(n**2-1)/(12*q**2))*(1/(2*(1+q**n)))))
C_claim = 2*sp.sqrt(2)*(n-1)/(3*sp.sqrt((n+1)*q*(1+q**n)))
print("C_n - claimed = ", sp.simplify(sp.powsimp(sp.simplify(C_expr - C_claim),
                                                 force=True)))

# ---- monotonicity of Phi_n and characterisation of q_c ----
Phi = n*(z-g)*z**(n-1)/(1+z**n)
dPhi = sp.simplify(sp.diff(Phi, z))
num = sp.simplify(sp.numer(sp.together(dPhi)))
print("numerator of Phi' =", sp.factor(sp.expand(num)))
print("Phi(inf) =", sp.limit(Phi, z, sp.oo))
