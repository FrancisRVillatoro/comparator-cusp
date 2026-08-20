"""Numerical check of the uniqueness proof for the cusp.

Along F parameterise by u = Psi_1(X).  With L_i(z) = m_i - (m_i-1)gamma_i/z,
    H(X) = 1 + L_1(X) - Psi_1(X)[1 + L_2(phi(X))]
vanishes exactly at cusp points, and

    d/du log[ (1+L_1)/(u(1+L_2)) ] = (1/u)[ C - 1 ],
    C = A_1/(1+L_1) + A_2/(1+L_2),   A_i = theta_i (L_i-1)/(L_i-Psi_i).

At a zero of H one has L_i - Psi_i = (L_1 L_2 - 1)/(1+L_{3-i}), hence C = t1 u1 + t2 u2
with u_i = (L_i-1)(1+L_j)/[(1+L_i)(L_1L_2-1)] and theta_i < 1, and
    (1+L_1)(1+L_2)(L_1L_2-1) [1 - u1 - u2] = (L_1L_2-1)^2 + (L_1-L_2)^2 >= 0.
"""
import numpy as np
from scipy.optimize import brentq

def setup(m1, m2, g1, g2):
    P1 = lambda z: m1*(z-g1)*z**(m1-1)/(1+z**m1)
    P2 = lambda z: m2*(z-g2)*z**(m2-1)/(1+z**m2)
    L1 = lambda z: m1 - (m1-1)*g1/z
    L2 = lambda z: m2 - (m2-1)*g2/z
    def phi(X):
        t = 1.0/P1(X); lo, hi = g2+1e-14, max(1.0, g2+1.0)
        while P2(hi) < t: hi *= 2
        return brentq(lambda Y: P2(Y)-t, lo, hi, xtol=1e-15, rtol=8.9e-16)
    return P1,P2,L1,L2,phi

# --- algebraic identity  (1+a)(1+b)(ab-1)(1-u1-u2) = (ab-1)^2+(a-b)^2 ------
import sympy as sp
a,b = sp.symbols('a b', positive=True)
u1 = (a-1)*(1+b)/((1+a)*(a*b-1)); u2 = (b-1)*(1+a)/((1+b)*(a*b-1))
lhs = sp.simplify((1+a)*(1+b)*(a*b-1)*(1-u1-u2))
print("identity check, lhs - [(ab-1)^2+(a-b)^2] =",
      sp.simplify(sp.expand(lhs - ((a*b-1)**2+(a-b)**2))))

# --- C < 1 at every cusp, over a wide grid --------------------------------
worst = -1e9; arg=None; nroots={}
grid = [(m1,m2,g1,g2)
        for m1 in [0.4,0.7,1.0,1.5,2,3,5,8,12,20]
        for m2 in [0.4,0.7,1.0,1.5,2,3,5,8,12,20]
        for g1 in [0.0,0.02,0.2,1.0,5.0,20.0]
        for g2 in [0.0,0.02,0.2,1.0,5.0,20.0] if m1*m2 > 1.0001]
for (m1,m2,g1,g2) in grid:
    try:
        P1,P2,L1,L2,phi = setup(m1,m2,g1,g2)
        Xs = brentq(lambda X: P1(X)-1.0/m2, g1+1e-13, 1e9, xtol=1e-15, rtol=8.9e-16)
        H = lambda X: 1+L1(X) - P1(X)*(1+L2(phi(X)))
        gg = Xs*np.exp(np.linspace(1e-9, np.log(1e8/Xs), 900))
        v = np.array([H(x) for x in gg])
        idx = [i for i in range(len(v)-1) if v[i]*v[i+1] < 0]
        nroots[len(idx)] = nroots.get(len(idx),0)+1
        for i in idx:
            Xc = brentq(H, gg[i], gg[i+1], xtol=1e-14, rtol=8.9e-16); Yc = phi(Xc)
            l1,l2 = L1(Xc), L2(Yc); t1,t2 = g1/Xc, g2/Yc
            A1 = t1*(l1-1)/(l1-P1(Xc)); A2 = t2*(l2-1)/(l2-P2(Yc))
            C = A1/(1+l1)+A2/(1+l2)
            if C > worst: worst, arg = C, (m1,m2,g1,g2,Xc)
    except Exception:
        pass
print("number of cusp roots found, histogram:", nroots)
print("max C over all cusp points = %.8f   at (m1,m2,g1,g2,Xc)=%s" % (worst, arg))
