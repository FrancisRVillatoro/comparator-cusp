import numpy as np
from scipy.optimize import brentq

def mk(m, g):
    Psi  = lambda z: m*(z-g)*z**(m-1)/(1+z**m)
    dPsi = lambda z: m*z**(m-2)*(g*z**m + m*z - g*(m-1))/(1+z**m)**2
    L    = lambda z: m - (m-1)*g/z
    return Psi, dPsi, L

def fold_partner(Y_of, val):
    return Y_of(val)

def setup(m1, m2, g1, g2):
    P1, dP1, L1 = mk(m1, g1)
    P2, dP2, L2 = mk(m2, g2)
    def phi(X):                      # Y on the fold curve
        target = 1.0/P1(X)
        lo, hi = g2+1e-14, max(1.0, g2+1.0)
        while P2(hi) < target: hi *= 2
        return brentq(lambda Y: P2(Y)-target, lo, hi, xtol=1e-15, rtol=8.9e-16)
    return P1, dP1, L1, P2, dP2, L2, phi

def cusp_new(m1, m2, g1, g2, Xs=None):
    """zeros of H(X) = 1+L1(X) - Psi1(X)*(1+L2(phi(X)))"""
    P1, dP1, L1, P2, dP2, L2, phi = setup(m1, m2, g1, g2)
    Xstar = brentq(lambda X: P1(X)-1.0/m2, g1+1e-14, 1e6, xtol=1e-15, rtol=8.9e-16)
    H = lambda X: 1+L1(X) - P1(X)*(1+L2(phi(X)))
    grid = Xstar*np.exp(np.linspace(1e-9, np.log(1e6/Xstar), 4000))
    v = np.array([H(x) for x in grid])
    roots = [brentq(H, grid[i], grid[i+1], xtol=1e-14, rtol=8.9e-16)
             for i in range(len(v)-1) if v[i]*v[i+1] < 0]
    return Xstar, roots, phi, P1, P2, L1, L2

# --- 1. reproduce Table 1 of the previous version -------------------------
print("check against previously computed cusps")
for (m1,m2,g1,g2,lab) in [(2,2,0.05,0.05,'sym'), (2,3,0.05,0.05,'p=3,q=2'),
                          (2,2,0.05,0.20,'asym leaks'), (4,1,0.02,0.02,'m1=4,m2=1')]:
    Xs, roots, phi, P1,P2,L1,L2 = cusp_new(m1,m2,g1,g2)
    for Xc in roots:
        Yc = phi(Xc)
        a1 = (Xc-g1)*(1+Yc**m2); a2 = (Yc-g2)*(1+Xc**m1)
        print("  %-12s X_c=%.8f Y_c=%.8f a1=%.6f a2=%.6f  (#roots=%d)"
              % (lab, Xc, Yc, a1, a2, len(roots)))

# --- 2. zero-leak closed form --------------------------------------------
print("\nzero-leak closed form  Psi1(Xc)=(1+m1)/(1+m2)")
for (n_, m_) in [(2,2),(2,3),(3,2),(4,1.5),(1,4),(0.6,3)]:
    m1, m2 = n_, m_
    if m1*m2 <= 1: continue
    Xs, roots, phi, P1,P2,L1,L2 = cusp_new(m1,m2,0.0,0.0)
    Xc_pred = ((m1+1)/(m1*m2-1))**(1.0/m1)
    Yc_pred = ((m2+1)/(m1*m2-1))**(1.0/m2)
    print("  m1=%-4s m2=%-4s  #roots=%d  X_c=%.10f (pred %.10f)  Y_c=%.10f (pred %.10f)"
          % (m1, m2, len(roots), roots[0], Xc_pred, phi(roots[0]), Yc_pred))
