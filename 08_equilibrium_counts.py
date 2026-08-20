import numpy as np
from scipy.optimize import brentq

def make(m1,m2,g1,g2,R=1.0):
    P1 = lambda z: m1*(z-g1)*z**(m1-1)/(1+z**m1)
    P2 = lambda z: m2*(z-g2/R)*z**(m2-1)/(1+z**m2)
    def phi(X):
        t = 1.0/P1(X); lo, hi = g2/R+1e-14, max(1.0, g2/R+1.0)
        while P2(hi) < t: hi *= 2
        return brentq(lambda Y: P2(Y)-t, lo, hi, xtol=1e-15, rtol=8.9e-16)
    A = lambda X,Y: ((X-g1)*(1+Y**m2), R*(Y-g2/R)*(1+X**m1))
    def eqs(a1,a2):
        F = lambda X: g1 + a1/(1+(g2/R + a2/(R*(1+X**m1)))**m2) - X
        gr = g1 + (a1+1e-9)*np.linspace(1e-13,1,60000)
        v = np.array([F(x) for x in gr])
        return [brentq(F, gr[i], gr[i+1], xtol=1e-15, rtol=8.9e-16)
                for i in range(len(v)-1) if v[i]*v[i+1] < 0]
    detJ = lambda X,Y: R*(1 - P1(X)*P2(Y))
    return P1,P2,phi,A,eqs,detJ

def cuspX(P1,phi,L1,L2,Xs):
    H = lambda X: 1+L1(X) - P1(X)*(1+L2(phi(X)))
    gr = Xs*np.exp(np.linspace(1e-9, np.log(1e7/Xs), 3000))
    v = [H(x) for x in gr]
    for i in range(len(v)-1):
        if v[i]*v[i+1] < 0:
            return brentq(H, gr[i], gr[i+1], xtol=1e-14, rtol=8.9e-16)

for (m1,m2,g1,g2,R) in [(2,2,0.018,0.018,1.0),(3,2,0.05,0.05,1.0),
                        (2,2,0.05,0.20,1.0),(4,1,0.02,0.02,1.0),
                        (2,3,0.0,0.0,1.0),(2,2,0.1,0.3,2.5),(0.5,4,0.1,0.1,1.0)]:
    P1,P2,phi,A,eqs,detJ = make(m1,m2,g1,g2,R)
    L1 = lambda z: m1-(m1-1)*g1/z
    L2 = lambda z: m2-(m2-1)*(g2/R)/z
    Xs = brentq(lambda X: P1(X)-1.0/m2, g1+1e-13, 1e9, xtol=1e-15, rtol=8.9e-16)
    Xc = cuspX(P1,phi,L1,L2,Xs)
    X0 = 2.5*Xc                       # well away from the cusp
    Y0 = phi(X0); a1f,a2f = A(X0,Y0)
    out=[]
    for eps in (-0.004, 0.004):
        E = eqs(a1f*(1+eps), a2f)
        nod = sum(1 for X in E if detJ(X, g2/R + a2f*(1+0)/(R*(1+X**m1))) > 0)
        out.append((len(E), nod, len(E)-nod))
    print("m=(%s,%s) g=(%s,%s) R=%s  Xc=%.5f | a1(1-0.4%%): %d eq (%d node,%d saddle) | a1(1+0.4%%): %d eq (%d node,%d saddle)"
          % (m1,m2,g1,g2,R,Xc,out[0][0],out[0][1],out[0][2],out[1][0],out[1][1],out[1][2]))
