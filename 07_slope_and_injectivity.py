"""Slope of the discriminant along F:  S = da2/da1 = R m1 X^{m1-1}(phi-g2R)/(1+phi^{m2}).

Claim:   (X-g1)(L2-Psi2) u  dlog S/dX  =  u^2 + D u + 1,   D = (L1-1)(L2-1) - 2,
so S is strictly increasing whenever (m1-1)(m2-1) >= 0.
"""
import mpmath as mp
mp.mp.dps = 35

def setup(m1,m2,g1,g2,R=1):
    m1,m2,g1,g2,R = map(mp.mpf,(m1,m2,g1,g2,R)); h2 = g2/R
    P1 = lambda z: m1*(z-g1)*z**(m1-1)/(1+z**m1)
    P2 = lambda z: m2*(z-h2)*z**(m2-1)/(1+z**m2)
    L1 = lambda z: m1-(m1-1)*g1/z
    L2 = lambda z: m2-(m2-1)*h2/z
    phi = lambda X: mp.findroot(lambda y: P2(y)-1/P1(X), mp.mpf('1.0'))
    a1 = lambda X: (X-g1)*(1+phi(X)**m2)
    a2 = lambda X: R*(phi(X)-h2)*(1+X**m1)
    S  = lambda X: R*m1*X**(m1-1)*(phi(X)-h2)/(1+phi(X)**m2)
    return m1,m2,g1,h2,R,P1,P2,L1,L2,phi,a1,a2,S

for (m1,m2,g1,g2,R,X) in [(2,3,0.05,0.20,1,0.7),(2,3,0.05,0.20,1,1.4),
                          (4,1,0.02,0.02,1,2.0),(3,2,0.05,0.05,1,0.8),
                          (2,2,0.1,0.3,2.5,1.3),(0.5,4,0.1,0.1,1,3.0)]:
    m1_,m2_,g1_,h2_,R_,P1,P2,L1,L2,phi,a1,a2,S = setup(m1,m2,g1,g2,R)
    Xv = mp.mpf(X); Y = phi(Xv); u = P1(Xv)
    S_direct = mp.diff(a2,Xv)/mp.diff(a1,Xv)
    D = (L1(Xv)-1)*(L2(Y)-1) - 2
    lhs = mp.diff(lambda t: mp.log(S(t)), Xv)*(Xv-g1_)*(L2(Y)-P2(Y))*u
    rhs = u**2 + D*u + 1
    print("m=(%s,%s) g=(%s,%s) R=%s X=%s | S formula err=%.1e | identity err=%.1e | u=%.4f D=%+.4f u^2+Du+1=%+.5f"
          % (m1,m2,g1,g2,R,X, float(abs(S(Xv)-S_direct)), float(abs(lhs-rhs)),
             float(u), float(D), float(rhs)))

print()
print("global monotonicity of S and separation of the two branches")
import numpy as np
from scipy.optimize import brentq
def scan(m1,m2,g1,g2,R=1.0):
    h2=g2/R
    P1=lambda z: m1*(z-g1)*z**(m1-1)/(1+z**m1)
    P2=lambda z: m2*(z-h2)*z**(m2-1)/(1+z**m2)
    def phi(X):
        t=1.0/P1(X); lo,hi=h2+1e-14,max(1.0,h2+1.0)
        for _ in range(200):
            if P2(hi)>=t: break
            hi*=1.6
        return brentq(lambda y:P2(y)-t, lo,hi, xtol=1e-15, rtol=8.9e-16)
    Xs=brentq(lambda X:P1(X)-1.0/m2, g1+1e-13,1e9, xtol=1e-15, rtol=8.9e-16)
    gr=np.exp(np.linspace(np.log(Xs*(1+1e-7)), np.log(Xs*1e6), 4000))
    S=np.array([R*m1*X**(m1-1)*(phi(X)-h2)/(1+phi(X)**m2) for X in gr])
    a1=np.array([(X-g1)*(1+phi(X)**m2) for X in gr])
    a2=np.array([R*(phi(X)-h2)*(1+X**m1) for X in gr])
    ic=int(np.argmin(a1))
    bL=(a1[:ic][::-1],a2[:ic][::-1]); bR=(a1[ic+1:],a2[ic+1:])
    lo=max(bL[0].min(),bR[0].min()); hi=min(bL[0].max(),bR[0].max())
    xs=np.exp(np.linspace(np.log(lo*1.0001),np.log(hi*0.9999),600))
    d=np.interp(xs,bL[0],bL[1])-np.interp(xs,bR[0],bR[1])
    return bool(np.all(np.diff(S)>0)), bool(np.all(d<0) or np.all(d>0))
for c in [(2,2,0.018,0.018,1),(2,3,0.05,0.20,1),(3,2,0.05,0.05,1),(4,1,0.02,0.02,1),
          (2,2,0.1,0.3,2.5),(6,0.5,0.5,0.5,1),(0.5,4,0.1,0.1,1),(1,8,2.0,0.0,1)]:
    mono,sep = scan(*c)
    print("  m=(%-4s,%-4s) g=(%-4s,%-4s) R=%-3s : (m1-1)(m2-1)>=0 %-5s | S increasing %-5s | branches separated %s"
          % (c[0],c[1],c[2],c[3],c[4], (c[0]-1)*(c[1]-1)>=0, mono, sep))
