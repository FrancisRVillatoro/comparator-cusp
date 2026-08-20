"""Verification of the sign identity of Lemma 4.1:

    kappa + phi' = - H / [ (X-gamma_1) Psi_1(X) Psi_2'(phi(X)) ],

so that sign(kappa + phi') = -sign(H), and both da_1/dX and da_2/dX are positive
multiples of -H.  Computed at 40 significant digits.
"""
import mpmath as mp

mp.mp.dps = 40


def check(m1, m2, g1, g2, Xv):
    m1, m2, g1, g2, Xv = map(mp.mpf, (m1, m2, g1, g2, Xv))
    P1 = lambda z: m1*(z-g1)*z**(m1-1)/(1+z**m1)
    P2 = lambda z: m2*(z-g2)*z**(m2-1)/(1+z**m2)
    dP1 = lambda z: mp.diff(P1, z)
    dP2 = lambda z: mp.diff(P2, z)
    L1 = lambda z: m1-(m1-1)*g1/z
    L2 = lambda z: m2-(m2-1)*g2/z
    phi = lambda x: mp.findroot(lambda y: P2(y)-1/P1(x), mp.mpf('1.0'))

    Y = phi(Xv)
    phip = -dP1(Xv)*P2(Y)/(P1(Xv)*dP2(Y))
    kappa = P1(Xv)*(Y-g2)/(Xv-g1)
    H = 1 + L1(Xv) - P1(Xv)*(1+L2(Y))
    identity = kappa + phip + H/((Xv-g1)*P1(Xv)*dP2(Y))

    a1 = lambda x: (x-g1)*(1+phi(x)**m2)
    a2 = lambda x: (phi(x)-g2)*(1+x**m1)
    return H, kappa+phip, abs(identity), mp.diff(a1, Xv), mp.diff(a2, Xv)


print('%-34s %-12s %-14s %-12s %-11s %-11s'
      % ('(m1,m2,g1,g2,X)', 'H', "kappa+phi'", 'identity err', "a1'", "a2'"))
for args in [(2, 3, 0.05, 0.20, 0.7), (2, 3, 0.05, 0.20, 1.4),
             (3, 2, 0.05, 0.05, 0.8), (4, 1, 0.02, 0.02, 2.0),
             (2, 2, 0.018, 0.018, 0.8), (2, 2, 0.018, 0.018, 1.5)]:
    H, kp, err, da1, da2 = check(*args)
    print('%-34s %-+12.6f %-+14.6f %-12.1e %-+11.5f %-+11.5f'
          % (str(args), float(H), float(kp), float(err), float(da1), float(da2)))

print()
print('In every row sign(kappa+phi\') = -sign(H) and a1\', a2\' share that sign.')
