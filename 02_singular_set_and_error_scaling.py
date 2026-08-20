"""Numerical verification for the short deterministic paper (fast version)."""
import mpmath as mp

mp.mp.dps = 30

n = 2
g = mp.mpf('0.1') / (mp.mpf('0.069') * 80)


def Phi(zv):
    return n * (zv - g) * zv ** (n - 1) / (1 + zv ** n)


def dPhi(zv):
    return n * zv ** (n - 2) * (g * zv ** n + n * zv - g * (n - 1)) \
        / (1 + zv ** n) ** 2


def Phi_inv(val, x0=None):
    """Newton inversion of the strictly increasing Phi on (g, inf)."""
    x = mp.mpf(1) if x0 is None else mp.mpf(x0)
    if x <= g:
        x = g + mp.mpf('0.5')
    for _ in range(200):
        fx = Phi(x) - val
        step = fx / dPhi(x)
        xn = x - step
        if xn <= g:
            xn = (x + g) / 2
        if abs(xn - x) < mp.mpf('1e-28') * max(1, abs(xn)):
            return xn
        x = xn
    return x


q_c = Phi_inv(mp.mpf(1))
q_closed = g + mp.sqrt(1 + g ** 2)
a_c = (q_c - g) * (1 + q_c ** n)
A_mu = (n - 1) / (2 * q_c * (1 + q_c ** n))
b = (n ** 2 - 1) / (12 * q_c ** 2)
h = 1 / (2 * (1 + q_c ** n))
C_n = 2 * mp.sqrt(2) * (n - 1) / (3 * mp.sqrt((n + 1) * q_c * (1 + q_c ** n)))
X_min = Phi_inv(mp.mpf(1) / n, x0='0.6')

print('=' * 70)
print('1.  Cusp point and reduced coefficients')
print('=' * 70)
print('g              =', mp.nstr(g, 12))
print('q_c   (Phi=1)  =', mp.nstr(q_c, 12))
print('q_c   (closed) =', mp.nstr(q_closed, 12))
print('|difference|   =', mp.nstr(abs(q_c - q_closed), 4))
print('a_c            =', mp.nstr(a_c, 12))
print('A_mu =', mp.nstr(A_mu, 10), '  b =', mp.nstr(b, 10),
      '  h =', mp.nstr(h, 10), '  C_n =', mp.nstr(C_n, 10))
print('X_min (Phi=1/n, asymptote) =', mp.nstr(X_min, 12))


def fold_partner(Xv):
    return Phi_inv(1 / Phi(Xv), x0=q_c ** 2 / Xv)


def params(Xv, Yv):
    return (Xv - g) * (1 + Yv ** n), (Yv - g) * (1 + Xv ** n)


def detJ(Xv, Yv, a1, a2):
    fpY = -n * Yv ** (n - 1) / (1 + Yv ** n) ** 2
    fpX = -n * Xv ** (n - 1) / (1 + Xv ** n) ** 2
    return 1 - a1 * a2 * fpX * fpY


print()
print('=' * 70)
print('2.  det J on the curve Phi(X)Phi(Y)=1 (independent check)')
print('=' * 70)
for Xv in ['1.2', '1.8', '3.0', '8.0']:
    Xv = mp.mpf(Xv)
    Yv = fold_partner(Xv)
    a1, a2 = params(Xv, Yv)
    print('X=%-6s Y=%-14s a1=%-12s a2=%-12s detJ=%s'
          % (mp.nstr(Xv, 4), mp.nstr(Yv, 8), mp.nstr(a1, 7),
             mp.nstr(a2, 7), mp.nstr(detJ(Xv, Yv, a1, a2), 4)))

print()
print('=' * 70)
print('3.  Exact fold locus versus cubic normal form near the cusp')
print('=' * 70)
print('%-13s %-18s %-18s %-12s' % ('rho', '|Delta a| exact', 'C_n rho^{3/2}',
                                   'rel err'))
rows = []
for k in range(2, 15):
    Xv = q_c + mp.mpf(10) ** (-mp.mpf(k) / 2)
    Yv = fold_partner(Xv)
    a1, a2 = params(Xv, Yv)
    abar, da = (a1 + a2) / 2, abs(a1 - a2)
    rho = abar - a_c
    nf = C_n * rho ** mp.mpf('1.5')
    rel = abs(da - nf) / da
    rows.append((rho, rel))
    print('%-13s %-18s %-18s %-12s'
          % (mp.nstr(rho, 5), mp.nstr(da, 9), mp.nstr(nf, 9), mp.nstr(rel, 4)))

print()
print('local slope d log(rel err) / d log(rho):')
for i in range(len(rows) - 1):
    (r0, e0), (r1, e1) = rows[i], rows[i + 1]
    print('   rho ~ %-10s slope = %s'
          % (mp.nstr(r1, 3), mp.nstr(mp.log(e1 / e0) / mp.log(r1 / r0), 6)))


print()
print('Sections 4 and 5 of an earlier draft (experimental slice, far field) are')
print('covered by 03_experimental_slice.py.')
