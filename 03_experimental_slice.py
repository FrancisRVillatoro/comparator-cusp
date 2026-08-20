"""Experimental slice, equilibrium counts and far field (fast, guarded)."""
import mpmath as mp

mp.mp.dps = 30
n = 2
g = mp.mpf('0.1') / (mp.mpf('0.069') * 80)
A, m, kap = mp.mpf(3), 2, mp.mpf('1.5')


def Phi(z):
    return n * (z - g) * z ** (n - 1) / (1 + z ** n)


def dPhi(z):
    return n * z ** (n - 2) * (g * z ** n + n * z - g * (n - 1)) / (1 + z ** n) ** 2


def Phi_inv(val, x0=1):
    x = mp.mpf(x0)
    for _ in range(200):
        step = (Phi(x) - val) / dPhi(x)
        xn = x - step
        if xn <= g:
            xn = (x + g) / 2
        if abs(xn - x) < mp.mpf('1e-27') * max(1, abs(xn)):
            return xn
        x = xn
    return x


q_c = Phi_inv(1)
a_c = (q_c - g) * (1 + q_c ** n)
C_n = 2 * mp.sqrt(2) * (n - 1) / (3 * mp.sqrt((n + 1) * q_c * (1 + q_c ** n)))
X_min = Phi_inv(mp.mpf(1) / n, x0='0.6')


def fold_partner(X):
    return Phi_inv(1 / Phi(X), x0=q_c ** 2 / X)


def params(X, Y):
    return (X - g) * (1 + Y ** n), (Y - g) * (1 + X ** n)


def sr_of_fold(X):
    Y = fold_partner(X)
    a1, a2 = params(X, Y)
    if a1 >= A:
        return None
    sF = (a1 / (A - a1)) ** (mp.mpf(1) / m)
    rF = (a2 / A) * (kap ** m + sF ** m) / sF ** m
    return sF, rF


def bisect(fun, lo, hi, it=80):
    flo = fun(lo)
    for _ in range(it):
        mid = (lo + hi) / 2
        fm = fun(mid)
        if fm == 0:
            return mid
        if (flo < 0) != (fm < 0):
            hi = mid
        else:
            lo, flo = mid, fm
    return (lo + hi) / 2


# scan the valid part of the fold curve
grid = [g + mp.mpf(i) / 1000 * 3 for i in range(1, 1200)]
scan = []
for X in grid:
    v = sr_of_fold(X)
    if v is not None:
        scan.append((X, v[0], v[1]))

print('valid fold-curve range: X in [%s, %s], r_F in [%s, %s]'
      % (mp.nstr(scan[0][0], 6), mp.nstr(scan[-1][0], 6),
         mp.nstr(min(t[2] for t in scan), 8), mp.nstr(max(t[2] for t in scan), 8)))
print('s_* =', mp.nstr((a_c / (A - a_c)) ** (mp.mpf(1) / m), 12),
      ' r_* =', mp.nstr(a_c / A + (1 - a_c / A) * kap ** m, 12))


def exact_folds_at_r(rt):
    out = []
    for i in range(len(scan) - 1):
        if (scan[i][2] - rt) * (scan[i + 1][2] - rt) < 0:
            X = bisect(lambda t: sr_of_fold(t)[1] - rt, scan[i][0], scan[i + 1][0])
            out.append(sr_of_fold(X)[0])
    return sorted(out)


def nf_folds_at_r(rt):
    def F(s):
        a1 = A * s ** m / (1 + s ** m)
        a2 = rt * A * s ** m / (kap ** m + s ** m)
        ab = (a1 + a2) / 2
        if ab <= a_c:
            return mp.mpf(1)
        return abs(a1 - a2) - C_n * (ab - a_c) ** mp.mpf('1.5')
    gr = [mp.mpf(1) + mp.mpf(i) / 200 for i in range(0, 1200)]
    vs = [F(s) for s in gr]
    return sorted(bisect(F, gr[i], gr[i + 1], 70)
                  for i in range(len(gr) - 1) if vs[i] * vs[i + 1] < 0)


print()
print('%-6s %-13s %-13s %-13s %-13s %-9s %-9s'
      % ('r', 'exact s-', 'exact s+', 'NF s-', 'NF s+', 'err s-', 'err s+'))
for rt in ['1.39', '1.35', '1.30', '1.25', '1.20']:
    rt = mp.mpf(rt)
    ex, nf = exact_folds_at_r(rt), nf_folds_at_r(rt)
    if len(ex) == 2 and len(nf) == 2:
        e0 = abs(nf[0] - ex[0]) / ex[0] * 100
        e1 = abs(nf[1] - ex[1]) / ex[1] * 100
        print('%-6s %-13s %-13s %-13s %-13s %-9s %-9s'
              % (mp.nstr(rt, 4), mp.nstr(ex[0], 7), mp.nstr(ex[1], 7),
                 mp.nstr(nf[0], 7), mp.nstr(nf[1], 7),
                 mp.nstr(e0, 3) + '%', mp.nstr(e1, 3) + '%'))
    else:
        print(mp.nstr(rt, 4), [mp.nstr(v, 7) for v in ex],
              [mp.nstr(v, 7) for v in nf])


# ---- equilibrium counts across a fold ----
def equilibria(a1, a2):
    def G(X):
        Y = a2 / (1 + X ** n) + g
        return a1 / (1 + Y ** n) + g - X
    gr = [g + mp.mpf(i) / 300 * 6 for i in range(1, 301)]
    vs = [G(x) for x in gr]
    return [bisect(G, gr[i], gr[i + 1], 60)
            for i in range(len(gr) - 1) if vs[i] * vs[i + 1] < 0]


print()
X0 = mp.mpf('1.8')
Y0 = fold_partner(X0)
a1f, a2f = params(X0, Y0)
print('fold point: X=%s Y=%s a1=%s a2=%s'
      % (mp.nstr(X0, 6), mp.nstr(Y0, 8), mp.nstr(a1f, 8), mp.nstr(a2f, 8)))
for dl in ['-0.02', '-0.001', '0.001', '0.02']:
    dl = mp.mpf(dl)
    print('  a1 -> a1(1%+.3f):  %d equilibria'
          % (float(dl), len(equilibria(a1f * (1 + dl), a2f))))

# ---- far field ----
print()
print('%-8s %-14s %-14s %-16s' % ('X', 'Y', 'a1', 'a2/a1^n'))
for X in [10, 30, 100, 300, 1000]:
    X = mp.mpf(X)
    Y = fold_partner(X)
    a1, a2 = params(X, Y)
    print('%-8s %-14s %-14s %-16s'
          % (mp.nstr(X, 5), mp.nstr(Y, 9), mp.nstr(a1, 7), mp.nstr(a2 / a1 ** n, 9)))
print('predicted limit (X_min-g)/(1+X_min^n)^n =',
      mp.nstr((X_min - g) / (1 + X_min ** n) ** n, 9))
