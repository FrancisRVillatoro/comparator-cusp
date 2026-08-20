"""Section 7: analytic expansion of the fold about the symmetric cusp.

Verifies, in this order,

  (a) Proposition 7.1 for symbolic Hill exponent m:
          eta_2 = ((m+1)q^m + 1) / (2 q (1+q^m)),
          r_2   = (m+1)(1+q^m) / (2q)      = 3b/A_mu,
          t_3   = -(m^2-1)(1+q^m) / (3q^2) = -2b/h,
      the last two obtained from the exact singular set alone, with no
      center-manifold input, and then compared with Theorem 6.3;

  (b) Proposition 7.2, the closed form of c_1 for m = 2 and arbitrary leak,
      together with its numerical values at the three leaks quoted in the paper
      and the small-leak limit c_1 -> g/12;

  (c) Corollary 7.3, the zero-leak quadratic expansion
          |Delta a| / (C_2 rho^{3/2}) = 1 - (5/288) rho^2 + (29/3456) rho^3 + ...;

  (d) Remark 7.4, the closed form of c_1 for m = 3 and the zero-leak values for
      m = 2, ..., 6, which show that the cancellation of (c) is specific to m = 2.
"""
import sympy as sp
import mpmath as mp

q, d, m, w = sp.symbols('q d m w', positive=True)
e2, e4 = sp.symbols('e2 e4')

print('=' * 74)
print('(a) leading coefficients for symbolic Hill exponent')
print('=' * 74)

N = 5
P = q ** m
g_of_q = q - (1 + P) / (m * q ** (m - 1))          # criticality Psi(q) = 1
Xm = sp.series((q + w) ** m, w, 0, N).removeO()
Xm1 = sp.series((q + w) ** (m - 1), w, 0, N).removeO()


def psi_shift(sh):
    """Psi(q + sh) truncated as a series in d."""
    return sp.series(sp.expand(m * (q + sh - g_of_q) * Xm1.subs(w, sh)
                               / (1 + Xm.subs(w, sh))), d, 0, N).removeO()


G = sp.expand(sp.series(psi_shift(e2 * d ** 2 + d) * psi_shift(e2 * d ** 2 - d) - 1,
                        d, 0, N).removeO())
E2 = sp.simplify(sp.solve(G.coeff(d, 2), e2)[0])
print('  eta_2      =', sp.factor(E2))
print('  claimed    = ((m+1)q^m + 1)/(2q(1+q^m))')
print('  difference =', sp.simplify(E2 - ((m + 1) * P + 1) / (2 * q * (1 + P))))

Xs = q + E2 * d ** 2 + d
Ys = q + E2 * d ** 2 - d
a1 = sp.series(sp.expand((Xs - g_of_q) * (1 + Xm.subs(w, Ys - q))), d, 0, N).removeO()
a2 = sp.series(sp.expand((Ys - g_of_q) * (1 + Xm.subs(w, Xs - q))), d, 0, N).removeO()
ac = (q - g_of_q) * (1 + P)
rho = sp.expand(sp.series((a1 + a2) / 2 - ac, d, 0, N).removeO())
da = sp.expand(sp.series(a1 - a2, d, 0, N).removeO())

r2 = sp.simplify(sp.factor(rho.coeff(d, 2)))
t3 = sp.simplify(sp.factor(da.coeff(d, 3)))
print('\n  rho(0) =', sp.simplify(rho.coeff(d, 0)),
      '   coefficient of d in Delta a =', sp.simplify(da.coeff(d, 1)))
print('  r_2 =', r2)
print('  t_3 =', t3)

Amu = (m - 1) / (2 * q * (1 + P))
b = (m ** 2 - 1) / (12 * q ** 2)
h = 1 / (2 * (1 + P))
print('  r_2 - 3b/A_mu =', sp.simplify(r2 - 3 * b / Amu))
print('  t_3 + 2b/h    =', sp.simplify(t3 + 2 * b / h))

print()
print('=' * 74)
print('(b) Proposition 7.2: c_1 for m = 2 and arbitrary leak')
print('=' * 74)

gq = (q ** 2 - 1) / (2 * q)                        # g in terms of q_c, m = 2
sig = e2 * d ** 2 + e4 * d ** 4
X2 = q + sig + d
Y2 = q + sig - d
F = sp.expand(4 * (X2 - gq) * X2 * (Y2 - gq) * Y2 - (1 + X2 ** 2) * (1 + Y2 ** 2))
S2 = sp.solve(sp.expand(F).coeff(d, 2).subs(e4, 0), e2)[0]
S4 = sp.solve(sp.expand(F.subs(e2, S2)).coeff(d, 4), e4)[0]
X2 = q + S2 * d ** 2 + S4 * d ** 4 + d
Y2 = q + S2 * d ** 2 + S4 * d ** 4 - d
A1 = sp.expand((X2 - gq) * (1 + Y2 ** 2))
A2 = sp.expand((Y2 - gq) * (1 + X2 ** 2))
rho2 = sp.expand(sp.simplify((A1 + A2) / 2 - (1 + q ** 2) ** 2 / (2 * q)))
da2 = sp.expand(sp.simplify(A1 - A2))
R2, R4 = sp.factor(sp.simplify(rho2.coeff(d, 2))), sp.factor(sp.simplify(rho2.coeff(d, 4)))
T3, T5 = sp.factor(sp.simplify(da2.coeff(d, 3))), sp.factor(sp.simplify(da2.coeff(d, 5)))
print('  r_2 =', R2)
print('  r_4 =', R4)
print('  t_3 =', T3)
print('  t_5 =', T5, '=', sp.expand(T5))
C1 = sp.simplify(sp.factor(sp.simplify((T5 / T3 - sp.Rational(3, 2) * R4 / R2) / R2)))
print('  c_1 =', C1)
print('  c_1 - (q^2-1)(3q^2+1)/(12q(1+q^2)^3) =',
      sp.simplify(C1 - (q ** 2 - 1) * (3 * q ** 2 + 1) / (12 * q * (1 + q ** 2) ** 3)))

mp.mp.dps = 30
gref = mp.mpf('0.1') / (mp.mpf('0.069') * 80)
print('\n  numerical values quoted in the paper:')
for gv, printed in [(mp.mpf('0.001'), '8.32083e-5'),
                    (gref, '1.468650e-3'),
                    (mp.mpf('0.2'), '1.182662e-2')]:
    qc = gv + mp.sqrt(1 + gv ** 2)
    print('    g = %-16s c_1 = %-18s paper: %s'
          % (mp.nstr(gv, 9), mp.nstr(gv * (3 * qc ** 2 + 1) / (6 * (1 + qc ** 2) ** 3), 10), printed))
tiny = mp.mpf('1e-10')
qt = tiny + mp.sqrt(1 + tiny ** 2)
print('    small-leak limit c_1/g = %s   vs 1/12 = %s'
      % (mp.nstr(tiny * (3 * qt ** 2 + 1) / (6 * (1 + qt ** 2) ** 3) / tiny, 10),
         mp.nstr(mp.mpf(1) / 12, 10)))

print()
print('=' * 74)
print('(c) Corollary 7.3: zero-leak quadratic expansion')
print('=' * 74)

e = sp.symbols('e')
Xz = 1 + e
a1z = 4 * Xz ** 3 / (3 * Xz ** 2 - 1)
a2z = (1 + Xz ** 2) ** sp.Rational(3, 2) / sp.sqrt(3 * Xz ** 2 - 1)
C2 = 2 / (3 * sp.sqrt(3))
NN = 12
rhoz = sp.series((a1z + a2z) / 2 - 2, e, 0, NN).removeO()
daz = sp.series(a1z - a2z, e, 0, NN).removeO()
print('  rho    =', sp.nsimplify(sp.expand(sp.series(rhoz, e, 0, 7).removeO())))
print('  Delta a=', sp.nsimplify(sp.expand(sp.series(daz, e, 0, 8).removeO())))
u = sp.symbols('u', positive=True)
Ac = sp.symbols('A1:8')
eexp = sum(Ac[i] * u ** (i + 1) for i in range(7))
eq = sp.series(sp.expand(rhoz.subs(e, eexp) - u ** 2), u, 0, 9).removeO()
sol = {}
for k in range(2, 9):
    c = sp.expand(eq.subs(sol)).coeff(u, k)
    s_ = sp.solve(c, Ac[k - 2], dict=True)
    if s_:
        sol[Ac[k - 2]] = sp.simplify(s_[0][Ac[k - 2]])
    eq = sp.expand(eq.subs(sol))
ratio = sp.series(sp.expand(daz / (C2 * rhoz ** sp.Rational(3, 2))), e, 0, 8).removeO()
ratio = sp.simplify(sp.expand(ratio.subs(e, sp.simplify(eexp.subs(sol)))))
ratio = sp.nsimplify(sp.expand(sp.series(ratio, u, 0, 7).removeO()))
print('  |Delta a|/(C_2 rho^{3/2}) - 1 in powers of u = sqrt(rho):')
print('   ', sp.expand(-ratio - 1), '   [expect  -5/288 u^4 + 29/3456 u^6]')
print('  5/288  =', float(sp.Rational(5, 288)), '   29/3456 =', float(sp.Rational(29, 3456)))

print()
print('=' * 74)
print('(d) Remark 7.4: the cancellation is specific to m = 2')
print('=' * 74)


def c1_integer_m(mint, zero_leak=False):
    mm = sp.Integer(mint)
    if zero_leak:
        gg = sp.Integer(0)
        Psi = lambda z: mm * z ** mm / (1 + z ** mm)
    else:
        gg = q - (1 + q ** mm) / (mm * q ** (mm - 1))
        Psi = lambda z: mm * (z - gg) * z ** (mm - 1) / (1 + z ** mm)
    s_ = e2 * d ** 2 + e4 * d ** 4
    XX, YY = q + s_ + d, q + s_ - d
    GG = sp.expand(sp.numer(sp.together(Psi(XX) * Psi(YY) - 1)))
    E2_ = sp.solve(sp.expand(GG.subs(e4, 0)).coeff(d, 2), e2)[0]
    E4_ = sp.solve(sp.expand(GG.subs(e2, E2_)).coeff(d, 4), e4)[0]
    XX = q + E2_ * d ** 2 + E4_ * d ** 4 + d
    YY = q + E2_ * d ** 2 + E4_ * d ** 4 - d
    aa1 = sp.expand((XX - gg) * (1 + YY ** mm))
    aa2 = sp.expand((YY - gg) * (1 + XX ** mm))
    acc = (q - gg) * (1 + q ** mm)
    rr = sp.expand(sp.simplify((aa1 + aa2) / 2 - acc))
    dd = sp.expand(sp.simplify(aa1 - aa2))
    return sp.simplify((dd.coeff(d, 5) / dd.coeff(d, 3)
                        - sp.Rational(3, 2) * rr.coeff(d, 4) / rr.coeff(d, 2)) / rr.coeff(d, 2))


c1_m3 = sp.simplify(sp.factor(c1_integer_m(3)))
print('  m = 3, arbitrary leak:  c_1 =', c1_m3)
print('  claimed 3(2q^6-2q^3-1)/(16 q (1+q^3)^3), difference =',
      sp.simplify(c1_m3 - 3 * (2 * q ** 6 - 2 * q ** 3 - 1) / (16 * q * (1 + q ** 3) ** 3)))
print('\n  zero-leak values, q_c = (m-1)^(-1/m):')
for mm in (2, 3, 4, 5, 6):
    val = sp.simplify(c1_integer_m(mm, zero_leak=True)
                      .subs(q, sp.Rational(1, mm - 1) ** sp.Rational(1, mm)))
    print('    m = %d :  c_1 = %-26s = %s' % (mm, val, sp.N(val, 10)))
print('\n  only m = 2 gives c_1 = 0, which is Corollary 7.3.')
