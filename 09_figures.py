"""Corrected Figures 1 and 3: regions shaded and verified by equilibrium counting."""
import os
import numpy as np
from scipy.optimize import brentq
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

OUT = 'figures'
os.makedirs(OUT, exist_ok=True)
plt.rcParams.update({
    'font.size': 9, 'axes.labelsize': 9, 'axes.titlesize': 9,
    'legend.fontsize': 7.6, 'xtick.labelsize': 8, 'ytick.labelsize': 8,
    'axes.linewidth': 0.7, 'lines.linewidth': 1.3,
    'figure.dpi': 160, 'savefig.bbox': 'tight',
})
BLUE, RED, GRY, FILL = '#1f4e79', '#c0392b', '#8a8a8a', '#dce7f2'


def model(m1, m2, g1, g2, R=1.0):
    h2 = g2 / R
    P1 = lambda z: m1*(z-g1)*z**(m1-1)/(1+z**m1)
    P2 = lambda z: m2*(z-h2)*z**(m2-1)/(1+z**m2)
    L1 = lambda z: m1-(m1-1)*g1/z
    L2 = lambda z: m2-(m2-1)*h2/z

    def phi(X):
        t = 1.0/P1(X); lo, hi = h2+1e-14, max(1.0, h2+1.0)
        for _ in range(400):
            if P2(hi) >= t: break
            hi *= 1.6
        return brentq(lambda Y: P2(Y)-t, lo, hi, xtol=1e-15, rtol=8.9e-16)

    A = lambda X, Y: ((X-g1)*(1+Y**m2), R*(Y-h2)*(1+X**m1))
    Xstar = brentq(lambda X: P1(X)-1.0/m2, g1+1e-13, 1e9, xtol=1e-15, rtol=8.9e-16)
    Ystar = brentq(lambda Y: P2(Y)-1.0/m1, h2+1e-13, 1e9, xtol=1e-15, rtol=8.9e-16)
    H = lambda X: 1+L1(X) - P1(X)*(1+L2(phi(X)))
    gr = Xstar*np.exp(np.linspace(1e-9, np.log(1e7/Xstar), 3000))
    v = [H(x) for x in gr]
    Xc = next(brentq(H, gr[i], gr[i+1], xtol=1e-14, rtol=8.9e-16)
              for i in range(len(v)-1) if v[i]*v[i+1] < 0)
    ker = lambda X, Y: (-(X-g1)*P2(Y)/(Y-h2), 1.0)

    def neq(a1, a2):
        F = lambda X: g1 + a1/(1+(h2 + a2/(R*(1+X**m1)))**m2) - X
        gg = g1 + (a1+1e-9)*np.linspace(1e-13, 1, 30000)
        w = np.array([F(x) for x in gg])
        return sum(1 for i in range(len(w)-1) if w[i]*w[i+1] < 0)

    return dict(P1=P1, P2=P2, phi=phi, A=A, Xstar=Xstar, Ystar=Ystar, Xc=Xc,
                ker=ker, h2=h2, m1=m1, m2=m2, g1=g1, R=R, neq=neq)


# =================== Figure 1 =====================
M = model(2, 3, 0.05, 0.20)
Xc, Yc = M['Xc'], M['phi'](M['Xc'])
a1c, a2c = M['A'](Xc, Yc)
lo = M['Xstar']*(1+1e-9)
Xs = np.concatenate([np.exp(np.linspace(np.log(lo), np.log(Xc), 900))[:-1],
                     np.exp(np.linspace(np.log(Xc), np.log(Xc*400), 900))])
Ys = np.array([M['phi'](x) for x in Xs])
a1s, a2s = M['A'](Xs, Ys)

fig, ax = plt.subplots(1, 2, figsize=(7.3, 3.1))

# --- panel (a): log axes so both asymptotes are visible
ax[0].plot(Xs, Ys, color=BLUE, lw=1.5, label=r'$\mathcal{F}:\ \Psi_1(X)\Psi_2(Y)=1$')
ax[0].axhline(M['Ystar'], color=GRY, ls=':', lw=0.9)
ax[0].axvline(M['Xstar'], color=GRY, ls=':', lw=0.9)
for xk in [Xc*0.55, Xc*0.78, Xc, Xc*1.45, Xc*2.6]:
    yk = M['phi'](xk)
    vx, vy = M['ker'](xk, yk)
    n = np.hypot(vx, vy); sc = 0.26*yk
    ax[0].annotate('', xy=(xk+sc*vx/n, yk+sc*vy/n), xytext=(xk-sc*vx/n, yk-sc*vy/n),
                   arrowprops=dict(arrowstyle='-',
                                   color=RED if abs(xk-Xc) < 1e-9 else '#444444',
                                   lw=1.7 if abs(xk-Xc) < 1e-9 else 0.9))
ax[0].plot([Xc], [Yc], 'o', ms=5.5, color=RED, zorder=6, label='cusp point')
ax[0].plot([], [], '-', color='#444444', lw=0.9, label=r'$\ker D\mathcal{A}$')
ax[0].set_xscale('log'); ax[0].set_yscale('log')
ax[0].set_xlim(0.42, 30); ax[0].set_ylim(0.55, 12)
ax[0].annotate(r'$X=X_\ast$', xy=(M['Xstar'], 6.5), xytext=(0.95, 7.6),
               fontsize=8, color='#444444',
               arrowprops=dict(arrowstyle='->', lw=0.7, color='#444444'))
ax[0].annotate(r'$Y=Y_\ast$', xy=(12, M['Ystar']), xytext=(3.0, 1.00),
               fontsize=8, color='#444444',
               arrowprops=dict(arrowstyle='->', lw=0.7, color='#444444'))
ax[0].set_xlabel(r'$X$'); ax[0].set_ylabel(r'$Y$')
ax[0].set_title(r'(a) singular set of $\mathcal{A}$ in the state plane')
ax[0].legend(frameon=False, loc='upper right')

# --- panel (b): shade the three-equilibrium region between the two branches
ic = int(np.argmin(a1s))
bL = (a1s[:ic][::-1], a2s[:ic][::-1])          # X < Xc : flat branch
bR = (a1s[ic+1:], a2s[ic+1:])                  # X > Xc : steep branch
grid = np.linspace(a1c, 9, 800)[1:]
lowb = np.interp(grid, bL[0], bL[1])
uppb = np.interp(grid, bR[0], bR[1], right=np.inf)
uppb = np.minimum(np.nan_to_num(uppb, posinf=9.0), 9.0)
ax[1].fill_between(grid, lowb, uppb, where=uppb > lowb, color=FILL, lw=0, zorder=0)
sel = (a1s < 9.2) & (a2s < 9.2)
ax[1].plot(a1s[sel], a2s[sel], color=BLUE, lw=1.5, zorder=3,
           label=r'discriminant $\mathcal{A}(\mathcal{F})$')
ax[1].plot([a1c], [a2c], 'o', ms=5.5, color=RED, zorder=6, label='cusp')
ax[1].text(7.1, 5.2, 'three equilibria\n(two nodes,\none saddle)', fontsize=7.6,
           color='#1a3a55', ha='center')
ax[1].text(2.4, 7.4, 'one equilibrium\n(globally attracting)', fontsize=7.6,
           color='#333333', ha='center')
ax[1].set_xlim(0, 9); ax[1].set_ylim(0, 9)
ax[1].set_xlabel(r'$a_1$'); ax[1].set_ylabel(r'$a_2$')
ax[1].set_title(r'(b) image under $\mathcal{A}$: bifurcation diagram')
ax[1].legend(frameon=False, loc='lower right')
fig.savefig(os.path.join(OUT, 'singular_set.pdf'))
fig.savefig(os.path.join(OUT, 'singular_set.png'))
plt.close(fig)

print('Figure 1 checks (m1=2, m2=3, g1=0.05, g2R=0.20):')
print('  cusp (a1,a2) = (%.6f, %.6f)' % (a1c, a2c))
print('  label "three equilibria" is centred at (7.10, 5.20); measured box')
print('  a1 in [5.68, 8.52], a2 in [5.12, 6.38].  Corner checks:')
for p in [(7.10, 5.20), (5.68, 5.12), (5.68, 6.38), (8.52, 5.12), (8.52, 6.38)]:
    print('      (%.2f, %.2f): %d equilibria' % (p[0], p[1], M['neq'](*p)))
print('  label "one equilibrium" is centred at (2.40, 7.40); its box spans about')
print('  measured box a1 in [0.51, 4.29], a2 in [7.32, 8.14].  Corner checks:')
for p in [(2.40, 7.40), (0.51, 7.32), (0.51, 8.14), (4.29, 7.32), (4.29, 8.14)]:
    print('      (%.2f, %.2f): %d equilibria' % (p[0], p[1], M['neq'](*p)))

# =================== Figure 3 =====================
g = 0.1/(0.069*80)
S = model(2, 2, g, g)
qc = S['Xc']; ac = (qc-g)*(1+qc**2)
A, ms, kap = 3.0, 2, 1.5
s_star = (ac/(A-ac))**(1.0/ms)
r_star = ac/A + (1-ac/A)*kap**ms
sU, rU, sL, rL = [], [], [], []
for X in np.linspace(0.60, qc, 900)[:-1]:      # branch with Y > qc
    Y = S['phi'](X); a1, a2 = S['A'](X, Y)
    if a1 < A:
        s = (a1/(A-a1))**(1.0/ms)
        sL.append(s); rL.append((a2/A)*(kap**ms+s**ms)/s**ms)
for X in np.linspace(qc, 2.2, 900):            # branch with Y < qc
    Y = S['phi'](X); a1, a2 = S['A'](X, Y)
    if a1 < A:
        s = (a1/(A-a1))**(1.0/ms)
        sU.append(s); rU.append((a2/A)*(kap**ms+s**ms)/s**ms)
sL, rL, sU, rU = map(np.array, (sL, rL, sU, rU))
o = np.argsort(sL); sL, rL = sL[o], rL[o]
o = np.argsort(sU); sU, rU = sU[o], rU[o]
top, bot = (sU, rU) if rU.mean() > rL.mean() else (sL, rL), \
           (sL, rL) if rU.mean() > rL.mean() else (sU, rU)

fig, ax = plt.subplots(figsize=(4.2, 3.2))
gs = np.linspace(s_star, 2.6, 700)[1:]
hi = np.interp(gs, top[0], top[1])
lo2 = np.interp(gs, bot[0], bot[1], right=1.10)
ax.fill_between(gs, np.maximum(lo2, 1.15), hi, where=hi > np.maximum(lo2, 1.15),
                color=FILL, lw=0, zorder=0)
ax.plot(top[0], top[1], color=BLUE, zorder=3)
ax.plot(bot[0], bot[1], color=BLUE, zorder=3)
ax.plot([s_star], [r_star], 'o', ms=5, color=RED, zorder=6, label='cusp')
ax.text(1.95, 1.252, 'bistable', color='#1a3a55')
ax.text(1.32, 1.21, 'monostable', color='#333333')
ax.set_xlabel(r'$s=p_s/k_1$'); ax.set_ylabel(r'$r=\bar\alpha_2/\bar\alpha_1$')
ax.set_xlim(1.2, 2.6); ax.set_ylim(1.15, 1.45)
ax.legend(frameon=False, loc='upper right')
fig.savefig(os.path.join(OUT, 'experimental_slice.pdf'))
fig.savefig(os.path.join(OUT, 'experimental_slice.png'))
plt.close(fig)


def neq_sr(s, r):
    a1 = A*s**ms/(1+s**ms); a2 = r*A*s**ms/(kap**ms+s**ms)
    return S['neq'](a1, a2)


print('\nFigure 3 checks (A=3, m_s=m=2, kappa_s=1.5, g=%.7f):' % g)
print('  cusp (s_*, r_*) = (%.9f, %.9f)' % (s_star, r_star))
print('  label "bistable" starts at (1.95, 1.252); its box spans about')
print('  s in [1.95, 2.24], r in [1.252, 1.269].  Corner checks:')
for p in [(1.95, 1.252), (1.95, 1.269), (2.24, 1.252), (2.24, 1.269)]:
    print('      (s,r)=(%.2f, %.3f): %d equilibria' % (p[0], p[1], neq_sr(*p)))
print('  label "monostable" starts at (1.32, 1.21); box spans about')
print('  s in [1.32, 1.70], r in [1.20, 1.228].  Corner checks:')
for p in [(1.32, 1.21), (1.70, 1.20), (1.70, 1.228)]:
    print('      (s,r)=(%.2f, %.3f): %d equilibria' % (p[0], p[1], neq_sr(*p)))
print('\nfigures written to', OUT)
