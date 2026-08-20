# Global fold and cusp geometry of an asymmetric two-gene cross-repressive module

Reproducibility package for

> F. R. Villatoro, *Global fold and cusp geometry of an asymmetric two-gene
> cross-repressive module*, submitted to **SIAM Journal on Applied Dynamical
> Systems** (2026). Contact: `frvillatoro@uma.es`

Every number, table and figure of the paper is recomputed here from scratch.
The scripts are independent of one another and each prints, next to every
computed quantity, the value claimed in the manuscript.

## The model

    dX/dt = a1/(1 + Y^m2) - X   + g1
    dY/dt = a2/(1 + X^m1) - R Y + g2

with Hill exponents `m1, m2 > 0`, basal leaks `g1, g2 >= 0`, degradation ratio
`R > 0` and synthesis strengths `a1, a2 > 0`. The equilibrium map is

    A(X,Y) = ( (X-g1)(1+Y^m2),  R (Y-g2/R)(1+X^m1) ),

its singular set is the level set `Psi1(X) Psi2(Y) = 1` with

    Psi_i(z) = m_i (z - gamma_i) z^(m_i - 1) / (1 + z^m_i),

and the cusp is the unique root of `1 + L1(X) = Psi1(X) [1 + L2(Y)]`, where
`L_i(z) = m_i - (m_i - 1) gamma_i / z`.

## Contents

| script | what it recomputes |
|---|---|
| `01_center_manifold_coefficients.py` | Theorem 6.3 and Appendix A: `A_mu`, `b = (m^2-1)/(12 q_c^2)`, `h`, `C_m`, symbolically in `m` and `q_c`; monotonicity of `Psi` |
| `02_singular_set_and_error_scaling.py` | Theorem 3.2 and Table 3: `det J = 0` on the singular set to 1e-31; cusp constants; slope-one decay of the normal-form error |
| `03_experimental_slice.py` | Table 4 and Figure 3: exact against normal-form bistable windows in `(s,r)`; far-field power law of Corollary 3.4 |
| `04_cusp_equation_and_zero_leak.py` | Lemma 4.1, Table 2 and Corollary 4.4: scalar cusp equation and the zero-leak closed form `X_c = ((m1+1)/(m1 m2 - 1))^(1/m1)` |
| `05_cusp_uniqueness.py` | Theorem 4.3: the sum-of-squares identity of Appendix B and the bound `Xi < 1` at every cusp, over 3132 parameter sets with `m_i` in [0.4, 20] and `gamma_i` in [0, 20] |
| `06_tangency_sign_identity.py` | Lemma 4.1, eq. (4.6): `kappa + phi' = -H / [(X-g1) Psi1 Psi2'(phi)]`, hence `sign(kappa + phi') = -sign(H)` |
| `07_slope_and_injectivity.py` | Lemma 5.1 and Theorem 5.2: the slope identity `w (X-g1)(L2-Psi2) dlogS/dX = w^2 + D w + 1`, monotonicity of the slope, and separation of the two branches of the discriminant |
| `08_equilibrium_counts.py` | Proposition 3.5 and Theorem 5.3: counts 1 and 3 across the discriminant, with index sum +1, for unequal exponents, unequal leaks, `R != 1` and `m_i < 1` |
| `09_figures.py` | Figures 1 and 3, with every region label verified by counting equilibria at the four corners of its bounding box |
| `10_section7_expansions.py` | Section 7: `r_2 = 3b/A_mu` and `t_3 = -2b/h` for symbolic `m`; the closed form of `c_1` for `m = 2`; the zero-leak expansion `1 - (5/288) rho^2 + (29/3456) rho^3`; and the zero-leak values of `c_1` for `m = 2,...,6` showing that the cancellation is specific to `m = 2` |

## Running

    python3 -m pip install -r requirements.txt
    ./run_all.sh

Logs are written to `outputs/`, one per script, and concatenated into
`outputs/ALL.log`. Figures are written to `figures/`. The logs of a reference
run are included in the archive, so the output can be compared without
rerunning anything.

Approximate run times on a laptop: scripts 01, 04, 06 a few seconds each;
02, 03, 07, 08, 09 under two minutes each; 05 about three minutes; 10 about
eight minutes, dominated by the symbolic expansions for exponents 5 and 6.

## Precision

Arbitrary-precision arithmetic is used throughout: 30 digits for the singular
set and the experimental slice, 50 digits for the error scaling of Figure 2.
Double precision is **not** sufficient near the cusp, because the
symmetry-breaking parameter `Delta a = a1 - a2` is a difference of two
quantities of order one whose gap is of order `rho^{3/2}`; at `rho = 1e-8` the
gap is of order `1e-12` and double precision retains no correct digit of the
relative error.

Two independent checks are built in. The singular set is located by safeguarded
Newton inversion of the monotone function `Psi`, and `det J` is then recomputed
from the Jacobian at every tabulated point, vanishing to the working precision.
The cusp locations are obtained from the scalar equation and cross-checked
against the minimum of the parameter-space speed along the singular set.

## Licence

Code: MIT (see `LICENSE`). Figures and the accompanying text: CC BY 4.0.

## Citing

See `CITATION.cff`. Please cite both the archived software and the paper.
