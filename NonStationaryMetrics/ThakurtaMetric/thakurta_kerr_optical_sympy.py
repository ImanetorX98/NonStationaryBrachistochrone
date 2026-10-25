# -*- coding: utf-8 -*-
"""
Optical metric depending on the conformal time, for Thakurta-Kerr.

From the indicatrix (equatorial, rail Ehat = -u_eta):

    (r^2/Dl) dr^2 + Pb (dphi - phi'_0 deta)^2 = R^2 deta^2 ,
    R^2 = Dl vb^2/Pb ,  phi'_0 = (2Ms/r) vb^2/Pb ,  vb^2 = 1 - A^2 f/Ehat^2

solving the quadratic in deta (future root):

    deta = n(eta, r) * alpha_K  +  beta_K                    [RANDERS]

    alpha_K^2 = r^2 dr^2/(f Dl) + Dl dphi^2/f^2      (Riemannian part)
    beta_K    = - (2Ms/(r f)) dphi                   (gravitomagnetic one-form)
    n(eta,r)  = 1/vb = Ehat / sqrt(Ehat^2 - A(eta)^2 f)

THEOREM (checked below): alpha_K and beta_K are EXACTLY the optical
Randers data of LIGHT in Kerr (null Fermat principle, Gibbons-Werner):
all the physics of mass/rail/expansion is in the scalar index
n(eta,r), which multiplies ONLY the Riemannian part. The one-form
(frame dragging) is rigid: conformally invariant and independent of Ehat.

Limits:
    Ehat -> oo : n -> 1: null optical metric of Kerr (conformal invariance ok)
    A = 1      : n = E/sqrt(w): stationary Perlick metric on Kerr
    s = 0      : beta = 0, alpha = Schwarzschild optics: Thakurta-Schw
    f = 1(s=0,M=0): n = Ehat/sqrt(Ehat^2-A^2) = 1/v FLRW

Singularities:
    f = 0 (ergosphere): alpha and beta diverge -- the optical REDUCTION
        dies at the ergosphere (as in static Kerr: intrinsic 1/F);
        the indicatrix/Hamiltonian stays regular down to the horizon.
    vb = 0 (conformal freezing): n -> oo, a moving "instant-surface".

Checks:
  V1  deta = n alpha + beta satisfies the quadric of the indicatrix
      (exact sympy identity, future root)
  V2  beta and alpha coincide with the Randers data of the NULL Fermat principle of
      Kerr (derived independently of them, from ds^2 = 0)
  V3  limits: s=0, A=1, Ehat->oo, and f->1 (FLRW)
  V4  numerical: dphi/dr from the EL equation of F = n alpha + beta (at fixed eta, plus
      a non-autonomous term that vanishes for A=1) vs the closed t-branch form of Paper II
"""

import numpy as np
import sympy as sp
from scipy.optimize import brentq

r, M, s_, E, A = sp.symbols('r M s Ehat A', positive=True)
dr_, dphi_, deta_ = sp.symbols('dr dphi deta', real=True)

f = 1 - 2 * M / r
Dl = r**2 - 2 * M * r + s_**2
P = r**2 + s_**2 + 2 * M * s_**2 / r
vb2 = 1 - A**2 * f / E**2
Pb = P + A**2 * (2 * M * s_ / r)**2 / E**2
php0 = (2 * M * s_ / r) * vb2 / Pb
R2 = Dl * vb2 / Pb

alpha2 = r**2 * dr_**2 / (f * Dl) + Dl * dphi_**2 / f**2
beta = -(2 * M * s_ / (r * f)) * dphi_
n_idx = 1 / sp.sqrt(vb2)
F = n_idx * sp.sqrt(alpha2) + beta

print("=" * 72)
print("[V1] F = n*alpha + beta satisfies the indicatrix")
print("=" * 72)
quadrica = (r**2 / Dl) * dr_**2 + Pb * (dphi_ - php0 * deta_)**2 \
    - R2 * deta_**2
res_expr = quadrica.subs(deta_, F)
# nested radicals: simplify blows up -- exact numerical check
rng = np.random.default_rng(7)
worst = 0.0
for _ in range(10):
    sub = {M: 1.0, s_: float(rng.uniform(0.1, 0.95)),
           E: float(rng.uniform(1.05, 2.5)), A: float(rng.uniform(0.4, 1.3)),
           r: float(rng.uniform(2.5, 12.0)),
           dr_: float(rng.uniform(-1, 1)), dphi_: float(rng.uniform(-1, 1))}
    if float(vb2.subs(sub)) <= 0:
        continue
    worst = max(worst, abs(float(res_expr.subs(sub))))
print(f"  quadric(deta = n*alpha + beta): max residual over 10 random"
      f" points = {worst:.2e}")

print()
print("=" * 72)
print("[V2] alpha, beta = Randers data of the NULL Fermat principle of Kerr")
print("=" * 72)
# ds^2 = 0 equatorial Kerr (seed; any A cancels):
# -f dt^2 - (4Ms/r) dt dphi + (r^2/Dl) dr^2 + P dphi^2 = 0
dt = sp.Symbol('dt', real=True)
null_cond = -f * dt**2 - (4 * M * s_ / r) * dt * dphi_ \
    + (r**2 / Dl) * dr_**2 + P * dphi_**2
# future root: dt = -(2Ms/(rf))dphi + sqrt(alpha_K^2)
F_null = beta + sp.sqrt(alpha2)
res2_expr = null_cond.subs(dt, F_null)
worst2 = 0.0
for _ in range(10):
    sub = {M: 1.0, s_: float(rng.uniform(0.1, 0.95)),
           r: float(rng.uniform(2.5, 12.0)),
           dr_: float(rng.uniform(-1, 1)), dphi_: float(rng.uniform(-1, 1))}
    worst2 = max(worst2, abs(float(res2_expr.subs(sub))))
print(f"  null condition( dt = beta + alpha ): max residual = {worst2:.2e}")
print("  => light: n=1; massive+expansion: ONLY n changes, the one-form")
print("     gravitomagnetic beta = -(2Ms/rf) dphi is RIGID")
print("     (conformally invariant and independent of Ehat).")

print()
print("=" * 72)
print("[V3] limits")
print("=" * 72)
F_s0 = F.subs(s_, 0)
target_s0 = (1 / sp.sqrt(vb2.subs(s_, 0))) \
    * sp.sqrt(dr_**2 / f**2 + r**2 * dphi_**2 / f)
worst3 = 0.0
for _ in range(8):
    sub = {M: 1.0, E: float(rng.uniform(1.05, 2.5)),
           A: float(rng.uniform(0.4, 1.3)), r: float(rng.uniform(2.5, 12.0)),
           dr_: float(rng.uniform(-1, 1)), dphi_: float(rng.uniform(-1, 1))}
    if float(vb2.subs(s_, 0).subs(sub)) <= 0:
        continue
    worst3 = max(worst3, abs(float((F_s0 - target_s0).subs(sub))))
print(f"  s=0:  F - n*alpha_Schw: max residual = {worst3:.2e}")
n_A1 = sp.simplify(n_idx.subs(A, 1))
print("  A=1:  n =", n_A1, " = E/sqrt(w)  (Perlick on Kerr)")
n_null = sp.limit(n_idx, E, sp.oo)
print("  Ehat->oo:  n =", n_null, " (null optics of Kerr)")
n_flrw = sp.simplify(n_idx.subs([(M, 0), (s_, 0)]))
print("  M=s=0:  n =", n_flrw, " = 1/v FLRW  (clock lam)")

print()
print("=" * 72)
print("[V4] numerical: EL equation of F (A=1) vs closed form doranT")
print("=" * 72)
M_n, s_n, E_n, J_n, A_n = 1.0, 0.9, 1.2, 1.3, 1.0
# EL with r as the parameter: F_r(phi', r) = n*sqrt(a_rr + a_pp phi'^2) + b phi'
php = sp.Symbol("php", real=True)
F_r = (n_idx * sp.sqrt(r**2 / (f * Dl) + Dl * php**2 / f**2)
       - (2 * M * s_ / (r * f)) * php)
F_rn = F_r.subs([(M, M_n), (s_, s_n), (E, E_n), (A, A_n)])
# conserved Fermat momentum: dF/dphi' = J
pJ = sp.diff(F_rn, php)
ok = True
for rv in (4.0, 6.0, 9.0):
    pJ_r = sp.lambdify(php, pJ.subs(r, rv), 'numpy')
    lo, hi = -3.0, 3.0
    php_star = brentq(lambda x: pJ_r(x) - J_n, lo, hi)
    fv = 1 - 2 * M_n / rv
    wv = E_n**2 - fv
    Dv = rv**2 - 2 * M_n * rv + s_n**2
    Kv = (fv * J_n + 2 * M_n * s_n / rv) / E_n
    cf = Kv * rv * np.sqrt(wv * fv) / (Dv * np.sqrt(Dv - Kv**2 * wv))
    ok = ok and abs(php_star - cf) < 1e-9 * cf
    print(f"  r={rv}:  dphi/dr from dF/dphi'=J: {php_star:.10f}"
          f"   closed form doranT: {cf:.10f}")
print("  match (the eta branch of the optical metric = the t branch of Kerr):", ok)
print("\nDONE.")
