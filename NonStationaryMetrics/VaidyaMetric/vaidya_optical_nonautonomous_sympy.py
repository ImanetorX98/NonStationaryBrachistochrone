# -*- coding: utf-8 -*-
"""
Non-autonomous "optical metrics" for generic Vaidya -- in a weak sense.

Inverting the indicatrix (inverse Legendre transform of the Hamiltonians of section 3d),
for ANY m(v) one obtains closed optical forms, with f = f(v,r):

    tau branch:  dtau = sqrt( dr^2 + f r^2 dphi^2 ) / sqrt(w)
                 [time-dependent Riemann-Jacobi]
    v branch:    dv   = dr/f + (E/(f sqrt(w))) sqrt( dr^2 + f r^2 dphi^2 )
                 [time-dependent Randers: drift dr/f + norm]

    w = E^2 - f,  f = 1 - 2m(v)/r.

With m' = 0 they are EXACTLY the two Perlick (1991) metrics in EF coordinates
(doranTau: L_tau = sqrt((B^2+Af)/w) with B=1, A = r^2 phi'^2, g_rr=0):
Perlick generalised verbatim with f -> f(v,r).

BUT: without a conformal Killing field nothing is rigid -- drift AND norm both depend on v.
There is no underlying static geometry (except for m = mu*v, the homothety,
R5-R6): "optical" only in a non-autonomous sense. The hierarchy:

    FLRW           separable index -> clock, trivial geometry
    Thakurta-Schw  static geometry x index n(eta,r)
    Thakurta-Kerr  static geometry + RIGID wind x index
    Vaidya mu*v    static geometry on the self-similar quotient
    Vaidya m(v)    everything breathes: weak sense only (= dual of the H)

Checks:
  V1  both forms satisfy the constraint/indicatrix for any m(v)
      (numerical residuals ~1e-14 at random points)
  V2  m'=0: dtau reproduces L_tau of doranTau (B=1, A=r^2 phi'^2, a=0)
  V3  Euler-Lagrange equation of the form dv (m'=0, r as parameter): dphi/dr = the
      closed form of doranT at a=0 (K = fJ/E)
"""

import numpy as np
import sympy as sp
from scipy.optimize import brentq

v, r, E, J = sp.symbols('v r E J', positive=True)
dr_, dphi_, dv_ = sp.symbols('dr dphi dv', real=True)
m = sp.Function('m', positive=True)

f = 1 - 2 * m(v) / r
w = E**2 - f

dl2 = dr_**2 + f * r**2 * dphi_**2
dtau_opt = sp.sqrt(dl2) / sp.sqrt(w)
dv_opt = dr_ / f + (E / (f * sp.sqrt(w))) * sp.sqrt(dl2)

print("=" * 72)
print("[V1] the optical forms satisfy the constraint (any m(v))")
print("=" * 72)
# worldline: -ds^2 = f dv^2 - 2 dv dr - r^2 dphi^2 = dtau^2  and
# constraint -u_v = E  <=>  dtau = (f dv - dr)/E
# => eliminating: f w dv^2 - 2 w dr dv - (dr^2 + E^2 r^2 dphi^2) = 0
quadrica = f * w * dv_**2 - 2 * w * dr_ * dv_ \
    - (dr_**2 + E**2 * r**2 * dphi_**2)
res_v = quadrica.subs(dv_, dv_opt)
res_tau = (f * dv_opt - dr_) / E - dtau_opt      # consistency dtau = (f dv-dr)/E

rng = np.random.default_rng(11)
w1 = w2 = 0.0
mv = sp.Symbol('mv', positive=True)
for _ in range(10):
    sub = {E: float(rng.uniform(1.05, 2.5)), r: float(rng.uniform(2.5, 12)),
           m(v): float(rng.uniform(0.5, 1.4)),
           dr_: float(rng.uniform(-1, 1)), dphi_: float(rng.uniform(-1, 1))}
    w1 = max(w1, abs(float(res_v.subs(sub))))
    w2 = max(w2, abs(float(res_tau.subs(sub))))
print(f"  quadric(dv = Randers):         max residual = {w1:.2e}")
print(f"  dtau_opt = (f dv_opt - dr)/E:  max residual = {w2:.2e}")
print("  => they hold for ANY m(v): f = f(v,r) inside the forms")

print()
print("=" * 72)
print("[V2] m'=0: dtau = L_tau di doranTau (B=1, A=r^2 phi'^2, a=0)")
print("=" * 72)
php = sp.Symbol("php", real=True)
L_tau_doran = sp.sqrt((1 + (r**2 * php**2) * f) / w)     # (B^2+Af)/w, B=1
L_tau_nostra = dtau_opt.subs([(dr_, 1), (dphi_, php)])   # per unit dr
diff = L_tau_doran - L_tau_nostra
w3 = 0.0
for _ in range(8):
    sub = {E: float(rng.uniform(1.05, 2.5)), r: float(rng.uniform(2.5, 12)),
           m(v): float(rng.uniform(0.5, 1.4)),
           php: float(rng.uniform(-1, 1))}
    w3 = max(w3, abs(float(diff.subs(sub))))
print(f"  L_tau(doranTau) - dtau_opt/dr: max residual = {w3:.2e}")

print()
print("=" * 72)
print("[V3] m'=0: Euler-Lagrange equation of dv_opt against the doranT closed form (a=0)")
print("=" * 72)
M_n, E_n, J_n = 1.0, 1.2, 1.3
L_v = dv_opt.subs([(dr_, 1), (dphi_, php), (m(v), M_n), (E, E_n)])
pJ = sp.diff(L_v, php)
ok = True
for rv in (4.0, 6.0, 9.0):
    pJ_r = sp.lambdify(php, pJ.subs(r, rv), 'numpy')
    php_star = brentq(lambda x: pJ_r(x) - J_n, -3, 3)
    fv = 1 - 2 * M_n / rv
    wv = E_n**2 - fv
    Dv = rv**2 * fv
    Kv = fv * J_n / E_n
    cf = Kv * rv * np.sqrt(wv * fv) / (Dv * np.sqrt(Dv - Kv**2 * wv))
    ok = ok and abs(php_star - cf) < 1e-9 * cf
    print(f"  r={rv}:  dphi/dr from dL/dphi'=J: {php_star:.10f}"
          f"   closed form: {cf:.10f}")
print("  match:", ok)
print("\n  => Perlick 1991 generalised verbatim with f -> f(v,r);")
print("     without a conformal Killing field, however, drift AND norm both breathe:")
print("     'optical' only in a non-autonomous sense (dual of the Hamiltonians).")
print("\nDONE.")
