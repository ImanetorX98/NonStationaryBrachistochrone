# -*- coding: utf-8 -*-
"""
t- and tau-brachistochrones in ingoing Vaidya -- NON-autonomous constrained worldline.

Metric (equatorial, c=G=1):
    ds^2 = -f dv^2 + 2 dv dr + r^2 dphi^2 ,   f = 1 - 2 m(v)/r

Rail constraint (direct generalisation of Kerr, -u_t = E):
    -u_v = E  held fixed       (d_v is NOT Killing: the rail imposes it,
                                 with a force f_v != 0, like f_eta = aH in FLRW)

Reduction (parameter lambda = r, on a stretch with increasing r; ' = d/dr):
    Lam = dtau/dr = sqrt(f v'^2 - 2 v' - r^2 phi'^2)
    constraint => Lam = (f v' - 1)/E
    eliminating Lam:   phi'^2 = W(v, v', r)
        W = [ f v'^2 - 2 v' - (f v' - 1)^2/E^2 ] / r^2

Reduced functionals (J = Fermat multiplier for Delta_phi, constant
because phi stays cyclic -- spherical symmetry survives):

    tau branch:  F_tau = (f v' - 1)/E - J*sqrt(W)     [ T_tau = int Lam dr ]
    t branch:    F_t   =  v'          - J*sqrt(W)     [ T_v  = int v'  dr ]

Non-autonomous structure:
    F depends on v ONLY through m(v)  =>  dF/dv = (dF/dm) m'(v):
    a forcing proportional to m'. In Schwarzschild (m'=0) p_v = dF/dv'
    is conserved; spatial arrival with free v gives the transversality
    p_v(r_arrival) = 0  =>  stationary case: p_v == 0 along the whole orbit.

Checks:
  V1  dF/dm printed: the explicit non-autonomous forcing (exact chain rule)
  V2  stationary (m'=0, p_v=0): dphi/dr in closed form =
        tau branch:  J  sqrt(w f) r / ( Dl sqrt(Dl - J^2 w) )
        t branch:    K  sqrt(w f) r / ( Dl sqrt(Dl - K^2 w) ),  K = f J / E
      with Dl = r^2 f, w = E^2 - f: the Kerr formulas at a=0 (doranTau.md,
      doranT.md with 2Ma/r -> 0). Numerical check at several points.
  V3  numerical: linear Vaidya m(v) = m0 + mu*v, tau branch, integration
      backwards from the arrival point with p_v(r1) = 0; comparison with
      Schwarzschild (mu=0) at equal (E, J, r0, r1, v_arrival).

Everything in plain symbols (no Function(r)): the Euler-Lagrange equation uses the
explicit total derivative  d/dr = d_r + v' d_v + v'' d_v'.
"""

import sympy as sp
import numpy as np
from scipy.optimize import brentq
from scipy.integrate import solve_ivp

r, vv, uu, q = sp.symbols('r vv uu q', real=True)     # r, v, v', v''
E, J, mm, mp, mu = sp.symbols('E J mm mp mu', positive=True)  # m(v), m'(v)

f = 1 - 2 * mm / r
W = (f * uu**2 - 2 * uu - (f * uu - 1)**2 / E**2) / r**2
phip = sp.sqrt(W)

F_tau = (f * uu - 1) / E - J * phip
F_t = uu - J * phip

def eq_moto(F):
    """Euler-Lagrange: d/dr(dF/du) = dF/dv with m=m(v): returns v'' = RHS and p_v."""
    p = sp.diff(F, uu)
    # total derivative along the curve: m -> m(v(r)), dm/dr = mp*uu
    dp_tot = (sp.diff(p, r) + sp.diff(p, vv) * uu + sp.diff(p, mm) * mp * uu
              + sp.diff(p, uu) * q)
    dFdv = sp.diff(F, mm) * mp                      # chain rule: F(v) via m(v)
    EL = dp_tot - dFdv
    vpp = sp.solve(sp.Eq(EL, 0), q)[0]              # linear in q: fast
    return vpp, p

print("=" * 72)
print("[V1] non-autonomous forcing: dF/dv = (dF/dm) * m'(v)  (chain rule)")
print("=" * 72)
for nome, F in (("tau", F_tau), ("t  ", F_t)):
    dFdm = sp.radsimp(sp.diff(F, mm))
    print(f"\n  branch {nome}:  dF/dm =", sp.factor_terms(dFdm))
print("""
  => m' = 0: p_v conserved; transversality p_v(r_arrival)=0 => p_v == 0.
     m' != 0: dp_v/dr ~ m'(v) != 0; p_v(r)!=0 along the path even with
     p_v(r1)=0: the brachistochrone 'anticipates' the future accretion.""")

print("=" * 72)
print("[V2] stationary limit: p_v = 0  =>  Kerr closed form at a=0")
print("=" * 72)

w = E**2 - f
Dl = r**2 * f

def vp_radiale(f_val, E_val):
    """Large root of W=0 in v' (radial motion, phi'=0)."""
    a2 = f_val - f_val**2 / E_val**2
    a1 = -2 + 2 * f_val / E_val**2
    a0 = -1 / E_val**2
    return (-a1 + np.sqrt(a1**2 - 4 * a2 * a0)) / (2 * a2)

vals = {mm: 1.0, E: 1.2, J: 1.3}
for nome, F, K in (("tau", F_tau, J), ("t  ", F_t, f * J / E)):
    pv = sp.diff(F, uu)
    target = K * sp.sqrt(w * f) * r / (Dl * sp.sqrt(Dl - K**2 * w))
    ok = True
    for rv in (4.0, 7.3, 12.0):
        sub = dict(vals);  sub[r] = rv
        pv_n = sp.lambdify(uu, pv.subs(sub), 'numpy')
        lo = vp_radiale(float(f.subs(sub)), vals[E]) * (1 + 1e-9)
        hi = lo * 2
        while pv_n(hi) < 0:
            hi *= 2
        up_star = brentq(pv_n, lo, hi)
        phi_num = float(phip.subs(sub).subs(uu, up_star))
        phi_cf = float(target.subs(sub))
        ok = ok and abs(phi_num - phi_cf) < 1e-9 * max(1, abs(phi_cf))
        print(f"  branch {nome} r={rv:5.1f}:  dphi/dr num = {phi_num:.12f}"
              f"   closed form = {phi_cf:.12f}")
    print(f"  branch {nome}: MATCH with the Kerr closed form at a=0:", ok)

print()
print("=" * 72)
print("[V3] numerical: tau branch in linear Vaidya m(v) = m0 + mu*v")
print("=" * 72)

m0v, E_n, J_n = 1.0, 1.2, 1.3
r1, r0 = 10.0, 3.5   # arrival at r1, start at r0 (periastron at r~3.01: above it)

vpp_tau, pv_tau = eq_moto(F_tau)
lin = [(mm, m0v + mu * vv), (mp, mu), (E, E_n), (J, J_n)]
args = (r, vv, uu, mu)
vpp_fn = sp.lambdify(args, vpp_tau.subs(lin), 'numpy')
pv_fn = sp.lambdify(args, pv_tau.subs(lin), 'numpy')
W_fn = sp.lambdify(args, W.subs(lin), 'numpy')
Lam_fn = sp.lambdify(args, ((f * uu - 1) / E).subs(lin), 'numpy')

def integra(mu_n, v_arrivo):
    # transversality: p_v(r1) = 0 fixes v'(r1)
    fl = 1 - 2 * (m0v + mu_n * v_arrivo) / r1
    lo = vp_radiale(fl, E_n) * (1 + 1e-9)
    hi = lo * 2
    while pv_fn(r1, v_arrivo, hi, mu_n) < 0:
        hi *= 2
    up1 = brentq(lambda z: pv_fn(r1, v_arrivo, z, mu_n), lo, hi)

    def rhs(rr, y):
        v_, u_ = y[0], y[1]
        return [u_, vpp_fn(rr, v_, u_, mu_n),
                np.sqrt(max(W_fn(rr, v_, u_, mu_n), 0.0)),
                Lam_fn(rr, v_, u_, mu_n)]

    sol = solve_ivp(rhs, [r1, r0], [v_arrivo, up1, 0.0, 0.0],
                    rtol=1e-10, atol=1e-12)
    y0 = sol.y[:, -1]
    return dict(v_part=y0[0], dphi=-y0[2], T_tau=-y0[3],
                T_v=v_arrivo - y0[0],
                pv_part=pv_fn(r0, y0[0], y0[1], mu_n))

v_arr = 40.0
res_v = integra(0.01, v_arr)       # accretion mu = 0.01
res_s = integra(0.0, v_arr)        # Schwarzschild m = m0

print(f"\n  ascent r0={r0} -> r1={r1}, arrival at v={v_arr}, E={E_n}, J={J_n}")
print(f"  {'':14s}{'Schwarzschild':>16s}{'Vaidya mu=0.01':>16s}")
for k, lab in (("v_part", "v start"), ("dphi", "Delta_phi"),
               ("T_v", "T_v = Dv"), ("T_tau", "T_tau"),
               ("pv_part", "p_v start")):
    print(f"  {lab:14s}{res_s[k]:16.6f}{res_v[k]:16.6f}")
print("\n  Vaidya: mass during the passage m(v) in [%.3f, %.3f]" %
      (m0v + 0.01 * res_v['v_part'], m0v + 0.01 * v_arr))
print("  p_v(start) != 0 in the Vaidya case (exactly 0 in Schwarzschild):")
print("  non-autonomous 'memory' of the accretion along the path.")
print("\nDONE.")
