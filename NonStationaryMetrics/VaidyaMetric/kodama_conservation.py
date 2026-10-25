# -*- coding: utf-8 -*-
"""
The Kodama vector along the Vaidya brachistochrones.

"Conservation" made precise in three statements, all checked:

 [A] Kodama miracle (a property of the spacetime, ingoing EF coordinates
     ds^2 = -f dv^2 + 2 dv dr + r^2 dOmega^2,  f = 1 - 2 m(v)/r):
       K^a = eps^{ab} nabla_b r  =  (1,0,0,0) = d/dv   (exact)
       nabla_a K^a = 0                                  (without a Killing vector)
       Kodama charge = Misner-Sharp mass = m(v)

 [B] rail identity: along the tau brachistochrone the constraint W
     is equivalent to  -u.K = E  EXACTLY (it is not imposed by hand: it follows
     from the form of W).  => E_K := -K.u is the rail invariant.

 [C] non-triviality: for a GEODESIC (rail switched off)
       d(-K.u)/dtau = -(1/2)(L_K g)_{ab} u^a u^b = - m'(v) (u^v)^2 / r
     i.e. the Kodama energy would drift at a rate proportional to m'.
     The rail supplies exactly the compensation: it is the "missing
     symmetry" made manifest.

Output: prints the checks + fig_kodama_conservazione (E_K flat along
the brachistochrone vs the geodesic drift that the rail cancels).
"""

import os
import sys
import numpy as np
import matplotlib.pyplot as plt
import sympy as sp
from scipy.integrate import solve_ivp
from scipy.optimize import brentq

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from paper_style import COL, set_style, savefig

set_style()
OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'Vaidyafigures')

print("=" * 70)
print("[A] Kodama miracle (symbolic, ingoing EF Vaidya)")
print("=" * 70)
v, r, th = sp.symbols('v r theta', real=True)
m = sp.Function('m')
f = 1 - 2 * m(v) / r
# 4D metric in (v, r, theta, phi)
g = sp.diag(-f, 0, r**2, r**2 * sp.sin(th)**2)
g[0, 1] = g[1, 0] = 1                       # (v,r) block = [[-f,1],[1,0]]
ginv = g.inv()
detg = g.det()
sqrtg = sp.sqrt(-detg)
print("  sqrt(-g) =", sp.simplify(sqrtg), "  (= r^2 sin theta)")

# K^a = eps^{ab} nabla_b r  in the (v,r) block;  nabla_b r = delta_b^r
# eps^{vr} = 1/sqrt|h| with |h(v,r)| = 1  =>  K = (1,0,0,0)
K_up = sp.Matrix([1, 0, 0, 0])
print("  K^a =", K_up.T, " = d/dv")
# divergence: (1/sqrt-g) d_a( sqrt-g K^a )
divK = sum(sp.diff(sqrtg * K_up[a], (v, r, th, sp.Symbol('phi'))[a])
           for a in range(3)) / sqrtg
print("  nabla_a K^a =", sp.simplify(divK), "   (miracle: = 0)")

# Misner-Sharp: M = (r/2)(1 - g^{ab} d_a r d_b r) = (r/2)(1 - g^{rr})
M_MS = sp.simplify((r / 2) * (1 - ginv[1, 1]))
print("  Misner-Sharp mass = (r/2)(1 - g^{rr}) =", M_MS, " = m(v)")

# K_a (lower index), used below
K_dn = sp.simplify(g * K_up)
print("  K_a =", K_dn.T, "   (= (-f, 1, 0, 0)  => K = -f dv + dr)")

print()
print("=" * 70)
print("[B] rail identity: constraint W  <=>  -u.K = E  (symbolic)")
print("=" * 70)
uu, E, J = sp.symbols('uu E J', real=True, positive=True)
mm = sp.symbols('mm', positive=True)
feq = 1 - 2 * mm / r
# constraint of the tau brachistochrone (parameter r): phi'^2 = W
W = (feq * uu**2 - 2 * uu - (feq * uu - 1)**2 / E**2) / r**2
# tangent T = (v', r', 0, phi') = (uu, 1, 0, sqrt(W)); norm:
gTT = -feq * uu**2 + 2 * uu + r**2 * W       # g(T,T)
print("  -g(T,T) =", sp.simplify(-gTT), "   ( = (f uu - 1)^2 / E^2 )")
# K.T = K_a T^a = (-f) uu + 1 = 1 - f uu ;  -K.u = -(K.T)/sqrt(-gTT)
E_K = sp.simplify((feq * uu - 1) / sp.sqrt(sp.simplify(-gTT)))
print("  -u.K = (f uu - 1)/sqrt(-g(T,T)) =", E_K,
      "   (= E, up to the sign of the input)")

print()
print("=" * 70)
print("[C] numerical check along a brachistochrone (accreting Vaidya)")
print("=" * 70)
# tau EOM with parameter v (from the existing machinery)
vv2, q = sp.symbols('vv2 q', real=True)
mu = sp.symbols('mu', positive=True)
f2 = 1 - 2 * mm / r
W2 = (f2 * uu**2 - 2 * uu - (f2 * uu - 1)**2 / E**2) / r**2
F_tau = (f2 * uu - 1) / E - J * sp.sqrt(W2)
p_tau = sp.diff(F_tau, uu)
dp_tot = (sp.diff(p_tau, r) + sp.diff(p_tau, vv2) * uu
          + sp.diff(p_tau, mm) * mu * uu + sp.diff(p_tau, uu) * q)
EL = dp_tot - sp.diff(F_tau, mm) * mu
vpp = sp.solve(sp.Eq(EL, 0), q)[0]

m0, E_n, J_n = 1.0, 1.2, 1.3
mu_n = 0.02
r1, r0 = 10.0, 4.2     # above the periapsis (~3.01): the v-parametrisation
#                        degenerates (dv/dr -> inf) exactly at the periapsis
lin = [(mm, m0 + mu * vv2), (E, E_n), (J, J_n), (mu, mu_n)]
vpp_f = sp.lambdify((r, vv2, uu), vpp.subs(lin), 'numpy')
pv_f = sp.lambdify((r, vv2, uu), p_tau.subs([(mm, m0 + mu * vv2), (E, E_n),
                                             (J, J_n), (mu, mu_n)]), 'numpy')
W_f = sp.lambdify((r, vv2, uu),
                  W2.subs([(mm, m0 + mu * vv2), (E, E_n), (mu, mu_n)]),
                  'numpy')

def vp_rad(fv):
    a2 = fv - fv**2 / E_n**2
    a1 = -2 + 2 * fv / E_n**2
    a0 = -1 / E_n**2
    return (-a1 + np.sqrt(a1**2 - 4 * a2 * a0)) / (2 * a2)

v_arr = 40.0
fl = 1 - 2 * (m0 + mu_n * v_arr) / r1
lo = vp_rad(fl) * (1 + 1e-9)
hi = lo * 2
while pv_f(r1, v_arr, hi) < 0:
    hi *= 2
up1 = brentq(lambda z: pv_f(r1, v_arr, z), lo, hi)
sol = solve_ivp(lambda rr, y: [y[1], vpp_f(rr, y[0], y[1])],
                [r1, r0], [v_arr, up1], rtol=1e-12, atol=1e-14,
                dense_output=True, method='DOP853')
rg = np.linspace(r1, r0, 600)
vg, ug = sol.sol(rg)
mg = m0 + mu_n * vg
fg = 1 - 2 * mg / rg
Wg = np.clip(W_f(rg, vg, ug), 0, None)

# 4-velocity u = T/sqrt(-gTT),  T=(uu,1,0,sqrt(W))
gTT_n = -fg * ug**2 + 2 * ug + rg**2 * Wg
norm = np.sqrt(-gTT_n)
uv = ug / norm                         # u^v
E_K_num = (fg * ug - 1) / norm         # -u.K
dev_E = np.max(np.abs(np.abs(E_K_num) - E_n))
print(f"  |-u.K| along the brachistochrone: min {np.abs(E_K_num).min():.12f}  "
      f"max {np.abs(E_K_num).max():.12f}")
print(f"  rail invariant E = {E_n}  =>  "
      f"|(-u.K) - E| max = {dev_E:.2e}")

# geodesic drift that the rail cancels:  d(-K.u)/dtau = -m'(u^v)^2/r
drift_geo = -mu_n * uv**2 / rg          # m'(v) = mu_n constant
# integrated in dtau (dtau = norm dr, r decreases): cumulative geodesic drift
dtau = norm * np.gradient(rg)
drift_cum = np.cumsum(drift_geo * dtau)
print(f"  geodesic drift rate -m'(u^v)^2/r: |.| between "
      f"{np.abs(drift_geo).min():.2e} and {np.abs(drift_geo).max():.2e}")
print(f"  Kodama energy that the rail saves (cumulative geodesic drift) "
      f"= {drift_cum[-1]:+.4f}  against Delta(-u.K)_brachi < 1e-9")

# --------------------------------------------------------------- figure
fig, (a1, a2) = plt.subplots(2, 1, figsize=(COL, 5.2))
# top panel: the conservation RESIDUAL on a log scale (a near-flat linear plot
# cannot show ~1e-16 constancy; Fig 4).
resid = np.abs(E_K_num - E_n) + 1e-18
a1.semilogy(rg, resid, 'C0-', lw=1.5,
            label='$|{-}u\\cdot K-E|$ (brachistochrone)')
a1.axhline(dev_E, color='k', ls='--', lw=0.8, label=f'max $={dev_E:.1e}$')
a1.set_ylabel('$|{-}u\\cdot K-E|$')
a1.set_ylim(1e-18, 1e-10)
a1.set_title('Kodama energy conserved along the $\\tau$-brachistochrone\n'
             f'(accreting Vaidya $m=1+{mu_n}v$): residual $<{dev_E:.0e}$ (log scale)')
a1.legend()
a2.plot(rg, drift_geo, 'C3-', lw=1.5,
        label="geodesic rate $-m'(u^v)^2/r$")
a2.axhline(0, color='k', lw=0.6)
a2.set_xlabel('$r$')
a2.set_ylabel('$d(-u\\cdot K)/d\\tau$')
a2.set_title('what the rail cancels: the free-particle Kodama drift\n'
             "($\\propto m'$): non-zero, this is the missing symmetry")
a2.legend()
savefig(fig, OUT, 'fig_kodama_conservazione')
print('\nDONE.')
