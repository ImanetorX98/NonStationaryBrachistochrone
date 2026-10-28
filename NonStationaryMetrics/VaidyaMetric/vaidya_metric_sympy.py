# -*- coding: utf-8 -*-
"""
Vaidya metric -- SymPy study + numerical demo.

INGOING Vaidya (advanced Eddington-Finkelstein coordinates, c=G=1):

    ds^2 = -(1 - 2 m(v)/r) dv^2 + 2 dv dr + r^2 (dth^2 + sin^2 th dph^2)

Checks:
  V1  Einstein: G_munu = (2 m'(v)/r^2) (dv)_mu (dv)_nu  -- the only component is G_vv
      => T_munu = [m'(v)/(4 pi r^2)] l_mu l_nu, with l_mu = -(dv)_mu:
      ingoing radial NULL DUST.
  V2  l^mu null and affinely geodesic (l^mu = -g^{mu v}: ingoing radial ray).
  V3  weak/null energy condition: rho ~ m'(v) >= 0
      => physical accretion; m'(v) < 0 violates the NEC (Hawking toy model).
  V4  m'=0 => Schwarzschild in advanced EF coordinates (Killing field dv recovered).
  V5  d_v is NOT Killing if m'(v) != 0:  (Lie_v g)_vv = +2 m'/r.
  V6  apparent horizon: outgoing expansion theta_+ ~ (1 - 2m(v)/r)/r
      => AH at r = 2 m(v), moving.
  V7  (numerical) m(v) = m0 + mu*v: EVENT horizon r_EH(v) > r_AH = 2m(v)
      (the EH anticipates the future accretion; the AH is local, the EH teleological).
"""

import sympy as sp

v, r, th, ph = sp.symbols('v r theta phi', real=True)
m = sp.Function('m', positive=True)(v)
X = [v, r, th, ph]

f = 1 - 2 * m / r
g = sp.Matrix([
    [-f, 1, 0, 0],
    [1, 0, 0, 0],
    [0, 0, r**2, 0],
    [0, 0, 0, r**2 * sp.sin(th)**2],
])
ginv = g.inv()

# --- Christoffel, Ricci, Einstein ---
Gam = [[[sum(ginv[a, d] * (g[d, b].diff(X[c]) + g[d, c].diff(X[b])
                           - g[b, c].diff(X[d])) for d in range(4)) / 2
         for c in range(4)] for b in range(4)] for a in range(4)]
Gam = [[[sp.simplify(Gam[a][b][c]) for c in range(4)]
        for b in range(4)] for a in range(4)]

def ricci(a, b):
    ex = sum(Gam[c][a][b].diff(X[c]) - Gam[c][a][c].diff(X[b])
             for c in range(4))
    ex += sum(Gam[c][c][d] * Gam[d][a][b] - Gam[c][b][d] * Gam[d][a][c]
              for c in range(4) for d in range(4))
    return sp.simplify(ex)

Ric = sp.Matrix(4, 4, lambda a, b: ricci(a, b))
Rs = sp.simplify(sum(ginv[a, b] * Ric[a, b] for a in range(4) for b in range(4)))
G = sp.simplify(Ric - g * Rs / 2)

print("=" * 70)
print("INGOING VAIDYA:  ds^2 = -(1-2m(v)/r) dv^2 + 2 dv dr + r^2 dOmega^2")
print("=" * 70)
print("\n[V1] Ricci scalar R =", Rs)
print("     Einstein G_munu: non-zero components:")
for a in range(4):
    for b in range(a, 4):
        if G[a, b] != 0:
            print(f"       G_{X[a]}{X[b]} =", G[a, b])

# null dust: l_mu = -(dv)_mu
l_dn = sp.Matrix([-1, 0, 0, 0])
T_expected = (m.diff(v) / (4 * sp.pi * r**2)) * (l_dn * l_dn.T)
resid = sp.simplify(G - 8 * sp.pi * T_expected)
print("     G - 8 pi [m'/(4 pi r^2)] l l  =", resid.norm(), " (expected 0)")
print("     => T_munu = m'(v)/(4 pi r^2) l_mu l_nu : radial NULL DUST")

print("\n[V2] l^mu = g^{mu nu} l_nu =", list(ginv * l_dn), " (radial, INGOING)")
l_up = ginv * l_dn
print("     l.l =", sp.simplify((l_up.T * g * l_up)[0]), " (null)")
# geodesic: l^b nabla_b l^a
geo = [sp.simplify(sum(l_up[b] * (sp.diff(l_up[a], X[b])
                                  + sum(Gam[a][b][c] * l_up[c] for c in range(4)))
                       for b in range(4))) for a in range(4)]
print("     l^b nabla_b l^a =", geo, " (affinely geodesic)")

print("\n[V3] density seen by an observer u:  T_uu = m'/(4 pi r^2) (l.u)^2")
print("     => NEC/WEC satisfied  <=>  m'(v) >= 0  (ACCRETION)")
print("        m'(v) < 0: negative energy flux (Hawking toy model, violates the NEC)")

print("\n[V4] static limit m'=0:")
G_static = G.subs(m.diff(v), 0)
print("     G_munu|_{m'=0} =", sp.simplify(G_static.norm()),
      " => vacuum: Schwarzschild in advanced EF coordinates")

print("\n[V5] Lie_{d_v} g: the only non-zero component:")
lie_vv = sp.simplify(g[0, 0].diff(v))
print("     (Lie_v g)_vv =", lie_vv, "  != 0 if m' != 0 => no timelike Killing field")

print("\n[V6] apparent horizon:")
print("     OUTGOING null congruence: k^v = 1, k^r = f/2 = (1-2m/r)/2")
k_up = sp.Matrix([1, f / 2, 0, 0])
print("     k.k =", sp.simplify((k_up.T * g * k_up)[0]))
# expansion of the spheres along k: theta_+ = (2/r) k^r
theta_p = sp.simplify(2 * k_up[1] / r)
print("     theta_+ = (2/r) dr/dlam =", theta_p)
print("     theta_+ = 0  <=>  r = 2 m(v):  moving APPARENT HORIZON")
print("     r_AH' = 2 m'(v): it grows under accretion, recedes under evaporation")

print("\n[note] OUTGOING Vaidya (retarded u):"
      " ds^2 = -(1-2m(u)/r) du^2 - 2 du dr + r^2 dOmega^2")
print("       T ~ -m'(u): physical for m'(u) <= 0 => RADIATING STAR"
      " (Vaidya 1951)")

# ---------------------------------------------------------------------------
# V7: numerical -- event horizon against apparent horizon, linear accretion
# ---------------------------------------------------------------------------
print()
print("=" * 70)
print("[V7] NUMERICAL: m(v) = m0 + mu*v, m0=1, mu=0.02, accretion for "
      "v in [0, 50], then constant m")
print("=" * 70)

import numpy as np
from scipy.integrate import solve_ivp

m0, mu_acc, v_stop = 1.0, 0.02, 50.0

def mass(vv):
    return m0 + mu_acc * min(vv, v_stop)

def rhs(vv, y):
    return [0.5 * (1.0 - 2.0 * mass(vv) / y[0])]

def escapes(r0, v_end=2000.0):
    sol = solve_ivp(rhs, [0, v_end], [r0], rtol=1e-10, atol=1e-12,
                    dense_output=True, max_step=1.0)
    return sol.y[0, -1] > 10.0 * mass(v_end), sol

# bisection on the initial radius of the critical generator (at v=0)
lo, hi = 2.0 * m0, 6.0 * m0
for _ in range(60):
    mid = 0.5 * (lo + hi)
    esc, _ = escapes(mid)
    if esc:
        hi = mid
    else:
        lo = mid
r_EH0 = 0.5 * (lo + hi)
print(f"\n  event horizon at v=0:     r_EH(0) = {r_EH0:.6f}")
print(f"  apparent horizon at v=0: r_AH(0) = 2*m(0) = {2*m0:.6f}")
print(f"  => EH OUTSIDE the AH by {r_EH0 - 2*m0:.4f} (teleological: it 'knows' "
      "about the future accretion)")

_, sol_c = escapes(r_EH0)
for vq in (0.0, 10.0, 25.0, 40.0, 50.0, 80.0):
    r_eh = float(sol_c.sol(vq)[0])
    print(f"  v={vq:5.1f}:  r_EH = {r_eh:8.4f}   r_AH = 2m(v) = "
          f"{2*mass(vq):8.4f}   gap = {r_eh - 2*mass(vq):+.4f}")
print("  after the stop (v>50) the gap -> 0: EH and AH converge to "
      f"r = {2*mass(1e9):.2f} (final Schwarzschild)")
print("\nDONE.")
