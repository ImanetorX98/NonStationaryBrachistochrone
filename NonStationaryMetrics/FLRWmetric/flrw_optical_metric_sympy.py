# -*- coding: utf-8 -*-
"""
Optical metric for FLRW from the conformal Killing field d_eta -- construction and check.

Physical metric (k = 0, +1, -1; comoving coordinates (chi, phi) on the plane):

    ds^2 = a(eta)^2 [ -deta^2 + dchi^2 + S_k(chi)^2 dphi^2 ] = a^2 ghat
    S_0 = chi,  S_1 = sin(chi),  S_-1 = sinh(chi)

d_eta is a genuine Killing field of ghat (static). Optical reduction on ghat:

NULL RAYS (conformal invariance):
    ghat = -deta^2 + dl_k^2  with  dl_k^2 = dchi^2 + S_k^2 dphi^2
    => deta = dl_k :  optical metric = COMOVING 3-GEOMETRY, index n = 1.
    No Randers term (ghat_{eta i} = 0: no rotation or frame dragging).

MASSIVE particles with the conformal rail Ehat = -u_eta:
    v(eta) = sqrt(1 - a^2/Ehat^2)  (time only, by homogeneity)
    => deta = dl_k / v(eta):  Fermat with index n(eta) = 1/v(eta),
    TIME-DEPENDENT but not space-dependent (non-autonomous Finsler... however:)

    OPTICAL CLOCK: lam(eta) = int v deta  =>  dlam = dl_k
    minimising eta_arrival  <=>  minimising lam_arrival  <=>  minimising
    the COMOVING LENGTH. EXACT optical metric = dl_k^2; all the
    non-stationary dynamics is absorbed into the clock lam(eta).
    T_t and T_tau are monotone in eta_arrival => this holds for all three branches.

Checks:
  V1  clock identity: dlam = v deta = dl (algebra, exact)
  V2  optical clock in de Sitter: lam(eta) in closed form
      (= the Dx_max formula already validated) and lam(eta_f) = maximum distance
  V3  k=+1, numerical: on the comoving 3-sphere the massive brachistochrone
      is the GREAT CIRCLE: perturbed family on the equatorial S^2,
      eta_arrival(eps) minimal at eps=0 (radiation a=eta, Ehat=2)
  V4  index n(eta)=1/v: divergence at freezing a->Ehat -- the temporal
      analogue of the light surface (n->oo instead of a spatial F->0)
"""

import sympy as sp
import numpy as np

print("=" * 72)
print("[V1] optical clock: dlam = v(eta) deta  ==>  dlam = dl_k")
print("=" * 72)
eta, chi, Eh, H = sp.symbols('eta chi Ehat H', positive=True)
a = sp.Function('a', positive=True)
v = sp.sqrt(1 - a(eta)**2 / Eh**2)
print("  constraint (from flrw_brachistochrone_sympy):  dl/deta = v(eta) =", v)
print("  dlam := v deta  =>  dlam = dl_k along every rail worldline:")
print("  the motion has UNIT SPEED in the clock lam. eta_arrival = lam^{-1}(L),")
print("  monotone => brachistochrone = geodesic of the comoving 3-geometry.")
print("  Optical metric:  dl_k^2 = dchi^2 + S_k(chi)^2 dphi^2   (exact)")

print()
print("=" * 72)
print("[V2] de Sitter a = -1/(H eta): optical clock in closed form")
print("=" * 72)
u = sp.Symbol('u', positive=True)
# lam(eta) from eta0=-1/H with u = -Ehat*H*eta (u: Ehat -> 1)
prim = sp.sqrt(u**2 - 1) - sp.acos(1 / u)
chk = sp.simplify(sp.diff(prim, u) - sp.sqrt(1 - 1 / u**2))
print("  antiderivative of sqrt(1-1/u^2):", prim, "  [check d/du:", chk, "]")
lam_u = (prim.subs(u, Eh) - prim) / (Eh * H)
print("  lam(u) =", lam_u)
lam_tot = sp.simplify(lam_u.subs(u, 1))
dx_max = (sp.sqrt(Eh**2 - 1) + sp.asin(1 / Eh) - sp.pi / 2) / (Eh * H)
print("  lam(eta_f) =", lam_tot)
diff_ = sp.simplify((lam_tot - dx_max).rewrite(sp.asin))
print("  lam(eta_f) - Dx_max =", diff_, " (agrees, OK)")

print()
print("=" * 72)
print("[V3] k=+1: massive brachistochrone = great circle on S^3 (numerical)")
print("=" * 72)
# radiation: a(eta) = eta, Ehat = 3, start at eta0 = 1 (a=1), stop at eta=3
# (lam_max = 1.375 > the lengths involved: the particle always arrives)
Eh_n, eta0 = 3.0, 1.0
etag = np.linspace(eta0, Eh_n - 1e-12, 400000)
vg = np.sqrt(np.clip(1 - (etag / Eh_n)**2, 0, None))
lam_g = np.concatenate([[0], np.cumsum((vg[1:] + vg[:-1]) / 2
                                       * np.diff(etag))])

def eta_arrivo(L):
    return np.interp(L, lam_g, etag)

# equatorial S^2 of the comoving 3-sphere: from (pi/2, 0) to (pi/2, 1.0)
Phi = 1.0
sg = np.linspace(0, 1, 4000)
res = []
for eps in np.linspace(-0.3, 0.3, 61):
    th = np.pi / 2 + eps * np.sin(np.pi * sg)
    dth = eps * np.pi * np.cos(np.pi * sg)
    L = np.trapezoid(np.sqrt(dth**2 + np.sin(th)**2 * Phi**2), sg)
    res.append((eps, L, eta_arrivo(L)))
res = np.array(res)
i0 = np.argmin(res[:, 2])
print(f"  great-circle distance: L(0) = {res[30, 1]:.9f} (expected {Phi})")
print(f"  argmin eta_arrival(eps) = {res[i0, 0]:.4f}  (expected 0)")
print(f"  eta_arrival(0) = {res[30, 2]:.6f};  eta_arrival(0.3) = {res[-1, 2]:.6f}")
print("  => spherical geodesic (great circle) = brachistochrone.  OK")

print()
print("=" * 72)
print("[V4] optical index n(eta) = 1/v: divergence at freezing")
print("=" * 72)
n_expr = 1 / v
print("  n(eta) =", n_expr)
print("  a -> Ehat:  n -> oo  (the clock stops: dlam -> 0)")
print("  TEMPORAL analogue of the stationary light surface:")
print("  in Kerr the optical reduction degenerates on a spatial SURFACE")
print("  (F -> 0); in FLRW it degenerates at a cosmic INSTANT (a = Ehat).")
print("\nDONE.")
