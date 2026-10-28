# -*- coding: utf-8 -*-
"""
Optical metric for Vaidya: is it possible?

- GENERIC Vaidya: NO. The optical/Randers reduction requires a timelike
  Killing field (Perlick 1991); d_v is not Killing if m'(v) != 0
  (checked in vaidya_metric_sympy.py), and there is not even a generic
  conformal Killing field.

- LINEAR Vaidya m(v) = mu*v: YES, through a HOMOTHETY. The field
      xi = v d_v + r d_r
  satisfies Lie_xi g = 2 g (homothetic Killing field, a special conformal one).
  Adapted coordinates:  s = ln v,  x = r/v:
      g = e^{2s} * ghat(x) ,   d_s a GENUINE Killing field of ghat
  => the same move as in FLRW (g = a^2 ghat), with e^s in place of a.
  On the static ghat the standard optical/Randers reduction is built.

Checks:
  V1  Lie_xi g = 2 g  (exact homothety, only for m = mu*v)
  V2  in (s, x, phi):  g = e^{2s} ghat(x),  d ghat/ds = 0
      ghat = (2x - f) ds^2 + 2 ds dx + x^2 dphi^2,  f = 1 - 2 mu/x
  V3  light surfaces of the reduction: A = -ghat_ss = f - 2x = 0
        <=>  2 x^2 - x + 2 mu = 0  =>  x_pm = (1 ± sqrt(1 - 16 mu))/4
      they EXIST only for mu <= 1/16:  critical mu_c = 1/16 (known in the
      literature on self-similar Vaidya / naked singularities)
  V4  x_pm are self-similar null hypersurfaces: outgoing ray
      dr/dv = f/2 with x = const  <=>  x = f(x)/2  <=>  2x^2 - x + 2mu = 0
      => x_- is the EVENT HORIZON (r_EH = x_- v), x_+ the
      outer escape surface. Numerical check: bisection on the rays.
  V5  optical (Randers) metric for null rays on ghat:
        ds_arrival = [ dx + sqrt( dx^2 + A x^2 dphi^2 ) ] / A
      that is  F = beta + alpha:
        alpha^2 = dx^2/A^2 + x^2 dphi^2 / A ,   beta = dx / A
      it minimises s_arrival = ln v_arrival: it is the null brachistochrone in v!
      Check: F satisfies the null condition of ghat. Singular at
      A = 0 (x_pm), like 1/F at the Kerr ergosphere.
"""

import sympy as sp

print("=" * 72)
print("[V1] homothety of linear Vaidya: Lie_xi g = 2g, xi = v d_v + r d_r")
print("=" * 72)

v, r, x, s, mu = sp.symbols('v r x s mu', positive=True)
mfun = sp.Function('m', positive=True)

# equatorial metric (v, r, phi)
def metrica(mass_expr):
    f = 1 - 2 * mass_expr / r
    return sp.Matrix([[-f, 1, 0], [1, 0, 0], [0, 0, r**2]])

X = [v, r, sp.Symbol('phi')]
xi = [v, r, 0]

def lie_g(g):
    L = sp.zeros(3, 3)
    for a in range(3):
        for b in range(3):
            L[a, b] = sum(xi[c] * g[a, b].diff(X[c])
                          + g[c, b] * sp.diff(xi[c], X[a])
                          + g[a, c] * sp.diff(xi[c], X[b]) for c in range(3))
    return sp.simplify(L)

g_lin = metrica(mu * v)
print("  m = mu*v:      Lie_xi g - 2g =", sp.simplify(lie_g(g_lin) - 2 * g_lin))
g_gen = metrica(mfun(v))
resid = sp.simplify(lie_g(g_gen) - 2 * g_gen)
print("  generic m(v): (Lie_xi g - 2g)_vv =", resid[0, 0])
print("  => homothety  <=>  v m' = m  <=>  m = mu*v  (linear ONLY)")

print()
print("=" * 72)
print("[V2] self-similar coordinates  s = ln v,  x = r/v")
print("=" * 72)

# transformation: v = e^s, r = x e^s; Jacobian d(v,r)/d(s,x)
Jm = sp.Matrix([[sp.exp(s), 0], [x * sp.exp(s), sp.exp(s)]])
Jfull = sp.eye(3)
Jfull[0:2, 0:2] = Jm
g_old = g_lin.subs([(v, sp.exp(s)), (r, x * sp.exp(s))])
g_new = sp.simplify(Jfull.T * g_old * Jfull)
ghat = sp.simplify(g_new / sp.exp(2 * s))
print("  g = e^{2s} * ghat,   ghat =")
sp.pprint(ghat)
print("  d ghat/ds =", sp.simplify(ghat.diff(s)), "  => d_s is a KILLING field of ghat")

print()
print("=" * 72)
print("[V3] light surfaces of the reduction: A = -ghat_ss = 0")
print("=" * 72)

A = sp.simplify(-ghat[0, 0])
print("  A = f - 2x =", sp.factor(A), "   [f = 1 - 2mu/x]")
sols = sp.solve(sp.Eq(A, 0), x)
print("  A = 0  <=>  2x^2 - x + 2mu = 0  =>  x_pm =", sols)
print("  real  <=>  mu <= 1/16:  CRITICAL ACCRETION mu_c = 1/16")

print()
print("=" * 72)
print("[V4] x_- = self-similar event horizon (numerical check)")
print("=" * 72)

# x = const is null iff dr/dv = x = f(x)/2
cond = sp.simplify(x - (1 - 2 * mu / x) / 2)
print("  outgoing ray with x=const:  x - f/2 =",
      sp.factor(sp.together(cond)), " = 0  <=> the same x_pm")

import numpy as np
from scipy.integrate import solve_ivp
mu_n = 0.02
xm = float(sols[0].subs(mu, mu_n))
xp = float(sols[1].subs(mu, mu_n))
xm_, xp_ = min(xm, xp), max(xm, xp)

def raggio(x0, v_end=1e6):
    # dx/dv = (f/2 - x)/v, da v=1
    sol = solve_ivp(lambda vv, y: [((1 - 2 * mu_n / y[0]) / 2 - y[0]) / vv],
                    [1.0, v_end], [x0], rtol=1e-10, atol=1e-13)
    return sol.y[0, -1]

lo, hi = 2 * mu_n * 1.01, 0.4     # between the AH (x=2mu) and x_+
for _ in range(60):
    mid = 0.5 * (lo + hi)
    if raggio(mid) > xm_:          # escapes towards x_+
        hi = mid
    else:
        lo = mid
print(f"  mu = {mu_n}:  x_- analytic = {xm_:.8f}   "
      f"x_EH bisection = {0.5*(lo+hi):.8f}")
print(f"             x_+ analytic = {xp_:.8f}   "
      f"ray from x=0.1 ends at x = {raggio(0.1):.8f}")
print(f"  apparent horizon: x_AH = 2mu = {2*mu_n:.4f}  <  x_- (EH outside)")

print()
print("=" * 72)
print("[V5] Randers optical metric on ghat (null brachistochrone in v)")
print("=" * 72)

dx, dphi = sp.symbols('dx dphi', real=True)
ds_opt = (dx + sp.sqrt(dx**2 + A * x**2 * dphi**2)) / A
null_res = sp.simplify(ghat[0, 0] * ds_opt**2 + 2 * ghat[0, 1] * ds_opt * dx
                       + ghat[2, 2] * dphi**2)
print("  ds_arrival = [dx + sqrt(dx^2 + A x^2 dphi^2)]/A")
print("  residual of the null condition ghat(ds,dx,dphi) =", null_res, " (expected 0)")
print("  => Randers:  F = alpha + beta,")
print("     alpha^2 = dx^2/A^2 + x^2 dphi^2/A ,   beta = dx/A")
print("     s_arrival = ln v_arrival: minimising it = null brachistochrone in v")
print("     Singular at A=0 (x_pm): exact analogue of 1/F at the ergosphere.")
print("\n  NOTE, massive case: -uhat_s = -u.xi is conserved by the Killing field of ghat;")
print("  for the physical g  -u.xi is NOT conserved (homothety: d(u.xi)/dtau = -1")
print("  on geodesics). Rail on Ehat_homot = -u.xi: a natural constraint,")
print("  the functional v_arrival = e^{s_arrival} is monotone in s_arrival =>")
print("  the massive t branch (arrival in v) reduces to a STATIC problem on")
print("  ghat, as in FLRW. The tau branch does NOT: int dtau = int e^s dtau_hat")
print("  is not conformally/homothetically invariant.")
print("\nDONE.")
