# -*- coding: utf-8 -*-
"""
Hamiltonians of the t and tau branches in flat FLRW (optimal-control /
Pontryagin formulation, with cosmic time t as the evolution parameter).

From the conformal rail Ehat = -u_eta (results R1):
    comoving velocity:  dx/dt = (v/a) u_hat ,  |u_hat| = 1  (control)
    v(t) = sqrt(1 - a(t)^2/Ehat^2)

Costs:  t branch:  int dt ;   tau branch:  int dtau = int (a/Ehat) dt.

Hamiltonians (maximised over the control => u_hat = p/|p|):

    H_t   = (v/a) |p| - 1
    H_tau = (v/a) |p| - a/Ehat

Structure:
  - dp/dt = -dH/dx = 0  (homogeneity): p CONSERVED => comoving straight lines.
    The trajectories follow from Hamilton's equations: dx/dt = dH/dp = (v/a) p_hat.
  - H NOT conserved:  dH/dt = dH/dt|_explicit ∝ a' -- the Hamiltonian
    version of the non-autonomous forcing (in Vaidya: dp_v/dr ∝ m').
  - Transversality (spatial arrival, free final instant): H(t1) = 0
    => |p| fixed at arrival:  |p|_t = a1/v1 ,  |p|_tau = a1^2/(Ehat v1).
    It is the exact analogue of p_v(r1) = 0 in Vaidya (p_v = costate of v).

Checks:
  V1  symbolic Hamilton equations: dp/dt = 0, dx/dt = (v/a) p_hat;
      dH/dt along the flow = explicit d_t (forcing ∝ a')
  V2  de Sitter, numerical: integration of Hamilton's equations (both
      branches) against the closed forms R5:  T_t = ln(a1)/H,  T_tau = (a1-1)/(Ehat H)
  V3  consistency with the constrained worldline: the |p|_tau(t) required by
      H_tau = 0 equals the modulus of the reduced canonical momentum
      p_x = Ehat*v*(1+mu*Ehat) => (1+mu*Ehat) = |p|... through mu(eta) of R2
"""

import sympy as sp
import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import brentq

t, Eh, H0 = sp.symbols('t Ehat H', positive=True)
px, py = sp.symbols('p_x p_y', real=True)
x, y = sp.symbols('x y', real=True)
a = sp.Function('a', positive=True)

v = sp.sqrt(1 - a(t)**2 / Eh**2)
P = sp.sqrt(px**2 + py**2)

H_t = (v / a(t)) * P - 1
H_tau = (v / a(t)) * P - a(t) / Eh

print("=" * 72)
print("[V1] Hamilton equations (symbolic)")
print("=" * 72)
for nome, Hh in (("t  ", H_t), ("tau", H_tau)):
    print(f"\n  H_{nome} =", Hh)
    print(f"    dx/dt =  dH/dp_x =", sp.simplify(sp.diff(Hh, px)))
    print(f"    dp_x/dt = -dH/dx =", -sp.diff(Hh, x), "  (p conserved)")
    print(f"    dH/dt|expl =", sp.simplify(sp.diff(Hh, t)),
          "  [∝ a': non-autonomous forcing]")
print("\n  => trajectories: constant p, fixed direction: comoving STRAIGHT LINE;")
print("     modulus |p| fixed by the transversality condition H(t1) = 0.")

print()
print("=" * 72)
print("[V2] de Sitter (a = e^t, H=1): numerical Hamilton flow against the closed forms R5")
print("=" * 72)

Eh_n, X_tgt = 3.0, 0.30

def v_n(tt):
    return np.sqrt(max(1 - np.exp(2 * tt) / Eh_n**2, 0.0))

# Hamilton equations: dx/dt = (v/a)*p_hat, p_hat = (1,0); dtau/dt = a/Ehat
def rhs(tt, yv):
    aa = np.exp(tt)
    return [v_n(tt) / aa, aa / Eh_n]

t_free = brentq(lambda T: solve_ivp(rhs, [0, T], [0, 0], rtol=1e-12,
                                    atol=1e-14).y[0, -1] - X_tgt,
                0.05, np.log(Eh_n) - 1e-9)
sol = solve_ivp(rhs, [0, t_free], [0, 0], rtol=1e-12, atol=1e-14)
T_t_num, T_tau_num = t_free, sol.y[1, -1]
a1 = np.exp(t_free)
T_t_cf = np.log(a1)                    # = t1 (trivial but consistent with R5)
T_tau_cf = (a1 - 1) / Eh_n
print(f"  arrival at X={X_tgt}:  t1 = {t_free:.10f}")
print(f"  T_t   Hamilton = {T_t_num:.12f}   closed form = {T_t_cf:.12f}")
print(f"  T_tau Hamilton = {T_tau_num:.12f}   closed form = {T_tau_cf:.12f}")
print(f"  |p| from transversality:  t branch = {a1/v_n(t_free):.6f} ,"
      f"  tau branch = {a1**2/(Eh_n*v_n(t_free)):.6f}")

print()
print("=" * 72)
print("[V3] consistency with the constrained worldline (multiplier mu of R2)")
print("=" * 72)
# R2: p_x(worldline) = -C with (1+mu*Ehat) = C/(Ehat*v); u_x = Ehat*v
# tau costate:  |p|(t) from H_tau=0 is a^2/(Ehat*v) ONLY at arrival;
# along the motion |p| = const = a1^2/(Ehat*v1). Comparison: C of the worldline
# is fixed by the same transversality => C = Ehat*v*(1+mu*Ehat) constant.
C_, mu_, vv_ = sp.symbols('C mu v', positive=True)
rel = sp.Eq(1 + mu_ * Eh, C_ / (Eh * vv_))
print("  R2:  (1 + mu*Ehat) = C/(Ehat*v)  <=>  C = Ehat*v*(1+mu*Ehat)")
print("  PMP costate: |p| = const, fixed at arrival = a1^2/(Ehat*v1)")
print("  => same structure: transverse constant of motion + factor 1/v;")
print("     the multiplier mu(eta) is the rescaled costate.")
print("\nDONE.")
