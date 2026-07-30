# -*- coding: utf-8 -*-
# Universal form of the off-shell source S_D for VAIDYA (mass-function, spherical, non-conformal).
#
# Result:  S_D(lambda) = [ r p_r ]_0^lambda  -  lambda      (a boundary term minus lambda)
#
# Proof structure (all checks pass exactly / to machine precision):
#   H = sqrt(f v^2) sqrt(f pr^2 + J^2/r^2) - 1,  f = 1 - 2m/r,  v^2 = 1 - f/E^2   (a=0 Randers)
#   (i)   Self-similarity of Schwarzschild: (r,m)->k(r,m) with E,pr invariant, J->kJ leaves H
#         invariant (weight 0)  =>  r H_r + J H_J + m H_m = 0  =>  m H_m = -(r H_r + J H_J).
#   (ii)  Finsler identity still holds:  (J d_J + pr d_pr) H = H + 1.
#   (iii) On the frozen shell H=0:  J H_J = 1 - pr H_pr, hence
#         m H_m = -(r H_r + J H_J) = (pr H_pr - r H_r) - 1.
#   (iv)  Hamilton (dr/dl=H_pr, dpr/dl=-H_r):  pr H_pr - r H_r = d(r pr)/dl.
#         Therefore  S_D = int m H_m dl = [ r pr ]_0^lambda - lambda.
import numpy as np, sympy as sp
from scipy.integrate import solve_ivp, cumulative_trapezoid as ct
from scipy.optimize import brentq

Ehat, m0, J0, r0 = 1.4, 1.0, 6.0, 12.0
rr, pr, E, J = sp.symbols('r pr E J_'); m = sp.symbols('m', positive=True)
f = 1 - 2*m/rr; v2 = 1 - f/E**2
H = sp.sqrt(f*v2)*sp.sqrt(f*pr**2 + J**2/rr**2) - 1

# ---- symbolic identities ----
mHm_s = m*sp.diff(H, m); rHr_s = rr*sp.diff(H, rr); JHJ_s = J*sp.diff(H, J)
print("Vaidya symbolic identities:")
print("  m H_m + (r H_r + J H_J)     =", sp.simplify(mHm_s + (rHr_s + JHJ_s)), " [expect 0]")
print("  m d_m f - (f-1)             =", sp.simplify(m*sp.diff(f, m) - (f-1)))
print("  (J d_J + pr d_pr)H - (H+1)  =", sp.simplify(J*sp.diff(H, J) + pr*sp.diff(H, pr) - (H+1)))

Hn = sp.lambdify((rr, pr, E, J, m), H, 'numpy')
Hp = sp.lambdify((rr, pr, E, J, m), sp.diff(H, pr), 'numpy')
Hr = sp.lambdify((rr, pr, E, J, m), sp.diff(H, rr), 'numpy')
HJ = sp.lambdify((rr, pr, E, J, m), sp.diff(H, J), 'numpy')
Hm = sp.lambdify((rr, pr, E, J, m), sp.diff(H, m), 'numpy')

def prof(rv, Jv):
    pg = np.linspace(-80, 80, 4001); Hv = Hn(rv, pg, Ehat, Jv, m0)
    rts = [brentq(lambda p: Hn(rv, p, Ehat, Jv, m0), pg[i], pg[i+1]) for i in range(len(pg)-1)
           if np.isfinite(Hv[i]) and np.isfinite(Hv[i+1]) and Hv[i]*Hv[i+1] < 0]
    ing = [p for p in rts if Hp(rv, p, Ehat, Jv, m0) < 0]
    return min(ing) if ing else np.nan

ev = lambda l, y: y[1]; ev.terminal = True; ev.direction = 1
def rhs(l, y):
    rv, pv, ph = y
    return [Hp(rv, pv, Ehat, J0, m0), -Hr(rv, pv, Ehat, J0, m0), HJ(rv, pv, Ehat, J0, m0)]
so = solve_ivp(rhs, [0, 300], [r0, prof(r0, J0), 0.0], rtol=1e-12, atol=1e-14,
               max_step=0.004, dense_output=True, events=ev)
lam = np.linspace(0, so.t[-1], 16000); Y = so.sol(lam); rF, prF = Y[0], Y[1]

Hval = Hn(rF, prF, Ehat, J0, m0)
mHm = m0*Hm(rF, prF, Ehat, J0, m0)
rHr = rF*Hr(rF, prF, Ehat, J0, m0); JHJ = J0*HJ(rF, prF, Ehat, J0, m0)
S_D_direct = ct(mHm, lam, initial=0)
S_D_boundary = (rF*prF - r0*prF[0]) - lam
print("\nVaidya numeric (Schwarzschild frozen orbit):")
print("  max|H| on frozen orbit          =", np.nanmax(np.abs(Hval)))
print("  max|m H_m + (r H_r + J H_J)|    =", np.nanmax(np.abs(mHm + (rHr + JHJ))))
print("  max|S_D - ([r p_r] - lambda)|   =", np.nanmax(np.abs(S_D_direct - S_D_boundary)))
print("  S_D end: direct=%.6f  boundary=%.6f" % (S_D_direct[-1], S_D_boundary[-1]))
