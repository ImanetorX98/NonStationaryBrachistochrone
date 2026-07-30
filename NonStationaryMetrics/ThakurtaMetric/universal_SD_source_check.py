# -*- coding: utf-8 -*-
# Universal form of the off-shell adiabatic source S_D (conformally stationary case).
#
# Result:  S_D(lambda) = lambda + int_0^lambda E_eff H_{E_eff} dlambda'
#
# Proof structure (all checks below pass exactly / to machine precision):
#   (i)  H is an arrival-time Randers form  H = beta*J + sqrt(A_rr pr^2 + A_ff J^2) - 1,
#        homogeneous degree 1 in (pr,J) apart from the -1 normalization.
#   (ii) Euler => (J d_J + pr d_pr) H = H + 1   [UNIVERSAL: holds for arbitrary beta,A_rr,A_ff].
#   (iii)Hence DH = E H_E + (H+1);  on the frozen shell H=0  =>  DH = 1 + E H_E.
#   (iv) Integrate: S_D = lambda + int E H_E dlambda.
# The lambda term is fully universal; E H_E is the single metric-specific "energy-Euler" piece,
# fixed by the same potentials via E d_E v^2 = 2(1-v^2), E d_E Pbar^2 = -2 b^2/E^2.
import numpy as np, sympy as sp
from scipy.integrate import solve_ivp, cumulative_trapezoid as ct
from scipy.optimize import brentq

# ---------- (ii) universality with ABSTRACT Randers functions ----------
pr_s, E_s, J_s, r_s = sp.symbols('pr E J_ r')
beta = sp.Function('beta')(r_s, E_s); Arr = sp.Function('A_rr')(r_s, E_s); Aff = sp.Function('A_ff')(r_s, E_s)
H_abs = beta*J_s + sp.sqrt(Arr*pr_s**2 + Aff*J_s**2) - 1
id_abs = sp.simplify((J_s*sp.diff(H_abs, J_s) + pr_s*sp.diff(H_abs, pr_s)) - (H_abs + 1))
DH_abs = sp.simplify((E_s*sp.diff(H_abs, E_s) + J_s*sp.diff(H_abs, J_s) + pr_s*sp.diff(H_abs, pr_s))
                     - (E_s*sp.diff(H_abs, E_s) + H_abs + 1))
print("UNIVERSAL (abstract Randers):")
print("  (J d_J + pr d_pr)H - (H+1) =", id_abs, " [expect 0]")
print("  DH - (E H_E + H + 1)       =", DH_abs, " [expect 0]")

# ---------- concrete Thakurta-Kerr Hamiltonian ----------
M, a, Ehat = 1.0, 0.9, 1.4
J0, r0 = 6.0, 12.0
rr, pr, E, J = sp.symbols('r pr E J_')
f2 = 1 - 2*M/rr; Dl2 = rr**2 - 2*M*rr + a**2; b2 = 2*M*a/rr
v2 = 1 - f2/E**2; P2 = rr**2 + a**2 + 2*M*a**2/rr; Pb2 = P2 + b2**2/E**2
H = J*b2*v2/Pb2 + sp.sqrt(Dl2*v2/Pb2)*sp.sqrt((Dl2/rr**2)*pr**2 + J**2/Pb2) - 1
print("E d_E v^2 - 2(1-v^2)   =", sp.simplify(E*sp.diff(v2, E) - 2*(1-v2)))
print("E d_E Pbar^2 + 2b^2/E^2 =", sp.simplify(E*sp.diff(Pb2, E) + 2*b2**2/E**2))

Hn = sp.lambdify((rr, pr, E, J), H, 'numpy')
Hp = sp.lambdify((rr, pr, E, J), sp.diff(H, pr), 'numpy')
Hr = sp.lambdify((rr, pr, E, J), sp.diff(H, rr), 'numpy')
HJ = sp.lambdify((rr, pr, E, J), sp.diff(H, J), 'numpy')
HE = sp.lambdify((rr, pr, E, J), sp.diff(H, E), 'numpy')

def prof(rv, Ev, Jv):
    pg = np.linspace(-80, 80, 3001); Hv = Hn(rv, pg, Ev, Jv)
    rts = [brentq(lambda p: Hn(rv, p, Ev, Jv), pg[i], pg[i+1]) for i in range(len(pg)-1)
           if np.isfinite(Hv[i]) and np.isfinite(Hv[i+1]) and Hv[i]*Hv[i+1] < 0]
    ing = [p for p in rts if Hp(rv, p, Ev, Jv) < 0]
    return min(ing) if ing else np.nan

ev = lambda l, y: y[1]; ev.terminal = True; ev.direction = 1
def rhs(l, y):
    rv, pv, ph = y
    return [Hp(rv, pv, Ehat, J0), -Hr(rv, pv, Ehat, J0), HJ(rv, pv, Ehat, J0)]
so = solve_ivp(rhs, [0, 300], [r0, prof(r0, Ehat, J0), 0.0], rtol=1e-12, atol=1e-14,
               max_step=0.004, dense_output=True, events=ev)
lam = np.linspace(0, so.t[-1], 16000); Y = so.sol(lam); rF, prF = Y[0], Y[1]

Hval = Hn(rF, prF, Ehat, J0)
DH = Ehat*HE(rF, prF, Ehat, J0) + J0*HJ(rF, prF, Ehat, J0) + prF*Hp(rF, prF, Ehat, J0)
EHE = Ehat*HE(rF, prF, Ehat, J0)
S_D_direct = ct(DH, lam, initial=0)
S_D_universal = lam + ct(EHE, lam, initial=0)
print("\nNUMERIC (Thakurta-Kerr frozen orbit):")
print("  max|H| on frozen orbit        =", np.nanmax(np.abs(Hval)))
print("  max|DH - (1 + E H_E)| (H=0)   =", np.nanmax(np.abs(DH - (1 + EHE))))
print("  max|S_D - (lam + int E H_E)|  =", np.nanmax(np.abs(S_D_direct - S_D_universal)))
print("  S_D end: direct=%.6f  universal=%.6f" % (S_D_direct[-1], S_D_universal[-1]))
