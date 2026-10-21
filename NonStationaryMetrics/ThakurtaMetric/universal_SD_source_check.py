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

# ---------- proper-time branch H_tau (paper eq:Htau-tk) ----------
# H_tau is NOT of the arrival-time form: its cost is f/E and its azimuthal momentum is shifted,
# pt = J - b/E.  Euler in (P_r, pt) then gives
#     (J d_J + pr d_pr) H_tau = H_tau + (f + b H_J)/E,
# and along the frozen flow (f + b dphi/dlambda)/E = dtau/dlambda, so on the shell
#     S_D = Delta tau + int_0^lambda E_eff H_{E_eff} dlambda'.
# The arrival-time formula S_D = lambda + int E H_E does NOT hold on this branch.
J0t = 2.5
Ht = (J - b2/E)*b2*v2/Pb2 + sp.sqrt(Dl2*v2/Pb2)*sp.sqrt((Dl2/rr**2)*pr**2 + (J - b2/E)**2/Pb2) - f2/E
# identity checked exactly in Mathematica; here at random points (sympy's simplify is slow)
idt_f = sp.lambdify((rr, pr, E, J), J*sp.diff(Ht, J) + pr*sp.diff(Ht, pr) - (Ht + (f2 + b2*sp.diff(Ht, J))/E), 'numpy')
rng = np.random.default_rng(1)
pts = [(rng.uniform(2.5, 20), rng.uniform(-3, 3), rng.uniform(1.05, 3), rng.uniform(-8, 8)) for _ in range(20)]
print("\nPROPER-TIME branch:")
print("  max|(J d_J + pr d_pr)H_tau - (H_tau + (f + b H_J)/E)| at 20 random points =",
      max(abs(idt_f(*q)) for q in pts), " [expect ~1e-15]")
Lt = lambda X: sp.lambdify((rr, pr, E, J), X, 'numpy')
tHn, tHp, tHr, tHJ, tHE = Lt(Ht), Lt(sp.diff(Ht, pr)), Lt(sp.diff(Ht, rr)), Lt(sp.diff(Ht, J)), Lt(sp.diff(Ht, E))
pg = np.linspace(-80, 80, 3001); hv = tHn(r0, pg, Ehat, J0t)
rts = [brentq(lambda p: tHn(r0, p, Ehat, J0t), pg[i], pg[i+1]) for i in range(len(pg)-1)
       if np.isfinite(hv[i]) and np.isfinite(hv[i+1]) and hv[i]*hv[i+1] < 0]
p0t = min(p for p in rts if tHp(r0, p, Ehat, J0t) < 0)
sot = solve_ivp(lambda l, y: [tHp(y[0], y[1], Ehat, J0t), -tHr(y[0], y[1], Ehat, J0t)],
                [0, 300], [r0, p0t], rtol=1e-12, atol=1e-14, max_step=0.004,
                dense_output=True, events=ev)
lt = np.linspace(0, sot.t[-1], 16000); rT, pT = sot.sol(lt)
DHt = Ehat*tHE(rT, pT, Ehat, J0t) + J0t*tHJ(rT, pT, Ehat, J0t) + pT*tHp(rT, pT, Ehat, J0t)
EHEt = Ehat*tHE(rT, pT, Ehat, J0t)
dtau = ((1 - 2*M/rT) + (2*M*a/rT)*tHJ(rT, pT, Ehat, J0t))/Ehat
S_dir_t = ct(DHt, lt, initial=0)
S_tau = ct(dtau, lt, initial=0) + ct(EHEt, lt, initial=0)
S_arr = lt + ct(EHEt, lt, initial=0)
print("  max|S_D - (Delta tau + int E H_E)| =", np.nanmax(np.abs(S_dir_t - S_tau)))
print("  S_D end: direct=%.6f  Delta tau + int E H_E=%.6f  (arrival formula lam + int E H_E=%.6f)"
      % (S_dir_t[-1], S_tau[-1], S_arr[-1]))
