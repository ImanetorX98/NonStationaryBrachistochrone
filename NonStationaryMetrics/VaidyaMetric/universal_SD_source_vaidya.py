# -*- coding: utf-8 -*-
# Universal form of the off-shell source S_D for VAIDYA (mass-function, spherical, non-conformal),
# Paper I eq:SD-vaidya, checked on the two branch Hamiltonians of the ingoing Vaidya rail
# (eq:Hv, eq:Htau-vaidya), whose flow parameter is the advanced time v:
#
#   S_D(lambda) = int_0^lambda m H_m dlambda' = [ r p_r ]_0^lambda - int_0^lambda ell dlambda',
#   ell = 1 on H_v (so the subtracted term is lambda), ell_tau = (f - d_pr H_tau)/Ehat on H_tau
#   (so the subtracted term is the elapsed proper time Delta tau).
#
# Proof structure: Euler identity of the support function, (J d_J + pr d_pr)H = H + ell;
# degree-zero self-similarity, r H_r + J H_J + m H_m = 0; on the shell H = 0 and along the
# Hamiltonian flow, d(r p_r)/dlambda = p_r H_pr - r H_r = ell + m H_m.
#
# An earlier version of this script integrated the static Schwarzschild Hamiltonian in the
# Killing time t instead of H_v; the identity it checked was the analogous one, not eq:SD-vaidya.
import numpy as np, sympy as sp
from scipy.integrate import solve_ivp, cumulative_trapezoid as ct
from scipy.optimize import brentq

Ehat, m0, r0 = 1.4, 1.0, 12.0
rr, pr, J = sp.symbols('r pr J_'); m = sp.symbols('m', positive=True)
f = 1 - 2*m/rr; w = Ehat**2 - f
Hv = pr*(f - Ehat**2) - 1 + sp.sqrt(w)*sp.sqrt(Ehat**2*pr**2 + J**2/rr**2)          # eq:Hv
Ht = pr*(f - Ehat**2) - Ehat + sp.sqrt(w)*sp.sqrt((Ehat*pr + 1)**2 + J**2/rr**2)    # eq:Htau-vaidya
ells = {'H_v': sp.Integer(1), 'H_tau': (f - sp.diff(Ht, pr))/Ehat}

ok = True
for name, H, J0 in [('H_v', Hv, 6.0), ('H_tau', Ht, 2.5)]:
    ell = ells[name]
    print(f"\n=== {name}  (Ehat={Ehat}, J={J0}, r0={r0}) ===")
    e1 = sp.simplify(m*sp.diff(H, m) + rr*sp.diff(H, rr) + J*sp.diff(H, J))
    e2 = sp.simplify(J*sp.diff(H, J) + pr*sp.diff(H, pr) - (H + ell))
    print("  m H_m + r H_r + J H_J          =", e1, " [expect 0]")
    print("  (J d_J + pr d_pr)H - (H + ell) =", e2, " [expect 0]")
    ok = ok and e1 == 0 and e2 == 0
    L = lambda X: sp.lambdify((rr, pr, J, m), X, 'numpy')
    Hn, Hp, Hr, HJ, Hm, elln = L(H), L(sp.diff(H, pr)), L(sp.diff(H, rr)), L(sp.diff(H, J)), L(sp.diff(H, m)), L(ell + 0*pr)
    pg = np.linspace(-80, 80, 4001); hv = Hn(r0, pg, J0, m0)
    rts = [brentq(lambda p: Hn(r0, p, J0, m0), pg[i], pg[i+1]) for i in range(len(pg)-1)
           if np.isfinite(hv[i]) and np.isfinite(hv[i+1]) and hv[i]*hv[i+1] < 0]
    p0 = min(p for p in rts if Hp(r0, p, J0, m0) < 0)                  # ingoing root
    ev = lambda l, y: y[1]; ev.terminal = True; ev.direction = 1
    so = solve_ivp(lambda l, y: [Hp(y[0], y[1], J0, m0), -Hr(y[0], y[1], J0, m0)],
                   [0, 300], [r0, p0], rtol=1e-12, atol=1e-14, max_step=0.004,
                   dense_output=True, events=ev)
    lam = np.linspace(0, so.t[-1], 16000); rF, prF = so.sol(lam)
    S_direct = ct(m0*Hm(rF, prF, J0, m0), lam, initial=0)
    S_bound = (rF*prF - r0*prF[0]) - ct(elln(rF, prF, J0, m0)*np.ones_like(lam), lam, initial=0)
    eH = np.nanmax(np.abs(Hn(rF, prF, J0, m0))); eS = np.nanmax(np.abs(S_direct - S_bound))
    print(f"  max|H| on the frozen orbit      = {eH:.1e}")
    print(f"  max|S_D - ([r p_r] - int ell)|  = {eS:.1e}   (S_D at the end: {S_direct[-1]:.6f})")
    ok = ok and eH < 1e-9 and eS < 1e-5
print("\nALL CHECKS PASSED" if ok else "\nFAILED")
raise SystemExit(0 if ok else 1)
