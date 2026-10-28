# -*- coding: utf-8 -*-
"""
Trichotomy of CAPTURE into the ergosphere for the equatorial brachistochrones
of Thakurta-Kerr, t (eta) and tau branches. We look for the critical J (in
particular the smallest NEGATIVE J) for which the ingoing orbit PENETRATES
the equatorial ergosphere r_e = 2M.

Equatorial Hamiltonians (theta=pi/2, Sigma=r^2, p_theta=0; conformal
factor A, Kodama/rail energy Ehat):
  f = 1-2M/r,  Dl = r^2-2Mr+a^2,  P = r^2+a^2+2Ma^2/r,  b = 2Ma/r
  vb2 = 1 - A^2 f/Ehat^2 ,  Pb = P + A^2 b^2/Ehat^2
  php0 = b vb2/Pb ,  R2 = vb2 Dl/Pb
  tau: ptp = J - A^2 b/Ehat ;
       H_tau = ptp php0 + sqrt(R2) sqrt((Dl/r^2)pr^2 + ptp^2/Pb) - A^2 f/Ehat
  t/eta: H_t  = J php0 + sqrt(R2) sqrt((Dl/r^2)pr^2 + J^2/Pb) - 1

Method: at r0 solve H=0 for pr (INGOING branch, dr/dlam<0);
integrate the Hamilton flow; events: r=2M (PENETRATES) or pr=0
(SCATTER, turning point above the ergosphere). Scan in J.
"""

import os
import sys
import numpy as np
import sympy as sp
from scipy.integrate import solve_ivp
from scipy.optimize import brentq

HERE = os.path.dirname(os.path.abspath(__file__))
M, a, Eh = 1.0, 0.9, 1.2
r_e = 2 * M                         # equatorial ergosphere
r_plus = M + np.sqrt(M**2 - a**2)
r0 = 8.0

r, pr, J = sp.symbols('r p_r J', real=True)
A = sp.symbols('A', positive=True)
f = 1 - 2 * M / r
Dl = r**2 - 2 * M * r + a**2
P = r**2 + a**2 + 2 * M * a**2 / r
b = 2 * M * a / r
vb2 = 1 - A**2 * f / Eh**2
Pb = P + A**2 * b**2 / Eh**2
php0 = b * vb2 / Pb
R2 = vb2 * Dl / Pb
rad = sp.sqrt((Dl / r**2) * pr**2)      # p_theta=0

def build(branch, A_val):
    if branch == 'tau':
        ptp = J - A**2 * b / Eh
        H = ptp * php0 + sp.sqrt(R2) * sp.sqrt((Dl / r**2) * pr**2
                                               + ptp**2 / Pb) - A**2 * f / Eh
    else:                               # t / eta
        H = J * php0 + sp.sqrt(R2) * sp.sqrt((Dl / r**2) * pr**2
                                             + J**2 / Pb) - 1
    sub = {A: A_val}
    Hf = sp.lambdify((r, pr, J), H.subs(sub), 'numpy')
    dHp = sp.lambdify((r, pr, J), sp.diff(H, pr).subs(sub), 'numpy')
    dHr = sp.lambdify((r, pr, J), sp.diff(H, r).subs(sub), 'numpy')
    return Hf, dHp, dHr

TOL = 1e-3          # r_min <= r_e + TOL  =>  capture (it reaches the ergosphere)

def r_min_orbita(Hf, dHp, dHr, J_val):
    """Minimum r reached by the INGOING orbit from r0 (None if it cannot start)."""
    pg = np.linspace(-400, 400, 200001)
    with np.errstate(invalid='ignore'):
        Hv = Hf(r0, pg, J_val)
    roots = [brentq(lambda p: Hf(r0, p, J_val), pg[i], pg[i + 1])
             for i in range(len(pg) - 1)
             if np.isfinite(Hv[i]) and np.isfinite(Hv[i + 1])
             and Hv[i] * Hv[i + 1] < 0]
    ing = [p for p in roots if dHp(r0, p, J_val) < 0]   # dr/dlam<0, ingoing
    if not ing:
        return None
    p0 = min(ing)                       # the most ingoing root
    ev_turn = lambda t_, y: y[1]
    ev_turn.terminal, ev_turn.direction = True, 1       # pr rises to 0 = turning point
    ev_erg = lambda t_, y: y[0] - (r_e - 1e-4)
    ev_erg.terminal, ev_erg.direction = True, -1        # touches r_e
    s = solve_ivp(lambda t_, y: [dHp(y[0], y[1], J_val),
                                 -dHr(y[0], y[1], J_val)],
                  [0, 5000], [r0, p0], rtol=1e-11, atol=1e-13,
                  method='DOP853', events=[ev_turn, ev_erg], max_step=0.5)
    return float(s.y[0].min())

def cattura(Hf, dHp, dHr, J_val):
    rm = r_min_orbita(Hf, dHp, dHr, J_val)
    return None if rm is None else (rm <= r_e + TOL)

def bordo(Hf, dHp, dHr, Jlo, Jhi):
    """Bisection of the capture/scatter boundary (capture towards the inside)."""
    for _ in range(60):
        Jm = 0.5 * (Jlo + Jhi)
        cm = cattura(Hf, dHp, dHr, Jm)
        if cm is None:
            return None
        if cm == cattura(Hf, dHp, dHr, Jlo):
            Jlo = Jm
        else:
            Jhi = Jm
    return 0.5 * (Jlo + Jhi)

print("=" * 70)
print(f"equatorial Thakurta-Kerr: M={M}, a={a}, Ehat={Eh}, r0={r0}")
print(f"ergosphere r_e={r_e}, horizon r_+={r_plus:.4f}")
print("Capture = the ingoing orbit REACHES r_e (grazing; by R12d the")
print("tau branch ends there). Capture band J in [J_neg, J_pos].")
print("=" * 70)

for A_val in (1.0, 1.1):
    Jc_cf = a * A_val**2 / Eh
    print(f"\n--- A = {A_val}  (closed tau separatrix: J_c = sA^2/Ehat = "
          f"{Jc_cf:.4f}) ---")
    for branch in ('tau', 't'):
        Hf, dHp, dHr = build(branch, A_val)
        Js = np.linspace(-14.0, 3.0, 171)
        st = [cattura(Hf, dHp, dHr, Jv) for Jv in Js]
        bordi = []
        for i in range(len(Js) - 1):
            if st[i] is not None and st[i + 1] is not None \
                    and st[i] != st[i + 1]:
                Jb = bordo(Hf, dHp, dHr, Js[i], Js[i + 1])
                if Jb is not None:
                    bordi.append(Jb)
        bordi.sort()
        Jneg = min([b for b in bordi if b < 0], default=None)
        Jpos = max([b for b in bordi if b >= 0], default=None)
        s_neg = f"{Jneg:+.4f}" if Jneg is not None else "  --  "
        s_pos = f"{Jpos:+.4f}" if Jpos is not None else "  --  "
        extra = ""
        if branch == 'tau' and Jpos is not None:
            extra = f"  (against the closed J_c {Jc_cf:.4f}, diff {abs(Jpos-Jc_cf):.1e})"
        print(f"  branch {branch:3s}: capture band J in [{s_neg}, {s_pos}]"
              f"{extra}")
        if Jneg is not None:
            print(f"           => smallest NEGATIVE J for capture = {Jneg:+.4f}")
print("\nDONE.")
