# -*- coding: utf-8 -*-
"""
Quasi-constants for NON-constant A(eta) -- dynamical 3D Thakurta-Kerr.

3D PMP Hamiltonian (new, off the equator, eta as the evolution parameter):
from the constraint -u_eta = Ehat on conformal Kerr (BL, Sigma = r^2+a^2cos^2th):

    f_S = 1 - 2Mr/Sigma ,   b = 2Mar sin^2th/Sigma ,
    G   = (r^2+a^2+2Ma^2 r sin^2th/Sigma) sin^2th
    vbS^2 = 1 - A^2 f_S/Ehat^2 ,   Gb = G + A^2 b^2/Ehat^2
    phi'_0 = b vbS^2/Gb ,          R^2 = vbS^2 * Dl * sin^2th / Gb
        [key identity: f_S G + b^2 = Dl sin^2th  (determinant of the t-phi block)]

    H_tau = pt_phi phi'_0 + R sqrt( (Dl/Sigma)p_r^2 + p_th^2/Sigma
                                    + pt_phi^2/Gb )  -  A^2 f_S/Ehat
    pt_phi = J - A^2 b/Ehat          (conformal gravitomagnetic shift)

Structural checks:
  V1  identity f_S G + b^2 = Dl sin^2th (exact, sympy) and equatorial
      reduction = H_tau of R7 (numerical, random points)
  V2  a = 0, dynamical A(eta): L^2 = p_th^2 + J^2/sin^2th conserved
      EXACTLY along the non-autonomous flow (spatial Noether theorem):
      expected drift ~ machine precision
  V3  hierarchy of the candidates along the dynamical flow (a=0.4):
        L^2 (reference), Q_std(E_eff(eta)), K_NLO(E_eff(eta)),
        K_NLO(E frozen at launch)
      + scaling of the new drift in eps = A'/A (H_c vs H_c/2)
  V4  source of the transport equation for the counterterm h:
        dK/deta|_expl = -eps * E * dK/dE = eps a^2 S1(r,th,p_th)
      S1 printed explicitly;  {h, H_0} = -S1  is the subsequent
      analytical computation (h = the first object with no analogue in Kerr).
"""

import numpy as np
import sympy as sp
from scipy.integrate import solve_ivp

M_n, Eh_n = 1.0, 1.2

# ------------------------------------------------------- symbolic 3D H
ee = sp.Symbol('ee', negative=True)          # eta
Hc = sp.Symbol('H_c', positive=True)
r, th = sp.symbols('r theta', positive=True)
pr, pth = sp.symbols('p_r p_theta', real=True)
a, J = sp.symbols('a J', positive=True)
A = sp.Symbol('A', positive=True)            # conformal factor (symbol)

Sig = r**2 + a**2 * sp.cos(th)**2
Dl = r**2 - 2 * M_n * r + a**2
f_S = 1 - 2 * M_n * r / Sig
b = 2 * M_n * a * r * sp.sin(th)**2 / Sig
G = (r**2 + a**2 + 2 * M_n * a**2 * r * sp.sin(th)**2 / Sig) * sp.sin(th)**2

vbS2 = 1 - A**2 * f_S / Eh_n**2
Gb = G + A**2 * b**2 / Eh_n**2
php0 = b * vbS2 / Gb
R2 = vbS2 * Dl * sp.sin(th)**2 / Gb
ptphi = J - A**2 * b / Eh_n
H3 = ptphi * php0 + sp.sqrt(R2) * sp.sqrt((Dl / Sig) * pr**2
                                          + pth**2 / Sig + ptphi**2 / Gb) \
    - A**2 * f_S / Eh_n

print("=" * 72)
print("[V1] key identity and equatorial reduction")
print("=" * 72)
idk = sp.simplify(f_S * G + b**2 - Dl * sp.sin(th)**2)
print("  f_S G + b^2 - Dl sin^2th =", idk)
# equatorial reduction: comparison with H_tau of R7 (P, Pb, etc.)
P_eq = r**2 + a**2 + 2 * M_n * a**2 / r
vb2_eq = 1 - A**2 * (1 - 2 * M_n / r) / Eh_n**2
Pb_eq = P_eq + A**2 * (2 * M_n * a / r)**2 / Eh_n**2
php0_eq = (2 * M_n * a / r) * vb2_eq / Pb_eq
R2_eq = Dl * vb2_eq / Pb_eq
ptp_eq = J - 2 * M_n * a * A**2 / (r * Eh_n)
H_eq = ptp_eq * php0_eq + sp.sqrt(R2_eq) * sp.sqrt(
    (Dl / r**2) * pr**2 + ptp_eq**2 / Pb_eq) - A**2 * (1 - 2 * M_n / r) / Eh_n
diff_eq = (H3.subs([(th, sp.pi / 2), (pth, 0)]) - H_eq)
rng = np.random.default_rng(3)
w_ = 0.0
for _ in range(8):
    sub = {r: rng.uniform(3, 10), a: rng.uniform(0.1, 0.95),
           A: rng.uniform(0.6, 1.3), J: rng.uniform(0.4, 2.0),
           pr: rng.uniform(-2, 2)}
    if float(vb2_eq.subs(sub)) <= 0:
        continue
    w_ = max(w_, abs(float(diff_eq.subs(sub))))
print(f"  max |H3(equator) - H_eq(R7)| over random points = {w_:.2e}")

# ------------------------------------------------ lambdify the dynamical flow
def make_flow(a_v, Hc_v):
    A_din = -1 / (Hc_v * ee)
    Hd = H3.subs([(A, A_din), (a, a_v)])
    args = (ee, r, th, pr, pth, J)
    return {q: sp.lambdify(args, sp.diff(Hd, v_), 'numpy')
            for q, v_ in (('dp_r', pr), ('dp_th', pth),
                          ('dr', r), ('dth', th))}

def integra(a_v, Hc_v, y0, Jv, eta0, eta1):
    F = make_flow(a_v, Hc_v)

    def rhs(e_, y):
        return [F['dp_r'](e_, *y, Jv), F['dp_th'](e_, *y, Jv),
                -F['dr'](e_, *y, Jv), -F['dth'](e_, *y, Jv)]
    ev = lambda e_, y: y[0] - 2.6
    ev.terminal, ev.direction = True, -1
    ev2 = lambda e_, y: y[0] - 25.0
    ev2.terminal, ev2.direction = True, 1
    s = solve_ivp(rhs, [eta0, eta1], y0, rtol=1e-11, atol=1e-13,
                  method='DOP853', events=[ev, ev2], dense_output=True)
    return s

print()
print("=" * 72)
print("[V2] a=0, dynamical A(eta): L^2 exact (spatial Noether theorem)")
print("=" * 72)
Hc_v = 0.005
eta0 = -1.0 / Hc_v                      # A(eta0) = 1
eta1 = eta0 / 1.25                      # A -> 1.25
s = integra(1e-12, Hc_v, [8.0, 1.1, -0.2, 0.8], 1.2, eta0, eta1)
lg = np.linspace(eta0, s.t[-1], 300)
Y = s.sol(lg)
L2 = Y[3]**2 + 1.2**2 / np.sin(Y[1])**2
print(f"  relative drift of L^2 over Delta_eta = {s.t[-1]-eta0:.0f}:"
      f"  {np.std(L2)/np.mean(L2):.2e}   (expected ~1e-12)")

print()
print("=" * 72)
print("[V3] dynamical hierarchy of the candidates (a=0.4)")
print("=" * 72)
a_v = 0.4
D0s = r**2 - 2 * M_n * r
Es = sp.Symbol('E_s', positive=True)
DE0s = (Es**2 - 1) * r**2 + 2 * M_n * r
N2s = (Es**2 - 1) * r**2 + 4 * M_n * r - 4 * M_n**2
f2s = N2s * sp.cos(2 * th) / (2 * r**2 * D0s * DE0s)
Q_std_s = pth**2 + J**2 / sp.tan(th)**2 - a**2 * Es**2 * sp.cos(th)**2
K_NLO_s = Q_std_s + a**2 * f2s * pth**2
Qf = sp.lambdify((r, th, pth, J, Es, a), Q_std_s, 'numpy')
Kf = sp.lambdify((r, th, pth, J, Es, a), K_NLO_s, 'numpy')

def drift_dinamico(Hc_v, n_traj=10):
    eta0 = -1.0 / Hc_v
    eta1 = eta0 / 1.25
    rng3 = np.random.default_rng(21)
    out = {k: [] for k in ('L2', 'Q_run', 'K_run', 'K_froz')}
    for _ in range(n_traj):
        y0 = [rng3.uniform(5, 9), rng3.uniform(0.8, 2.0),
              -0.2, rng3.uniform(0.5, 1.5)]
        Jv = rng3.uniform(0.5, 1.8)
        s = integra(a_v, Hc_v, y0, Jv, eta0, eta1)
        if s.t[-1] - eta0 < 0.3 * (eta1 - eta0):
            continue
        lg = np.linspace(eta0, s.t[-1], 400)
        Y = s.sol(lg)
        Aq = -1 / (Hc_v * lg)
        Eq = Eh_n / Aq
        sd = lambda q: np.std(q) / max(abs(np.mean(q)), 1e-14)
        out['L2'].append(sd(Y[3]**2 + Jv**2 / np.sin(Y[1])**2))
        out['Q_run'].append(sd(Qf(Y[0], Y[1], Y[3], Jv, Eq, a_v)))
        out['K_run'].append(sd(Kf(Y[0], Y[1], Y[3], Jv, Eq, a_v)))
        out['K_froz'].append(sd(Kf(Y[0], Y[1], Y[3], Jv, Eh_n, a_v)))
    return {k: np.median(v) for k, v in out.items()}

for Hc_v in (0.005, 0.0025):
    d = drift_dinamico(Hc_v)
    print(f"  H_c = {Hc_v}:  median drift:")
    print(f"    L^2 (a=0 reference)     : {d['L2']:.3e}")
    print(f"    Q_std(E_eff(eta))       : {d['Q_run']:.3e}")
    print(f"    K_NLO(E_eff(eta))       : {d['K_run']:.3e}")
    print(f"    K_NLO(E frozen)         : {d['K_froz']:.3e}")

print()
print("=" * 72)
print("[V4] source of the counterterm h (transport equation)")
print("=" * 72)
# dK/deta|_expl = (dK/dE)*Edot,  Edot = -E*(A'/A) => S1 = -E dK/dE / a^2
S1 = -Es * sp.diff(K_NLO_s, Es) / a**2
# compact form (using DE0 - N2 = -2M D0/r  =>  d_E f2 = -2ME cos2th/(r DE0^2)):
S1c = 2 * Es**2 * (sp.cos(th)**2
                   + M_n * sp.cos(2 * th) * pth**2 / (r * DE0s**2))
print("  dK/deta|_expl = eps * a^2 * S1,   eps = A'/A")
print("  S1 = 2 E_eff^2 [ cos^2(th) + M cos(2th) p_th^2 / (r DE0^2) ]")
print("  check of the compact form:", sp.simplify(S1 - S1c))
print("  NOTE: DE0 = r^2 w_eff  =>  DOUBLE POLE on the freezing wall")
print("\n  TRANSPORT EQUATION for the counterterm (order a^2 eps):")
print("     {h, H_0(Thakurta-Schw 3D)} = -S1(r, theta, p_theta)")
print("  with H_0 = vbar sqrt(f) sqrt(f p_r^2 + (p_th^2+J^2/sin^2th)/r^2)")
print("           - A^2 f/Ehat  (non-autonomous: A frozen at leading order).")
print("  h is the first object of the programme WITHOUT an analogue in Kerr.")
print("\nDONE.")
