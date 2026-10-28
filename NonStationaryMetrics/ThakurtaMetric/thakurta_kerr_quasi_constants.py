# -*- coding: utf-8 -*-
"""
Quasi-constants in Thakurta-Kerr: adiabatic invariant of the expansion
of the Hamiltonian (conformal analogue of the Kerr quasi-constants programme).

In Kerr the programme looked for Carter-like quasi-constants by expanding
the Hamiltonian and testing the drift along the trajectories. In Thakurta-Kerr
the small parameter is ALREADY there: eps = A'/A per orbital period (conformal
drift). The expansion of the non-autonomous Hamiltonian

    K = I_r + eps K_1 + eps^2 K_2 + ...

has at first order the ADIABATIC INVARIANT of the radial action:

    I_r = (1/2pi) \oint p_r dr

Adiabatic theorem: for a slowly varying H(eta), dI_r/deta = O(eps^2)
while dH/deta = O(eps): I_r is the quasi-constant.

Setting: bound regime E_eff = Ehat/A < 1 (A > Ehat): the tau branch librates
in a BOX between the inner centrifugal turning point and the outer freezing
wall (vbar = 0); with a growing A(eta) = -1/(H_c eta) the
box shrinks (capture by expansion): H drifts by O(1), I_r stays
quasi-constant.

Check: dynamical PMP flow of H_tau (s=0.9, Ehat=1.2, J=1.5,
H_c = 2e-4, A: 1.30 -> squeeze); cycles delimited by the inner turning points
(p_r = 0, direction +1); I_k = int p_r dr over cycle k.
Figure: Thakurtafigures/fig_thakurta_kerr_quasicostante
"""

import os
import sys
import numpy as np
import matplotlib.pyplot as plt
import sympy as sp
from scipy.integrate import solve_ivp

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(HERE))
from paper_style import COL, set_style, savefig

set_style()
OUT = os.path.join(HERE, 'Thakurtafigures')

M_n, s_n, E_n, J_n = 1.0, 0.9, 1.2, 1.5
H_c = 2e-4
A0 = 1.30
eta0 = -1.0 / (H_c * A0)
r_start = 4.0
r_plus = M_n + np.sqrt(M_n**2 - s_n**2)

ee = sp.Symbol('ee', negative=True)
r, pr = sp.symbols('r p_r', positive=True), sp.Symbol('p_r', real=True)
r = sp.Symbol('r', positive=True)
A_e = -1 / (H_c * ee)
f = 1 - 2 * M_n / r
Dl = r**2 - 2 * M_n * r + s_n**2
P = r**2 + s_n**2 + 2 * M_n * s_n**2 / r
vb2 = 1 - A_e**2 * f / E_n**2
Pb = P + A_e**2 * (2 * M_n * s_n / r)**2 / E_n**2
php0 = (2 * M_n * s_n / r) * vb2 / Pb
R_ = sp.sqrt(Dl * vb2 / Pb)
ptp = J_n - 2 * M_n * s_n * A_e**2 / (r * E_n)
H_tau = ptp * php0 + R_ * sp.sqrt((Dl / r**2) * pr**2 + ptp**2 / Pb) \
    - A_e**2 * f / E_n

Hf = sp.lambdify((ee, r, pr), H_tau, 'numpy')
dHp = sp.lambdify((ee, r, pr), sp.diff(H_tau, pr), 'numpy')
dHr = sp.lambdify((ee, r, pr), sp.diff(H_tau, r), 'numpy')
vb2f = sp.lambdify((ee, r), vb2, 'numpy')

def a_dyn(e_):
    return -1.0 / (H_c * e_)

def r_freeze(e_):
    th = E_n**2 / a_dyn(e_)**2
    return 2 * M_n / (1 - th) if th < 1 else np.inf

def rhs(e_, y):
    dr_ = dHp(e_, y[0], y[1])
    return [dr_, -dHr(e_, y[0], y[1]), y[1] * dr_]     # S' = p_r r'

ev_in = lambda e_, y: y[1]
ev_in.terminal, ev_in.direction = False, 1              # inner turning points
ev_fr = lambda e_, y: vb2f(e_, y[0]) - 1e-6             # wall: reflect
ev_fr.terminal, ev_fr.direction = True, -1
ev_h = lambda e_, y: y[0] - r_plus * 1.01
ev_h.terminal, ev_h.direction = True, -1

# piecewise integration: at the freezing wall the indicatrix
# degenerates (soft turn): regularised reflection p_r -> -p_r
# (H is even in p_r: its value is preserved)
eta_end = -1.0 / (H_c * 1.75)
segs = []
t_in = []
y_cur = [r_start, 0.0, 0.0]
e_cur = eta0
esito = 'end of window'
for _ in range(60):
    s = solve_ivp(rhs, [e_cur, eta_end], y_cur, rtol=1e-11, atol=1e-13,
                  method='DOP853', events=[ev_in, ev_fr, ev_h],
                  dense_output=True, max_step=5.0)
    segs.append(s)
    t_in.extend(list(s.t_events[0]))
    if len(s.t_events[2]):
        esito = 'horizon'
        break
    if len(s.t_events[1]):                # contact with the wall: reflect
        e_cur = s.t_events[1][0]
        y_c = s.y_events[1][0].copy()
        y_c[1] = -y_c[1]
        y_cur = y_c
        if e_cur >= eta_end - 1.0:
            esito = 'final freezing'
            break
        continue
    break

def sol_glob(e_):
    for s_ in segs:
        if s_.t[0] <= e_ <= s_.t[-1]:
            return s_.sol(e_)
    return segs[-1].sol(e_)

t_in = np.array(sorted(t_in))
e_last = segs[-1].t[-1]
print("=" * 72)
print(f"[quasi-constant] libration cycles: {max(len(t_in) - 1, 0)}   "
      f"(A: {A0} -> {a_dyn(e_last):.3f}, outcome: {esito}, "
      f"segments: {len(segs)})")
print("=" * 72)
Iks, Hks, Aks = [], [], []
for k in range(len(t_in) - 1):
    S0 = sol_glob(t_in[k])[2]
    S1 = sol_glob(t_in[k + 1])[2]
    Ik = abs(S1 - S0) / (2 * np.pi)
    yk = sol_glob(t_in[k])
    Hk = Hf(t_in[k], yk[0], yk[1])
    Iks.append(Ik)
    Hks.append(Hk)
    Aks.append(a_dyn(t_in[k]))
    print(f"  cycle {k+1:2d}:  A = {Aks[-1]:.4f}   E_eff = "
          f"{E_n/Aks[-1]:.4f}   H = {Hk:+.6f}   I_r = {Ik:.8f}")
Iks = np.array(Iks)
Hks = np.array(Hks)
if len(Iks) > 2:
    dI = (Iks.max() - Iks.min()) / abs(Iks.mean())
    dH = (Hks.max() - Hks.min()) / max(abs(Hks).max(), 1e-12)
    print(f"\n  relative spread:  I_r = {dI:.2e}   H = {dH:.2e}")
else:
    # RESULT: universal capture at the wall -- no radial libration
    e_c = segs[0].t_events[1][0] if len(segs[0].t_events[1]) else e_last
    y_c = segs[0].y_events[1][0] if len(segs[0].t_events[1]) \
        else segs[0].y[:, -1]
    v_orb = abs(dHp(e_c, y_c[0], y_c[1]))
    Ac = a_dyn(e_c)
    th = E_n**2 / Ac**2
    drfz = abs((2 * M_n / (1 - th)**2) * (2 * E_n**2 / Ac**3)
               * H_c * Ac**2)
    print("\n  RESULT (negative but exact): equatorial radial libration")
    print("  does NOT exist in Thakurta-Kerr: at the freezing wall")
    print("  the orbital speed vanishes while the wall")
    print("  moves down at finite speed -- whatever touches it is captured.")
    print(f"  at the first contact (eta={e_c:.1f}, r={y_c[0]:.3f}):")
    print(f"    |dr/deta| orbit = {v_orb:.2e}")
    print(f"    |dr/deta| wall  = {drfz:.2e}   (ratio "
          f"{drfz/max(v_orb,1e-30):.0f}x)")
    print("  => the radial action I_r cannot be defined; the adiabatic")
    print("     quasi-constant of the programme must be sought in the POLAR motion")
    print("     (I_theta, off the equator) or in the conformal quasi-Carter constant.")

# --------------------------------------------------------------- figure
fig, (axa, axb) = plt.subplots(2, 1, figsize=(COL, 6.0))
eg = np.linspace(eta0, e_last, 4000)
axa.plot(eg, [sol_glob(e_)[0] for e_ in eg], 'C0-', lw=0.9,
         label='orbit $r(\\eta)$')
axa.plot(eg, [r_freeze(e_) for e_ in eg], 'C3--', lw=1.3,
         label='freezing wall')
for tk in t_in:
    axa.axvline(tk, color='gray', lw=0.3)
axa.set_xlabel('$\\eta$')
axa.set_ylabel('$r$')
axa.set_title(f'libration in the shrinking box\n($s={s_n}$, '
              f'$\\hat E={E_n}$, $J={J_n}$, $H_c={H_c}$, $A_0={A0}$)')
axa.legend()

if len(Iks) > 2:
    kk = np.arange(1, len(Iks) + 1)
    axb.plot(kk, Iks / Iks[0], 'C0o-', ms=4, label='$I_r^{(k)}/I_r^{(1)}$')
    axb.plot(kk, np.array(Hks) / Hks[0], 'C3s--', ms=4,
             label='$H^{(k)}/H^{(1)}$')
    axb.set_xlabel('cycle $k$')
    axb.legend()
else:
    e_c = segs[0].t[-1]
    eg2 = np.linspace(eta0 + 0.7 * (e_c - eta0), e_c, 800)
    axb.plot(eg2, [sol_glob(e_)[0] for e_ in eg2], 'C0-', lw=1.2,
             label='orbit')
    axb.plot(eg2, [r_freeze(e_) for e_ in eg2], 'C3--', lw=1.2,
             label='freezing wall')
    axb.set_xlabel('$\\eta$')
    axb.set_ylabel('$r$')
    axb.set_title('universal capture: $v_{orbit}\\to 0$ at the wall,\n'
                  'while the wall descends at finite speed')
    axb.legend()
axb.grid(alpha=0.3)
savefig(fig, OUT, 'fig_thakurta_kerr_quasicostante')
print('\nDONE.')
