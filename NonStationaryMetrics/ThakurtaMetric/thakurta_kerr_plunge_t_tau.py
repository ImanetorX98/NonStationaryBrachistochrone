# -*- coding: utf-8 -*-
"""
Same-launch depth comparison t vs tau in Thakurta-Kerr (equatorial), protocol (P1) of Paper II:
same launch event, same numerical axial costate J, each radial momentum on ITS OWN branch shell
(H_br = 0), free arrival.  The quantity compared is the outer turning radius r_min reached inward
from r0; if a branch has no turning point in (2M, r0) it REACHES the conformal stationary limit
r = 2M, which is the end of the control domain: a limiting contact (type C), outside the
normality lemma, not an optimal crossing.

r_min of each branch is the outer root of the PMP turning condition H_br(r, p_r = 0; J) = 0 below
r0 (same functions as inversione_conformale_AJ.py), checked against the PMP Hamiltonian flow
launched on the shell (A = 1, where the branch Hamiltonians are those of fig_phi_validation_corrected.py).

An earlier version launched both QUADRATIC optical Hamiltonians with the same radial momentum
p_r = -0.5: that puts the two branches on DIFFERENT energy levels, and produced a spurious
"rotational inversion" in the (a, J) plane.  A later version used the cost -1 (H_eta) for the t-branch
at every A, i.e. compared eta with tau at J_eta = J_tau (J_t = A J_tau).  With the declared protocol,
H_t at J_t = J_tau, no inversion is found on the scanned grids, rotational or conformal: the t-branch
turns deeper or reaches r = 2M.

Panels: (a) Delta r = r_min^t - r_min^tau on (J, a) at A = 1, E = 1.2, r0 = 10, with the region
where the t-branch reaches r = 2M marked; (b) Delta r along J at a = 0.9 for A = 1, 1.3, 1.4, launched
at r0 = 20 or just inside the freezing surface.
"""
import os
import sys
import numpy as np
import sympy as sp
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp
from scipy.optimize import brentq

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(HERE))
from paper_style import COL, set_style, savefig

set_style()
OUT = os.path.join(HERE, 'Thakurtafigures')
M, Eh = 1.0, 1.2
R_E = 2 * M


def geom(r, A, a):
    f = 1 - 2 * M / r
    Dl = r**2 - 2 * M * r + a**2
    P = r**2 + a**2 + 2 * M * a**2 / r
    b = 2 * M * a / r
    vb2 = 1 - A**2 * f / Eh**2
    Pb = P + A**2 * b**2 / Eh**2
    return f, Dl, b, vb2, Pb


def H0(r, A, J, a, br):
    """PMP branch Hamiltonian at p_r = 0.  The coordinate-time branch is H_t = h - A (paper:
    H_eta = h - 1, H_t = h - A, h the support function); at A = 1 the two coincide.  An earlier
    version used the cost -1 (H_eta) for all A, which compares eta with tau at J_eta = J_tau,
    i.e. J_t = A J_tau, not the declared J_t = J_tau."""
    f, Dl, b, vb2, Pb = geom(r, A, a)
    if vb2 <= 0:
        return np.nan
    php0 = b * vb2 / Pb
    S = np.sqrt(vb2 * Dl) / Pb
    if br == 'tau':
        ptp = J - A**2 * b / Eh
        return ptp * php0 + S * abs(ptp) - A**2 * f / Eh
    return J * php0 + S * abs(J) - A


def r_min(J, a, A, r0, br):
    """Outer turning radius below r0; R_E if the branch reaches r = 2M; nan if r0 is off-shell."""
    if not H0(r0, A, J, a, br) <= 0:
        return np.nan                       # no real inward momentum at launch
    rs = np.linspace(r0, R_E + 1e-7, 6000)
    v = np.array([H0(x, A, J, a, br) for x in rs])
    for i in range(len(rs) - 1):
        if np.isfinite(v[i]) and np.isfinite(v[i + 1]) and v[i] <= 0 < v[i + 1]:
            return brentq(lambda x: H0(x, A, J, a, br), rs[i + 1], rs[i])
    return R_E


# ---- independent check at A = 1: PMP Hamiltonian flow launched on the shell ----
rr, pr, Js = sp.symbols('r p_r J')
def H_full(a, br):
    f = 1 - 2 * M / rr; Dl = rr**2 - 2 * M * rr + a**2; b = 2 * M * a / rr
    v = 1 - f / Eh**2; P = rr**2 + a**2 + 2 * M * a**2 / rr; Pb = P + b**2 / Eh**2
    if br == 'tau':
        pt = Js - b / Eh
        return pt * b * v / Pb + sp.sqrt(Dl * v / Pb) * sp.sqrt((Dl / rr**2) * pr**2 + pt**2 / Pb) - f / Eh
    return Js * b * v / Pb + sp.sqrt(Dl * v / Pb) * sp.sqrt((Dl / rr**2) * pr**2 + Js**2 / Pb) - 1


def r_min_flow(J, a, r0, br):
    H = H_full(a, br)
    Hn = sp.lambdify((rr, pr, Js), H, 'numpy')
    Hp = sp.lambdify((rr, pr, Js), sp.diff(H, pr), 'numpy')
    Hr = sp.lambdify((rr, pr, Js), sp.diff(H, rr), 'numpy')
    p0 = brentq(lambda p: Hn(r0, p, J), -60.0, 0.0)          # ingoing root (H even in p_r)
    ev_t = lambda s, y: y[1]; ev_t.terminal, ev_t.direction = True, 1
    ev_e = lambda s, y: y[0] - (R_E + 1e-4); ev_e.terminal, ev_e.direction = True, -1
    so = solve_ivp(lambda s, y: [Hp(y[0], y[1], J), -Hr(y[0], y[1], J)], [0, 5000], [r0, p0],
                   rtol=1e-11, atol=1e-13, events=[ev_t, ev_e])
    if len(so.t_events[0]):
        return so.y_events[0][0][0]
    return R_E


print("=" * 72)
print("check: turning radius from H(r,0)=0 vs PMP flow launched on the shell (A=1, r0=10)")
worst = 0.0
for a in (0.2, 0.5, 0.9):
    for J in (0.5, 2.5, 4.5, -2.0):
        for br in ('t', 'tau'):
            ra, rf = r_min(J, a, 1.0, 10.0, br), r_min_flow(J, a, 10.0, br)
            if np.isfinite(ra):
                worst = max(worst, abs(ra - rf) if ra > R_E else abs(rf - R_E) * (rf > R_E + 1e-3))
print(f"  max |r_min(turning condition) - r_min(flow)| = {worst:.1e}")
assert worst < 1e-5, "turning-condition and flow disagree"

# ---- (a) map on (J, a) at A = 1 ----
J_g = np.linspace(-6.0, 6.0, 97)
a_g = np.linspace(0.05, 0.98, 48)
D = np.full((len(a_g), len(J_g)), np.nan)
reach = np.zeros_like(D, dtype=bool)
for i, a in enumerate(a_g):
    for j, J in enumerate(J_g):
        rt, ru = r_min(J, a, 1.0, 10.0, 't'), r_min(J, a, 1.0, 10.0, 'tau')
        if np.isfinite(rt) and np.isfinite(ru):
            D[i, j] = rt - ru
            reach[i, j] = rt <= R_E
n_tau = int(np.nansum(D > 1e-9)); n_tot = int(np.sum(np.isfinite(D)))
print(f"(a) A=1: tau deeper in {n_tau}/{n_tot} cells; t reaches r=2M in {int(reach.sum())}")

# ---- (b) Delta r along J at a = 0.9; launch at r0 = 20 or just inside the freezing
#      surface A^2 f = Ehat^2 when that lies closer (r_w = 2M/(1 - Ehat^2/A^2)) ----
Jb = np.linspace(0.5, 30.0, 300)
curves = {}
A_b = (1.0, 1.3, 1.4)
def r0_of(A):
    return 20.0 if A <= Eh else min(20.0, 0.999 * 2 * M / (1 - Eh**2 / A**2))
for A in A_b:
    d = []
    r0b = r0_of(A)
    for J in Jb:
        rt, ru = r_min(J, 0.9, A, r0b, 't'), r_min(J, 0.9, A, r0b, 'tau')
        d.append(rt - ru if np.isfinite(rt) and np.isfinite(ru) else np.nan)
    curves[A] = np.array(d)
    pos = Jb[np.array(curves[A]) > 1e-9]
    print(f"(b) A={A} (r0={r0b:.3f}): tau deeper for J in "
          + (f"[{pos.min():.2f}, {pos.max():.2f}]" if len(pos) else "none"))

fig, (axa, axb) = plt.subplots(2, 1, figsize=(COL, 6.6))
vmax = np.nanmax(np.abs(D))
im = axa.pcolormesh(J_g, a_g, D, cmap='RdBu_r', vmin=-vmax, vmax=vmax, shading='auto')
plt.colorbar(im, ax=axa, label=r'$\Delta r = r_{\min}^{t}-r_{\min}^{\tau}$')
axa.contourf(J_g, a_g, reach.astype(float), levels=[0.5, 1.5], colors='none', hatches=['///'])
axa.plot([], [], ' ', label=r'hatched: $t$ reaches $r=2M$ (limiting contact)')
axa.set_xlabel(r'$J$ (axial costate)'); axa.set_ylabel(r'$a$')
axa.set_title(f'same launch, each branch on its own shell ($A=1$, $E={Eh}$, $r_0=10$):\n'
              f'$\\tau$ deeper in {n_tau} of {n_tot} cells', fontsize=7)
axa.legend(loc='lower right', fontsize=6)
for A, st in zip(A_b, ('k-', 'C0--', 'C3-.')):
    axb.plot(Jb, curves[A], st, lw=1.2, label=f'$A={A}$')
axb.axhline(0, color='0.5', lw=0.6)
axb.set_xlabel(r'$J$'); axb.set_ylabel(r'$\Delta r$')
axb.set_title(r'$a=0.9$, $t$ vs $\tau$ at $J_t=J_\tau$, launch inside the freezing surface', fontsize=7)
axb.legend(fontsize=6); axb.grid(alpha=0.3)
savefig(fig, OUT, 'fig_thakurta_kerr_plunge_t_tau')
print('\nFATTO.')
