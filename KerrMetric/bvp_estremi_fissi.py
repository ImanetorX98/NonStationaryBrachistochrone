# -*- coding: utf-8 -*-
"""
Brachistocrone t e tau a ESTREMI FISSI (problema ai limiti, BVP).

Fissati due punti A e B nel piano equatoriale (Schwarzschild statico,
ottica Riemanniana), per ciascun ramo si cerca la geodetica del PROPRIO
indice ottico che connette A e B (fastest-time curve), e se ne estrae il
raggio minimo r_min. t e tau connettono gli STESSI due punti ma sono
curve diverse -> confronto pulito di profondita' (chi affonda di piu').

Indici ottici (rotaia con invariante E, f=1-2M/r):
  n_tau = sqrt(f/(E^2-f)) ,   n_t = E/sqrt(f(E^2-f))
metrica spaziale base dl^2 = dr^2/f + r^2 dphi^2.
Geodetica (Beltrami/Clairaut):  N r^2/sqrt(r^2+r'^2/f) = J  (r'=dr/dphi)
=> dphi/dr = J /( r sqrt(f) sqrt(N^2 r^2 - J^2) ) ,  r_min: N(r_min) r_min = J.

Casi:
  (a) SIMMETRICO:  A=(r0,-Phi), B=(r0,+Phi)  (stesso raggio)
  (b) ASIMMETRICO: A=(rA,0),    B=(rB,dphi)  (raggi diversi)
Per ciascun ramo si spara J per far combaciare Delta_phi = phi_B - phi_A.
"""

import os
import sys
import numpy as np
from scipy.integrate import quad
from scipy.optimize import brentq

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from paper_style import COL, set_style, savefig
import matplotlib.pyplot as plt

set_style()
HERE = os.path.dirname(os.path.abspath(__file__))
M, E = 1.0, 1.2

def f(r):
    return 1 - 2 * M / r

def n_tau(r):
    return np.sqrt(f(r) / (E**2 - f(r)))

def n_t(r):
    return E / np.sqrt(f(r) * (E**2 - f(r)))

def rmin_of_J(N, J, r_hi):
    """r_min = radice ESTERNA di N(r) r = J sotto r_hi (scan verso l'interno)."""
    rs = np.linspace(r_hi * (1 - 1e-6), 2.0001 * M, 5000)
    g = N(rs) * rs - J
    for i in range(len(rs) - 1):
        if g[i] * g[i + 1] < 0:
            return brentq(lambda r: N(r) * r - J, rs[i + 1], rs[i])
    return None

def dphi_seg(N, J, r_lo, r_hi):
    """angolo da r_lo(>r_min) a r_hi lungo la geodetica."""
    integ = lambda r: J / (r * np.sqrt(f(r)) * np.sqrt(N(r)**2 * r**2 - J**2))
    return quad(integ, r_lo, r_hi, limit=200, points=[r_lo])[0]

def solve_bvp(N, rA, rB, dphi_target):
    """spara J: Delta_phi(A->rmin->B) = dphi_target. Ritorna J, r_min."""
    r_hi = min(rA, rB)
    def gap(J):
        rm = rmin_of_J(N, J, r_hi)
        if rm is None:
            return np.nan
        return dphi_seg(N, J, rm + 1e-9, rA) + dphi_seg(N, J, rm + 1e-9, rB) \
            - dphi_target
    Jmax = N(r_hi) * r_hi * (1 - 1e-6)
    Js = np.linspace(0.2 * Jmax, Jmax, 500)
    gv = np.array([gap(J) for J in Js])
    for i in range(len(Js) - 1):
        if np.isfinite(gv[i]) and np.isfinite(gv[i + 1]) \
                and gv[i] * gv[i + 1] < 0:
            J = brentq(gap, Js[i], Js[i + 1])
            return J, rmin_of_J(N, J, r_hi)
    return None, None

def curva(N, J, rA, rB, phiA):
    """(r,phi) della geodetica da A(rA,phiA) a B, via periasse r_min."""
    rm = rmin_of_J(N, J, min(rA, rB))
    r_in = np.linspace(rA, rm + 1e-7, 300)
    phi_in = phiA + np.array([dphi_seg(N, J, rm + 1e-9, rA)
                              - dphi_seg(N, J, rm + 1e-9, rr)
                              for rr in r_in])
    phi_peri = phiA + dphi_seg(N, J, rm + 1e-9, rA)
    r_out = np.linspace(rm + 1e-7, rB, 300)
    phi_out = phi_peri + np.array([dphi_seg(N, J, rm + 1e-9, rr)
                                   for rr in r_out])
    return (np.concatenate([r_in, r_out]),
            np.concatenate([phi_in, phi_out]), rm)

print("=" * 66)
print("Schwarzschild statico, M=1, E=1.2: brachistocrone t/tau a estremi fissi")
print("=" * 66)

# ---- (a) simmetrico ----
r0, Phi = 6.0, 0.9
print(f"\n(a) SIMMETRICO: A=({r0},-{Phi}), B=({r0},+{Phi})")
res_a = {}
for nome, N in (('tau', n_tau), ('t', n_t)):
    J, rm = solve_bvp(N, r0, r0, 2 * Phi)
    res_a[nome] = (J, rm)
    print(f"    ramo {nome}: J_opt={J:.4f}, r_min={rm:.4f}")
print(f"    => r_min^t - r_min^tau = "
      f"{res_a['t'][1]-res_a['tau'][1]:+.4f} (t piu' fondo se <0)")

# ---- (b) asimmetrico ----
rA, rB, dphi = 10.0, 6.0, 1.9
print(f"\n(b) ASIMMETRICO: A=({rA},0), B=({rB},{dphi})")
res_b = {}
for nome, N in (('tau', n_tau), ('t', n_t)):
    J, rm = solve_bvp(N, rA, rB, dphi)
    res_b[nome] = (J, rm)
    print(f"    ramo {nome}: J_opt={J:.4f}, r_min={rm:.4f}")
print(f"    => r_min^t - r_min^tau = "
      f"{res_b['t'][1]-res_b['tau'][1]:+.4f}")

# --------------------------------------------------------------- figura
fig, (a1, a2) = plt.subplots(2, 1, figsize=(COL, 6.4))
th = np.linspace(0, 2 * np.pi, 200)
for ax in (a1, a2):
    ax.plot(2 * M * np.cos(th), 2 * M * np.sin(th), 'k-', lw=0.8)  # r=2M
    ax.set_aspect('equal')
    ax.set_xlabel('$x$')
    ax.set_ylabel('$y$')

# (a) simmetrico
for nome, col in (('tau', 'C3'), ('t', 'C0')):
    J, rm = res_a[nome]
    r, ph, _ = curva(n_tau if nome == 'tau' else n_t, J, r0, r0, -Phi)
    a1.plot(r * np.cos(ph), r * np.sin(ph), col, lw=1.6,
            label=rf'${nome}$: $r_{{\min}}={rm:.3f}$')
a1.plot(r0 * np.cos(-Phi), r0 * np.sin(-Phi), 'ks', ms=6)
a1.plot(r0 * np.cos(Phi), r0 * np.sin(Phi), 'ks', ms=6, label='fixed A, B')
a1.set_title(f'(a) symmetric endpoints (same $r_0={r0}$, $\\pm\\Phi$):\n'
             '$t$ and $\\tau$ join the SAME A,B, differ in depth')
a1.legend(fontsize=6, loc='lower left')

# (b) asimmetrico
for nome, col in (('tau', 'C3'), ('t', 'C0')):
    J, rm = res_b[nome]
    r, ph, _ = curva(n_tau if nome == 'tau' else n_t, J, rA, rB, 0.0)
    a2.plot(r * np.cos(ph), r * np.sin(ph), col, lw=1.6,
            label=rf'${nome}$: $r_{{\min}}={rm:.3f}$')
a2.plot(rA, 0, 'ks', ms=6)
a2.plot(rB * np.cos(dphi), rB * np.sin(dphi), 'ks', ms=6, label='fixed A, B')
a2.set_title(f'(b) asymmetric endpoints ($r_A={rA}$, $r_B={rB}$):\n'
             'genuine two-point brachistochrones')
a2.legend(fontsize=6, loc='lower left')
savefig(fig, HERE, 'fig_bvp_estremi_fissi')
print("\nFATTO.")
