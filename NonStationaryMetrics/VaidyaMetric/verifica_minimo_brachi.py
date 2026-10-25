# -*- coding: utf-8 -*-
"""
Check that the curves found ARE the brachistochrones: a perturbation with
fixed endpoints and a check that the travel time is the minimum,
RESPECTIVELY (t minimises the coordinate time, tau the proper time),
within the perturbation family tested.

Static equatorial Schwarzschild (M=1), rail with invariant E:
  local speed  v = sqrt(1 - f/E^2),  f = 1 - 2M/r
  proper length element  dl = sqrt(dr^2/f + r^2 dphi^2)
  coordinate time   dt   = n_t  dl ,  n_t   = E / sqrt(f (E^2 - f))
  proper time       dtau = n_tau dl ,  n_tau = sqrt(f / (E^2 - f))

Each branch is the geodesic of ITS OWN optical index (Fermat). The curve
is built from the Beltrami ODE  dr/dphi = (r/J) sqrt(f (n^2 r^2 - J^2))
with J = n(r_min) r_min (INDEPENDENT of the time functionals below).
Then r(phi) is perturbed with fixed endpoints and T_t, T_tau are evaluated by
direct quadrature.

Expected:  T_tau minimal (eps=0) on the tau curve but NOT on the t curve;
           T_t   minimal (eps=0) on the t curve   but NOT on the tau curve.
"""

import os
import sys
import numpy as np
from scipy.integrate import solve_ivp
from scipy.interpolate import CubicSpline

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from paper_style import COL, set_style, savefig
import matplotlib.pyplot as plt

set_style()
OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'Vaidyafigures')
M, E = 1.0, 1.2

def f(r):
    return 1 - 2 * M / r

def n_tau(r):
    return np.sqrt(f(r) / (E**2 - f(r)))

def n_t(r):
    return E / np.sqrt(f(r) * (E**2 - f(r)))

def costruisci(indice, r_min, r0):
    """Brachistochrone curve (half branch) from the Beltrami ODE."""
    J = indice(r_min) * r_min
    def drdphi(phi, y):
        r = y[0]
        val = f(r) * (indice(r)**2 * r**2 - J**2)
        return [(r / J) * np.sqrt(max(val, 0.0))]
    ev = lambda phi, y: y[0] - r0
    ev.terminal, ev.direction = True, 1
    s = solve_ivp(drdphi, [0, 10], [r_min + 1e-9], rtol=1e-11, atol=1e-13,
                  events=[ev], dense_output=True, max_step=1e-3)
    Phi = s.t_events[0][0]
    ph = np.linspace(0, Phi, 400)
    r_half = s.sol(ph)[0]
    # complete symmetric branch phi in [-Phi, Phi]
    phi_full = np.concatenate([-ph[::-1], ph[1:]])
    r_full = np.concatenate([r_half[::-1], r_half[1:]])
    return phi_full, r_full, Phi, J

def tempo(phi, r, indice):
    """T = int indice(r) sqrt(r^2 + (dr/dphi)^2/f) dphi."""
    rp = np.gradient(r, phi)
    integ = indice(r) * np.sqrt(r**2 + rp**2 / f(r))
    return np.trapezoid(integ, phi)

def test_ramo(nome, indice, r_min, r0):
    phi, r, Phi, J = costruisci(indice, r_min, r0)
    bump = 1 - (phi / Phi)**2                       # 0 at the endpoints
    eps = np.linspace(-0.30, 0.30, 61)
    Tt, Ttau = [], []
    for e in eps:
        rr = r + e * bump
        Tt.append(tempo(phi, rr, n_t))
        Ttau.append(tempo(phi, rr, n_tau))
    Tt, Ttau = np.array(Tt), np.array(Ttau)
    i0 = len(eps) // 2
    print(f"[{nome}]  r_min={r_min}, r0={r0}, span 2Phi={2*Phi:.3f} rad, "
          f"J={J:.4f}")
    print(f"    argmin T_t   at eps = {eps[np.argmin(Tt)]:+.4f}"
          f"   (curve {nome})")
    print(f"    argmin T_tau at eps = {eps[np.argmin(Ttau)]:+.4f}"
          f"   (curve {nome})")
    return dict(nome=nome, phi=phi, r=r, eps=eps, Tt=Tt, Ttau=Ttau, i0=i0,
                r_min=r_min)

print("=" * 66)
print("building the curves from the ODEs + minimum test by perturbation")
print("=" * 66)
Ctau = test_ramo('tau', n_tau, 3.0, 7.0)
Ct = test_ramo('t', n_t, 3.0, 7.0)

# --------------------------------------------------------------- figure
fig, (a1, a2) = plt.subplots(2, 1, figsize=(COL, 5.6))
for C, ax, titolo in ((Ctau, a1, r'$\tau$-brachistochrone'),
                      (Ct, a2, r'$t$-brachistochrone')):
    i0 = C['i0']
    ax.plot(C['eps'], C['Ttau'] / C['Ttau'][i0], 'C3-',
            label=r'$T_\tau(\varepsilon)/T_\tau(0)$ (proper time)')
    ax.plot(C['eps'], C['Tt'] / C['Tt'][i0], 'C0--',
            label=r'$T_t(\varepsilon)/T_t(0)$ (coord. time)')
    ax.axvline(0, color='k', lw=0.6)
    ax.set_ylabel('travel time (norm.)')
    ax.set_title(f'{titolo}: perturbing $r(\\varphi)$ at fixed endpoints')
    ax.legend()
a2.set_xlabel(r'perturbation amplitude $\varepsilon$ '
              r'($\delta r = \varepsilon(1-(\varphi/\Phi)^2)$)')
savefig(fig, OUT, 'fig_verifica_minimo_brachi')
print("\nexpected reading: on the tau curve -> minimum of T_tau at eps=0 (not T_t);")
print("                  on the t curve   -> minimum of T_t   at eps=0 (not T_tau).")
print('DONE.')
