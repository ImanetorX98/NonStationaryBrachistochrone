# -*- coding: utf-8 -*-
"""
Verifica che le curve trovate SONO le brachistocrone: perturbazione a
estremi fissi e controllo che il tempo di percorrenza sia il minimo,
RISPETTIVAMENTE (t minimizza il tempo coordinato, tau il tempo proprio).

Statico Schwarzschild equatoriale (M=1), rotaia con invariante E:
  velocita' locale  v = sqrt(1 - f/E^2),  f = 1 - 2M/r
  elemento di lunghezza propria  dl = sqrt(dr^2/f + r^2 dphi^2)
  tempo coordinato   dt   = n_t  dl ,  n_t   = E / sqrt(f (E^2 - f))
  tempo proprio      dtau = n_tau dl ,  n_tau = sqrt(f / (E^2 - f))

Ogni ramo e' la geodetica del PROPRIO indice ottico (Fermat). La curva
si costruisce dalla ODE di Beltrami  dr/dphi = (r/J) sqrt(f (n^2 r^2 - J^2))
con J = n(r_min) r_min (INDIPENDENTE dai funzionali di tempo qui sotto).
Poi si perturba r(phi) a estremi fissi e si valutano T_t, T_tau per
quadratura diretta.

Atteso:  T_tau minimo (eps=0) sulla curva tau ma NON sulla curva t;
         T_t   minimo (eps=0) sulla curva t   ma NON sulla curva tau.
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
    """Curva brachistocrona (mezzo ramo) dalla ODE di Beltrami."""
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
    # ramo completo simmetrico phi in [-Phi, Phi]
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
    bump = 1 - (phi / Phi)**2                       # 0 agli estremi
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
    print(f"    argmin T_t   a eps = {eps[np.argmin(Tt)]:+.4f}"
          f"   (curva {nome})")
    print(f"    argmin T_tau a eps = {eps[np.argmin(Ttau)]:+.4f}"
          f"   (curva {nome})")
    return dict(nome=nome, phi=phi, r=r, eps=eps, Tt=Tt, Ttau=Ttau, i0=i0,
                r_min=r_min)

print("=" * 66)
print("costruzione curve dalle ODE + test di minimo per perturbazione")
print("=" * 66)
Ctau = test_ramo('tau', n_tau, 3.0, 7.0)
Ct = test_ramo('t', n_t, 3.0, 7.0)

# --------------------------------------------------------------- figura
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
print("\nletto atteso: sulla curva tau -> min di T_tau a eps=0 (T_t no);")
print("              sulla curva t   -> min di T_t   a eps=0 (T_tau no).")
print('FATTO.')
