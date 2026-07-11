# -*- coding: utf-8 -*-
"""
Il FATTORE CONFORME A genera inversione a ESTREMI FISSI (simmetrici)?

Transfer theorem: Thakurta-Schwarzschild(A, Ehat) == Schwarzschild con
E_eff = Ehat/A (r coordinata invariante). Quindi variare A a estremi
fissi <=> variare l'energia di rotaia E nel BVP di Schwarzschild.
Studio Delta_r(E) = r_min^t - r_min^tau; se cambia segno, A INVERTE.

Ottica (Beltrami, niente ODE): f=1-2M/r,
  n_tau = sqrt(f/(E^2-f)),   n_t = E/sqrt(f(E^2-f))
  dl^2 = dr^2/f + r^2 dphi^2,   J = N r^2/sqrt(r^2+r'^2/f) conservato
  => r_min: N(r_min) r_min = J ,  Phi = int_{rmin}^{r0} J dr/(r sqrt(f)
                                        sqrt(N^2 r^2 - J^2)).

Analitico: r_min = radice esterna di N^2 r^2 = J^2:
  tau: f r^2/(E^2-f) = J^2 ;   t: E^2 r^2/(f(E^2-f)) = J^2.
J e' fissato dalla BVP ellittica. Delta_r(E) e la sua eventuale radice.
"""

import os
import sys
import numpy as np
from scipy.integrate import quad
from scipy.optimize import brentq

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from paper_style import COL, set_style, savefig
import matplotlib.pyplot as plt
try:
    from tqdm import tqdm
except ImportError:
    def tqdm(x, **k):
        return x

set_style()
HERE = os.path.dirname(os.path.abspath(__file__))
M, r0, Phi = 1.0, 6.0, 0.9

def f(r):
    return 1 - 2 * M / r

def n_tau(r, E):
    return np.sqrt(f(r) / (E**2 - f(r)))

def n_t(r, E):
    return E / np.sqrt(f(r) * (E**2 - f(r)))

def rmin_of_J(N, E, J, r_hi):
    rs = np.linspace(r_hi * (1 - 1e-6), 2.0001 * M, 6000)
    g = N(rs, E) * rs - J
    for i in range(len(rs) - 1):
        if g[i] * g[i + 1] < 0:
            return brentq(lambda r: N(r, E) * r - J, rs[i + 1], rs[i])
    return None

def dphi_half(N, E, J, rm):
    integ = lambda r: J / (r * np.sqrt(f(r))
                           * np.sqrt(N(r, E)**2 * r**2 - J**2))
    return quad(integ, rm + 1e-10, r0, limit=200)[0]

def rmin_bvp(N, E):
    """spara J: dphi_half = Phi. Ritorna r_min shallow (brachistocrona min)."""
    Jmax = N(r0, E) * r0 * (1 - 1e-7)
    Js = np.linspace(0.02 * Jmax, Jmax, 600)
    sols = []
    prev = None
    for J in Js:
        rm = rmin_of_J(N, E, J, r0)
        if rm is None:
            prev = None
            continue
        d = dphi_half(N, E, J, rm) - Phi
        if prev is not None and prev[1] * d < 0:
            try:
                Jstar = brentq(lambda x: dphi_half(N, E, x,
                               rmin_of_J(N, E, x, r0)) - Phi, prev[0], J)
                sols.append(rmin_of_J(N, E, Jstar, r0))
            except (ValueError, TypeError):
                pass
        prev = (J, d)
    return max(sols) if sols else None      # shallow = brachistocrona minima

print("=" * 60)
print(f"Conforme a estremi fissi simmetrici (Schw, r0={r0}, +/-{Phi})")
print(f"E_eff = Ehat/A  =>  variare E == variare A (inversa)")
print("=" * 60)
Egrid = np.linspace(1.03, 2.6, 24)
res = []
print(f"{'E':>6} {'A(Eh=1.2)':>9} {'r_min^tau':>10} {'r_min^t':>9} "
      f"{'Delta_r':>9}", flush=True)
for E in tqdm(Egrid, desc='scan E (=1/A)', ncols=70):
    rt = rmin_bvp(n_tau, E)
    rr = rmin_bvp(n_t, E)
    if rt is None or rr is None:
        continue
    res.append((E, rt, rr, rr - rt))
    print(f"{E:6.2f} {1.2/E:9.3f} {rt:10.4f} {rr:9.4f} {rr-rt:+9.4f}",
          flush=True)

res = np.array(res)
print("\n--- tabella completa (E, A, r_tau, r_t, Delta_r) ---", flush=True)
for row in res:
    print(f"  E={row[0]:.3f}  A={1.2/row[0]:.3f}  r_tau={row[1]:.4f}  "
          f"r_t={row[2]:.4f}  Dr={row[3]:+.4f}", flush=True)
print(f"\nDelta_r range: [{res[:,3].min():+.4f}, {res[:,3].max():+.4f}]")
if res[:, 3].min() * res[:, 3].max() < 0:
    print(">>> INVERSIONE: Delta_r cambia segno!")
    i = np.where(np.diff(np.sign(res[:, 3])))[0][0]
    E_inv = np.interp(0, [res[i, 3], res[i+1, 3]], [res[i, 0], res[i+1, 0]])
    print(f"    E_inv ~ {E_inv:.3f}  (A_inv = {1.2/E_inv:.3f} per Ehat=1.2)")
else:
    print(">>> NESSUNA inversione: Delta_r sempre dello stesso segno.")

# --------------------------------------------------------------- figura
fig, ax = plt.subplots(figsize=(COL, COL * 0.75))
ax.plot(res[:, 0], res[:, 3], 'C0o-', ms=3)
ax.axhline(0, color='k', lw=0.6)
ax.set_xlabel('$E_{eff}=\\hat E/A$  (conformal factor via transfer)')
ax.set_ylabel(r'$\Delta r = r_{\min}^{t}-r_{\min}^{\tau}$')
ax.set_title('Conformal inversion at fixed symmetric endpoints?\n'
             f'(Schwarzschild, $r_0={r0}$, $\\pm{Phi}$)')
savefig(fig, HERE, 'fig_bvp_conforme_inversione')
print("FATTO.")
