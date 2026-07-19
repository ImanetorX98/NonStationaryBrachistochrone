# -*- coding: utf-8 -*-
"""
Famiglie "che respirano" (ordine adiabatico dominante WKB) per ENTRAMBI i rami
t e tau della brachistocrona Thakurta-Kerr.

Rail energy Ehat FISSO; fattore conforme A che cresce => E_eff = Ehat/A decresce.
Le forme chiuse Kerr(E_eff) sono l'ordine dominante WKB (famiglia congelata);
la traiettoria vera le infila a O(A'/A) (correzione quasi-costante, Sec.V).

Entrambi i rami respirano via LO STESSO E_eff (il ramo tau ha in piu' solo il
prefattore A^-2 sul timing, non sulla forma). Mostro anche la soglia
J_c = a/E_eff = a A/Ehat che CRESCE con A (la finestra di regimi si sposta).

Uso le forme chiuse d(phi)/dr (BL):
  tau: J r sqrt(wf)/(Δ sqrt(Δ - J^2 w))
  t:   K(r)/sqrt(R6),  K=r[(E^2-1)r+2M](J(r-2M)+2Ma)/Δ
integrando l'orbita di scattering (turning esterno) per una scala di A.
"""
import os
import sys
import numpy as np
from scipy.integrate import quad
from scipy.optimize import brentq

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from paper_style import COL, set_style, savefig
import matplotlib.pyplot as plt
from matplotlib import cm

set_style()
HERE = os.path.dirname(os.path.abspath(__file__))
M, a = 1.0, 0.9
Ehat = 1.4
r0 = 12.0


def f(r): return 1 - 2*M/r
def Dl(r): return r**2 - 2*M*r + a**2
def w(r, E): return E**2 - f(r)


def dphidr_tau(r, E, J):
    rad = Dl(r) - J**2*w(r, E)
    return J*r*np.sqrt(w(r, E)*f(r))/(Dl(r)*np.sqrt(rad))


def R6(r, E, J):
    Q2 = (2*E**2*J**2*M*r - E**2*J**2*r**2 - 4*E**2*J*M*a*r + 2*E**2*M*a**2*r
          + E**2*a**2*r**2 + E**2*r**4 + 4*J**2*M**2 - 4*J**2*M*r + J**2*r**2
          - 8*J*M**2*a + 4*J*M*a*r + 4*M**2*a**2)
    return r*Q2*((E**2 - 1)*r + 2*M)


def dphidr_t(r, E, J):
    K = r*((E**2-1)*r+2*M)*(J*(r-2*M)+2*M*a)/Dl(r)
    return K/np.sqrt(R6(r, E, J))


def rmin_tau(E, J):
    # turning esterno: Δ - J^2 w = 0 (se sopra r_e)
    try:
        return brentq(lambda r: Dl(r)-J**2*w(r, E), 2.0+1e-6, r0)
    except Exception:
        return np.nan


def rmin_t(E, J):
    rs = np.linspace(1.5, r0, 40000); v = R6(rs, E, J)
    idx = np.where(np.diff(np.sign(v)))[0]
    roots = [brentq(lambda r: R6(r, E, J), rs[i], rs[i+1]) for i in idx]
    roots = [r for r in roots if r > 2.0]
    return max(roots) if roots else np.nan


def orbit(dphidr, rmin, E, J):
    rst = rmin(E, J)
    if not np.isfinite(rst):
        return None
    rr = np.linspace(r0, rst+1e-4, 500)
    ph = np.array([quad(lambda x: dphidr(x, E, J), r, r0, limit=200)[0]
                   for r in rr])
    return (np.concatenate([rr, rr[::-1]]),
            np.concatenate([ph, 2*ph[-1]-ph[::-1]]), rst)


As = np.array([1.0, 1.1, 1.2, 1.3, 1.4])
Eeffs = Ehat/As
Jt, Jtau = 6.0, 2.5   # scattering per entrambi (|J|>J_c nel range)

fig, axes = plt.subplots(1, 2, figsize=(2*COL, COL*1.02))
th = np.linspace(0, 2*np.pi, 200)
colors = cm.viridis(np.linspace(0, 0.85, len(As)))
for ax, (dphidr, rmin, Jv, name) in zip(
        axes, [(dphidr_t, rmin_t, Jt, 't'),
               (dphidr_tau, rmin_tau, Jtau, r'\tau')]):
    rmins = []
    for E, A, c in zip(Eeffs, As, colors):
        o = orbit(dphidr, rmin, E, Jv)
        if o is None:
            continue
        rr, ph, rst = o
        rmins.append((A, rst))
        ax.plot(rr*np.cos(ph), rr*np.sin(ph), color=c, lw=1.4,
                label=rf'$A={A:.1f}$')
    ax.plot(2*M*np.cos(th), 2*M*np.sin(th), 'b--', lw=0.7)
    ax.set_aspect('equal'); ax.set_xlabel('$x$'); ax.set_ylabel('$y$')
    ax.set_title(rf'${name}$-branch breathing family ($J={Jv}$)'
                 '\n' r'$\hat E$ fixed, $A\!\uparrow\Rightarrow E_{\rm eff}\!'
                 r'\downarrow$', fontsize=6.6)
    ax.legend(fontsize=5.2, loc='upper right', framealpha=0.9,
              title=r'$E_{\rm eff}=\hat E/A$', title_fontsize=5.2)
    print(f"{name}: r_min(A) = {[(round(A,1),round(r,3)) for A,r in rmins]}")
    print(f"     J_c=a/E_eff=aA/Ehat: "
          f"{[(round(A,1),round(a*A/Ehat,3)) for A in As]}  (cresce con A)")
savefig(fig, HERE, 'fig_breathing_families')
print('FATTO.')
