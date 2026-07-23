# -*- coding: utf-8 -*-
"""
Orbite PENETRANTI ramo t (J=3 < J_c^+): la brachistocrona congelata phi_0(r;E_eff)
ATTRAVERSA l'ergosfera (r_e=2M) e gira dentro, vicino all'orizzonte r_+.
Famiglia che respira: A cresce -> E_eff=Ehat/A decresce -> l'orbita evolve.
BL: dphi/dr=K/sqrt(R6); il turning interno (R6=0) e' DENTRO l'ergosfera.
Pannello (a): famiglia penetrante x-y. (b): raggio di turning interno vs A.
"""
import os, sys
import numpy as np
from scipy.integrate import quad
from scipy.optimize import brentq
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from paper_style import COL, set_style, savefig
import matplotlib.pyplot as plt
from matplotlib import cm

set_style()
HERE = os.path.dirname(os.path.abspath(__file__))
M, a, Ehat, J = 1.0, 0.9, 1.4, 3.0
r0 = 12.0
r_plus = M + np.sqrt(M**2 - a**2); r_e = 2*M
print(f"r_+={r_plus:.4f}, r_e={r_e}")


def R6(r, E):
    Q2 = (2*E**2*J**2*M*r-E**2*J**2*r**2-4*E**2*J*M*a*r+2*E**2*M*a**2*r
          + E**2*a**2*r**2+E**2*r**4+4*J**2*M**2-4*J**2*M*r+J**2*r**2
          - 8*J*M**2*a+4*J*M*a*r+4*M**2*a**2)
    return r*Q2*((E**2-1)*r+2*M)


def Kf(r, E): return r*((E**2-1)*r+2*M)*(J*(r-2*M)+2*M*a)/(r**2-2*M*r+a**2)


def inner_turning(E):
    """radice piu' esterna di R6 SOPRA r_+ (turning interno all'ergosfera)."""
    rs = np.linspace(r_plus+0.005, r0, 30000); v = R6(rs, E)
    idx = np.where(np.diff(np.sign(v)))[0]
    roots = [brentq(lambda r: R6(r, E), rs[i], rs[i+1]) for i in idx]
    roots = [r for r in roots if r > r_plus]
    return max(roots) if roots else np.nan


def orbit(E):
    rt = inner_turning(E)
    if not np.isfinite(rt):
        return None
    # griglia GRADUATA (passo ridotto vicino al turning, integrando ~1/sqrt(r-rt))
    u = np.linspace(0, 1, 2500)
    rr = r0 - (r0-(rt+5e-4))*(1-np.cos(u*np.pi/2))
    # phi(r) con substituzione che addolcisce la singolarita' al turning
    ph = np.empty_like(rr)
    for i, rc in enumerate(rr):
        ph[i] = quad(lambda x: Kf(x, E)/np.sqrt(max(R6(x, E), 1e-300)),
                     rc, r0, limit=400, points=[rt] if rc < rt+0.05 else None)[0]
    # incoming + outgoing (mirror alla svolta)
    return (np.concatenate([rr, rr[::-1]]),
            np.concatenate([ph, 2*ph[-1]-ph[::-1]]), rt)


As = np.array([1.0, 1.15, 1.3, 1.45, 1.6])
Eeffs = Ehat/As
fig, (ax, axb) = plt.subplots(1, 2, figsize=(2*COL, COL*1.02))
th = np.linspace(0, 2*np.pi, 300)
cols = cm.plasma(np.linspace(0.05, 0.8, len(As)))
rts = []
for E, A, c in zip(Eeffs, As, cols):
    o = orbit(E)
    if o is None:
        print(f"  A={A}: nessuna orbita"); continue
    rr, ph, rt = o; rts.append((A, rt))
    if rt < r_e:      # pannello (a): solo orbite che PENETRANO l'ergosfera
        ax.plot(rr*np.cos(ph), rr*np.sin(ph), color=c, lw=1.5,
                label=rf'$A={A:.2f}$ ($r_{{\rm t}}={rt:.2f}$)')
        # turning: marca che sta FUORI dall'orizzonte (dentro l'ergosfera)
        i_t = len(rr)//2; ax.plot(rr[i_t]*np.cos(ph[i_t]), rr[i_t]*np.sin(ph[i_t]),
                                  'o', color=c, ms=4, mec='k', mew=0.4)
    print(f"  A={A:.2f}: E_eff={E:.3f}, turning interno r={rt:.3f} "
          f"({'DENTRO' if rt < r_e else 'fuori'} ergosfera)")
ax.fill(r_e*np.cos(th), r_e*np.sin(th), color='b', alpha=0.06)
ax.plot(r_e*np.cos(th), r_e*np.sin(th), 'b:', lw=0.9, label='$r_e$ (ergosphere)')
ax.plot(r_plus*np.cos(th), r_plus*np.sin(th), 'k-', lw=0.9, label='$r_+$ (horizon)')
ax.plot(r0, 0, 'ks', ms=3)
ax.set_aspect('equal'); ax.set_xlabel('$x$'); ax.set_ylabel('$y$')
ax.set_xlim(-3.5, 12.5); ax.set_ylim(-3.2, 3.2)
ax.set_title(rf'$t$-branch PENETRATING ($J={J}$): enters ergosphere,'
             '\n' r'turns ($\circ$) OUTSIDE horizon, then exits', fontsize=6.6)
ax.annotate(r'turnings $r_+<r_{\rm t}<r_e$', xy=(0.3, -1.5),
            fontsize=5.6, ha='center', color='0.2')
ax.legend(fontsize=5.2, loc='lower right', framealpha=0.9,
          title=r'$E_{\rm eff}=\hat E/A$ (penetrating)', title_fontsize=5.2)
# pannello b: turning interno vs A (evoluzione adiabatica)
Ap, rtp = zip(*rts)
axb.plot(Ap, rtp, 'o-', color='C3', lw=1.4)
axb.axhline(r_e, color='b', ls=':', lw=0.9); axb.axhline(r_plus, color='k', lw=0.9)
axb.text(As[0], r_e+0.02, '$r_e$ (ergosphere)', color='b', fontsize=6)
axb.text(As[0], r_plus-0.05, '$r_+$ (horizon)', fontsize=6)
axb.set_xlabel('$A$'); axb.set_ylabel('inner turning $r_{\\rm turn}$')
axb.set_title('penetration depth vs conformal factor\n'
              r'($A\uparrow$: orbit turns deeper toward $r_+$)', fontsize=6.6)
savefig(fig, HERE, 'fig_phi_penetrating')
print('FATTO.')
