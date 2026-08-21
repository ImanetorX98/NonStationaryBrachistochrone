# -*- coding: utf-8 -*-
"""
Orbita numerica Thakurta-Kerr (ramo t, J=3.2) che PENETRA l'ergosfera, RIMBALZA
al turning interno (dentro r_e, fuori r_+) e RIESCE. Marca i due attraversamenti
di r_e (entrata e uscita) e il turning.
"""
import os, sys
import numpy as np
from scipy.integrate import quad
from scipy.optimize import brentq
_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, _ROOT)
# paper_style.py lives in NonStationaryMetrics/; a fresh clone must find it
# without relying on an unarchived sibling directory (CQG-116884, major 10)
sys.path.insert(0, os.path.join(_ROOT, "NonStationaryMetrics"))
from paper_style import COL, set_style, savefig
import matplotlib.pyplot as plt

set_style()
HERE = os.path.dirname(os.path.abspath(__file__))
M, a, E, J = 1.0, 0.9, 1.4, 3.2
r0 = 12.0
r_plus = M+np.sqrt(M**2-a**2); r_e = 2*M


def R6(r):
    Q2 = (2*E**2*J**2*M*r-E**2*J**2*r**2-4*E**2*J*M*a*r+2*E**2*M*a**2*r
          + E**2*a**2*r**2+E**2*r**4+4*J**2*M**2-4*J**2*M*r+J**2*r**2
          - 8*J*M**2*a+4*J*M*a*r+4*M**2*a**2)
    return r*Q2*((E**2-1)*r+2*M)


def Kf(r): return r*((E**2-1)*r+2*M)*(J*(r-2*M)+2*M*a)/(r**2-2*M*r+a**2)


rs = np.linspace(r_plus+1e-3, r0, 30000); v = R6(rs)
idx = np.where(np.diff(np.sign(v)))[0]
r_turn = max(brentq(R6, rs[i], rs[i+1]) for i in idx if rs[i] > r_plus)
print(f"r_+={r_plus:.4f}, r_e={r_e}, turning={r_turn:.4f} "
      f"({'DENTRO' if r_turn < r_e else 'fuori'} ergosfera)")

# orbita: incoming r0->turning, poi mirror (outgoing)
u = np.linspace(0, 1, 3000)
rr = r0 - (r0-(r_turn+3e-4))*(1-np.cos(u*np.pi/2))
ph = np.array([quad(lambda x: Kf(x)/np.sqrt(max(R6(x), 1e-300)), r, r0,
                    limit=400, points=[r_turn] if r < r_turn+0.05 else None)[0]
               for r in rr])
rr_full = np.concatenate([rr, rr[::-1]])
ph_full = np.concatenate([ph, 2*ph[-1]-ph[::-1]])
x, y = rr_full*np.cos(ph_full), rr_full*np.sin(ph_full)

fig, ax = plt.subplots(figsize=(COL*1.15, COL*1.1))
th = np.linspace(0, 2*np.pi, 400)
# colora il tratto dentro l'ergosfera diverso
inside = rr_full < r_e
ax.plot(x, y, 'C0-', lw=1.6, label='orbit (outside $r_e$)')
xi = np.where(inside, x, np.nan); yi = np.where(inside, y, np.nan)
ax.plot(xi, yi, 'C3-', lw=2.2, label='inside ergosphere (co-rotating)')
ax.fill(r_e*np.cos(th), r_e*np.sin(th), color='b', alpha=0.06)
ax.plot(r_e*np.cos(th), r_e*np.sin(th), 'b:', lw=1.0, label='$r_e$ (ergosphere)')
ax.plot(r_plus*np.cos(th), r_plus*np.sin(th), 'k-', lw=1.0, label='$r_+$ (horizon)')
# turning (rimbalzo) e attraversamenti di r_e
it = len(rr_full)//2
ax.plot(x[it], y[it], 'r^', ms=7, mec='k', label='turning (bounce)')
cross = np.where(np.diff(np.sign(rr_full-r_e)))[0]
ax.plot(x[cross], y[cross], 'ko', ms=4)
ax.plot(r0, 0, 'ks', ms=4)
ax.annotate('enters', xy=(x[cross[0]], y[cross[0]]), fontsize=6,
            xytext=(x[cross[0]]+1.2, y[cross[0]]+0.8),
            arrowprops=dict(arrowstyle='->', lw=0.5))
ax.annotate('exits', xy=(x[cross[-1]], y[cross[-1]]), fontsize=6,
            xytext=(x[cross[-1]]+0.5, y[cross[-1]]-1.4),
            arrowprops=dict(arrowstyle='->', lw=0.5))
ax.set_aspect('equal'); ax.set_xlabel('$x$'); ax.set_ylabel('$y$')
ax.set_xlim(-3.5, 12.5); ax.set_ylim(-4.5, 4)
ax.set_title(rf'Thakurta--Kerr $t$-brachistochrone ($J={J}$): penetrates $r_e$,'
             '\n' rf'bounces at $r_{{\rm t}}={r_turn:.2f}$ (inside ergosphere, '
             r'outside $r_+$), exits', fontsize=6.6)
ax.legend(fontsize=5.4, loc='lower right', framealpha=0.95)
savefig(fig, HERE, 'fig_penetrate_bounce_orbit')
print('FATTO.')
