# -*- coding: utf-8 -*-
"""
Ergosphere-penetration threshold (t-branch) in the (A, J) plane.
Two mechanisms:
  (1) freezing wall: E_eff<1 -> r_w=2M/(1-E_eff^2); launch forbidden if r_w<r0
      => A_c^wall = Ehat/sqrt(1-2M/r0)  (vertical, independent of J).
  (2) ergosphere: the orbit penetrates if the outer turning point < r_e; boundary
      J_c^+(A) = 2M^2/a + a + a A^2/(2 Ehat^2).
Each (A,J) is classified numerically as forbidden / scattering / penetrating,
and the analytic boundaries are overlaid.
"""
import numpy as np
import os, sys
from scipy.optimize import brentq
_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, _ROOT)
# paper_style.py lives in NonStationaryMetrics/; a fresh clone must find it
# without relying on an unarchived sibling directory
sys.path.insert(0, os.path.join(_ROOT, "NonStationaryMetrics"))
from paper_style import COL, set_style, savefig
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap

set_style()
HERE = os.path.dirname(os.path.abspath(__file__))
M, a, Ehat, r0 = 1.0, 0.9, 1.4, 12.0
r_plus = M + np.sqrt(M**2 - a**2); r_e = 2*M


def R6(r, E, J):
    Q2 = (2*E**2*J**2*M*r-E**2*J**2*r**2-4*E**2*J*M*a*r+2*E**2*M*a**2*r
          + E**2*a**2*r**2+E**2*r**4+4*J**2*M**2-4*J**2*M*r+J**2*r**2
          - 8*J*M**2*a+4*J*M*a*r+4*M**2*a**2)
    return r*Q2*((E**2-1)*r+2*M)


def classify(A, J):
    E = Ehat/A
    if R6(r0, E, J) <= 0:                      # launch forbidden (beyond the wall)
        return 0                               # forbidden
    rs = np.linspace(r_plus+1e-3, r0-1e-3, 4000); v = R6(rs, E, J)
    idx = np.where(np.diff(np.sign(v)))[0]
    roots = [brentq(lambda r: R6(r, E, J), rs[i], rs[i+1]) for i in idx if rs[i] > r_plus]
    if not roots:
        return 3                               # PLUNGE: no turning point -> horizon
    r_turn = max(roots)                         # outer turning point (where it stops)
    if r_turn < r_e:
        return 2                               # penetrates and bounces inside the ergosphere
    return 1                                    # scattering (turning point outside r_e)

# --- grid (A, J) ---
As = np.linspace(0.5, 2.4, 150); Js = np.linspace(-10.0, 5.0, 220)
Z = np.zeros((len(Js), len(As)))
for i, J in enumerate(Js):
    for j, A in enumerate(As):
        Z[i, j] = classify(A, J)

# --- analytic boundaries ---
Ac_wall = Ehat/np.sqrt(1-2*M/r0)               # the freezing wall reaches r0
Jcp = lambda A: 2*M**2/a + a + a*A**2/(2*Ehat**2)   # ergosphere threshold J_c^+
print(f"A_c^wall(r0={r0}) = {Ac_wall:.4f}")
print(f"J_c^+(A=1) = {Jcp(1.0):.4f},  J_c^+(A=1.5) = {Jcp(1.5):.4f}")

fig, ax = plt.subplots(figsize=(1.35*COL, COL*1.05))
# 0 forbidden, 1 scattering, 2 penetrate-bounce, 3 plunge
cmap = ListedColormap(['0.85', '#f4a582', '#4393c3', '#08519c'])
ax.pcolormesh(As, Js, Z, cmap=cmap, shading='auto', vmin=-0.5, vmax=3.5)
# analytic boundaries
ax.axvline(Ac_wall, color='k', ls='--', lw=1.2)
ax.text(Ac_wall+0.03, -8.5, r'$A_c^{\rm wall}=\hat E/\sqrt{1-2M/r_0}$',
        rotation=90, fontsize=6, va='bottom')
Aline = np.linspace(As[0], As[-1], 200)
ax.plot(Aline, Jcp(Aline), 'k-', lw=1.4)
ax.text(0.6, Jcp(0.6)+0.2, r'$J_c^+(A)$', fontsize=6.5)
from matplotlib.patches import Patch
leg = [Patch(fc='#08519c', label='plunge (continuation, optimality not claimed)'),
       Patch(fc='#4393c3', label='penetrate + bounce (continuation, optimality not claimed)'),
       Patch(fc='#f4a582', label='scattering'),
       Patch(fc='0.85', label='forbidden (beyond wall)')]
ax.legend(handles=leg, fontsize=5.4, loc='lower right', framealpha=0.95)
ax.set_xlabel('$A$ (conformal factor)'); ax.set_ylabel('$J$')
ax.set_title(r'$t$-branch ergosphere-penetration phase diagram'
             '\n' rf'($M=1,a=0.9,\hat E=1.4,r_0=12$); $A_c^{{\rm wall}}={Ac_wall:.3f}$',
             fontsize=6.8)
savefig(fig, HERE, 'fig_penetration_threshold')
# --- check of the boundaries against the numerics ---
print("\ncheck of J_c^+ (penetration/scattering boundary) at a few A:")
for A in [0.8, 1.0, 1.3]:
    E = Ehat/A
    # find J where the outer turning point = r_e (numerically)
    def turn_at_re(J):
        rs = np.linspace(r_plus+1e-3, r0, 4000); v = R6(rs, E, J)
        idx = np.where(np.diff(np.sign(v)))[0]
        roots = [brentq(lambda r: R6(r, E, J), rs[i], rs[i+1]) for i in idx if rs[i] > r_plus]
        return (max(roots) if roots else r_plus) - r_e
    try:
        Jc_num = brentq(turn_at_re, 2.0, 5.0)
        print(f"  A={A}: J_c^+ num={Jc_num:.4f}, analytic={Jcp(A):.4f}, "
              f"diff={abs(Jc_num-Jcp(A)):.2e}")
    except Exception as e:
        print(f"  A={A}: {e}")
print('DONE.')
