# -*- coding: utf-8 -*-
"""
Master plot: behaviour of the tau- and t-branch controlled-rail extremals at the
conformal stationary limit r_e = 2M, as a function of the control costate J.

REVISED for CQG-116884, referee major comment 4.  The previous version painted
the whole closed interval [-Jc, +Jc] as "cusp bounce", which is wrong at BOTH
endpoints.  The correct tau-branch law has FOUR cases (Proposition 3.1 of the
manuscript), established in
    NonStationaryMetrics/paper2/verification/verify_cusp_corner.wls
    NonStationaryMetrics/paper2/verification/verify_marginal_hamilton.wls

  |J| > Jc      smooth periapsis at r_min > r_e
  |J| < Jc      semicubical cusp AT r_e, reached in finite parameter
  J   = -Jc     double root of p_r^2, approach rate does NOT degenerate:
                r_e is an ASYMPTOTIC endpoint (log-divergent parameter),
                with a finite one-sided tangent.  Neither cusp nor reflection.
  J   = +Jc     double root of p_r^2 AND simple root of the approach rate;
                the 0/0 resolves to a finite ingoing rate: crosses r_e.

The shape radical S(r) is even in J, so it cannot distinguish +Jc from -Jc.
What does is the conformal gravitomagnetic momentum at the stationary limit,
ptilde_phi(r_e) = J - Jc, which is linear in J.

t-branch (dichotomy): crosses over the whole window (Jc-, Jc+), reflects
smoothly outside.  Unchanged.

General BL forms (no Doran shift):
  dphi/dr|_tau = J r sqrt(wf)/(D sqrt(D - J^2 w))
  dphi/dr|_t   = K(r) r sqrt(wf)/(D sqrt(D - K^2 w)),  K(r) = (fJ + 2Ma/r)/E
"""
import os
import sys
import numpy as np

_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, _ROOT)
# paper_style.py lives in NonStationaryMetrics/; a fresh clone must find it
# without relying on an unarchived sibling directory (CQG-116884, major 10)
sys.path.insert(0, os.path.join(_ROOT, "NonStationaryMetrics"))
from paper_style import COL, set_style, savefig
import matplotlib.pyplot as plt
from matplotlib.patches import Patch

set_style()
HERE = os.path.dirname(os.path.abspath(__file__))
M, a, E = 1.0, 0.9, 1.2
re = 2*M
Jc = a/E
f = lambda r: 1 - 2*M/r
w = lambda r: E**2 - f(r)
Dl = lambda r: r**2 - 2*M*r + a**2
Kt = lambda r, J: (f(r)*J + 2*M*a/r)/E
rr = np.linspace(re + 1e-6, 20, 500000)

# ------------------------------------------------------------------ t window
pen_t = np.array([np.sum((Dl(rr) - Kt(rr, J)**2*w(rr)) < -1e-9) == 0
                  for J in np.linspace(-12, 6, 36001)])
Jscan = np.linspace(-12, 6, 36001)
edges = Jscan[np.where(np.diff(pen_t.astype(int)))[0]]
Jcm, Jcp = edges.min(), edges.max()

# ------------------------------------------------- self-check of the tau law
# p_r^2 = S(r) / (D_E^2 Delta^2) with S = r (r-2M) D_E (r Delta - J^2 D_E);
# classify the root of S at r_e by the sign of the bracket B(2M).
DE = lambda r: (E**2 - 1)*r + 2*M
B = lambda r, J: r*Dl(r) - J**2*DE(r)
assert abs(B(re, Jc)) < 1e-12, "B(2M) must vanish at |J| = Jc"
assert abs(B(re, -Jc)) < 1e-12, "B(2M) is even in J"
assert B(re, 0.6*Jc) > 0, "sub-marginal: simple root, cusp"
assert B(re, 1.4*Jc) < 0, "super-marginal: turning point outside r_e"
# the approach rate T(r_e) = -a (J - Jc)/Pbar_e vanishes only at J = +Jc
Pbar_e = 4*M**2 + 2*a**2 + a**2/E**2
T_re = lambda J: -a*(J - Jc)/Pbar_e
assert abs(T_re(Jc)) < 1e-15, "rate must degenerate at the prograde marginal"
assert abs(T_re(-Jc)) > 1e-3, "rate must NOT degenerate at the retrograde one"
# closed-form retrograde decay rate, cross-checked in the .wls scripts
C_retro = -np.sqrt(4*M**2 + Jc**2)/(8*M**2)

print(f"tau: Jc = {Jc:+.3f}")
print(f"     J = +Jc  crosses  (finite ingoing rate)")
print(f"     J = -Jc  asymptotic, dr/ds -> C (r-r_e), C = {C_retro:.9f}")
print(f"t:   window ({Jcm:.3f}, {Jcp:.3f})")

# ------------------------------------------------------------------- figure
fig, ax = plt.subplots(figsize=(COL, 3.1))
Jmin, Jmax = -10, 6

# Okabe-Ito: safe under all three common CVD types; hatching carries the same
# distinction in greyscale and in print.
c_cross  = '#009E73'   # crosses r_e
c_smooth = '#0072B2'   # smooth periapsis outside r_e
c_cusp   = '#E69F00'   # semicubical cusp at r_e
c_asym   = '#CC79A7'   # asymptotic endpoint

# --- t branch (y = 1) ---
yt = 1.0
ax.add_patch(plt.Rectangle((Jmin, yt-0.18), Jcm-Jmin, 0.36,
                           facecolor=c_smooth, edgecolor='white', lw=0.6))
ax.add_patch(plt.Rectangle((Jcm, yt-0.18), Jcp-Jcm, 0.36,
                           facecolor=c_cross, edgecolor='white', lw=0.6,
                           hatch='///'))
ax.add_patch(plt.Rectangle((Jcp, yt-0.18), Jmax-Jcp, 0.36,
                           facecolor=c_smooth, edgecolor='white', lw=0.6))
ax.plot([Jcm, Jcm], [yt-0.18, yt+0.18], 'k', lw=1.0)
ax.plot([Jcp, Jcp], [yt-0.18, yt+0.18], 'k', lw=1.0)
ax.text(Jcm, yt+0.26, r'$J_c^-=%.2f$' % Jcm, ha='center', fontsize=6)
ax.text(Jcp, yt+0.26, r'$J_c^+=%.2f$' % Jcp, ha='center', fontsize=6)
ax.text(Jmin-0.4, yt, r'$t$', ha='right', va='center', fontsize=9)

# --- tau branch (y = 0) ---
ytau = 0.0
ax.add_patch(plt.Rectangle((Jmin, ytau-0.18), (-Jc)-Jmin, 0.36,
                           facecolor=c_smooth, edgecolor='white', lw=0.6))
ax.add_patch(plt.Rectangle((-Jc, ytau-0.18), 2*Jc, 0.36,
                           facecolor=c_cusp, edgecolor='white', lw=0.6,
                           hatch='...'))
ax.add_patch(plt.Rectangle((Jc, ytau-0.18), Jmax-Jc, 0.36,
                           facecolor=c_smooth, edgecolor='white', lw=0.6))
# the two marginal values are measure-zero cases and are NOT part of the
# cusp interval: mark them individually
ax.plot([Jc], [ytau], marker='*', color=c_cross, ms=13, mec='k', mew=0.5,
        zorder=6, clip_on=False)
ax.plot([-Jc], [ytau], marker='o', color='white', ms=7, mec=c_asym, mew=1.6,
        zorder=6, clip_on=False)
ax.plot([-Jc, -Jc], [ytau-0.18, ytau+0.18], 'k', lw=1.0)
ax.plot([Jc, Jc], [ytau-0.18, ytau+0.18], 'k', lw=1.0)
# +-Jc are only 1.5 apart, so the two labels are led away from the bar to
# opposite sides instead of being centred under it
ax.annotate(r'$-J_c=-%.2f$' % Jc, xy=(-Jc, ytau-0.19),
            xytext=(-Jc-2.6, ytau-0.52), fontsize=6, ha='center',
            arrowprops=dict(arrowstyle='-', lw=0.5, color='0.35'))
ax.annotate(r'$+J_c=+%.2f$' % Jc, xy=(Jc, ytau-0.19),
            xytext=(Jc+2.2, ytau-0.52), fontsize=6, ha='center',
            arrowprops=dict(arrowstyle='-', lw=0.5, color='0.35'))
ax.text(Jmin-0.4, ytau, r'$\tau$', ha='right', va='center', fontsize=9)

ax.set_xlim(Jmin-1.2, Jmax+0.5)
ax.set_ylim(-0.8, 1.7)
ax.set_yticks([])
ax.set_xlabel(r'control costate $J$')
ax.set_title(r'behaviour at $r_e=2M$: $\tau$ four cases, $t$ dichotomy')

legend = [Patch(facecolor=c_cross, hatch='///', label=r'crosses $r_e$'),
          Patch(facecolor=c_smooth, label=r'smooth periapsis, $r_{\min}>r_e$'),
          Patch(facecolor=c_cusp, hatch='...',
                label=r'cusp at $r_e$ ($|J|<J_c$, open)'),
          plt.Line2D([], [], marker='*', color=c_cross, mec='k', ms=9, ls='',
                     label=r'$J=+J_c$: crosses, finite rate'),
          plt.Line2D([], [], marker='o', color='white', mec=c_asym, mew=1.6,
                     ms=7, ls='', label=r'$J=-J_c$: asymptotic, not attained')]
ax.legend(handles=legend, fontsize=5.5, loc='upper left',
          bbox_to_anchor=(0.0, -0.24), ncol=2, framealpha=0.9)
savefig(fig, HERE, 'fig_master_penetration_taut')
print('DONE.')
