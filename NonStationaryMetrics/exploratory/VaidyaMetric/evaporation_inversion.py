# -*- coding: utf-8 -*-
"""
############################################################################
## DEPRECATED -- DO NOT USE. WRONG result (negative-mass artefact).       ##
## The inversion found here (mu_inv ~ -0.058) falls where m=1+mu*v is     ##
## already NEGATIVE at the anchor v1=40 (naked singularity): unphysical.  ##
## Refuted by VaidyaMetric/physical_inversion.py (2D scan, linear and      ##
## exponential, m>0). Correct conclusion: PHYSICAL evaporation does not   ##
## invert the plunge; r_min^t < r_min^tau stays robust (R13). See         ##
## VaidyaResults.md R16. Correct figure: fig_vaidya_no_inversion_evaporation. ##
## Kept only as the historical record of the false alarm.                 ##
############################################################################

Plunge inversion from EVAPORATION in Vaidya (non-rotating). [WRONG]

R13 (static) and R15 (small mu): r_min^t < r_min^tau, t goes deeper,
the sign is locked -- it SEEMED not invertible without rotation. FALSE for
strong evaporation: the dynamical correction ~m' (teleological p_v memory,
R15) NARROWS the gap for m' < 0 and INVERTS it beyond a threshold mu_inv < 0.

[NB: the statement above is WRONG -- see the banner. The inversion
exists only because m becomes negative, an unphysical regime.]

Finds mu_inv (gap = r_tau - r_t changes sign) and makes a figure of gap(mu).
"""

import os
import sys
import numpy as np
import sympy as sp
from scipy.integrate import solve_ivp
from scipy.optimize import brentq

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))   # NonStationaryMetrics/, for paper_style
from paper_style import COL, set_style, savefig
import matplotlib.pyplot as plt

set_style()
OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'Vaidyafigures')
m0, E_n, J_n, r1 = 1.0, 1.2, 8.0, 15.0

vv, rr, rp, q = sp.symbols('vv rr rp q', real=True)
E, J, mm, mp, mu = sp.symbols('E J mm mp mu', real=True)
f = 1 - 2 * mm / rr
Wv = (f - 2 * rp - (f - rp)**2 / E**2) / rr**2
sub = [(mm, m0 + mu * vv), (mp, mu), (E, E_n), (J, J_n)]

def mk(F):
    p_r = sp.diff(F, rp)
    dp = (sp.diff(p_r, vv) + sp.diff(p_r, mm) * mp
          + sp.diff(p_r, rr) * rp + sp.diff(p_r, rp) * q)
    rpp = sp.solve(sp.Eq(dp - sp.diff(F, rr), 0), q)[0]
    return sp.lambdify((vv, rr, rp, mu), rpp.subs(sub), 'numpy')

rpp_tau = mk((f - rp) / E - J * sp.sqrt(Wv))
rpp_t = mk(1 - J * sp.sqrt(Wv))
u_ = sp.Symbol('u_', real=True)
Wr = (f * u_**2 - 2 * u_ - (f * u_ - 1)**2 / E**2) / rr**2
pvT = sp.lambdify((rr, vv, u_, mu),
                  sp.diff((f * u_ - 1) / E - J * sp.sqrt(Wr), u_).subs(sub),
                  'numpy')
pvt = sp.lambdify((rr, vv, u_, mu),
                  sp.diff(u_ - J * sp.sqrt(Wr), u_).subs(sub), 'numpy')

def vpr(fv):
    a2 = fv - fv**2 / E_n**2
    a1 = -2 + 2 * fv / E_n**2
    a0 = -1 / E_n**2
    return (-a1 + np.sqrt(a1**2 - 4 * a2 * a0)) / (2 * a2)

def orb(rpp, pv, mu_n, v1=40.0):
    fl = 1 - 2 * (m0 + mu_n * v1) / r1
    lo = vpr(fl) * (1 + 1e-9)
    hi = lo * 2
    while pv(r1, v1, hi, mu_n) < 0:
        hi *= 2
    rp1 = 1.0 / brentq(lambda z: pv(r1, v1, z, mu_n), lo, hi)
    ev = lambda v_, y: y[1]
    ev.terminal, ev.direction = True, 0
    s = solve_ivp(lambda v_, y: [y[1], rpp(v_, y[0], y[1], mu_n)],
                  [v1, -1e6], [r1, rp1], rtol=1e-12, atol=1e-14, events=[ev])
    return s.y_events[0][0][0], s.t_events[0][0]

def gap(mu_n):
    return orb(rpp_tau, pvT, mu_n)[0] - orb(rpp_t, pvt, mu_n)[0]

# inversion threshold
lo, hi = -0.06, -0.04
for _ in range(60):
    m = 0.5 * (lo + hi)
    if gap(m) > 0:
        hi = m
    else:
        lo = m
mu_inv = 0.5 * (lo + hi)
r_inv = orb(rpp_tau, pvT, mu_inv)[0]
print(f"mu_inv = {mu_inv:.5f}   (E={E_n}, J={J_n:.0f}, r1={r1})")
print(f"  at mu_inv: r_tau = r_t = {r_inv:.5f}")
print("  mu > mu_inv (static/accretion):    r_t < r_tau (t deeper)")
print("  mu < mu_inv (strong evaporation):  r_tau < r_t (INVERTED)")

# gap(mu) curve
mus = np.linspace(0.02, -0.115, 40)
gaps, rtau_a, rt_a = [], [], []
for mv in mus:
    try:
        rt = orb(rpp_tau, pvT, mv)[0]
        rtt = orb(rpp_t, pvt, mv)[0]
    except Exception:
        rt = rtt = np.nan
    rtau_a.append(rt)
    rt_a.append(rtt)
    gaps.append(rt - rtt)
gaps = np.array(gaps)
rtau_a, rt_a = np.array(rtau_a), np.array(rt_a)

fig, (a1, a2) = plt.subplots(2, 1, figsize=(COL, 5.4), sharex=True)
a1.plot(mus, rtau_a, 'C3-', label=r'$r_{\min}^{\tau}$')
a1.plot(mus, rt_a, 'C0-', label=r'$r_{\min}^{t}$')
a1.axvline(mu_inv, color='k', ls=':', lw=1.0)
a1.axvline(0, color='k', lw=0.5)
a1.set_ylabel('$r_{\\min}$')
a1.set_title(f'Evaporative plunge inversion (Vaidya, non-rotating)\n'
             f'$E={E_n}$, $J={J_n:.0f}$: curves cross at '
             f'$\\mu_{{inv}}={mu_inv:.4f}$')
a1.legend()
a2.plot(mus, gaps, 'C2-')
a2.axhline(0, color='k', lw=0.6)
a2.axvline(mu_inv, color='k', ls=':', lw=1.0)
a2.axvline(0, color='k', lw=0.5)
a2.annotate(f'$\\mu_{{inv}}={mu_inv:.4f}$', (mu_inv, 0),
            xytext=(mu_inv - 0.005, 1.2), fontsize=6.5,
            arrowprops=dict(arrowstyle='->', lw=0.6))
a2.text(-0.10, -0.7, r'$\tau$ deeper', fontsize=6.5, color='C2')
a2.text(-0.01, 2.2, r'$t$ deeper (static)', fontsize=6.5, color='C2')
a2.set_xlabel(r"accretion/evaporation rate $\mu = m'$")
a2.set_ylabel(r'gap $r_{\min}^{\tau}-r_{\min}^{t}$')
a2.set_title('gap changes sign: a second, EVAPORATIVE inversion mechanism')
savefig(fig, OUT, 'fig_vaidya_inversion_evaporation')
print("DONE.")
