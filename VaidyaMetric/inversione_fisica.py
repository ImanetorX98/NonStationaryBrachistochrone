# -*- coding: utf-8 -*-
"""
L'inversione evaporativa di R16 e' fisica o artefatto del regime a massa
negativa? Il modello lineare m=1+mu*v va NEGATIVO per v>1/|mu|, e a
mu_inv~-0.058 con v1=40 si ha m(v1)<0 (singolarita' nuda).

Due controlli, tracciando la massa MINIMA lungo l'orbita:
 A) LINEARE con orbita corta (r1 ridotto) -> m puo' restare >0;
 B) ESPONENZIALE m=m0*exp(-lam*v) -> m>0 SEMPRE per costruzione.

Per ciascuno: scan del tasso di evaporazione, gap = r_min^tau - r_min^t,
e m_min lungo l'orbita. L'inversione e' GENUINA se gap<0 mentre m_min>0.
"""

import os
import sys
import numpy as np
import sympy as sp
from scipy.integrate import solve_ivp
from scipy.optimize import brentq

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from paper_style import COL, set_style, savefig
import matplotlib.pyplot as plt

set_style()
OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'Vaidyafigures')
E_n, J_n = 1.2, 8.0

# ---- derivazione simbolica UNA VOLTA (mm, mp indipendenti) -----------
vv, rr, rp, q = sp.symbols('vv rr rp q', real=True)
E, J, mm, mp = sp.symbols('E J mm mp', real=True)
f = 1 - 2 * mm / rr
Wv = (f - 2 * rp - (f - rp)**2 / E**2) / rr**2

def rpp_sym(F):
    p_r = sp.diff(F, rp)
    dp = (sp.diff(p_r, mm) * mp + sp.diff(p_r, rr) * rp + sp.diff(p_r, rp) * q)
    return sp.solve(sp.Eq(dp - sp.diff(F, rr), 0), q)[0]

rpp_tau_s = rpp_sym((f - rp) / E - J * sp.sqrt(Wv))
rpp_t_s = rpp_sym(1 - J * sp.sqrt(Wv))
u_ = sp.Symbol('u_', real=True)
Wr = (f * u_**2 - 2 * u_ - (f * u_ - 1)**2 / E**2) / rr**2
pv_tau_s = sp.diff((f * u_ - 1) / E - J * sp.sqrt(Wr), u_)
pv_t_s = sp.diff(u_ - J * sp.sqrt(Wr), u_)

def lambdify_model(mfun):
    """mfun = m(vv) simbolico. Ritorna funzioni numeriche."""
    dm = sp.diff(mfun, vv)
    sub = [(mm, mfun), (mp, dm), (E, E_n), (J, J_n)]
    L = lambda expr: sp.lambdify((vv, rr, rp), expr.subs(sub), 'numpy')
    Lp = lambda expr: sp.lambdify((rr, vv, u_), expr.subs(sub), 'numpy')
    mfn = sp.lambdify(vv, mfun, 'numpy')
    return (L(rpp_tau_s), L(rpp_t_s), Lp(pv_tau_s), Lp(pv_t_s), mfn)

def vp_radiale(fv):
    a2 = fv - fv**2 / E_n**2
    a1 = -2 + 2 * fv / E_n**2
    a0 = -1 / E_n**2
    return (-a1 + np.sqrt(a1**2 - 4 * a2 * a0)) / (2 * a2)

def orbita(rpp, pv, mfn, r1, v1):
    fl = 1 - 2 * mfn(v1) / r1
    if fl <= 0 or not np.isfinite(fl):
        return None
    z_rad = vp_radiale(fl)                # u = dv/dr radiale (>0)
    zs = np.linspace(z_rad * (1 + 1e-9), z_rad * 60, 4000)
    with np.errstate(invalid='ignore'):
        pvv = np.array([pv(r1, v1, z) for z in zs])
    br = None
    for i in range(len(zs) - 1):
        if np.isfinite(pvv[i]) and np.isfinite(pvv[i + 1]) \
                and pvv[i] * pvv[i + 1] < 0:
            br = (zs[i], zs[i + 1])
            break
    if br is None:
        return None
    rp1 = 1.0 / brentq(lambda z: pv(r1, v1, z), *br)
    ev = lambda v_, y: y[1]
    ev.terminal, ev.direction = True, 0
    s = solve_ivp(lambda v_, y: [y[1], rpp(v_, y[0], y[1])],
                  [v1, -1e6], [r1, rp1], rtol=1e-11, atol=1e-13,
                  events=[ev], dense_output=True)
    if s.status != 1:
        return None
    v_peri = s.t_events[0][0]
    r_min = s.y_events[0][0][0]
    vg = np.linspace(v_peri, v1, 400)
    m_min = float(np.min(mfn(vg)))       # massa minima lungo l'orbita
    return r_min, v_peri, m_min

def studio(nome, r1, v1, param_grid, mfun_of, pname):
    print(f"\n{'='*66}\n{nome}\n{'='*66}")
    print(f"{pname:>8} {'r_tau':>8} {'r_t':>8} {'gap':>9} {'m_min':>8} "
          f"{'stato':>14}")
    prev_gap = None
    for p in param_grid:
        f_tau, f_t, p_tau, p_t, mfn = lambdify_model(mfun_of(p))
        oa = orbita(f_tau, p_tau, mfn, r1, v1)
        ob = orbita(f_t, p_t, mfn, r1, v1)
        if oa is None or ob is None:
            print(f"{p:8.4f}  {'--':>8} {'--':>8} {'--':>9} {'--':>8} "
                  f"{'no orbita':>14}")
            continue
        gap = oa[0] - ob[0]
        mmin = min(oa[2], ob[2])
        stato = 'INVERSO' if gap < 0 else 'normale'
        fis = 'fisico' if mmin > 0 else 'm<0 NON FIS'
        print(f"{p:8.4f}  {oa[0]:8.4f} {ob[0]:8.4f} {gap:+9.4f} "
              f"{mmin:8.4f}  {stato:>7} {fis}")

# ============ STUDIO A: lineare, orbita corta =========================
lam = sp.Symbol('lam', real=True)
studio("A) LINEARE m=1+mu*v, orbita corta (r1=6, v1=8)",
       6.0, 8.0, np.linspace(0.0, -0.40, 17),
       lambda mu: 1 + mu * vv, "mu")

# ============ STUDIO B: esponenziale (m>0 sempre) =====================
studio("B) ESPONENZIALE m=exp(-lam*v), r1=15, v1=40 (m sempre>0)",
       15.0, 40.0, np.linspace(0.0, 0.20, 17),
       lambda l_: sp.exp(-l_ * vv), "lam")

print("\nLettura: l'inversione e' GENUINA solo dove gap<0 E m_min>0.")

# --------------------------------------------------------------- figura
def curva(mfun_of, r1, v1, grid):
    g, mm_ = [], []
    for p in grid:
        f_tau, f_t, p_tau, p_t, mfn = lambdify_model(mfun_of(p))
        oa = orbita(f_tau, p_tau, mfn, r1, v1)
        ob = orbita(f_t, p_t, mfn, r1, v1)
        if oa is None or ob is None:
            g.append(np.nan); mm_.append(np.nan)
        else:
            g.append(oa[0] - ob[0]); mm_.append(min(oa[2], ob[2]))
    return np.array(g), np.array(mm_)

print("\ngenerazione figura...")
mu_g = np.linspace(0.02, -0.11, 27)
gL, mL = curva(lambda mu: 1 + mu * vv, 15.0, 40.0, mu_g)
lam_g = np.linspace(0.0, 0.25, 26)
gE, mE = curva(lambda l_: sp.exp(-l_ * vv), 15.0, 40.0, lam_g)

fig, (a1, a2) = plt.subplots(2, 1, figsize=(COL, 5.6))
# (a) lineare: gap con zona m<0 ombreggiata
mu_zero = -1.0 / 40.0            # m(v1=40)=1+40mu=0
a1.axvspan(mu_g.min(), mu_zero, color='0.85', zorder=0)
a1.text(-0.09, 2.2, 'm < 0 at anchor\n(unphysical)', fontsize=6, color='0.35')
a1.plot(mu_g, gL, 'C3-', label=r'gap $r_{\min}^{\tau}-r_{\min}^{t}$')
a1.axhline(0, color='k', lw=0.6)
a1.axvline(mu_zero, color='0.4', ls=':', lw=0.9)
a1.axvline(0, color='k', lw=0.4)
a1.set_ylabel('gap')
a1.set_xlabel(r"evaporation rate $\mu=m'$")
a1.set_title('Linear $m=1+\\mu v$: gap crosses 0 ONLY inside $m<0$\n'
             '(the "inversion" is a negative-mass artifact)')
a1.legend(loc='upper left')
# (b) esponenziale: m>0 sempre, gap non inverte
a2.plot(lam_g, gE, 'C0-', label=r'gap $r_{\min}^{\tau}-r_{\min}^{t}$')
a2.axhline(0, color='k', lw=0.6)
a2.set_ylim(0, 3.3)
a2.set_xlabel(r'evaporation rate $\lambda$ ($m=e^{-\lambda v}>0$)')
a2.set_ylabel('gap')
a2.set_title('Exponential $m=e^{-\\lambda v}$ ($m>0$ always):\n'
             'gap shrinks toward a plateau but NEVER inverts')
a2.legend(loc='upper right')
savefig(fig, OUT, 'fig_vaidya_no_inversione_evaporazione')
print("FATTO.")
