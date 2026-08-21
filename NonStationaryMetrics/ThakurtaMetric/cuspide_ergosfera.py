# -*- coding: utf-8 -*-
"""
Dimostrazione: le brachistocrone tau equatoriali di Thakurta-Kerr con
|J| != J_c che RAGGIUNGONO l'ergosfera r_e = 2M vi si riflettono con una
CUSPIDE (non un periasse liscio).

Struttura (via transfer theorem, in variabili efficaci E=Ehat/A, J->J/A;
r_e=2M conforme-invariante, J_c = a/E):

   dphi/dr = J sqrt(w f) r / ( Dl sqrt(Dl - J^2 w) ),   w = E^2 - f

Vicino a r_e (f -> 0):
   f ~ (r - r_e)/(2M),   w -> E^2,   Dl -> a^2,
   Dl - J^2 w -> a^2 - J^2 E^2
 => dphi/dr ~ K sqrt(r - r_e),   K = J E r_e / (a^2 sqrt(a^2 - J^2 E^2) sqrt(2M))
 => phi - phi_e ~ (2K/3)(r - r_e)^{3/2}   (CUSPIDE: potenza 3/2)

Classificazione locale a r_e -- QUATTRO casi, non tre (rivisto per CQG-116884,
commento maggiore 4).  Il radicale di forma S(r) dipende da J solo via J^2 e
quindi NON puo' distinguere +J_c da -J_c; le equazioni di Hamilton si', perche'
ptilde_phi(r_e) = J - J_c e' lineare in J.

  |J| > J_c :  a^2 - J^2 E^2 < 0  -> svolta PRIMA di r_e (periasse liscio,
               dphi/dr -> infinito; esponente 2 nel piano (phi, r))
  |J| < J_c :  a^2 - J^2 E^2 > 0 finito, f->0 -> dphi/dr ~ sqrt(r-r_e)
               -> CUSPIDE semicubica (esponente 2/3), tangente radiale,
               raggiunta in parametro FINITO
  J   = -J_c: radice DOPPIA di p_r^2 a r_e, ma il rate di avvicinamento
               T(r_e) = -a(J-J_c)/Pbar_e NON degenera: dr/ds ~ C (r-r_e) con
               C = -sqrt(4M^2+J_c^2)/(8M^2), parametro log-divergente.
               r_e e' un ENDPOINT ASINTOTICO: ne' cuspide ne' riflessione.
               La tangente resta finita -> esponente 1.
  J   = +J_c: radice doppia di p_r^2 E radice semplice di T; il 2 0/0 si risolve
               a un rate ingoing finito -> raggiunge r_e e lo attraversa.

Verificato in:
  NonStationaryMetrics/paper2/verification/verify_cusp_corner.wls
  NonStationaryMetrics/paper2/verification/verify_marginal_hamilton.wls
Gli esponenti 2/3, 1, 2 sono misurati dal pannello inferiore di questa figura e
stampati a fine esecuzione.
"""

import os
import sys
import numpy as np
import sympy as sp
from scipy.integrate import quad
from scipy.optimize import brentq

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from paper_style import COL, set_style, savefig
import matplotlib.pyplot as plt

set_style()
OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'Thakurtafigures')
M, a, E = 1.0, 0.9, 1.2
r_e = 2 * M
Jc = a / E

print("=" * 70)
print(f"Kerr efficace: M={M}, a={a}, E(eff)={E}, r_e={r_e}, J_c=a/E={Jc:.4f}")
print("(Thakurta-Kerr: E=Ehat/A, J->J/A; r_e e la struttura sono le stesse)")
print("=" * 70)

# ---- [1] espansione simbolica di dphi/dr vicino a r_e -----------------
r, J = sp.symbols('r J', positive=True)
f = 1 - 2 * M / r
w = E**2 - f
Dl = r**2 - 2 * M * r + a**2
dphidr = J * sp.sqrt(w * f) * r / (Dl * sp.sqrt(Dl - J**2 * w))

eps = sp.symbols('epsilon', positive=True)     # r = r_e + eps
ser = sp.series(dphidr.subs(r, r_e + eps), eps, 0, 1).removeO()
print("\n[1] dphi/dr vicino a r_e (r = r_e + eps), termine dominante:")
lead = sp.simplify(sp.limit(dphidr.subs(r, r_e + eps) / sp.sqrt(eps),
                            eps, 0, '+'))
print(f"    dphi/dr ~ K*sqrt(eps),  K = {sp.nsimplify(lead, rational=False)}")
K_expr = J * E * r_e / (a**2 * sp.sqrt(a**2 - J**2 * E**2) * sp.sqrt(2 * M))
print(f"    formula:  K = J E r_e/(a^2 sqrt(a^2-J^2E^2) sqrt(2M))")
print(f"    coincidono? {sp.simplify(lead - K_expr)==0}")

# ---- [2] esponente numerico per |J|<J_c, |J|=J_c, |J|>J_c -------------
print("\n[2] esponente numerico di dphi/dr ~ (r-r_e)^p vicino a r_e:")
def dphidr_n(rv, Jv):
    fv = 1 - 2 * M / rv
    wv = E**2 - fv
    Dlv = rv**2 - 2 * M * rv + a**2
    rad = Dlv - Jv**2 * wv
    if rad <= 0 or wv * fv < 0:
        return np.nan
    return Jv * np.sqrt(wv * fv) * rv / (Dlv * np.sqrt(rad))

for lab, Jv in [('|J|<J_c', 0.5), ('|J|<J_c', 0.7),
                ('|J|=J_c', Jc), ('|J|>J_c', 0.9)]:
    ds = r_e + np.array([1e-4, 1e-5, 1e-6])
    vals = np.array([dphidr_n(rv, Jv) for rv in ds])
    if np.any(~np.isfinite(vals)):
        # trova dove diventa complesso: svolta prima di r_e
        rr = np.linspace(r_e + 1e-6, r_e + 0.5, 20000)
        vv = np.array([dphidr_n(x, Jv) for x in rr])
        rturn = rr[np.isfinite(vv)][0] if np.any(np.isfinite(vv)) else np.nan
        print(f"    J={Jv:.4f} ({lab}): svolta a r={rturn:.5f} > r_e "
              f"(periasse LISCIO, dphi/dr->inf)")
        continue
    p = np.polyfit(np.log(ds - r_e), np.log(vals), 1)[0]
    fin = dphidr_n(r_e + 1e-8, Jv)
    tag = ('CUSPIDE (p~1/2)' if abs(p - 0.5) < 0.05
           else ('LISCIO/attraversa (p~0, dphi/dr finito)'
                 if abs(p) < 0.05 else f'p={p:.3f}'))
    print(f"    J={Jv:.4f} ({lab}): esponente p = {p:+.4f}  -> {tag}")

# ---- [3] forma della traiettoria: phi-phi_e ~ (r-r_e)^{3/2} -----------
print("\n[3] Delta_phi da r a r_e (|J|<J_c): deve ~ (r-r_e)^{3/2}")
Jv = 0.5
def integ(rv):
    return dphidr_n(rv, Jv)
for dr in [1e-2, 1e-3, 1e-4]:
    val = quad(integ, r_e + 1e-12, r_e + dr, limit=200)[0]
    print(f"    r-r_e={dr:.0e}: Dphi={val:.3e},  Dphi/(r-r_e)^1.5="
          f"{val/dr**1.5:.4f}  (deve tendere a 2K/3)")
K = Jv * E * r_e / (a**2 * np.sqrt(a**2 - Jv**2 * E**2) * np.sqrt(2 * M))
print(f"    2K/3 = {2*K/3:.4f}")

print("\nCONCLUSIONE:")
print("  |J|>J_c: periasse liscio sopra r_e (scatter).")
print("  |J|=J_c: dphi/dr finito a r_e -> attraversa liscio (separatrice).")
print("  |J|<J_c: dphi/dr ~ sqrt(r-r_e) -> phi-phi_e ~ (r-r_e)^{3/2}:")
print("           CUSPIDE, tangente radiale (dphi/dr->0). QED.")

# --------------------------------------------------------------- figura
print("\ngenerazione figura...")
r0 = 5.0
r_plus = M + np.sqrt(M**2 - a**2)
c_sep = a / E

def dphiBL_sep(rv):
    """separatrice J=J_c: sqrt(f) si CANCELLA -> regolare all'ergosfera."""
    Dlv = rv**2 - 2 * M * rv + a**2
    return Jc * np.sqrt(E**2 - (1 - 2 * M / rv)) * rv \
        / (Dlv * np.sqrt(rv**2 + c_sep**2))

def dshift_doran(rv):
    """shift di Doran a*beta/Dl; phi_D = phi_BL - shift regolare a r_+."""
    Dlv = rv**2 - 2 * M * rv + a**2
    return a * np.sqrt(2 * M * rv / (rv**2 + a**2)) / Dlv

def traiettoria_crossing():
    """J=J_c: PLUNGE monotono r0 -> r_e -> r_+, in coord. di Doran."""
    rg = np.linspace(r_plus + 0.02, r0, 800)
    # phi_D(r) = int_r^{r0} (dphi_BL/dr - dshift/dr) dr  (Doran, regolare)
    phiD = np.array([quad(lambda x: dphiBL_sep(x) - dshift_doran(x),
                          rv, r0, limit=400)[0] for rv in rg])
    phi_at_re = np.interp(r_e, rg, phiD)   # allinea l'attraversamento a phi=0
    return rg, phiD - phi_at_re

def r_turn(Jv):
    """Turning radius: r_e if |J| <= J_c, else the root of Dl - J^2 w.

    The previous version returned the grid point at which the sign change was
    detected, which lies just INSIDE the forbidden region (Dl - J^2 w < 0).
    Quadratures started there returned NaN, silently dropping the |J| > J_c
    curve from the log-log panel.  Bracket the sign change and solve.
    """
    if Jv <= Jc + 1e-9:
        return r_e
    g = lambda x: (x**2 - 2*M*x + a**2) - Jv**2*(E**2 - (1 - 2*M/x))
    rr = np.linspace(r_e + 1e-6, r0, 40000)
    gv = np.array([g(x) for x in rr])
    idx = np.where(gv[:-1]*gv[1:] < 0)[0]
    if not len(idx):
        return r_e
    return brentq(g, rr[idx[0]], rr[idx[0] + 1], xtol=1e-15, rtol=8.9e-16)

def traiettoria(Jv):
    rt = r_turn(Jv)
    rg = np.linspace(rt, r0, 500)
    # phin(r) = int_r^{r0} |dphi/dr'| dr'  (angolo da r a r0); phin[0]=Phi_half
    phin = np.array([quad(lambda x: dphidr_n(x, Jv), rv, r0,
                          limit=200)[0] for rv in rg])
    Phi_half = phin[0]
    # entrante r0->rt: phi 0->Phi_half ; uscente rt->r0: phi Phi_half->2Phi_half
    r_full = np.concatenate([rg[::-1], rg])
    phi_full = np.concatenate([phin[::-1], 2 * Phi_half - phin])
    return r_full, phi_full, rt, Phi_half

def traiettoria_retro_marginale():
    """J = -J_c: ingoing arc only.

    At the marginal value the shape integrand is |dphi/dr| = dphiBL_sep, which
    is FINITE at r_e; and the radial Hamilton equation gives dr/ds ~ C (r-r_e)
    with C = -sqrt(4M^2+J_c^2)/(8M^2), so the parameter diverges logarithmically
    and r_e is an asymptotic endpoint.  The extremal therefore neither cusps nor
    reflects, and only the ingoing branch exists.  (CQG-116884, major 4.)
    """
    rg = np.linspace(r_e, r0, 600)
    ph = np.array([-quad(dphiBL_sep, r_e, rv, limit=300)[0] for rv in rg])
    return rg, ph


fig, (ax, axb) = plt.subplots(2, 1, figsize=(COL, 6.6))
# cusp e periapsis: rimbalzo (BL); crossing (J_c): PLUNGE in coord. Doran
for Jv, col, lab, ls in [(0.5, '#E69F00', r'$|J|<J_c$: cusp at $r_e$', '-'),
                         (0.9, '#0072B2', r'$|J|>J_c$: smooth periapsis', ':')]:
    rf, phf, rt, phe = traiettoria(Jv)
    ax.plot(phf - phe, rf, col, ls=ls, lw=1.8, label=lab)
rg_m, ph_m = traiettoria_retro_marginale()
ax.plot(ph_m, rg_m, color='#CC79A7', ls='-.', lw=1.8,
        label=r'$J{=}{-}J_c$: asymptotic, not attained')
ax.plot(ph_m[0], rg_m[0], marker='o', ms=5, mfc='white', mec='#CC79A7',
        mew=1.4, zorder=6)
rg_x, ph_x = traiettoria_crossing()
inside = rg_x < r_e
ax.plot(ph_x[~inside], rg_x[~inside], color='#009E73', ls='-', lw=1.8,
        label=r'$J{=}{+}J_c$: crosses $r_e$')
ax.plot(ph_x[inside], rg_x[inside], color='#009E73', ls=(0, (3, 2)), lw=1.3,
        alpha=0.85,
        label=r'continuation, $r<2M$ (no optimality claimed)')
ax.axhline(r_e, color='k', ls='--', lw=0.9)
ax.axhline(r_plus, color='k', ls='-', lw=0.8)
ax.text(0.02, r_e + 0.03, 'ergosphere $r_e=2M$', fontsize=6)
ax.text(0.02, r_plus - 0.07, 'horizon $r_+$', fontsize=6)
ax.set_xlim(-0.3, 0.3)
ax.set_ylim(1.3, 2.55)
ax.set_xlabel('$\\varphi - \\varphi_e$  (Doran for the crossing)')
ax.set_ylabel('$r$')
ax.set_title('At the conformal stationary limit (Thakurta-Kerr, $s=0.9$, '
             f'$E_{{eff}}={E}$, $J_c={Jc:.2f}$):\nfour distinct local '
             r'behaviours; $\pm J_c$ are NOT equivalent')
ax.legend(fontsize=5.0, loc='lower left', framealpha=0.92)

# pannello b: r-r_turn vs |phi-phi_e| log-log (pendenza = potenza)
# the three exponents are the decisive discriminant:
#   2/3 cusp | 1 marginal (finite tangent) | 2 smooth periapsis
def _int_sing(g, rt, d):
    """Integral of g from rt to rt+d.

    At a SIMPLE root of the turning polynomial (the |J|>J_c periapsis) the
    integrand carries an inverse-square-root endpoint singularity, on which a
    plain adaptive quadrature returns NaN.  The substitution x = rt + u^2
    removes it exactly and is harmless in the regular cases.
    """
    return quad(lambda u: 2.0*u*g(rt + u*u), 0.0, np.sqrt(d), limit=300)[0]


dr = np.logspace(-4.5, -1.0, 45)
for Jv, col, ls, lab in [(0.5, '#E69F00', '-', 'cusp ($|J|<J_c$): slope 2/3'),
                         (0.9, '#0072B2', ':',
                          'periapsis ($|J|>J_c$): slope 2')]:
    rt = r_turn(Jv)
    dphi = np.array([_int_sing(lambda x: dphidr_n(x, Jv), rt, d) for d in dr])
    axb.loglog(dphi, dr, col, ls=ls, lw=2.0, label=lab)
# marginal |J| = J_c: dphi/dr is finite at r_e, hence slope exactly 1
dr = np.logspace(-4.5, -1.0, 45)
dphi_m = np.array([quad(dphiBL_sep, r_e, r_e + d, limit=300)[0] for d in dr])
axb.loglog(dphi_m, dr, color='#CC79A7', ls='-.', lw=2.0,
           label=r'marginal ($|J|{=}J_c$): slope 1')
xr = np.array([2e-4, 3e-2])
axb.loglog(xr, 1.35 * xr**(2 / 3), 'k:', lw=0.7)
axb.loglog(xr, 0.9 * xr**2, 'k:', lw=0.7)
axb.text(1e-2, 1.35 * 1e-2**(2/3) * 1.3, 'slope 2/3', fontsize=6)
axb.text(3e-3, 0.9 * 3e-3**2 * 0.3, 'slope 2', fontsize=6)
axb.set_xlabel('$|\\varphi-\\varphi_e|$')
axb.set_ylabel('$r-r_{turn}$')
axb.set_title('exponent of the turn: $2/3$ cusp (radial tangent),\n'
              '$1$ marginal (finite tangent), $2$ smooth parabola')
axb.legend(fontsize=6, loc='lower right')

savefig(fig, OUT, 'fig_thakurta_cuspide_ergosfera')

# --- measured exponents, printed so the figure is self-certifying ----------
# r - r_turn ~ |phi - phi_e|^p :  p = 2/3 cusp, 1 marginal, 2 periapsis.
def _slope(x, y, n=12):
    x, y = np.asarray(x, float), np.asarray(y, float)
    ok = np.isfinite(x) & np.isfinite(y) & (x > 0) & (y > 0)
    x, y = x[ok][:n], y[ok][:n]
    if len(x) < 3:
        return float('nan')
    return np.polyfit(np.log(x), np.log(y), 1)[0]


print("\nmeasured log-log exponents (smallest decade):")
for Jv, nm in [(0.5, '|J|<Jc  (cusp)'), (0.9, '|J|>Jc  (periapsis)')]:
    rt = r_turn(Jv)
    dphi_c = np.array([_int_sing(lambda x: dphidr_n(x, Jv), rt, d)
                       for d in dr])
    print(f"  {nm:22s} p = {_slope(dphi_c, dr):.4f}")
print(f"  {'|J|=Jc  (marginal)':22s} p = {_slope(dphi_m, dr):.4f}")
print("FATTO.")
