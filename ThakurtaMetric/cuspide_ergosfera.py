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

Tricotomia analitica all'ergosfera:
  |J| > J_c :  a^2 - J^2 E^2 < 0  -> svolta PRIMA di r_e (periasse liscio,
               dphi/dr -> infinito, potenza 1/2 classica)
  |J| = J_c :  a^2 - J_c^2 E^2 = 0 e f->0 insieme -> dphi/dr -> FINITO
               (attraversamento liscio: separatrice)
  |J| < J_c :  a^2 - J^2 E^2 > 0 finito, f->0 -> dphi/dr ~ sqrt(r-r_e)
               -> CUSPIDE (potenza 3/2), tangente radiale
"""

import os
import sys
import numpy as np
import sympy as sp
from scipy.integrate import quad

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
    """raggio di svolta: r_e se |J|<=J_c, altrimenti radice di Dl-J^2 w."""
    if Jv <= Jc + 1e-9:
        return r_e
    rr = np.linspace(r_e + 1e-6, r0, 40000)
    g = np.array([(x**2 - 2 * M * x + a**2)
                  - Jv**2 * (E**2 - (1 - 2 * M / x)) for x in rr])
    idx = np.where(g[:-1] * g[1:] < 0)[0]
    return rr[idx[0]] if len(idx) else r_e

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

fig, (ax, axb) = plt.subplots(2, 1, figsize=(COL, 6.6))
# cusp e periapsis: rimbalzo (BL); crossing (J_c): PLUNGE in coord. Doran
for Jv, col, lab, ls in [(0.5, 'C0', r'$|J|<J_c$: CUSP (bounce)', '-'),
                         (0.9, 'C3', r'$|J|>J_c$: smooth periapsis', ':')]:
    rf, phf, rt, phe = traiettoria(Jv)
    ax.plot(phf - phe, rf, col, ls=ls, lw=1.8, label=lab)
rg_x, ph_x = traiettoria_crossing()
ax.plot(ph_x, rg_x, 'C2--', lw=1.8,
        label=r'$|J|=J_c$: crossing (Doran) $\to r_+$')
ax.axhline(r_e, color='k', ls='--', lw=0.9)
ax.axhline(r_plus, color='k', ls='-', lw=0.8)
ax.text(0.02, r_e + 0.03, 'ergosphere $r_e=2M$', fontsize=6)
ax.text(0.02, r_plus - 0.07, 'horizon $r_+$', fontsize=6)
ax.set_xlim(-0.3, 0.3)
ax.set_ylim(1.3, 2.55)
ax.set_xlabel('$\\varphi - \\varphi_e$  (Doran for the crossing)')
ax.set_ylabel('$r$')
ax.set_title('At the ergosphere (Thakurta-Kerr, $s=0.9$, $E_{eff}$='
             f'{E}, $J_c={Jc:.2f}$):\nCUSP bounce / smooth periapsis / '
             '$J_c$ CROSSES to $r_+$ (Doran)')
ax.legend(fontsize=5.5, loc='upper right')

# pannello b: r-r_turn vs |phi-phi_e| log-log (pendenza = potenza)
for Jv, col, ls, lab in [(0.5, 'C0', '-', 'cusp ($|J|<J_c$): slope 2/3'),
                         (0.9, 'C3', ':', 'periapsis ($|J|>J_c$): slope 2')]:
    rt = r_turn(Jv)
    dr = np.logspace(-4.5, -1.0, 45)
    dphi = np.array([quad(lambda x: dphidr_n(x, Jv), rt + 1e-13, rt + d,
                          limit=300)[0] for d in dr])
    axb.loglog(dphi, dr, col, ls=ls, lw=2.0, label=lab)
xr = np.array([2e-4, 3e-2])
axb.loglog(xr, 1.35 * xr**(2 / 3), 'k:', lw=0.7)
axb.loglog(xr, 0.9 * xr**2, 'k:', lw=0.7)
axb.text(1e-2, 1.35 * 1e-2**(2/3) * 1.3, 'slope 2/3', fontsize=6)
axb.text(3e-3, 0.9 * 3e-3**2 * 0.3, 'slope 2', fontsize=6)
axb.set_xlabel('$|\\varphi-\\varphi_e|$')
axb.set_ylabel('$r-r_{turn}$')
axb.set_title('power law of the turn: $2/3$ = cusp (radial tangent),\n'
              '$2$ = smooth parabola')
axb.legend(fontsize=6, loc='lower right')
savefig(fig, OUT, 'fig_thakurta_cuspide_ergosfera')
print("FATTO.")
