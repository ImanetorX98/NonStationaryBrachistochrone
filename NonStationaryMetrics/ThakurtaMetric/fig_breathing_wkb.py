# -*- coding: utf-8 -*-
"""
WKB adiabatico: la famiglia "che respira" delle brachistocrone Thakurta-Kerr e
la traiettoria non-autonoma che la infila.

Rail energy Ehat FISSO; fattore conforme A che cresce lungo l'orbita (universo
in espansione, eta aumenta) => E_eff = Ehat/A decresce. La forma chiusa
Kerr(E_eff) e' l'ordine adiabatico DOMINANTE (famiglia congelata). La
traiettoria non-autonoma vera (E che scorre come A(lambda)) HUGGA la famiglia,
con errore adiabatico O(A'/A).

Pannello (a): famiglia congelata (grigia) + traiettoria non-autonoma (colorata).
Pannello (b): errore adiabatico (scostamento della r_min vera dalla predizione
  congelata istantanea) vs tasso A'/A -> lineare = O(A'/A) confermato.
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
from matplotlib import cm

set_style()
HERE = os.path.dirname(os.path.abspath(__file__))
M, a = 1.0, 0.9
Ehat = 1.4
J = 6.0
r0 = 12.0

r, pr, Es, Js = sp.symbols('r pr E J_')
f = 1 - 2*M/r
Dl = r**2 - 2*M*r + a**2
b = 2*M*a/r
v2 = 1 - f/Es**2
P = r**2 + a**2 + 2*M*a**2/r
Pb = P + b**2/Es**2
H = Js*b*v2/Pb + sp.sqrt(Dl*v2/Pb)*sp.sqrt((Dl/r**2)*pr**2 + Js**2/Pb) - 1
Hn = sp.lambdify((r, pr, Es, Js), H, 'numpy')
dHp = sp.lambdify((r, pr, Es, Js), sp.diff(H, pr), 'numpy')
dHr = sp.lambdify((r, pr, Es, Js), sp.diff(H, r), 'numpy')
dHj = sp.lambdify((r, pr, Es, Js), sp.diff(H, Js), 'numpy')


def prof(rv, E):
    pg = np.linspace(-80, 80, 1601)
    Hv = Hn(rv, pg, E, J)
    rts = [brentq(lambda p: Hn(rv, p, E, J), pg[i], pg[i+1])
           for i in range(len(pg)-1)
           if np.isfinite(Hv[i]) and np.isfinite(Hv[i+1]) and Hv[i]*Hv[i+1] < 0]
    ing = [p for p in rts if dHp(rv, p, E, J) < 0]
    return min(ing) if ing else np.nan


def _turn_event(t, y):
    return y[1]           # p_r = 0
_turn_event.terminal = True
_turn_event.direction = 1


def frozen_rmin(E):
    p0 = prof(r0, E)
    if not np.isfinite(p0):
        return np.nan
    s = solve_ivp(lambda t, y: [dHp(y[0], y[1], E, J), -dHr(y[0], y[1], E, J)],
                  [0, 120], [r0, p0], rtol=1e-9, atol=1e-11,
                  max_step=0.1, dense_output=True, events=_turn_event)
    return s.y[0][-1] if len(s.t_events[0]) else np.nan


def frozen_orbit(E):
    """incoming a turning, poi mirror."""
    p0 = prof(r0, E)
    s = solve_ivp(lambda t, y: [dHp(y[0], y[1], E, J), -dHr(y[0], y[1], E, J),
                                dHj(y[0], y[1], E, J)],
                  [0, 120], [r0, p0, 0.0], rtol=1e-9, atol=1e-11, max_step=0.1,
                  dense_output=True, events=_turn_event)
    te = s.t_events[0][0] if len(s.t_events[0]) else s.t[-1]
    tt = np.linspace(0, te, 400)
    rr = s.sol(tt)[0]; ph = s.sol(tt)[2]
    # mirror (outgoing)
    rr2 = rr[::-1]; ph2 = 2*ph[-1] - ph[::-1]
    return np.concatenate([rr, rr2]), np.concatenate([ph, ph2])


def nonauto(eps, want_orbit=False):
    """A(lambda)=1+eps*lambda; ferma alla svolta. Ritorna (r_turn, A_turn)
       o l'orbita completa (mirror)."""
    Ael = lambda lam: 1.0 + eps*lam
    p0 = prof(r0, Ehat/Ael(0.0))
    def rhs(lam, y):
        e = Ehat/Ael(lam)
        return [dHp(y[0], y[1], e, J), -dHr(y[0], y[1], e, J),
                dHj(y[0], y[1], e, J)]
    s = solve_ivp(rhs, [0, 120], [r0, p0, 0.0], rtol=1e-9, atol=1e-11,
                  max_step=0.1, dense_output=True, events=_turn_event)
    lam_t = s.t_events[0][0] if len(s.t_events[0]) else s.t[-1]
    r_turn = s.y[0][-1]; A_turn = Ael(lam_t)
    if not want_orbit:
        return r_turn, A_turn
    tt = np.linspace(0, lam_t, 400)
    rr = s.sol(tt)[0]; ph = s.sol(tt)[2]
    rr2 = rr[::-1]; ph2 = 2*ph[-1] - ph[::-1]
    return np.concatenate([rr, rr2]), np.concatenate([ph, ph2]), r_turn, A_turn


import time
try:
    from tqdm import tqdm
except ImportError:
    def tqdm(x, **k): return x
_t0 = time.time()
def log(m): print(f"[{time.time()-_t0:6.1f}s] {m}", flush=True)

# --- famiglia congelata (sfondo) ---
As = np.array([1.0, 1.1, 1.2, 1.3, 1.35])
log("famiglia congelata (orbite t per scala di A)...")
fam = [frozen_orbit(Ehat/A) for A in tqdm(As, desc="frozen family")]
log("  fatta.")

# --- traiettoria non-autonoma adiabatica (eps piccolo) ---
log("traiettoria non-autonoma eps=0.02...")
rr_na, ph_na, rturn_na, Aturn_na = nonauto(0.02, want_orbit=True)
log(f"  turning r={rturn_na:.4f} a A={Aturn_na:.4f}; "
    f"congelato r_min(E_eff={Ehat/Aturn_na:.3f})={frozen_rmin(Ehat/Aturn_na):.4f}")

# --- errore adiabatico vs eps (A'/A) ---
log("scan errore adiabatico vs eps...")
epss_all = [0.004, 0.008, 0.015, 0.03, 0.05]
epss, errs = [], []
for e in tqdm(epss_all, desc="eps scan"):
    rt, At = nonauto(e)
    pred = frozen_rmin(Ehat/At)
    if np.isfinite(rt) and np.isfinite(pred):
        epss.append(e); errs.append(abs(rt - pred))
        log(f"  eps={e:.3f}: r_turn={rt:.4f}, congelato={pred:.4f}, "
            f"err={errs[-1]:.2e}")
    else:
        log(f"  eps={e:.3f}: skip (nan, E_eff fuori range)")
epss = np.array(epss); errs = np.array(errs)

# ---------------------------------------------------------------- figura
fig, (ax, axb) = plt.subplots(1, 2, figsize=(2*COL, COL*1.0))
th = np.linspace(0, 2*np.pi, 200)
for (rr, ph), A in zip(fam, As):
    ax.plot(rr*np.cos(ph), rr*np.sin(ph), color='0.7', lw=1.0)
ax.plot([], [], color='0.7', lw=1.0, label='frozen family (ladder of $A$)')
sc = ax.scatter(rr_na*np.cos(ph_na), rr_na*np.sin(ph_na),
                c=1.0+0.02*np.linspace(0, len(rr_na), len(rr_na)),
                cmap='viridis', s=2, zorder=5)
ax.plot(rr_na*np.cos(ph_na), rr_na*np.sin(ph_na), 'C2-', lw=0.8, alpha=0.5,
        label='non-autonomous ($A$ grows)')
ax.plot(2*M*np.cos(th), 2*M*np.sin(th), 'b--', lw=0.7, label='$r_e$')
ax.set_aspect('equal'); ax.set_xlabel('$x$'); ax.set_ylabel('$y$')
ax.set_title('breathing family + non-autonomous orbit\n'
             '(threads through the frozen family)', fontsize=6.8)
ax.legend(fontsize=5.4, loc='upper right', framealpha=0.9)

axb.loglog(epss, errs, 'o-', color='C3', label='adiabatic error')
k = len(epss)//2
axb.loglog(epss, epss*errs[k]/epss[k], 'k--', lw=0.8,
           label=r'$\propto \epsilon$ (i.e. $A^\prime/A$)')
axb.set_xlabel(r"rate $\epsilon \sim A'/A$")
axb.set_ylabel(r'$|r_{\min}^{\rm actual}-r_{\min}^{\rm frozen}|$')
axb.set_title('WKB error is $O(A^\\prime/A)$\n(adiabatic theorem)', fontsize=6.8)
axb.legend(fontsize=6, loc='upper left')
savefig(fig, HERE, 'fig_breathing_wkb')
print('FATTO.')
