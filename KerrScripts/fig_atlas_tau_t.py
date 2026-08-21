# -*- coding: utf-8 -*-
"""
Atlante delle brachistocrone all'ergosfera: piccoli multipli, 5 pannelli per
ramo, uno per regime rispetto ai valori critici.

Ramo tau (critici +-J_c, J_c=a/E):
  J>+J_c  smooth bounce | J=+J_c penetrates | -J_c<J<+J_c cusp bounce |
  J=-J_c cusp (retro) | J<-J_c smooth bounce
Ramo t (critici J_c^-, J_c^+ = -8.05, 3.43):
  J>J_c^+ smooth bounce | J=J_c^+ tangent | J_c^-<J<J_c^+ penetrates |
  J=J_c^- spirals onto r_* | J<J_c^- smooth bounce
"""
import os
import sys
import numpy as np
import sympy as sp
from scipy.integrate import quad, solve_ivp
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
M, a, E = 1.0, 0.9, 1.2
re = 2*M
rp = M + np.sqrt(M**2 - a**2)
rstar = 3.5139
Jc = a/E
f = lambda r: 1 - 2*M/r
w = lambda r: E**2 - f(r)
Dl = lambda r: r**2 - 2*M*r + a**2
r0 = 4.0


def Q2(r, J):
    return (2*E**2*J**2*M*r - E**2*J**2*r**2 - 4*E**2*J*M*a*r + 2*E**2*M*a**2*r
            + E**2*a**2*r**2 + E**2*r**4 + 4*J**2*M**2 - 4*J**2*M*r + J**2*r**2
            - 8*J*M**2*a + 4*J*M*a*r + 4*M**2*a**2)


def R6(r, J):
    return r*Q2(r, J)*((E**2 - 1)*r + 2*M)


def clip_wind(rr, ph, phimax=1.9*np.pi):
    """taglia l'arco quando |phi-phi0| supera phimax (leggibilita')."""
    d = np.abs(ph - ph[0])
    k = np.argmax(d > phimax) if np.any(d > phimax) else len(ph)
    return rr[:k], ph[:k]


def orbit_reflect(dphidr, r_stop, retro=False):
    rin = np.linspace(r0, r_stop, 500)
    s = -1.0 if retro else 1.0
    phin = np.array([s*quad(dphidr, r, r0, limit=400)[0] for r in rin])
    rout = np.linspace(r_stop, r0, 500)
    phout = 2*phin[-1] - np.array([s*quad(dphidr, r, r0, limit=400)[0]
                                   for r in rout])
    return np.concatenate([rin, rout]), np.concatenate([phin, phout])


def orbit_asymptote(dphidr, r_end, retro=False):
    """Ingoing arc only.

    At J = -Jc the radial momentum has a DOUBLE root at r_e while the approach
    rate T(r_e) = -a(J-Jc)/Pbar_e does not degenerate, so dr/ds ~ C (r - r_e)
    with C = -sqrt(4M^2+Jc^2)/(8M^2): the parameter diverges logarithmically and
    r_e is an asymptotic endpoint, never attained.  The extremal therefore does
    NOT reflect, and drawing a mirrored outgoing branch (as the previous version
    did) is wrong.  dphi/dr stays finite there, so the arc terminates at a
    definite angle with a definite tangent.  See CQG-116884 major comment 4 and
    NonStationaryMetrics/paper2/verification/verify_marginal_hamilton.wls.
    """
    rr = np.linspace(r0, r_end, 700)
    s = -1.0 if retro else 1.0
    ph = np.array([s*quad(dphidr, r, r0, limit=400)[0] for r in rr])
    return rr, ph


def orbit_through(dphidr, r_end, retro=False):
    rr = np.concatenate([np.linspace(r0, re, 250), np.linspace(re, r_end, 800)])
    s = -1.0 if retro else 1.0
    ph = np.array([s*quad(dphidr, r, r0, limit=500)[0] for r in rr])
    return rr, ph


# --- tau ---
def dphidr_tau(r, J):
    return abs(J)*r*np.sqrt(w(r)*f(r))/(Dl(r)*np.sqrt(Dl(r) - J**2*w(r)))


g_cancel = lambda r: a*np.sqrt(r**2*w(r))/(Dl(r)*np.sqrt(E**2*r**2 + a**2))


def tau_turn(J):
    if Dl(re) - J**2*w(re) < 0:
        return brentq(lambda r: Dl(r) - J**2*w(r), re + 1e-9, r0)
    return re


# --- t: flusso di Hamilton (adattivo, solve_ivp) ---
_r, _pr, _pp = sp.symbols('r pr pp')
_f = 1 - 2*M/_r
_Dl = _r**2 - 2*M*_r + a**2
_b = 2*M*a/_r
_v2 = 1 - _f/E**2
_P = _r**2 + a**2 + 2*M*a**2/_r
_Pb = _P + _b**2/E**2
_H = (_pp*_b*_v2/_Pb + sp.sqrt(_Dl*_v2/_Pb)
      * sp.sqrt((_Dl/_r**2)*_pr**2 + _pp**2/_Pb) - 1)
_Hn = sp.lambdify((_r, _pr, _pp), _H, 'numpy')
_dHp = sp.lambdify((_r, _pr, _pp), sp.diff(_H, _pr), 'numpy')
_dHr = sp.lambdify((_r, _pr, _pp), sp.diff(_H, _r), 'numpy')
_dHj = sp.lambdify((_r, _pr, _pp), sp.diff(_H, _pp), 'numpy')


def _prof(rv, J):
    pg = np.linspace(-200, 200, 20001)
    Hv = _Hn(rv, pg, J)
    rts = [brentq(lambda p: _Hn(rv, p, J), pg[i], pg[i+1])
           for i in range(len(pg)-1)
           if np.isfinite(Hv[i]) and np.isfinite(Hv[i+1]) and Hv[i]*Hv[i+1] < 0]
    ing = [p for p in rts if _dHp(rv, p, J) < 0]
    return min(ing) if ing else np.nan


rt0 = 6.0   # start radius per il ramo t (fuori da tutti i turning, max 4.77)


def t_orbit_ivp(J, lam_max=400.0, r_floor=None):
    """orbita t per flusso di Hamilton adattivo; ferma a r_+ o a lam_max."""
    if r_floor is None:
        r_floor = rp + 0.02
    p0 = _prof(rt0, J)
    ev_in = lambda t, y: y[0] - r_floor        # stop vicino orizzonte/r_*
    ev_in.terminal = True
    ev_in.direction = -1
    ev_out = lambda t, y: y[0] - (rt0 + 0.01)   # stop se riesce oltre rt0
    ev_out.terminal = True
    ev_out.direction = 1
    s = solve_ivp(lambda t, y: [_dHp(y[0], y[1], J), -_dHr(y[0], y[1], J),
                                _dHj(y[0], y[1], J)],
                  [0, lam_max], [rt0, p0, 0.0], rtol=1e-10, atol=1e-12,
                  max_step=0.02, dense_output=True, events=[ev_in, ev_out])
    return s.y[0], s.y[2]


Jtp = 2*M**2/a + a + a/(2*E**2)
Jtm = -8.0535


def tau_orbit(J):
    if abs(J - Jc) < 1e-9:
        return orbit_through(g_cancel, rp + 0.01)
    if abs(J + Jc) < 1e-9:
        return orbit_asymptote(g_cancel, re, retro=True)
    return orbit_reflect(lambda r: dphidr_tau(r, J), tau_turn(J), retro=(J < 0))


def t_orbit(J):
    if abs(J - 0.99*Jtm) < 1e-3:                # cattura retrograda: plunge a r_+
        rr, ph = t_orbit_ivp(J, lam_max=6000, r_floor=rp + 0.02)
        return clip_wind(rr, ph, phimax=4.0*np.pi)
    rr, ph = t_orbit_ivp(J, lam_max=600)
    return clip_wind(rr, ph, phimax=2.2*np.pi)


tau_cases = [(1.50, r'$J{=}1.5>J_c$', 'smooth periapsis'),
             (Jc,   r'$J{=}{+}J_c$', 'CROSSES'),
             (0.30, r'$-J_c{<}J{<}{+}J_c$', 'cusp at $r_e$'),
             (-Jc,  r'$J{=}{-}J_c$', 'asymptotic: not attained'),
             (-1.50, r'$J{=}{-}1.5<{-}J_c$', 'smooth periapsis')]
t_cases = [(5.0, r'$J{=}5>J_c^+$', 'smooth bounce'),
           (Jtp, r'$J{=}J_c^+$', 'tangent'),
           (3.0, r'$J_c^-{<}J{<}J_c^+$', 'PENETRATES'),
           (0.99*Jtm, r'$J{=}0.99\,J_c^-$', 'captured (to $r_+$)'),
           (-9.0, r'$J{=}{-}9<J_c^-$', 'smooth bounce')]


def panel(ax, rr, ph, col, ttl, sub, L, rstart, kind='plain'):
    """kind: 'plain' | 'asymptote' | 'cross'.

    'cross'     -- the segment at r < 2M is an ANALYTIC CONTINUATION of the
                   frozen equations, not a controlled-rail extremal: there
                   g(W,W) > 0 for W = d/d(eta), the admissible-velocity set is
                   not compact, and no optimality is claimed.  Drawn dashed
                   (Protocol 2.1 type (A), referee minor comment M6).
    'asymptote' -- r_e is approached but never attained; mark the open endpoint.
    """
    th = np.linspace(0, 2*np.pi, 300)
    ax.fill(re*np.cos(th), re*np.sin(th), color='0.92', zorder=0)
    ax.plot(re*np.cos(th), re*np.sin(th), 'b--', lw=0.7)
    ax.plot(rp*np.cos(th), rp*np.sin(th), 'k--', lw=0.6)
    if rr is not None:
        x, y = rr*np.cos(ph), rr*np.sin(ph)
        if kind == 'cross':
            out = rr >= re
            ax.plot(x[out], y[out], col, lw=1.5)
            ax.plot(x[~out], y[~out], col, lw=1.1, ls=(0, (3, 2)), alpha=0.85)
        else:
            ax.plot(x, y, col, lw=1.5)
        if kind == 'asymptote':
            ax.plot(x[-1], y[-1], marker='o', ms=3.6, mfc='white',
                    mec=col, mew=0.9, zorder=6)
        ax.plot(rstart, 0, 'ks', ms=2.5)
    ax.set_aspect('equal'); ax.set_xticks([]); ax.set_yticks([])
    ax.set_title(ttl + '\n' + sub, fontsize=5.6, pad=1.5)
    ax.set_xlim(-L, L); ax.set_ylim(-L, L)


def make_fig(branch):
    cases = tau_cases if branch == 'tau' else t_cases
    orbit = tau_orbit if branch == 'tau' else t_orbit
    L = 4.6 if branch == 'tau' else 6.5
    rstart = r0 if branch == 'tau' else rt0
    cols = ['C0', 'C2', 'C1', 'C3', 'C4']
    with plt.rc_context({'figure.constrained_layout.use': False}):
        fig, axes = plt.subplots(2, 3, figsize=(COL, COL*0.80))
    axes = axes.ravel()
    for i, (J, ttl, sub) in enumerate(cases):
        rr, ph = orbit(J)
        if branch == 'tau' and abs(J + Jc) < 1e-9:
            kind = 'asymptote'
        elif rr is not None and np.min(rr) < re - 1e-6:
            kind = 'cross'
        else:
            kind = 'plain'
        panel(axes[i], rr, ph, cols[i], ttl, sub, L, rstart, kind=kind)
    # ultimo pannello: didascalia riassuntiva del ramo
    axes[5].axis('off')
    lab = (r'$\bf \tau$-branch''\n'r'(grey: $r<2M$, $W$ spacelike)''\n\n'
           'FOUR CASES:\ncrosses only\n'r'at $J{=}{+}J_c{=}0.75$'
           '\n'r'($-J_c$: asymptotic)''\n\n'
           'dashed inside $2M$:\ncontinuation,\nno optimality claimed'
           if branch == 'tau' else
           r'$\bf t$-branch''\n'r'(grey: $r<2M$, $W$ spacelike)''\n\n'
           'DICHOTOMY:\ncrosses all\n'r'$(J_c^-,J_c^+)$'
           '\n'r'$=(-8.05,3.43)$''\n\n'
           'dashed inside $2M$:\ncontinuation,\nno optimality claimed')
    axes[5].text(0.5, 0.5, lab, ha='center', va='center', fontsize=6.4,
                 transform=axes[5].transAxes)
    fig.subplots_adjust(top=0.90, bottom=0.02, left=0.02, right=0.98,
                        hspace=0.40, wspace=0.12)
    savefig(fig, HERE, 'fig_atlas_' + branch)
    print(f'  {branch} done')


make_fig('tau')
make_fig('t')
print('FATTO.')
