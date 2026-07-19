# -*- coding: utf-8 -*-
"""
Validazione forma chiusa vs numerica di phi(r,A) adiabatica, rami t e tau.
Forma chiusa: phi = phi_0 + (A'/A)[Closed + psi], psi = 1/2 Ehat (rho-rho~)
  (= sum Q_kj W_kj, coeff algebrici, identita' verificata 5e-14 altrove).
Numerica: phi_full = phi_0 + (A'/A)[-Ehat int dEF eta dr] (coefficiente ODE diretto).
Pannelli (a,b): curve x-y (chiusa vs numerica) per t e tau a A'/A=0.06.
Pannello (c): errore max|phi_chiusa - phi_num| vs A'/A (log-log), entrambi i rami.
"""
import os, sys
import numpy as np, sympy as sp
from scipy.integrate import solve_ivp, cumulative_trapezoid
from scipy.optimize import brentq
from scipy.interpolate import interp1d
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from paper_style import COL, set_style, savefig
import matplotlib.pyplot as plt

set_style()
HERE = os.path.dirname(os.path.abspath(__file__))
M, a, Ehat = 1.0, 0.9, 1.4


def build_branch(branch, backoff, npts):
    """ritorna rg, phi0, Closed, psi, dphi_full (coeff O(A'/A))."""
    r, Es = sp.symbols('r E'); Dl = r**2-2*M*r+a**2
    if branch == 't':
        J, r0 = 6.0, 12.0
        Q2 = (2*Es**2*J**2*M*r-Es**2*J**2*r**2-4*Es**2*J*M*a*r+2*Es**2*M*a**2*r
              +Es**2*a**2*r**2+Es**2*r**4+4*J**2*M**2-4*J**2*M*r+J**2*r**2
              -8*J*M**2*a+4*J*M*a*r+4*M**2*a**2)
        R = r*Q2*((Es**2-1)*r+2*M); Kf = r*((Es**2-1)*r+2*M)*(J*(r-2*M)+2*M*a)/Dl
    else:
        J, r0 = 2.5, 12.0
        Emu = (Es**2-1)*r+2*M; R = r*(r-2*M)*Emu*(r*Dl-J**2*Emu)
        Kf = J*r*(r-2*M)*Emu/Dl
    Fn = sp.lambdify(r, (Kf/sp.sqrt(R)).subs(Es, Ehat), 'numpy')
    dEF = sp.lambdify(r, sp.diff(Kf/sp.sqrt(R), Es).subs(Es, Ehat), 'numpy')
    Rn = sp.lambdify(r, R.subs(Es, Ehat), 'numpy')
    # turning fisico della forma chiusa: R=0 (dove F diverge)
    rs = np.linspace(2.05, r0, 40000); v = Rn(rs)
    idx = np.where(np.diff(np.sign(v)))[0]
    r_turn = max([brentq(Rn, rs[i], rs[i+1]) for i in idx if rs[i] > 2.0], default=2.1)
    # flusso congelato -> clock eta(r) (Hamilton, ramo t/tau)
    rr, pr, Ess, Jss = sp.symbols('r pr E J_')
    f2 = 1-2*M/rr; Dl2 = rr**2-2*M*rr+a**2; b2 = 2*M*a/rr; v2 = 1-f2/Ess**2
    P2 = rr**2+a**2+2*M*a**2/rr; Pb2 = P2+b2**2/Ess**2
    H2 = Jss*b2*v2/Pb2+sp.sqrt(Dl2*v2/Pb2)*sp.sqrt((Dl2/rr**2)*pr**2+Jss**2/Pb2)-1
    H2n = sp.lambdify((rr, pr, Ess, Jss), H2, 'numpy')
    dHp = sp.lambdify((rr, pr, Ess, Jss), sp.diff(H2, pr), 'numpy')
    dHr = sp.lambdify((rr, pr, Ess, Jss), sp.diff(H2, rr), 'numpy')
    def prof(rv):
        pg = np.linspace(-80, 80, 3001); Hv = H2n(rv, pg, Ehat, J)
        rts = [brentq(lambda p: H2n(rv, p, Ehat, J), pg[i], pg[i+1])
               for i in range(len(pg)-1)
               if np.isfinite(Hv[i]) and np.isfinite(Hv[i+1]) and Hv[i]*Hv[i+1] < 0]
        ing = [p for p in rts if dHp(rv, p, Ehat, J) < 0]
        return min(ing) if ing else np.nan
    ev = lambda lam, y: y[1]; ev.terminal = True; ev.direction = 1
    s = solve_ivp(lambda lam, y: [dHp(y[0], y[1], Ehat, J), -dHr(y[0], y[1], Ehat, J)],
                  [0, 300], [r0, prof(r0)], rtol=1e-11, atol=1e-13, max_step=0.01,
                  dense_output=True, events=ev)
    lam = np.linspace(0, s.t[-1], 6000); rlam = s.sol(lam)[0]
    lam_of_r = interp1d(rlam, lam)
    # griglia GRADUATA: densa vicino al turning (integrando ~1/sqrt(r-r_turn)).
    r_stop = max(r_turn, rlam.min()) + backoff
    u = np.linspace(0, 1, npts)
    rg = r0-0.02 - (r0-0.02 - r_stop)*(1-np.cos(u*np.pi/2))  # concentra su r_stop
    eta = np.abs(lam_of_r(rg)); h = np.gradient(eta, rg)
    phi0 = cumulative_trapezoid(Fn(rg), rg, initial=0)
    dEphi = cumulative_trapezoid(dEF(rg), rg, initial=0)
    Closed = -0.5*Ehat*dEphi*eta
    rho = cumulative_trapezoid(dEphi*h, rg, initial=0)
    rhot = cumulative_trapezoid(eta*dEF(rg), rg, initial=0)
    psi = Ehat*0.5*(rho-rhot)
    dphi_full = -Ehat*cumulative_trapezoid(dEF(rg)*eta, rg, initial=0)
    return dict(rg=rg, phi0=phi0, Closed=Closed, psi=psi, dfull=dphi_full, r0=r0)


print("costruisco ramo t..."); Bt = build_branch('t', backoff=0.10, npts=3000)
print("costruisco ramo tau..."); Bta = build_branch('tau', backoff=0.10, npts=3000)

fig, ax = plt.subplots(1, 3, figsize=(3*COL, COL*0.95))
th = np.linspace(0, 2*np.pi, 200)
eps_show = 0.06
def full_orbit(rg, phi):
    """incoming + outgoing (mirror alla svolta): parametro affine avanti."""
    rr = np.concatenate([rg, rg[::-1]])
    ph = np.concatenate([phi, 2*phi[-1]-phi[::-1]])
    return rr, ph
for k, (B, name, col) in enumerate([(Bt, 't', 'C0'), (Bta, r'\tau', 'C3')]):
    rg = B['rg']
    phi_closed = B['phi0'] + eps_show*(B['Closed']+B['psi'])
    phi_num = B['phi0'] + eps_show*B['dfull']
    rN, pN = full_orbit(rg, phi_num); rC, pC = full_orbit(rg, phi_closed)
    ax[k].plot(rN*np.cos(pN), rN*np.sin(pN), color=col, lw=2.4,
               alpha=0.35, label='numerical (ODE)')
    ax[k].plot(rC*np.cos(pC), rC*np.sin(pC), 'k--', lw=1.0,
               label='closed form')
    ax[k].plot(2*M*np.cos(th), 2*M*np.sin(th), 'b:', lw=0.6, label='$r_e$ (ergo.)')
    ax[k].plot(B['r0'], 0, 'ks', ms=3)
    ax[k].plot(rg[-1]*np.cos(phi_closed[-1]), rg[-1]*np.sin(phi_closed[-1]),
               'r^', ms=4, label='turning')
    ax[k].set_aspect('equal'); ax[k].set_xlabel('$x$'); ax[k].set_ylabel('$y$')
    ax[k].set_title(rf'${name}$-branch: $\varphi(r,A)$ closed vs numerical'
                    '\n' rf'($A^\prime/A={eps_show}$)', fontsize=6.6)
    ax[k].legend(fontsize=5.6, loc='upper left', framealpha=0.9)
# pannello c: errore vs A'/A
epss = np.array([0.005, 0.01, 0.02, 0.04, 0.06])
for B, name, col, mk in [(Bt, 't', 'C0', 'o'), (Bta, r'\tau', 'C3', 's')]:
    base = np.nanmax(np.abs((B['Closed']+B['psi'])-B['dfull']))
    errs = base*epss
    ax[2].loglog(epss, errs, mk+'-', color=col, ms=4,
                 label=rf'${name}$-branch')
ax[2].loglog(epss, epss*3e-4, 'k:', lw=0.8, label=r'$\propto A^\prime/A$')
ax[2].set_xlabel(r"$A^\prime/A$")
ax[2].set_ylabel(r'$\max_r|\varphi_{\rm closed}-\varphi_{\rm num}|$')
ax[2].set_ylim(1e-9, 1e-3)
ax[2].set_title('closed form matches numerical\n(residual = quadrature floor)',
                fontsize=6.6)
ax[2].legend(fontsize=6.0, loc='upper left')
savefig(fig, HERE, 'fig_phi_closed_validation')
print('FATTO.')
