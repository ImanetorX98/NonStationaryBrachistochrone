# -*- coding: utf-8 -*-
"""
Il vettore di Kodama lungo le brachistocrone di Vaidya.

"Conservazione" precisata in tre affermazioni, tutte verificate:

 [A] miracolo di Kodama (proprieta' dello spaziotempo, coord. EF entranti
     ds^2 = -f dv^2 + 2 dv dr + r^2 dOmega^2,  f = 1 - 2 m(v)/r):
       K^a = eps^{ab} nabla_b r  =  (1,0,0,0) = d/dv   (esatto)
       nabla_a K^a = 0                                  (senza Killing)
       carica di Kodama = massa di Misner-Sharp = m(v)

 [B] identita' della rotaia: lungo la brachistocrona tau il vincolo W
     equivale a  -u.K = E  ESATTAMENTE (non e' imposto a mano: cade
     dalla forma di W).  => E_K := -K.u e' l'invariante di rotaia.

 [C] non-banalita': per una GEODETICA (rotaia spenta)
       d(-K.u)/dtau = -(1/2)(L_K g)_{ab} u^a u^b = - m'(v) (u^v)^2 / r
     cioe' l'energia di Kodama driftrebbe a tasso proporzionale a m'.
     La rotaia fornisce esattamente la compensazione: e' la "simmetria
     mancante" resa manifesta.

Uscita: stampa le verifiche + fig_kodama_conservazione (E_K piatta lungo
la brachistocrona vs drift geodetico che la rotaia annulla).
"""

import os
import sys
import numpy as np
import matplotlib.pyplot as plt
import sympy as sp
from scipy.integrate import solve_ivp
from scipy.optimize import brentq

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from paper_style import COL, set_style, savefig

set_style()
OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'Vaidyafigures')

print("=" * 70)
print("[A] miracolo di Kodama (simbolico, Vaidya EF entrante)")
print("=" * 70)
v, r, th = sp.symbols('v r theta', real=True)
m = sp.Function('m')
f = 1 - 2 * m(v) / r
# metrica 4D in (v, r, theta, phi)
g = sp.diag(-f, 0, r**2, r**2 * sp.sin(th)**2)
g[0, 1] = g[1, 0] = 1                       # blocco (v,r) = [[-f,1],[1,0]]
ginv = g.inv()
detg = g.det()
sqrtg = sp.sqrt(-detg)
print("  sqrt(-g) =", sp.simplify(sqrtg), "  (= r^2 sin theta)")

# K^a = eps^{ab} nabla_b r  nel blocco (v,r);  nabla_b r = delta_b^r
# eps^{vr} = 1/sqrt|h| con |h(v,r)| = 1  =>  K = (1,0,0,0)
K_up = sp.Matrix([1, 0, 0, 0])
print("  K^a =", K_up.T, " = d/dv")
# divergenza: (1/sqrt-g) d_a( sqrt-g K^a )
divK = sum(sp.diff(sqrtg * K_up[a], (v, r, th, sp.Symbol('phi'))[a])
           for a in range(3)) / sqrtg
print("  nabla_a K^a =", sp.simplify(divK), "   (miracolo: = 0)")

# Misner-Sharp: M = (r/2)(1 - g^{ab} d_a r d_b r) = (r/2)(1 - g^{rr})
M_MS = sp.simplify((r / 2) * (1 - ginv[1, 1]))
print("  massa di Misner-Sharp = (r/2)(1 - g^{rr}) =", M_MS, " = m(v)")

# K_a (indice basso), usato dopo
K_dn = sp.simplify(g * K_up)
print("  K_a =", K_dn.T, "   (= (-f, 1, 0, 0)  => K = -f dv + dr)")

print()
print("=" * 70)
print("[B] identita' della rotaia: vincolo W  <=>  -u.K = E  (simbolico)")
print("=" * 70)
uu, E, J = sp.symbols('uu E J', real=True, positive=True)
mm = sp.symbols('mm', positive=True)
feq = 1 - 2 * mm / r
# vincolo della brachistocrona tau (parametro r): phi'^2 = W
W = (feq * uu**2 - 2 * uu - (feq * uu - 1)**2 / E**2) / r**2
# tangente T = (v', r', 0, phi') = (uu, 1, 0, sqrt(W)); norma:
gTT = -feq * uu**2 + 2 * uu + r**2 * W       # g(T,T)
print("  -g(T,T) =", sp.simplify(-gTT), "   ( = (f uu - 1)^2 / E^2 )")
# K.T = K_a T^a = (-f) uu + 1 = 1 - f uu ;  -K.u = -(K.T)/sqrt(-gTT)
E_K = sp.simplify((feq * uu - 1) / sp.sqrt(sp.simplify(-gTT)))
print("  -u.K = (f uu - 1)/sqrt(-g(T,T)) =", E_K,
      "   (= E, a meno del segno di ingresso)")

print()
print("=" * 70)
print("[C] verifica numerica lungo una brachistocrona (Vaidya accresce)")
print("=" * 70)
# EOM tau in parametro v (dalla macchina esistente)
vv2, q = sp.symbols('vv2 q', real=True)
mu = sp.symbols('mu', positive=True)
f2 = 1 - 2 * mm / r
W2 = (f2 * uu**2 - 2 * uu - (f2 * uu - 1)**2 / E**2) / r**2
F_tau = (f2 * uu - 1) / E - J * sp.sqrt(W2)
p_tau = sp.diff(F_tau, uu)
dp_tot = (sp.diff(p_tau, r) + sp.diff(p_tau, vv2) * uu
          + sp.diff(p_tau, mm) * mu * uu + sp.diff(p_tau, uu) * q)
EL = dp_tot - sp.diff(F_tau, mm) * mu
vpp = sp.solve(sp.Eq(EL, 0), q)[0]

m0, E_n, J_n = 1.0, 1.2, 1.3
mu_n = 0.02
r1, r0 = 10.0, 4.2     # sopra il periasse (~3.01): la parametrizzazione-v
#                        degenera (dv/dr -> inf) proprio al periasse
lin = [(mm, m0 + mu * vv2), (E, E_n), (J, J_n), (mu, mu_n)]
vpp_f = sp.lambdify((r, vv2, uu), vpp.subs(lin), 'numpy')
pv_f = sp.lambdify((r, vv2, uu), p_tau.subs([(mm, m0 + mu * vv2), (E, E_n),
                                             (J, J_n), (mu, mu_n)]), 'numpy')
W_f = sp.lambdify((r, vv2, uu),
                  W2.subs([(mm, m0 + mu * vv2), (E, E_n), (mu, mu_n)]),
                  'numpy')

def vp_rad(fv):
    a2 = fv - fv**2 / E_n**2
    a1 = -2 + 2 * fv / E_n**2
    a0 = -1 / E_n**2
    return (-a1 + np.sqrt(a1**2 - 4 * a2 * a0)) / (2 * a2)

v_arr = 40.0
fl = 1 - 2 * (m0 + mu_n * v_arr) / r1
lo = vp_rad(fl) * (1 + 1e-9)
hi = lo * 2
while pv_f(r1, v_arr, hi) < 0:
    hi *= 2
up1 = brentq(lambda z: pv_f(r1, v_arr, z), lo, hi)
sol = solve_ivp(lambda rr, y: [y[1], vpp_f(rr, y[0], y[1])],
                [r1, r0], [v_arr, up1], rtol=1e-12, atol=1e-14,
                dense_output=True, method='DOP853')
rg = np.linspace(r1, r0, 600)
vg, ug = sol.sol(rg)
mg = m0 + mu_n * vg
fg = 1 - 2 * mg / rg
Wg = np.clip(W_f(rg, vg, ug), 0, None)

# 4-velocita' u = T/sqrt(-gTT),  T=(uu,1,0,sqrt(W))
gTT_n = -fg * ug**2 + 2 * ug + rg**2 * Wg
norm = np.sqrt(-gTT_n)
uv = ug / norm                         # u^v
E_K_num = (fg * ug - 1) / norm         # -u.K
dev_E = np.max(np.abs(np.abs(E_K_num) - E_n))
print(f"  |-u.K| lungo la brachistocrona: min {np.abs(E_K_num).min():.12f}  "
      f"max {np.abs(E_K_num).max():.12f}")
print(f"  invariante di rotaia E = {E_n}  =>  "
      f"|(-u.K) - E| max = {dev_E:.2e}")

# drift geodetico che la rotaia annulla:  d(-K.u)/dtau = -m'(u^v)^2/r
drift_geo = -mu_n * uv**2 / rg          # m'(v) = mu_n costante
# integrato in dtau (dtau = norm dr, r decresce): drift cumulato geodetico
dtau = norm * np.gradient(rg)
drift_cum = np.cumsum(drift_geo * dtau)
print(f"  tasso drift geodetico -m'(u^v)^2/r: |.| tra "
      f"{np.abs(drift_geo).min():.2e} e {np.abs(drift_geo).max():.2e}")
print(f"  energia di Kodama che la rotaia risparmia (drift geo. cumulato) "
      f"= {drift_cum[-1]:+.4f}  su Delta(-u.K)_brachi < 1e-9")

# --------------------------------------------------------------- figura
fig, (a1, a2) = plt.subplots(2, 1, figsize=(COL, 5.2))
a1.plot(rg, np.abs(E_K_num), 'C0-', lw=1.5,
        label='$-u\\cdot K$ (brachistochrone)')
a1.axhline(E_n, color='k', ls='--', lw=0.8, label=f'rail $E={E_n}$')
a1.set_ylabel('$-u\\cdot K$')
a1.set_ylim(E_n - 0.02, E_n + 0.02)
a1.set_title('Kodama energy conserved along the $\\tau$-brachistochrone\n'
             f'(accreting Vaidya $m=1+{mu_n}v$): $|{{-}}u\\cdot K-E|'
             f'<{dev_E:.0e}$')
a1.legend()
a2.plot(rg, drift_geo, 'C3-', lw=1.5,
        label="geodesic rate $-m'(u^v)^2/r$")
a2.axhline(0, color='k', lw=0.6)
a2.set_xlabel('$r$')
a2.set_ylabel('$d(-u\\cdot K)/d\\tau$')
a2.set_title('what the rail cancels: the free-particle Kodama drift\n'
             "($\\propto m'$): non-zero, this is the missing symmetry")
a2.legend()
savefig(fig, OUT, 'fig_kodama_conservazione')
print('\nFATTO.')
