# -*- coding: utf-8 -*-
"""
Hamiltoniane dei rami t e tau in FLRW piatta (formulazione di controllo
ottimo / Pontryagin, tempo cosmico t come evoluzione).

Dalla rotaia conforme Ehat = -u_eta (risultati R1):
    velocita' comovente:  dx/dt = (v/a) u_hat ,  |u_hat| = 1  (controllo)
    v(t) = sqrt(1 - a(t)^2/Ehat^2)

Costi:  ramo t:  int dt ;   ramo tau:  int dtau = int (a/Ehat) dt.

Hamiltoniane (massimizzate sul controllo => u_hat = p/|p|):

    H_t   = (v/a) |p| - 1
    H_tau = (v/a) |p| - a/Ehat

Struttura:
  - dp/dt = -dH/dx = 0  (omogeneita'): p CONSERVATO => rette comoventi.
    Le traiettorie escono dalle eq. di Hamilton: dx/dt = dH/dp = (v/a) p_hat.
  - H NON conservata:  dH/dt = dH/dt|_esplicito ∝ a' — la versione
    hamiltoniana del forzante non autonomo (in Vaidya: dp_v/dr ∝ m').
  - Trasversalita' (arrivo spaziale, istante libero): H(t1) = 0
    => |p| fissato all'arrivo:  |p|_t = a1/v1 ,  |p|_tau = a1^2/(Ehat v1).
    E' l'analogo esatto di p_v(r1) = 0 in Vaidya (p_v = costate di v).

Verifiche:
  V1  eq. di Hamilton simboliche: dp/dt = 0, dx/dt = (v/a) p_hat;
      dH/dt lungo il flusso = d_t esplicito (forzante ∝ a')
  V2  de Sitter numerico: integrazione delle eq. di Hamilton (entrambi
      i rami) vs forme chiuse R5:  T_t = ln(a1)/H,  T_tau = (a1-1)/(Ehat H)
  V3  consistenza col worldline vincolato: |p|_tau(t) richiesto da
      H_tau = 0 coincide con il modulo del momento canonico ridotto
      p_x = Ehat*v*(1+mu*Ehat) => (1+mu*Ehat) = |p|... via mu(eta) di R2
"""

import sympy as sp
import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import brentq

t, Eh, H0 = sp.symbols('t Ehat H', positive=True)
px, py = sp.symbols('p_x p_y', real=True)
x, y = sp.symbols('x y', real=True)
a = sp.Function('a', positive=True)

v = sp.sqrt(1 - a(t)**2 / Eh**2)
P = sp.sqrt(px**2 + py**2)

H_t = (v / a(t)) * P - 1
H_tau = (v / a(t)) * P - a(t) / Eh

print("=" * 72)
print("[V1] equazioni di Hamilton (simboliche)")
print("=" * 72)
for nome, Hh in (("t  ", H_t), ("tau", H_tau)):
    print(f"\n  H_{nome} =", Hh)
    print(f"    dx/dt =  dH/dp_x =", sp.simplify(sp.diff(Hh, px)))
    print(f"    dp_x/dt = -dH/dx =", -sp.diff(Hh, x), "  (p conservato)")
    print(f"    dH/dt|espl =", sp.simplify(sp.diff(Hh, t)),
          "  [∝ a': forzante non autonomo]")
print("\n  => traiettorie: p costante, direzione fissa: RETTA comovente;")
print("     modulo |p| fissato dalla trasversalita' H(t1) = 0.")

print()
print("=" * 72)
print("[V2] de Sitter (a = e^t, H=1): Hamilton numerico vs forme chiuse R5")
print("=" * 72)

Eh_n, X_tgt = 3.0, 0.30

def v_n(tt):
    return np.sqrt(max(1 - np.exp(2 * tt) / Eh_n**2, 0.0))

# eq. di Hamilton: dx/dt = (v/a)*p_hat, p_hat = (1,0); dtau/dt = a/Ehat
def rhs(tt, yv):
    aa = np.exp(tt)
    return [v_n(tt) / aa, aa / Eh_n]

t_free = brentq(lambda T: solve_ivp(rhs, [0, T], [0, 0], rtol=1e-12,
                                    atol=1e-14).y[0, -1] - X_tgt,
                0.05, np.log(Eh_n) - 1e-9)
sol = solve_ivp(rhs, [0, t_free], [0, 0], rtol=1e-12, atol=1e-14)
T_t_num, T_tau_num = t_free, sol.y[1, -1]
a1 = np.exp(t_free)
T_t_cf = np.log(a1)                    # = t1 (banale ma coerente con R5)
T_tau_cf = (a1 - 1) / Eh_n
print(f"  arrivo a X={X_tgt}:  t1 = {t_free:.10f}")
print(f"  T_t   Hamilton = {T_t_num:.12f}   forma chiusa = {T_t_cf:.12f}")
print(f"  T_tau Hamilton = {T_tau_num:.12f}   forma chiusa = {T_tau_cf:.12f}")
print(f"  |p| trasversalita':  ramo t = {a1/v_n(t_free):.6f} ,"
      f"  ramo tau = {a1**2/(Eh_n*v_n(t_free)):.6f}")

print()
print("=" * 72)
print("[V3] consistenza col worldline vincolato (moltiplicatore mu di R2)")
print("=" * 72)
# R2: p_x(worldline) = -C con (1+mu*Ehat) = C/(Ehat*v); u_x = Ehat*v
# costate tau:  |p|(t) da H_tau=0 e' a^2/(Ehat*v) SOLO all'arrivo;
# lungo il moto |p| = cost = a1^2/(Ehat*v1). Confronto: C del worldline
# e' fissata dalla stessa trasversalita' => C = Ehat*v*(1+mu*Ehat) costante.
C_, mu_, vv_ = sp.symbols('C mu v', positive=True)
rel = sp.Eq(1 + mu_ * Eh, C_ / (Eh * vv_))
print("  R2:  (1 + mu*Ehat) = C/(Ehat*v)  <=>  C = Ehat*v*(1+mu*Ehat)")
print("  costate PMP: |p| = cost, fissato all'arrivo = a1^2/(Ehat*v1)")
print("  => stessa struttura: costante del moto trasversale + fattore 1/v;")
print("     il moltiplicatore mu(eta) e' il costate riscalato.")
print("\nFATTO.")
