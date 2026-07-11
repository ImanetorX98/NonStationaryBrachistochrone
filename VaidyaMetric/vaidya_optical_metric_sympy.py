# -*- coding: utf-8 -*-
"""
Metrica ottica per Vaidya: si puo'?

- Vaidya GENERICA: NO. La riduzione ottica/Randers richiede un Killing
  di tipo tempo (Perlick 1991); d_v non e' Killing se m'(v) != 0
  (verificato in vaidya_metric_sympy.py) e non c'e' nemmeno un Killing
  conforme generico.

- Vaidya LINEARE m(v) = mu*v: SI', via OMOTETIA. Il campo
      xi = v d_v + r d_r
  soddisfa Lie_xi g = 2 g (Killing omotetico, caso speciale di conforme).
  Coordinate adattate:  s = ln v,  x = r/v:
      g = e^{2s} * ghat(x) ,   d_s Killing GENUINO di ghat
  => stessa mossa di FLRW (g = a^2 ghat), con e^s al posto di a.
  Su ghat statica si costruisce la riduzione ottica/Randers standard.

Verifiche:
  V1  Lie_xi g = 2 g  (omotetia esatta, solo per m = mu*v)
  V2  in (s, x, phi):  g = e^{2s} ghat(x),  d ghat/ds = 0
      ghat = (2x - f) ds^2 + 2 ds dx + x^2 dphi^2,  f = 1 - 2 mu/x
  V3  superfici di luce della riduzione: A = -ghat_ss = f - 2x = 0
        <=>  2 x^2 - x + 2 mu = 0  =>  x_pm = (1 ± sqrt(1 - 16 mu))/4
      ESISTONO solo per mu <= 1/16:  mu_c = 1/16 critico (noto in
      letteratura per Vaidya autosimilare / singolarita' nude)
  V4  x_pm sono ipersuperfici nulle autosimilari: raggio uscente
      dr/dv = f/2 con x = cost  <=>  x = f(x)/2  <=>  2x^2 - x + 2mu = 0
      => x_- e' l'ORIZZONTE DEGLI EVENTI (r_EH = x_- v), x_+ la
      superficie esterna di fuga. Check numerico: bisezione sui raggi.
  V5  metrica ottica (Randers) per i raggi nulli su ghat:
        ds_arrivo = [ dx + sqrt( dx^2 + A x^2 dphi^2 ) ] / A
      cioe'  F = beta + alpha:
        alpha^2 = dx^2/A^2 + x^2 dphi^2 / A ,   beta = dx / A
      minimizza s_arrivo = ln v_arrivo: e' la brachistocrona nulla in v!
      Verifica: F soddisfa la condizione nulla di ghat. Singolare a
      A = 0 (x_pm), come 1/F all'ergosfera di Kerr.
"""

import sympy as sp

print("=" * 72)
print("[V1] omotetia di Vaidya lineare: Lie_xi g = 2g, xi = v d_v + r d_r")
print("=" * 72)

v, r, x, s, mu = sp.symbols('v r x s mu', positive=True)
mfun = sp.Function('m', positive=True)

# metrica equatoriale (v, r, phi)
def metrica(mass_expr):
    f = 1 - 2 * mass_expr / r
    return sp.Matrix([[-f, 1, 0], [1, 0, 0], [0, 0, r**2]])

X = [v, r, sp.Symbol('phi')]
xi = [v, r, 0]

def lie_g(g):
    L = sp.zeros(3, 3)
    for a in range(3):
        for b in range(3):
            L[a, b] = sum(xi[c] * g[a, b].diff(X[c])
                          + g[c, b] * sp.diff(xi[c], X[a])
                          + g[a, c] * sp.diff(xi[c], X[b]) for c in range(3))
    return sp.simplify(L)

g_lin = metrica(mu * v)
print("  m = mu*v:      Lie_xi g - 2g =", sp.simplify(lie_g(g_lin) - 2 * g_lin))
g_gen = metrica(mfun(v))
resid = sp.simplify(lie_g(g_gen) - 2 * g_gen)
print("  m(v) generica: (Lie_xi g - 2g)_vv =", resid[0, 0])
print("  => omotetia  <=>  v m' = m  <=>  m = mu*v  (SOLO lineare)")

print()
print("=" * 72)
print("[V2] coordinate autosimilari  s = ln v,  x = r/v")
print("=" * 72)

# trasformazione: v = e^s, r = x e^s; Jacobiano d(v,r)/d(s,x)
Jm = sp.Matrix([[sp.exp(s), 0], [x * sp.exp(s), sp.exp(s)]])
Jfull = sp.eye(3)
Jfull[0:2, 0:2] = Jm
g_old = g_lin.subs([(v, sp.exp(s)), (r, x * sp.exp(s))])
g_new = sp.simplify(Jfull.T * g_old * Jfull)
ghat = sp.simplify(g_new / sp.exp(2 * s))
print("  g = e^{2s} * ghat,   ghat =")
sp.pprint(ghat)
print("  d ghat/ds =", sp.simplify(ghat.diff(s)), "  => d_s KILLING di ghat")

print()
print("=" * 72)
print("[V3] superfici di luce della riduzione: A = -ghat_ss = 0")
print("=" * 72)

A = sp.simplify(-ghat[0, 0])
print("  A = f - 2x =", sp.factor(A), "   [f = 1 - 2mu/x]")
sols = sp.solve(sp.Eq(A, 0), x)
print("  A = 0  <=>  2x^2 - x + 2mu = 0  =>  x_pm =", sols)
print("  reali  <=>  mu <= 1/16:  ACCRESCIMENTO CRITICO mu_c = 1/16")

print()
print("=" * 72)
print("[V4] x_- = orizzonte eventi autosimilare (check numerico)")
print("=" * 72)

# x = cost e' nulla sse dr/dv = x = f(x)/2
cond = sp.simplify(x - (1 - 2 * mu / x) / 2)
print("  raggio uscente con x=cost:  x - f/2 =",
      sp.factor(sp.together(cond)), " = 0  <=> stessi x_pm")

import numpy as np
from scipy.integrate import solve_ivp
mu_n = 0.02
xm = float(sols[0].subs(mu, mu_n))
xp = float(sols[1].subs(mu, mu_n))
xm_, xp_ = min(xm, xp), max(xm, xp)

def raggio(x0, v_end=1e6):
    # dx/dv = (f/2 - x)/v, da v=1
    sol = solve_ivp(lambda vv, y: [((1 - 2 * mu_n / y[0]) / 2 - y[0]) / vv],
                    [1.0, v_end], [x0], rtol=1e-10, atol=1e-13)
    return sol.y[0, -1]

lo, hi = 2 * mu_n * 1.01, 0.4     # tra AH (x=2mu) e x_+
for _ in range(60):
    mid = 0.5 * (lo + hi)
    if raggio(mid) > xm_:          # sfugge verso x_+
        hi = mid
    else:
        lo = mid
print(f"  mu = {mu_n}:  x_- analitico = {xm_:.8f}   "
      f"x_EH bisezione = {0.5*(lo+hi):.8f}")
print(f"             x_+ analitico = {xp_:.8f}   "
      f"raggio da x=0.1 finisce a x = {raggio(0.1):.8f}")
print(f"  orizzonte apparente: x_AH = 2mu = {2*mu_n:.4f}  <  x_- (EH fuori)")

print()
print("=" * 72)
print("[V5] metrica ottica di Randers su ghat (brachistocrona nulla in v)")
print("=" * 72)

dx, dphi = sp.symbols('dx dphi', real=True)
ds_opt = (dx + sp.sqrt(dx**2 + A * x**2 * dphi**2)) / A
null_res = sp.simplify(ghat[0, 0] * ds_opt**2 + 2 * ghat[0, 1] * ds_opt * dx
                       + ghat[2, 2] * dphi**2)
print("  ds_arrivo = [dx + sqrt(dx^2 + A x^2 dphi^2)]/A")
print("  residuo condizione nulla ghat(ds,dx,dphi) =", null_res, " (atteso 0)")
print("  => Randers:  F = alpha + beta,")
print("     alpha^2 = dx^2/A^2 + x^2 dphi^2/A ,   beta = dx/A")
print("     s_arrivo = ln v_arrivo: minimizzarlo = brachistocrona nulla in v")
print("     Singolare ad A=0 (x_pm): analogo esatto di 1/F all'ergosfera.")
print("\n  NOTA massivo: -uhat_s = -u.xi e' conservata dal Killing di ghat;")
print("  per g fisica  -u.xi NON e' conservata (omotetia: d(u.xi)/dtau = -1")
print("  sulle geodetiche). Rotaia su Ehat_omot = -u.xi: vincolo naturale,")
print("  il funzionale v_arrivo = e^{s_arrivo} e' monotono in s_arrivo =>")
print("  il ramo t (arrivo in v) massivo si riduce a problema STATICO su")
print("  ghat, come in FLRW. Il ramo tau NO: int dtau = int e^s dtau_hat")
print("  non e' conformemente/omoteticamente invariante.")
print("\nFATTO.")
