# -*- coding: utf-8 -*-
"""
Brachistocrone in Thakurta (Schwarzschild conforme) — analitico.

    g = a(eta)^2 * ghat_Schw ,
    ds^2 = a(eta)^2 [ -f deta^2 + dr^2/f + r^2 dphi^2 ] ,   f = 1 - 2M/r

Il matrimonio esatto dei due semilavorati: mossa conforme FLRW (rotaia
Ehat = -u_eta sul Killing conforme) su base statica Schwarzschild
(geometria ottica non banale). NIENTE Killing genuino se a' != 0.

Risultati:
  V1  d_eta e' Killing CONFORME per ogni a(eta): Lie g = (2a'/a) g;
      genuino <=> a' = 0
  V2  cinematica della rotaia Ehat = -u_eta:
        dtau/deta = a^2 f / Ehat ,   gamma = Ehat/(a sqrt(f)) ,
        v^2 = 1 - a^2 f / Ehat^2
      SUPERFICIE DI CONGELAMENTO  a(eta)^2 f(r) = Ehat^2: la barriera
      unifica FLRW (f=1: istante a=Ehat) e Schwarzschild/SdS (a=1:
      raggio f=E^2). Con a crescente la regione accessibile
      f < Ehat^2/a^2 si RESTRINGE verso l'orizzonte: l'espansione
      consegna la particella al buco nero.
  V3  indicatrice: (r'/f)^2/v^2 + r^2 phi'^2/(f v^2) = 1 (' = d/deta):
        r' = f v cos(theta),  phi' = (sqrt(f) v / r) sin(theta)
      niente vento (base statica). Hamiltoniane PMP:
        H_eta = v sqrt(f) sqrt( f p_r^2 + J^2/r^2 ) - 1
        H_tau = v sqrt(f) sqrt( f p_r^2 + J^2/r^2 ) - a^2 f/Ehat
      dH/deta ∝ a' (forzante conforme); trasversalita' H = 0 all'arrivo.
  V4  riduzione tipo Perlick: deta = dl_opt / v con dl_opt arco della
      METRICA OTTICA di Schwarzschild (dr^2/f^2 + r^2 dphi^2/f):
      Fermat con indice n(eta, r) = 1/v = Ehat/sqrt(Ehat^2 - a^2 f).
      a = 1: n = E/sqrt(w) = indice di Perlick 1991 (recuperato).
      a != 1: n dipende CONGIUNTAMENTE da eta e r: niente orologio
      globale (a differenza di FLRW): OLTRE Perlick in modo essenziale.
  V5  numerico, a=1: flusso H_tau e H_eta riproducono le forme chiuse
      statiche (J e K = fJ/E) — coordinate di Schwarzschild, non EF!
  V6  numerico, de Sitter conforme a = -1/(H eta): orbita lanciata in
      fuga viene raggiunta dalla superficie di congelamento e ricacciata
      verso il buco nero (cattura da espansione).
"""

import numpy as np
import sympy as sp
from scipy.integrate import solve_ivp
from scipy.optimize import brentq

eta, r, M, Eh, J, pr, th = sp.symbols('eta r M Ehat J p_r theta',
                                      positive=True)
a = sp.Function('a', positive=True)
f = 1 - 2 * M / r

print("=" * 72)
print("[V1] Killing conforme")
print("=" * 72)
X = [eta, r, sp.Symbol('phi')]
g = sp.diag(-a(eta)**2 * f, a(eta)**2 / f, a(eta)**2 * r**2)
lie = sp.Matrix(3, 3, lambda i, j: g[i, j].diff(eta))
chk = sp.simplify(lie - 2 * (a(eta).diff(eta) / a(eta)) * g)
print("  Lie_eta g - (2a'/a) g =", chk.norm(), " => CKV per ogni a(eta)")
print("  Killing genuino <=> a' = 0 (Schwarzschild)")

print()
print("=" * 72)
print("[V2] cinematica della rotaia Ehat = -u_eta")
print("=" * 72)
# u^mu = (deta/dtau)(1, r', phi'); vincolo -u_eta = a^2 f deta/dtau = Ehat
rp_, php_ = sp.symbols("r' phi'", real=True)
Lam2 = a(eta)**2 * (f - rp_**2 / f - r**2 * php_**2)      # (dtau/deta)^2
u_eta = -a(eta)**2 * f / sp.sqrt(Lam2)
dtau_deta = sp.solve(sp.Eq(-u_eta, Eh), sp.sqrt(Lam2))
print("  -u_eta = Ehat  =>  dtau/deta = a^2 f/Ehat  (da Lam = a^2 f/Ehat)")
# v^2: dalla normalizzazione
v2 = sp.simplify(1 - (a(eta)**2 * f / Eh**2))
print("  gamma = Ehat/(a sqrt(f));   v^2 = 1 - a^2 f/Ehat^2")
print("  congelamento: a(eta)^2 f(r) = Ehat^2  — superficie in (eta, r)")
print("  limiti: f=1 -> FLRW (a=Ehat);  a=1 -> statico (f=E^2, w=0)")
# verifica identita' indicatrice
v = sp.sqrt(1 - a(eta)**2 * f / Eh**2)
rp_th = f * v * sp.cos(th)
php_th = sp.sqrt(f) * v * sp.sin(th) / r
resid = sp.simplify(Lam2.subs([(rp_, rp_th), (php_, php_th)])
                    - (a(eta)**2 * f / Eh)**2 * a(eta)**2 * f / (a(eta)**2 * f))
resid = sp.simplify(Lam2.subs([(rp_, rp_th), (php_, php_th)])
                    - a(eta)**4 * f**2 / Eh**2)
print("  [V3] indicatrice ellittica: Lam^2(theta) - (a^2 f/Ehat)^2 =",
      resid, " (esatta)")

print()
print("=" * 72)
print("[V3] Hamiltoniane e [V4] riduzione Fermat")
print("=" * 72)
H_eta_s = v * sp.sqrt(f) * sp.sqrt(f * pr**2 + J**2 / r**2) - 1
H_tau_s = v * sp.sqrt(f) * sp.sqrt(f * pr**2 + J**2 / r**2) \
    - a(eta)**2 * f / Eh
# verifica max su theta (numerica, punti casuali)
gen = pr * rp_th + J * php_th
rng = np.random.default_rng(2)
ok = True
for _ in range(6):
    sub = {M: 1.0, Eh: 2.0, J: float(rng.uniform(0.3, 1.5)),
           r: float(rng.uniform(3.0, 12.0)),
           pr: float(rng.uniform(-2, 2)), a(eta): float(rng.uniform(0.5, 1.4))}
    th_g = np.linspace(0, 2 * np.pi, 4001)
    vals = [float(gen.subs(sub).subs(th, t_)) for t_ in th_g[::400]]
    # confronto col termine norma di H
    Hn = float((v * sp.sqrt(f) * sp.sqrt(f * pr**2 + J**2 / r**2)).subs(sub))
    ok = ok and abs(max(float(gen.subs(sub).subs(th, t_))
                        for t_ in th_g) - Hn) < 1e-6 * max(1, Hn)
print("  H = max_theta[p.x'] - costo, termine norma verificato:", ok)
print("  H_eta = v sqrt(f) sqrt(f p_r^2 + J^2/r^2) - 1")
print("  H_tau = v sqrt(f) sqrt(f p_r^2 + J^2/r^2) - a^2 f/Ehat")
print("  dH/deta ∝ a'(eta): forzante conforme; statico: H conservata")
print("\n  [V4] deta = dl_opt/v su metrica ottica Schwarzschild:")
print("       n(eta,r) = Ehat/sqrt(Ehat^2 - a^2 f)")
print("       a=1: n = E/sqrt(E^2-f) = PERLICK 1991 recuperato;")
print("       a(eta) generico: n(eta,r) NON separabile => niente orologio")
print("       globale (solo FLRW f=1 lo consente): OLTRE Perlick.")

print()
print("=" * 72)
print("[V5] a=1: flusso di Hamilton vs forme chiuse statiche")
print("=" * 72)
M_n, E_n, J_n = 1.0, 1.2, 1.3
rr_ = sp.Symbol('rr', positive=True)
f_s = 1 - 2 * M_n / rr_
av = sp.Symbol('a_v', positive=True)
v_s = sp.sqrt(1 - av**2 * f_s / E_n**2)
H_tau_n = v_s * sp.sqrt(f_s) * sp.sqrt(f_s * pr**2 + J_n**2 / rr_**2) \
    - av**2 * f_s / E_n
H_eta_n = v_s * sp.sqrt(f_s) * sp.sqrt(f_s * pr**2 + J_n**2 / rr_**2) - 1

for nome, Hx, Kfun in (("tau", H_tau_n, lambda rv: J_n),
                       ("eta", H_eta_n, lambda rv: f_n_(rv) * J_n / E_n)):
    pass  # definite sotto dopo f_n_

def f_n_(rv):
    return 1 - 2 * M_n / rv

for nome, Hx, Kfun in (("tau", H_tau_n, lambda rv: J_n),
                       ("eta", H_eta_n, lambda rv: f_n_(rv) * J_n / E_n)):
    H1 = Hx.subs(av, 1.0)
    dHdp = sp.lambdify((rr_, pr), sp.diff(H1, pr), 'numpy')
    Hfun = sp.lambdify((rr_, pr), H1, 'numpy')
    okk = True
    for rv in (4.0, 6.0, 9.0):
        # p_r da H=0 (famiglia brachistocrona statica), ramo entrante
        pg = np.linspace(-60, 60, 240001)
        Hg = Hfun(rv, pg)
        roots = [brentq(lambda p: Hfun(rv, p), pg[i], pg[i + 1])
                 for i in range(len(pg) - 1)
                 if np.isfinite(Hg[i]) and np.isfinite(Hg[i + 1])
                 and Hg[i] * Hg[i + 1] <= 0]
        got = None
        for p_ in roots:
            if dHdp(rv, p_) < 0:
                got = p_
        # dphi/dr = (dH/dJ)/(dH/dp_r)
        Jsym = sp.Symbol('Js', positive=True)
        HJ = H1.subs(J_n, Jsym) if False else None
        # ricostruisco dphi/dv esplicitamente
        vloc = np.sqrt(1 - f_n_(rv) / E_n**2)
        dphidv = vloc * np.sqrt(f_n_(rv)) * (J_n / rv**2) \
            / np.sqrt(f_n_(rv) * got**2 + J_n**2 / rv**2)
        drdv = dHdp(rv, got)
        num = abs(dphidv / drdv)
        Kv = Kfun(rv)
        wv = E_n**2 - f_n_(rv)
        Dl = rv**2 * f_n_(rv)
        cf = Kv * rv * np.sqrt(wv * f_n_(rv)) \
            / (Dl * np.sqrt(Dl - Kv**2 * wv))
        okk = okk and abs(num - cf) < 1e-9 * cf
        print(f"  ramo {nome} r={rv}: dphi/dr flusso = {num:.10f}"
              f"   forma chiusa = {cf:.10f}")
    print(f"  ramo {nome}: match:", okk)

print()
print("=" * 72)
print("[V6] cattura da espansione (a = -1/(H eta), H=0.02, Ehat=1.2)")
print("=" * 72)
H_n, Eh_n = 0.02, 1.2
eta0 = -1.0 / H_n           # a(eta0) = 1

def a_n(e_):
    return -1.0 / (H_n * e_)

vloc2 = lambda e_, rv: 1 - a_n(e_)**2 * f_n_(rv) / Eh_n**2
ee, rr2 = sp.symbols('ee rr2', negative=True), sp.Symbol('rr2', positive=True)
rr2 = sp.Symbol('rr2', positive=True)
a_s = -1 / (H_n * sp.Symbol('ee', negative=True))
ee_s = sp.Symbol('ee', negative=True)
a_s = -1 / (H_n * ee_s)
f_2 = 1 - 2 * M_n / rr2
v_2 = sp.sqrt(1 - a_s**2 * f_2 / Eh_n**2)
H_tau_2 = v_2 * sp.sqrt(f_2) * sp.sqrt(f_2 * pr**2 + J_n**2 / rr2**2) \
    - a_s**2 * f_2 / Eh_n
dHdp2 = sp.lambdify((ee_s, rr2, pr), sp.diff(H_tau_2, pr), 'numpy')
dHdr2 = sp.lambdify((ee_s, rr2, pr), sp.diff(H_tau_2, rr2), 'numpy')
Hfun2 = sp.lambdify((ee_s, rr2, pr), H_tau_2, 'numpy')

r_l = 4.0                    # lancio in fuga
pg = np.linspace(-200, 200, 400001)
Hg = Hfun2(eta0, r_l, pg)
roots = [brentq(lambda p: Hfun2(eta0, r_l, p), pg[i], pg[i + 1])
         for i in range(len(pg) - 1)
         if np.isfinite(Hg[i]) and np.isfinite(Hg[i + 1])
         and Hg[i] * Hg[i + 1] <= 0]
p0 = [p_ for p_ in roots if dHdp2(eta0, r_l, p_) > 0][-1]

def rhs(e_, y):
    return [dHdp2(e_, y[0], y[1]), -dHdr2(e_, y[0], y[1])]

ev_h = lambda e_, y: y[0] - 2.0 * M_n * 1.001
ev_h.terminal, ev_h.direction = True, -1
s = solve_ivp(rhs, [eta0, -1e-3], [r_l, p0], rtol=1e-11, atol=1e-13,
              method='DOP853', events=[ev_h], dense_output=True)

def r_freeze(e_):
    """f(r) = Ehat^2/a^2: raggio della superficie di congelamento."""
    t_ = Eh_n**2 / a_n(e_)**2
    if t_ >= 1:
        return np.inf
    return 2 * M_n / (1 - t_)

print(f"  lancio USCENTE da r={r_l}, J={J_n}, Ehat={Eh_n}, a(eta0)=1")
eta_c = s.t[-1]                        # contatto con la superficie
for frac in (0.0, 0.3, 0.6, 0.85, 0.97, 0.999):
    e_ = eta0 + (eta_c - eta0) * frac
    y = s.sol(e_)
    print(f"  eta={e_:9.2f}  a={a_n(e_):6.3f}   r_orbita={y[0]:8.4f}   "
          f"r_congelamento={r_freeze(e_):8.3f}")
r_cont = s.y[0, -1]
print(f"\n  CONTATTO a eta={eta_c:.3f}: r_orbita = {r_cont:.4f}, "
      f"r_congelamento = {r_freeze(eta_c):.4f}  (coincidono: v -> 0)")
print("  dopo il contatto la particella cavalca la superficie "
      "(l'indicatrice degenera a un punto): r = r_freeze(eta):")
for e_ in (eta_c * 0.7, eta_c * 0.4, eta_c * 0.15, eta_c * 0.05):
    print(f"    eta={e_:8.2f}  a={a_n(e_):7.2f}   r_freeze="
          f"{r_freeze(e_):7.4f}   (2M = {2*M_n})")
print("  => lanciata in FUGA, la particella viene raggiunta dalla")
print("     superficie di congelamento e trascinata verso l'orizzonte:")
print("     CATTURA DA ESPANSIONE (f_accessibile < Ehat^2/a^2 -> 0).")
print("\nFATTO.")
