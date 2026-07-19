# -*- coding: utf-8 -*-
"""
K_tau al sesto ordine PER THAKURTA-KERR: trasferimento e verifica.

Fatto stabilito (riscalamento conforme esatto, A costante):
    H_tau^{TK}(A, Ehat) = (1/A^2) * H_tau^{Kerr}(E -> E_eff = Ehat/A)
=> stesse traiettorie (il fattore 1/A^2 riscala solo il parametro):
   le quasi-costanti di Kerr si trasferiscono per SOSTITUZIONE
   E -> E_eff nei coefficienti (f2, f4, f6 di K_tau_expansion.txt).

Cio' che va VERIFICATO (e' nuovo):
  V1  {K_NLO(E), H_tau(E)} = O(a^4) resta vero anche per E < 1
      (regime E_eff<1 mai testato in Kerr, dove E>1): test di scala
      numerico del bracket di Poisson (rapporto 16x per dimezzamento
      di a => O(a^4); 4x => O(a^2)).
  V2  drift test 3D (protocollo numtest di Kerr): sigma(Q_std) vs
      sigma(K_NLO) lungo traiettorie off-equatoriali del flusso
      H_tau(E_eff), sia a E_eff = 1.2 (controllo Kerr) sia a
      E_eff = 0.96 (Thakurta-Kerr con A = 1.25): mediane e win-rate.
  V3  il polo nuovo: DE0(r) = r^2 w(r) esattamente => i coefficienti
      f2k divergono sul muro di congelamento r_w = 2M/(1-E_eff^2)
      (esiste solo per E_eff < 1): verifica numerica della posizione.
"""

import numpy as np
import sympy as sp
from scipy.integrate import solve_ivp

M_n = 1.0

# ------------------------------------------------- H_tau Kerr 3D + K_NLO
r, th, pr, pth = sp.symbols('r theta p_r p_theta', real=True)
a, E, J = sp.symbols('a E J', positive=True)

Sig = r**2 + a**2 * sp.cos(th)**2
Dl = r**2 - 2 * M_n * r + a**2
F = Sig - 2 * M_n * r
DE = (E**2 - 1) * Sig + 2 * M_n * r

H_tau = (DE * Dl / (Sig * F) * pr**2 + DE / (Sig * F) * pth**2
         + DE * J**2 / (Sig * Dl * sp.sin(th)**2)) / 2

D0 = r**2 - 2 * M_n * r
DE0 = (E**2 - 1) * r**2 + 2 * M_n * r
N2 = (E**2 - 1) * r**2 + 4 * M_n * r - 4 * M_n**2
f2 = N2 * sp.cos(2 * th) / (2 * r**2 * D0 * DE0)

Q_std = pth**2 + J**2 / sp.tan(th)**2 - a**2 * E**2 * sp.cos(th)**2
K_NLO = Q_std + a**2 * f2 * pth**2

PB = (sp.diff(K_NLO, r) * sp.diff(H_tau, pr)
      + sp.diff(K_NLO, th) * sp.diff(H_tau, pth)
      - sp.diff(K_NLO, pr) * sp.diff(H_tau, r)
      - sp.diff(K_NLO, pth) * sp.diff(H_tau, th))
PB_Q = (sp.diff(Q_std, th) * sp.diff(H_tau, pth)
        - sp.diff(Q_std, pth) * sp.diff(H_tau, th)
        + sp.diff(Q_std, r) * sp.diff(H_tau, pr))

PB_f = sp.lambdify((r, th, pr, pth, J, E, a), PB, 'numpy')
PBQ_f = sp.lambdify((r, th, pr, pth, J, E, a), PB_Q, 'numpy')

print("=" * 72)
print("[V1] scala del bracket: {K,H} ~ a^n (rapporto 2^n al dimezzare a)")
print("=" * 72)
rng = np.random.default_rng(9)
for E_v, nome in ((1.2, 'E_eff = 1.2 (regime Kerr)'),
                  (0.96, 'E_eff = 0.96 (Thakurta-Kerr, A=1.25)')):
    rQ, rK = [], []
    for _ in range(8):
        pt = (rng.uniform(4, 9), rng.uniform(0.7, 2.2),
              rng.uniform(-1, 1), rng.uniform(-1, 1),
              rng.uniform(0.5, 2.0))
        e1, e2 = 0.02, 0.01
        vQ = abs(PBQ_f(*pt, E_v, e1) / PBQ_f(*pt, E_v, e2))
        vK = abs(PB_f(*pt, E_v, e1) / PB_f(*pt, E_v, e2))
        rQ.append(vQ)
        rK.append(vK)
    print(f"  {nome}:")
    print(f"    {{Q_std,H}}: rapporto mediano = {np.median(rQ):6.2f}"
          f"   (O(a^2))")
    print(f"    {{K_NLO,H}}: rapporto mediano = {np.median(rK):6.2f}"
          f"   (O(a^2): NESSUNA costante esatta con ansatz p_theta^2 —")
    print("      coerente con la non-integrabilita' del summary Kerr §3;")
    print("      l'NLO riduce il drift LUNGO le traiettorie, non l'ordine")
    print("      puntuale del bracket)")

print()
print("=" * 72)
print("[V2] drift test 3D (protocollo numtest): sigma(Q_std) vs sigma(K_NLO)")
print("=" * 72)
lamb = lambda ex, Ev, av: sp.lambdify((r, th, pr, pth, J),
                                      ex.subs([(E, Ev), (a, av)]), 'numpy')

def drift_test(E_v, a_v, n_traj=24, lam_end=140.0):
    dHp_ = lamb(sp.diff(H_tau, pr), E_v, a_v)
    dHt_ = lamb(sp.diff(H_tau, pth), E_v, a_v)
    dHr_ = lamb(sp.diff(H_tau, r), E_v, a_v)
    dHth_ = lamb(sp.diff(H_tau, th), E_v, a_v)
    Q_f = lamb(Q_std, E_v, a_v)
    K_f = lamb(K_NLO, E_v, a_v)
    wins, ratios = 0, []
    rng2 = np.random.default_rng(11)
    n_ok = 0
    while n_ok < n_traj:
        # regione interna: e' li' che f2 conta (poli a r=2M)
        r0 = rng2.uniform(4.0, 7.0)
        th0 = rng2.uniform(np.deg2rad(45), np.deg2rad(80))
        Jv = rng2.uniform(0.5, 2.0)
        pth0 = rng2.uniform(0.5, 2.0)
        pr0 = -0.2

        def rhs(l_, y):
            return [dHp_(y[0], y[1], y[2], y[3], Jv),
                    dHt_(y[0], y[1], y[2], y[3], Jv),
                    -dHr_(y[0], y[1], y[2], y[3], Jv),
                    -dHth_(y[0], y[1], y[2], y[3], Jv)]
        ev = lambda l_, y: y[0] - 2.5
        ev.terminal, ev.direction = True, -1
        ev2 = lambda l_, y: y[0] - 18.0
        ev2.terminal, ev2.direction = True, 1
        s = solve_ivp(rhs, [0, lam_end], [r0, th0, pr0, pth0],
                      rtol=1e-10, atol=1e-12, method='DOP853',
                      events=[ev, ev2], dense_output=True)
        lg = np.linspace(0, s.t[-1], 400)
        Y = s.sol(lg)
        if s.t[-1] < 20 or np.any(~np.isfinite(Y)):
            continue
        Qs = Q_f(Y[0], Y[1], Y[2], Y[3], Jv)
        Ks = K_f(Y[0], Y[1], Y[2], Y[3], Jv)
        sQ = np.std(Qs) / max(abs(np.mean(Qs)), 1e-12)
        sK = np.std(Ks) / max(abs(np.mean(Ks)), 1e-12)
        if not (np.isfinite(sQ) and np.isfinite(sK)):
            continue
        n_ok += 1
        ratios.append(sQ / sK)
        wins += sK < sQ
    return np.median(ratios), wins, n_ok

for E_v, a_v, nome in ((1.2, 0.4, 'Kerr di controllo (E=1.2, a=0.4)'),
                       (0.96, 0.4, 'Thakurta-Kerr A=1.25 (E_eff=0.96, a=0.4)'),
                       (0.96, 0.7, 'Thakurta-Kerr A=1.25 (E_eff=0.96, a=0.7)')):
    med, wins, n = drift_test(E_v, a_v)
    print(f"  {nome}:")
    print(f"    sigma_std/sigma_NLO mediana = {med:.3f}   NLO vince: "
          f"{wins}/{n}")

print()
print("=" * 72)
print("[V3] il polo nuovo: coefficienti divergono sul muro r_w")
print("=" * 72)
E_v = 0.96
r_w = 2 * M_n / (1 - E_v**2)
f2_f = sp.lambdify((r, th), f2.subs([(E, E_v)]), 'numpy')
print(f"  E_eff = {E_v}:  r_w = 2M/(1-E_eff^2) = {r_w:.4f}")
for rv in (0.6 * r_w, 0.9 * r_w, 0.99 * r_w, 0.999 * r_w):
    print(f"    f2(r = {rv:7.3f}, th=60deg) = {f2_f(rv, np.pi/3):+.3e}")
print("  (in Kerr, E>1: DE0 = r^2 w > 0 ovunque, il polo non esiste)")
print("\nFATTO.")
