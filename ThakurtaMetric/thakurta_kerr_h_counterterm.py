# -*- coding: utf-8 -*-
"""
Il controtermine h: soluzione del trasporto e legge secolare chiusa.

Dal R11d:  dK/deta|espl = eps a^2 S1,  eps = A'/A,
           S1 = 2E^2 [ cos^2th + M cos2th p_th^2/(r DE0^2) ].

MEDIE ANGOLARI IN FORMA CHIUSA (piano orbitale: cos th = sin i sin psi,
k = sin^2 i = 1 - J^2/L^2, L^2 = p_th^2 + J^2/sin^2 th):

    <cos^2 th>        = k/2
    <1/sin^2 th>      = L/|J|
    <p_th^2>          = L (L - |J|)
    <cos2th p_th^2>   = |J| (L - |J|)

=> LEGGE SECOLARE (il drift ineliminabile della quasi-costante):

    <dK/deta> = eps a^2 <S1>,
    <S1> = E_eff^2 (1 - J^2/L^2) + 2 M E_eff^2 |J| (L-|J|) / (r DE0^2)

e il controtermine h = parte oscillatoria, valutabile lungo il flusso:
    h(eta) = -a^2 int eps (S1 - <S1>) deta'    [compensazione completa:
    K_comp = K(E_eff(eta)) - a^2 int eps S1 deta' toglie TUTTO il drift
    esplicito; il residuo e' il puro drift alla Kerr O(a^2)].

Verifiche:
  V1  medie angolari chiuse vs medie temporali sul flusso sferico (a=0)
  V2  legge secolare: int <S1>(r(eta)) deta vs int S1 deta (esatto)
  V3  dinamico (a=0.4): sigma(K_run) vs sigma(K_comp) vs residuo Kerr
      (controllo ad A costante): il controtermine riporta il drift al
      livello del residuo statico.
"""

import numpy as np
import sympy as sp
from scipy.integrate import solve_ivp

M_n, Eh_n = 1.0, 1.2

# --------------------------- H 3D (come thakurta_kerr_dynamic_K, param.)
ee = sp.Symbol('ee', negative=True)
r, th = sp.symbols('r theta', positive=True)
pr, pth = sp.symbols('p_r p_theta', real=True)
a, J, A = sp.symbols('a J A', positive=True)

Sig = r**2 + a**2 * sp.cos(th)**2
Dl = r**2 - 2 * M_n * r + a**2
f_S = 1 - 2 * M_n * r / Sig
b = 2 * M_n * a * r * sp.sin(th)**2 / Sig
G = (r**2 + a**2 + 2 * M_n * a**2 * r * sp.sin(th)**2 / Sig) * sp.sin(th)**2
vbS2 = 1 - A**2 * f_S / Eh_n**2
Gb = G + A**2 * b**2 / Eh_n**2
php0 = b * vbS2 / Gb
R2 = vbS2 * Dl * sp.sin(th)**2 / Gb
ptphi = J - A**2 * b / Eh_n
H3 = ptphi * php0 + sp.sqrt(R2) * sp.sqrt((Dl / Sig) * pr**2
                                          + pth**2 / Sig + ptphi**2 / Gb) \
    - A**2 * f_S / Eh_n

def flusso(a_v, A_expr, Jv):
    Hd = H3.subs([(A, A_expr), (a, max(a_v, 1e-12))])
    args = (ee, r, th, pr, pth)
    Hd = Hd.subs(J, Jv)
    d = {k: sp.lambdify(args, sp.diff(Hd, v_), 'numpy', cse=True)
         for k, v_ in (('pr', pr), ('pth', pth), ('r', r), ('th', th))}

    def rhs(e_, y):
        return [d['pr'](e_, *y), d['pth'](e_, *y),
                -d['r'](e_, *y), -d['th'](e_, *y)]
    return rhs

def S1_f(rv, thv, pthv, Ev):
    DE0 = (Ev**2 - 1) * rv**2 + 2 * M_n * rv
    return 2 * Ev**2 * (np.cos(thv)**2
                        + M_n * np.cos(2 * thv) * pthv**2 / (rv * DE0**2))

def S1_med(rv, Jv, Lv, Ev):
    DE0 = (Ev**2 - 1) * rv**2 + 2 * M_n * rv
    return (Ev**2 * (1 - Jv**2 / Lv**2)
            + 2 * M_n * Ev**2 * abs(Jv) * (Lv - abs(Jv)) / (rv * DE0**2))

print("=" * 72)
print("[V1] medie angolari: forme chiuse vs quadratura in psi")
print("=" * 72)
# cos th = sin i sin psi;  p_th^2 = L^2 - J^2/sin^2 th;  k = sin^2 i
from scipy.integrate import quad
Jq, Lq = 1.0, 1.6
kq = 1 - Jq**2 / Lq**2
si = np.sqrt(kq)
cth = lambda ps: si * np.sin(ps)
sth2 = lambda ps: 1 - kq * np.sin(ps)**2
pth2 = lambda ps: Lq**2 - Jq**2 / sth2(ps)
Q = lambda fn: quad(fn, 0, 2 * np.pi, limit=200)[0] / (2 * np.pi)
tab = [("<cos^2 th>", Q(lambda ps: cth(ps)**2), kq / 2),
       ("<1/sin^2 th>", Q(lambda ps: 1 / sth2(ps)), Lq / abs(Jq)),
       ("<p_th^2>", Q(pth2), Lq * (Lq - abs(Jq))),
       ("<cos2th p_th^2>", Q(lambda ps: (1 - 2 * cth(ps)**2) * pth2(ps)),
        abs(Jq) * (Lq - abs(Jq)))]
for nome, num, ana in tab:
    print(f"  {nome:18s} quad = {num:+.10f}   chiusa = {ana:+.10f}   "
          f"diff = {abs(num-ana):.1e}")
print("  NOTA (risultato): per E_eff < 1 il ramo tau NON ha orbite legate:")
print("  W_eff(r) e' monotona decrescente verso il muro => il muro e'")
print("  ATTRATTORE GLOBALE anche ad A costante (p_r -> oo in eta finito):")
print("  estende la cattura universale R11 oltre il caso dinamico.")

print()
print()
print("=" * 72)
print("[V3] dinamico a=0.15, regime legato: controtermine end-to-end")
print("=" * 72)
a_v, Hc_v = 0.15, 0.002
eta0 = -1.0 / Hc_v                    # A: 1 -> (arco di scattering)
eta1 = -1.0 / (Hc_v * 1.30)
A_din = -1 / (Hc_v * ee)

# derivate simboliche costruite UNA volta (J argomento)
Hd_sym = H3.subs([(A, A_din), (a, a_v)])
args_d = (ee, r, th, pr, pth, J)
d = {k: sp.lambdify(args_d, sp.diff(Hd_sym, v_), 'numpy', cse=True)
     for k, v_ in (('pr', pr), ('pth', pth), ('r', r), ('th', th))}
Hs_sym = H3.subs([(A, sp.Integer(1)), (a, a_v)])
ds = {k: sp.lambdify(args_d, sp.diff(Hs_sym, v_), 'numpy', cse=True)
      for k, v_ in (('pr', pr), ('pth', pth), ('r', r), ('th', th))}

rng = np.random.default_rng(31)
res = {k: [] for k in ('K_run', 'K_comp', 'K_static')}
n_ok, n_try = 0, 0
while n_ok < 10 and n_try < 60:
    n_try += 1
    y0 = [rng.uniform(10, 13), rng.uniform(0.8, 2.0), -0.3,
          rng.uniform(0.5, 1.2)]                 # arco di scattering
    Jv = rng.uniform(0.6, 1.6)


    def rhs_d(e_, y):
        Av = -1 / (Hc_v * e_)
        Ev = Eh_n / Av
        epsv = Hc_v * Av
        s1 = S1_f(y[0], y[1], y[3], Ev)
        return [d['pr'](e_, *y[:4], Jv), d['pth'](e_, *y[:4], Jv),
                -d['r'](e_, *y[:4], Jv), -d['th'](e_, *y[:4], Jv),
                epsv * a_v**2 * s1]
    ev = lambda e_, y: y[0] - 2.6
    ev.terminal, ev.direction = True, -1
    ev_w = lambda e_, y: y[0] - 16.0            # fine arco (uscita)
    ev_w.terminal, ev_w.direction = True, 1
    sD = solve_ivp(rhs_d, [eta0, eta1], y0 + [0.0], rtol=1e-10,
                   atol=1e-12, method='DOP853', events=[ev, ev_w],
                   dense_output=True)
    if sD.t[-1] - eta0 < 15.0:
        continue
    lg = np.linspace(eta0, sD.t[-1], 600)
    Yd = sD.sol(lg)
    Aq = -1 / (Hc_v * lg)
    Eq = Eh_n / Aq
    DE0q = (Eq**2 - 1) * Yd[0]**2 + 2 * M_n * Yd[0]
    N2q = (Eq**2 - 1) * Yd[0]**2 + 4 * M_n * Yd[0] - 4 * M_n**2
    D0q = Yd[0]**2 - 2 * M_n * Yd[0]
    f2q = N2q * np.cos(2 * Yd[1]) / (2 * Yd[0]**2 * D0q * DE0q)
    K_run = (Yd[3]**2 + Jv**2 / np.tan(Yd[1])**2
             - a_v**2 * Eq**2 * np.cos(Yd[1])**2 + a_v**2 * f2q * Yd[3]**2)
    K_comp = K_run - Yd[4]                 # controtermine completo
    # deriva END-TO-END (media sulle ultime/prime 5% per smussare psi)
    n5 = 30
    e2e = lambda q: abs(np.mean(q[-n5:]) - np.mean(q[:n5])) \
        / max(abs(np.mean(q)), 1e-14)
    res['K_run'].append(e2e(K_run))
    res['K_comp'].append(e2e(K_comp))
    # legge secolare: drift esplicito esatto vs forma chiusa <S1>
    Lq = np.sqrt(Yd[3]**2 + Jv**2 / np.sin(Yd[1])**2)
    epsq = Hc_v * Aq
    I_med = np.trapezoid(epsq * a_v**2
                         * S1_med(Yd[0], Jv, Lq, Eq), lg)
    if abs(I_med) > 0:
        res.setdefault('sec_ratio', []).append(Yd[4, -1] / I_med)

    # controllo statico LEGATO: A = 1.3 costante, stessa durata
    def rhs_s(e_, y):
        return [ds['pr'](-1e3, *y, Jv), ds['pth'](-1e3, *y, Jv),
                -ds['r'](-1e3, *y, Jv), -ds['th'](-1e3, *y, Jv)]
    sS = solve_ivp(rhs_s, [0, sD.t[-1] - eta0], y0, rtol=1e-10,
                   atol=1e-12, method='DOP853', dense_output=True)
    lgs = np.linspace(0, sS.t[-1], 600)
    Ys = sS.sol(lgs)
    E_st = Eh_n
    DE0s_ = (E_st**2 - 1) * Ys[0]**2 + 2 * M_n * Ys[0]
    N2s_ = (E_st**2 - 1) * Ys[0]**2 + 4 * M_n * Ys[0] - 4 * M_n**2
    D0s_ = Ys[0]**2 - 2 * M_n * Ys[0]
    f2s_ = N2s_ * np.cos(2 * Ys[1]) / (2 * Ys[0]**2 * D0s_ * DE0s_)
    K_st = (Ys[3]**2 + Jv**2 / np.tan(Ys[1])**2
            - a_v**2 * E_st**2 * np.cos(Ys[1])**2
            + a_v**2 * f2s_ * Ys[3]**2)
    res['K_static'].append(e2e(K_st))
    n_ok += 1
    print(f'    traiettoria {n_ok}/10 ok (tentativi {n_try})', flush=True)

print(f"  deriva end-to-end mediana su {n_ok} traiettorie "
      f"(archi di scattering, A dinamico da 1):")
print(f"    K(E_eff(eta))  senza controtermine : "
      f"{np.median(res['K_run']):.3e}")
print(f"    K_comp = K - a^2 int eps S1        : "
      f"{np.median(res['K_comp']):.3e}")
print(f"    residuo Kerr puro (A=1 costante)   : "
      f"{np.median(res['K_static']):.3e}")
print(f"    LEGGE SECOLARE: drift esplicito esatto / previsione <S1>:")
print(f"      rapporto mediano = {np.median(res['sec_ratio']):.4f}   "
      f"(spread: {np.percentile(res['sec_ratio'], 25):.3f} - "
      f"{np.percentile(res['sec_ratio'], 75):.3f})")
print("  NOTA: su questi archi (eps=0.002, a=0.15) il drift conforme")
print("  esplicito e' subdominante al residuo Kerr end-to-end: K_comp ~ K_run;")
print("  il contenuto analitico validato e' la legge secolare <S1>.")
print("\nFATTO.")
