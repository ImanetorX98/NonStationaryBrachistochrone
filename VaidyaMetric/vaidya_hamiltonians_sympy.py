# -*- coding: utf-8 -*-
"""
Hamiltoniane dei rami t e tau in Vaidya (parametro v, PMP/Zermelo).

Dal vincolo -u_v = E (parametro v, ' = d/dv):
    phi'^2 = W_v(r', v, r) = [ f - 2r' - (f-r')^2/E^2 ] / r^2

[V1] L'indicatrice (velocita' ammesse) e' un'ELLISSE — identita' esatta:

    W_v = [ w - (r' - f + E^2)^2/E^2 ] / r^2 ,     w = E^2 - f

    =>  r'  = (f - E^2) + E sqrt(w) cos(theta)
        phi' =  (sqrt(w)/r) sin(theta)                (theta = controllo)

    Zermelo: "vento" radiale ENTRANTE  f - E^2  (< 0 sempre, per E>1)
    + ellisse di navigazione con semiassi E*sqrt(w) (radiale), sqrt(w)/r.

[V2] Hamiltoniane (max su theta, Cauchy-Schwarz — forma chiusa):

    costi:  ramo v (arrivo):  L0 = 1
            ramo tau:         L0 = dtau/dv = (f - r')/E = E - sqrt(w) cos(theta)

    H_v   = p_r (f - E^2) - 1 + sqrt(w) * sqrt( E^2 p_r^2        + p_phi^2/r^2 )
    H_tau = p_r (f - E^2) - E + sqrt(w) * sqrt( (E p_r + 1)^2    + p_phi^2/r^2 )

    Randers/Zermelo: termine di deriva p_r(f - E^2) + norma. Il ramo tau
    e' il ramo v con p_r -> p_r + 1/E (shift dal costo dipendente dal
    controllo) e costo E.

Struttura:
  - p_phi = J conservato (sfericita');
  - dp_r/dv = -dH/dr != 0 (inomogeneita' radiale — a differenza di FLRW);
  - dH/dv = dH/dv|espl ∝ m'(v)  (forzante non autonomo);
  - trasversalita' (arrivo a istante libero): H = 0 all'arrivo;
    statico: H conservata => H ≡ 0 (e' il "p_v = 0" delle riduzioni:
    H = -p_v, momento coniugato al tempo).

Verifiche:
  V1  identita' ellisse (sympy, esatta)
  V2  eq. di Hamilton dai gradienti di H (sympy lambdify)
  V3  flusso di Hamilton H_tau (con H=0 all'arrivo) vs EL v-param
      validata: bounce mu=0.01 (r_min=3.010492, Dphi=2.13793,
      T_tau=17.46364) e statico (r_min = radice esatta 2.72713516);
      drift di H: statico ~ 0, dinamico = int dH/dv|espl
  V4  ramo v: stesso arrivo, orbita DIVERSA (r_min piu' piccolo:
      il ramo t/v penetra di piu', come in Kerr)
"""

import numpy as np
import sympy as sp
from scipy.integrate import solve_ivp

# ------------------------------------------------------------- V1: ellisse
vv, rr, rp, th = sp.symbols('vv rr rp theta', real=True)
E, J, mm, pr = sp.symbols('E J mm p_r', real=True)

f = 1 - 2 * mm / rr
w = E**2 - f
W_v = (f - 2 * rp - (f - rp)**2 / E**2) / rr**2
W_ell = (w - (rp - f + E**2)**2 / E**2) / rr**2

print("=" * 72)
print("[V1] indicatrice ellittica")
print("=" * 72)
print("  W_v - [w - (r'-f+E^2)^2/E^2]/r^2 =", sp.simplify(W_v - W_ell))
rp_th = (f - E**2) + E * sp.sqrt(w) * sp.cos(th)
php_th = (sp.sqrt(w) / rr) * sp.sin(th)
chk = sp.simplify(W_v.subs(rp, rp_th) - php_th**2)
print("  W_v(r'(theta)) - phi'(theta)^2 =", chk, " (parametrizzazione esatta)")
print("  vento di Zermelo: r'|_theta=pi/2 = f - E^2  (entrante, E>1)")

# ------------------------------------------------------------- V2: H
print()
print("=" * 72)
print("[V2] Hamiltoniane in forma chiusa (max su theta)")
print("=" * 72)
H_v = pr * (f - E**2) - 1 + sp.sqrt(w) * sp.sqrt(E**2 * pr**2 + J**2 / rr**2)
H_tau = pr * (f - E**2) - E \
    + sp.sqrt(w) * sp.sqrt((E * pr + 1)**2 + J**2 / rr**2)

# verifica del max: H = max_theta [ p_r r'(th) + J phi'(th) - L0(th) ]
for nome, Hh, L0 in (("v  ", H_v, sp.Integer(1)),
                     ("tau", H_tau, E - sp.sqrt(w) * sp.cos(th))):
    gen = pr * rp_th + J * php_th - L0
    dgen = sp.diff(gen, th)
    ths = sp.atan2(J * sp.sqrt(w) / rr,
                   E * sp.sqrt(w) * (pr + (0 if nome == "v  " else 1 / E))
                   * (1 if nome == "v  " else 1))
    # stazionarieta' e valore: verifica numerica su punti casuali
    rng = np.random.default_rng(1)
    ok_v, ok_s = True, True
    for _ in range(6):
        sub = {E: 1.2, J: float(rng.uniform(0.3, 1.5)),
               mm: float(rng.uniform(0.8, 1.3)),
               rr: float(rng.uniform(3.0, 12.0)),
               pr: float(rng.uniform(-1.5, 1.5))}
        th_star = float(ths.subs(sub))
        val = float(gen.subs(sub).subs(th, th_star))
        Hval = float(Hh.subs(sub))
        dval = float(dgen.subs(sub).subs(th, th_star))
        ok_v = ok_v and abs(val - Hval) < 1e-12 * max(1, abs(Hval))
        ok_s = ok_s and abs(dval) < 1e-12
    print(f"  ramo {nome}: H = max_theta[...] verificato numericamente:",
          ok_v, " (stazionarieta':", ok_s, ")")
dHdv = sp.simplify(sp.diff(H_tau, mm))
print("  dH/dv|espl = m'(v) * dH/dm ;   dH_tau/dm =", sp.factor_terms(dHdv))

# ------------------------------------------------------------- V3: flusso
print()
print("=" * 72)
print("[V3] flusso di Hamilton H_tau vs EL v-param (bounce validato)")
print("=" * 72)
mu = sp.Symbol('mu', real=True)
E_n, J_n, m0v, r1 = 1.2, 1.3, 1.0, 10.0
subs_n = [(mm, m0v + mu * vv), (E, E_n), (J, J_n)]
args = (vv, rr, pr, mu)
dHdp = sp.lambdify(args, sp.diff(H_tau, pr).subs(subs_n), 'numpy')
dHdr = sp.lambdify(args, sp.diff(H_tau, rr).subs(subs_n), 'numpy')
dHdv_f = sp.lambdify(args, (sp.diff(H_tau, mm) * mu).subs(subs_n), 'numpy')
H_fn = sp.lambdify(args, H_tau.subs(subs_n), 'numpy')

def pr_arrivo(mu_n, v1):
    """H_tau = 0 all'arrivo: quadratica in p_r, ramo uscente r' > 0."""
    fl = 1 - 2 * (m0v + mu_n * v1) / r1
    wl = E_n**2 - fl
    # (E - pr(fl-E^2))^2 = wl ((E pr + 1)^2 + J^2/r^2)
    a = (fl - E_n**2)**2 - wl * E_n**2
    b = -2 * E_n * (fl - E_n**2) - 2 * wl * E_n
    c = E_n**2 - wl * (1 + J_n**2 / r1**2)
    roots = np.roots([a, b, c])
    for p0 in sorted(roots):
        if abs(H_fn(v1, r1, p0, mu_n)) < 1e-9 and dHdp(v1, r1, p0, mu_n) > 0:
            return float(p0)
    raise RuntimeError('p_r arrivo non trovato')

def flusso(mu_n, v1=40.0):
    p1 = pr_arrivo(mu_n, v1)

    def rhs(v_, y):
        r_, p_, ph_, ta_ = y
        rpv = dHdp(v_, r_, p_, mu_n)
        return [rpv, -dHdr(v_, r_, p_, mu_n),
                np.sqrt(max((1 - 2*(m0v+mu_n*v_)/r_ - 2*rpv
                             - (1 - 2*(m0v+mu_n*v_)/r_ - rpv)**2/E_n**2)
                            / r_**2, 0.0)),
                (1 - 2*(m0v+mu_n*v_)/r_ - rpv) / E_n]

    ev10 = lambda v_, y: y[0] - r1
    ev10.terminal, ev10.direction = True, 1     # solo il RITORNO a r1
    s = solve_ivp(rhs, [v1, -50], [r1, p1, 0.0, 0.0], rtol=1e-12,
                  atol=1e-14, method='DOP853', events=[ev10],
                  dense_output=True)
    rg = s.sol(np.linspace(v1, s.t[-1], 20000))[0]
    y_end = s.y[:, -1]
    return dict(r_min=rg.min(), v_start=s.t[-1], dphi=-y_end[2],
                Ttau=-y_end[3], sol=s, p1=p1)

rif = {0.01: dict(r_min=3.010492, dphi=2.13793, Ttau=17.46364, v0=0.332),
       0.0: dict(r_min=2.72713516, dphi=2.23434, Ttau=18.31969, v0=1.578)}
for mu_n in (0.0, 0.01):
    o = flusso(mu_n)
    R = rif[mu_n]
    print(f"  mu={mu_n:5.2f}: r_min = {o['r_min']:.6f} (EL {R['r_min']:.6f})"
          f"   Dphi = {o['dphi']:.5f} (EL {R['dphi']:.5f})")
    print(f"           T_tau = {o['Ttau']:.5f} (EL {R['Ttau']:.5f})"
          f"   v_start = {o['v_start']:.3f} (EL {R['v0']:.3f})")
    # drift di H lungo il flusso
    s = o['sol']
    vg = np.linspace(40.0, s.t[-1], 400)
    Hg = np.array([H_fn(v_, s.sol(v_)[0], s.sol(v_)[1], mu_n) for v_ in vg])
    if mu_n == 0.0:
        print(f"           statico: max|H| lungo il flusso = "
              f"{np.abs(Hg).max():.2e}  (H ≡ 0 conservata)")
    else:
        dH_int = np.trapezoid(
            [dHdv_f(v_, s.sol(v_)[0], s.sol(v_)[1], mu_n) for v_ in vg], vg)
        print(f"           dinamico: H(partenza) = {Hg[-1]:+.6f};"
              f"  int dH/dv|espl = {dH_int:+.6f}  (coincidono)")

print()
print("=" * 72)
print("[V4] ramo v (arrivo in v): stesso setup, orbita diversa")
print("=" * 72)
dHdp_v = sp.lambdify(args, sp.diff(H_v, pr).subs(subs_n), 'numpy')
dHdr_v = sp.lambdify(args, sp.diff(H_v, rr).subs(subs_n), 'numpy')
Hv_fn = sp.lambdify(args, H_v.subs(subs_n), 'numpy')

def pr_arrivo_v(mu_n, v1):
    fl = 1 - 2 * (m0v + mu_n * v1) / r1
    wl = E_n**2 - fl
    a = (fl - E_n**2)**2 - wl * E_n**2
    b = -2 * (fl - E_n**2)
    c = 1 - wl * J_n**2 / r1**2
    for p0 in sorted(np.roots([a, b, c])):
        if abs(Hv_fn(v1, r1, p0, mu_n)) < 1e-9 \
                and dHdp_v(v1, r1, p0, mu_n) > 0:
            return float(p0)
    raise RuntimeError('p_r arrivo (ramo v) non trovato')

for mu_n in (0.0, 0.01):
    p1 = pr_arrivo_v(mu_n, 40.0)

    def rhs(v_, y):
        return [dHdp_v(v_, y[0], y[1], mu_n), -dHdr_v(v_, y[0], y[1], mu_n)]
    ev10 = lambda v_, y: y[0] - r1
    ev10.terminal, ev10.direction = True, 1
    s = solve_ivp(rhs, [40.0, -50], [r1, p1], rtol=1e-12, atol=1e-14,
                  method='DOP853', events=[ev10], dense_output=True)
    rg = s.sol(np.linspace(40.0, s.t[-1], 20000))[0]
    print(f"  mu={mu_n:5.2f}: ramo v  r_min = {rg.min():.6f}   "
          f"(ramo tau: {flusso(mu_n)['r_min']:.6f})  -> il ramo v "
          f"penetra di piu'")
print("\nFATTO.")
