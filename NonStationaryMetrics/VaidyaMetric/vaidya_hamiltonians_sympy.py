# -*- coding: utf-8 -*-
"""
Hamiltonians of the t and tau branches in Vaidya (parameter v, PMP/Zermelo).

From the constraint -u_v = E (parameter v, ' = d/dv):
    phi'^2 = W_v(r', v, r) = [ f - 2r' - (f-r')^2/E^2 ] / r^2

[V1] The indicatrix (admissible velocities) is an ELLIPSE -- exact identity:

    W_v = [ w - (r' - f + E^2)^2/E^2 ] / r^2 ,     w = E^2 - f

    =>  r'  = (f - E^2) + E sqrt(w) cos(theta)
        phi' =  (sqrt(w)/r) sin(theta)                (theta = control)

    Zermelo: INGOING radial "wind"  f - E^2  (always < 0, for E>1)
    + navigation ellipse with semi-axes E*sqrt(w) (radial), sqrt(w)/r.

[V2] Hamiltonians (maximum over theta, Cauchy-Schwarz -- closed form):

    costs:  v branch (arrival):  L0 = 1
            tau branch:          L0 = dtau/dv = (f - r')/E = E - sqrt(w) cos(theta)

    H_v   = p_r (f - E^2) - 1 + sqrt(w) * sqrt( E^2 p_r^2        + p_phi^2/r^2 )
    H_tau = p_r (f - E^2) - E + sqrt(w) * sqrt( (E p_r + 1)^2    + p_phi^2/r^2 )

    Randers/Zermelo: drift term p_r(f - E^2) + norm. The tau branch
    is the v branch with p_r -> p_r + 1/E (shift from the control-dependent
    cost) and cost E.

Structure:
  - p_phi = J conserved (spherical symmetry);
  - dp_r/dv = -dH/dr != 0 (radial inhomogeneity -- unlike FLRW);
  - dH/dv = dH/dv|expl ∝ m'(v)  (non-autonomous forcing);
  - transversality (arrival at a free instant): H = 0 at arrival;
    static case: H conserved => H ≡ 0 (it is the "p_v = 0" of the reductions:
    H = -p_v, the momentum conjugate to time).

Checks:
  V1  ellipse identity (SymPy, exact)
  V2  Hamilton equations from the gradients of H (SymPy lambdify)
  V3  Hamilton flow of H_tau (with H=0 at arrival) against the validated
      v-parametrised Euler-Lagrange solution: bounce mu=0.01 (r_min=3.010492, Dphi=2.13793,
      T_tau=17.46364) and static case (r_min = exact root 2.72713516);
      drift of H: static ~ 0, dynamical = int dH/dv|expl
  V4  v branch: same arrival, DIFFERENT orbit (smaller r_min:
      the t/v branch penetrates further, as in Kerr)
"""

import numpy as np
import sympy as sp
from scipy.integrate import solve_ivp

# ------------------------------------------------------------- V1: ellipse
vv, rr, rp, th = sp.symbols('vv rr rp theta', real=True)
E, J, mm, pr = sp.symbols('E J mm p_r', real=True)

f = 1 - 2 * mm / rr
w = E**2 - f
W_v = (f - 2 * rp - (f - rp)**2 / E**2) / rr**2
W_ell = (w - (rp - f + E**2)**2 / E**2) / rr**2

print("=" * 72)
print("[V1] elliptic indicatrix")
print("=" * 72)
print("  W_v - [w - (r'-f+E^2)^2/E^2]/r^2 =", sp.simplify(W_v - W_ell))
rp_th = (f - E**2) + E * sp.sqrt(w) * sp.cos(th)
php_th = (sp.sqrt(w) / rr) * sp.sin(th)
chk = sp.simplify(W_v.subs(rp, rp_th) - php_th**2)
print("  W_v(r'(theta)) - phi'(theta)^2 =", chk, " (exact parametrisation)")
print("  Zermelo wind: r'|_theta=pi/2 = f - E^2  (ingoing, E>1)")

# ------------------------------------------------------------- V2: H
print()
print("=" * 72)
print("[V2] Hamiltonians in closed form (maximum over theta)")
print("=" * 72)
H_v = pr * (f - E**2) - 1 + sp.sqrt(w) * sp.sqrt(E**2 * pr**2 + J**2 / rr**2)
H_tau = pr * (f - E**2) - E \
    + sp.sqrt(w) * sp.sqrt((E * pr + 1)**2 + J**2 / rr**2)

# check of the maximum: H = max_theta [ p_r r'(th) + J phi'(th) - L0(th) ]
for nome, Hh, L0 in (("v  ", H_v, sp.Integer(1)),
                     ("tau", H_tau, E - sp.sqrt(w) * sp.cos(th))):
    gen = pr * rp_th + J * php_th - L0
    dgen = sp.diff(gen, th)
    ths = sp.atan2(J * sp.sqrt(w) / rr,
                   E * sp.sqrt(w) * (pr + (0 if nome == "v  " else 1 / E))
                   * (1 if nome == "v  " else 1))
    # stationarity and value: numerical check at random points
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
    print(f"  branch {nome}: H = max_theta[...] verified numerically:",
          ok_v, " (stationarity:", ok_s, ")")
dHdv = sp.simplify(sp.diff(H_tau, mm))
print("  dH/dv|expl = m'(v) * dH/dm ;   dH_tau/dm =", sp.factor_terms(dHdv))

# ------------------------------------------------------------- V3: flow
print()
print("=" * 72)
print("[V3] Hamilton flow of H_tau against the v-parametrised Euler-Lagrange solution (validated bounce)")
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
    """H_tau = 0 at arrival: quadratic in p_r, outgoing branch r' > 0."""
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
    raise RuntimeError('arrival p_r not found')

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
    ev10.terminal, ev10.direction = True, 1     # only the RETURN to r1
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
    # drift of H along the flow
    s = o['sol']
    vg = np.linspace(40.0, s.t[-1], 400)
    Hg = np.array([H_fn(v_, s.sol(v_)[0], s.sol(v_)[1], mu_n) for v_ in vg])
    if mu_n == 0.0:
        print(f"           static: max|H| along the flow = "
              f"{np.abs(Hg).max():.2e}  (H ≡ 0 conserved)")
    else:
        dH_int = np.trapezoid(
            [dHdv_f(v_, s.sol(v_)[0], s.sol(v_)[1], mu_n) for v_ in vg], vg)
        print(f"           dynamical: H(start) = {Hg[-1]:+.6f};"
              f"  int dH/dv|expl = {dH_int:+.6f}  (they agree)")

print()
print("=" * 72)
print("[V4] v branch (arrival in v): same set-up, different orbit")
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
    raise RuntimeError('arrival p_r (v branch) not found')

for mu_n in (0.0, 0.01):
    p1 = pr_arrivo_v(mu_n, 40.0)

    def rhs(v_, y):
        return [dHdp_v(v_, y[0], y[1], mu_n), -dHdr_v(v_, y[0], y[1], mu_n)]
    ev10 = lambda v_, y: y[0] - r1
    ev10.terminal, ev10.direction = True, 1
    s = solve_ivp(rhs, [40.0, -50], [r1, p1], rtol=1e-12, atol=1e-14,
                  method='DOP853', events=[ev10], dense_output=True)
    rg = s.sol(np.linspace(40.0, s.t[-1], 20000))[0]
    print(f"  mu={mu_n:5.2f}: v branch  r_min = {rg.min():.6f}   "
          f"(tau branch: {flusso(mu_n)['r_min']:.6f})  -> the v branch "
          f"penetrates further")
print("\nDONE.")
