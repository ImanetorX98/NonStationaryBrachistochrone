"""Kerr thought experiment scaled to M87* (section `sec:observational-context`).

SUPERSEDED: Paper II now places the thought experiment around Gaia BH3
(probe_gaia_bh3.py, paragraph `par:obs-kerr-thought`).  Kept for the record.

Frozen equatorial Kerr, A = M = 1, Ehat = 1, arrival-time (t-branch) rail with endpoints
(r0, 0) and (r0, 1.8), zero winding and one radial turn; r0 = 120M and 38M, a/M = -0.8, 0, 0.8.
Proper acceleration b of the controlled worldline from the exact Christoffel symbols, written in a
form regular at the turning point. Reports turning radius, peak |b| in g for M = 6.5e9 Msun,
scalar budget I = int |b| dtau (units of c), ideal photon-rocket mass ratio exp(I), coordinate
duration and the spin shift of the arrival. Checks: geodesics have b = 0, the rail keeps the
charge and the mass shell, and g(b, d_t) = 0 (no Killing power at A = 1).
Writes probe_kerr_m87.json.
"""
from pathlib import Path
import json
import math
import numpy as np
import sympy as sp
from scipy.optimize import brentq, minimize_scalar
from numpy.polynomial.legendre import leggauss

HERE = Path(__file__).resolve().parent
C = 299792458.; GM_SUN = 1.32712440018e20; G0 = 9.80665; MASS = 6.5e9
checks = []


def check(name, value, atol):
    checks.append(dict(name=name, value=float(value), atol=atol, status='PASS' if abs(value) <= atol else 'FAIL'))


def integrate(fun, lo, hi, n):
    x, w = leggauss(n); pts = (lo + hi)/2 + (hi - lo)*x/2
    return float((hi - lo)/2*np.sum(w*np.array([fun(p) for p in pts])))


def acceleration_functions():
    r, a, j = sp.symbols('r a j', real=True)
    f = 1 - 2/r; d = r*r - 2*r + a*a
    g = sp.Matrix([[-f, 0, -2*a/r], [0, r*r/d, 0], [-2*a/r, 0, r*r + a*a + 2*a*a/r]])  # (t, r, phi)
    gi = sp.simplify(g.inv())
    dg = lambda i, k, l: sp.diff(g[i, k], r) if l == 1 else sp.S.Zero
    Gam = [[[sp.simplify(sum(gi[i, l]*(dg(l, k, h) + dg(l, h, k) - dg(k, h, l))/2 for l in range(3)))
             for h in range(3)] for k in range(3)] for i in range(3)]
    up = 2*(j*(r - 2) + 2*a)/(r*r*d)          # u^phi on the Ehat = 1 arrival rail, costate j
    ut = (1 - 2*a/r*up)/f
    V = sp.factor((-1 - g[0, 0]*ut*ut - 2*g[0, 2]*ut*up - g[2, 2]*up*up)/g[1, 1])   # (u^r)^2
    Ft = sp.diff(ut, r) + 2*Gam[0][0][1]*ut + 2*Gam[0][1][2]*up   # b^t / u^r
    Fp = sp.diff(up, r) + 2*Gam[2][0][1]*ut + 2*Gam[2][1][2]*up   # b^phi / u^r
    Br = sp.diff(V, r)/2 + Gam[1][0][0]*ut*ut + 2*Gam[1][0][2]*ut*up + Gam[1][2][2]*up*up + Gam[1][1][1]*V
    check('zero Killing power g(b, d_t) symbolic', 0 if sp.simplify(g[0, 0]*Ft + g[0, 2]*Fp) == 0 else 1, 0)
    Ug = gi*sp.Matrix([-1, 0, 2])            # geodesic E = 1, L = 2
    Vg = sp.factor((-1 - g[0, 0]*Ug[0]**2 - 2*g[0, 2]*Ug[0]*Ug[2] - g[2, 2]*Ug[2]**2)/g[1, 1])
    res = [sp.diff(Ug[0], r) + 2*Gam[0][0][1]*Ug[0] + 2*Gam[0][1][2]*Ug[2],
           sp.diff(Ug[2], r) + 2*Gam[2][0][1]*Ug[0] + 2*Gam[2][1][2]*Ug[2],
           sp.diff(Vg, r)/2 + Gam[1][0][0]*Ug[0]**2 + 2*Gam[1][0][2]*Ug[0]*Ug[2] + Gam[1][2][2]*Ug[2]**2 + Gam[1][1][1]*Vg]
    check('geodesic acceleration identically zero', 0 if all(sp.simplify(z) == 0 for z in res) else 1, 0)
    return sp.lambdify((r, a, j), (ut, up, V, Ft, Fp, Br), 'numpy', cse=True), sp.lambdify((r, a), g, 'numpy', cse=True)


def coefficients(r, a):
    """Arrival functional sqrt(P dr^2 + Q dphi^2) + b dphi at Ehat = 1, A = 1."""
    f = 1 - 2/r; d = r*r - 2*r + a*a
    return r**2/(f*d*(1 - f)), d/(f*f*(1 - f)), -2*a/(r*f)


def rail(a, r0, fn, metric, n):
    def acc(r, j):
        ut, up, V, Ft, Fp, Br = fn(r, a, j); gg = metric(r, a)
        return math.sqrt(max(0., V*(gg[0, 0]*Ft*Ft + 2*gg[0, 2]*Ft*Fp + gg[2, 2]*Fp*Fp) + gg[1, 1]*Br*Br))

    def kernels(rm, full=False):
        m = np.longdouble(rm); aa = np.longdouble(a)
        Pm, Qm, bm = coefficients(m, aa); rootQ = np.sqrt(Qm); j = float(bm + rootQ)

        def at(z):
            rr = m + np.longdouble(z)**2; P, Q, b = coefficients(rr, aa); db = b - bm
            gap = Q - Qm + 2*rootQ*db - db*db
            kp = 2*z*(rootQ - db)*np.sqrt(P/Q/gap); kal = 2*z*np.sqrt(P*Q/gap)
            if not full:
                return float(kp)
            ktau = (1 - 2/rr)*kal
            return np.array([kp, kal + b*kp, ktau, acc(float(rr), j)*ktau], float)
        end = math.sqrt(r0 - rm)
        if not full:
            return 2*integrate(at, 0, end, n)
        return j, np.array([2*integrate(lambda z: at(z)[i], 0, end, n) for i in range(4)])
    rm = brentq(lambda x: kernels(x) - 1.8, max(2.05, .28*r0), .75*r0, xtol=2e-12)
    j, vals = kernels(rm, True)
    opt = minimize_scalar(lambda x: -acc(x, j), bounds=(rm, r0), method='bounded')
    peak_r = max([rm, r0, float(opt.x)], key=lambda x: acc(x, j))
    shell = charge = 0.
    for x in np.linspace(rm, r0, 151):
        ut, up, V, Ft, Fp, Br = fn(x, a, j); gg = metric(x, a)
        U = np.array([ut, math.sqrt(max(0., V)), up])
        shell = max(shell, abs(U@gg@U + 1)); charge = max(charge, abs(-(gg@U)[0] - 1))
    check(f'a={a} r0={r0} endpoint angle', vals[0] - 1.8, 2e-9)
    check(f'a={a} r0={r0} mass shell', shell, 1e-10)
    check(f'a={a} r0={r0} charge Ehat=1', charge, 1e-10)
    return dict(a=a, r0_M=r0, rmin_M=rm, t_M=vals[1], tau_M=vals[2], impulse_over_c=vals[3],
                peak_accel_over_M=acc(peak_r, j), peak_r_M=peak_r, beta_static_peak=math.sqrt(2/rm))


def main():
    fn, metric = acceleration_functions()
    rg = MASS*GM_SUN/C**2; tg = rg/C
    rows = []
    for r0 in [120., 38.]:
        for a in [-.8, 0., .8]:
            row = rail(a, r0, fn, metric, 224)
            ref = rail(a, r0, fn, metric, 160)
            check(f'quadrature refinement a={a} r0={r0} t', (row['t_M'] - ref['t_M'])/row['t_M'], 3e-8)
            row.update(peak_g=row['peak_accel_over_M']*C*C/rg/G0, t_days=row['t_M']*tg/86400,
                       photon_mass_ratio=math.exp(row['impulse_over_c']))
            rows.append(row)
    for r0 in [120., 38.]:
        grp = {z['a']: z for z in rows if z['r0_M'] == r0}
        for a in [-.8, .8]:
            grp[a]['spin_shift_minutes'] = (grp[a]['t_M'] - grp[0.]['t_M'])*tg/60
    out = dict(model='frozen equatorial Kerr, Ehat = 1, arrival-time rail, endpoints (r0,0)-(r0,1.8), one turn, '
                     'zero winding; spin illustrative; M = 6.5e9 Msun', rows=rows,
               photon_rocket_power_W_per_tonne_at_1g=1000*G0*C, checks=checks)
    (HERE/'probe_kerr_m87.json').write_text(json.dumps(out, indent=2))
    print('checks', len(checks), 'PASS', sum(c['status'] == 'PASS' for c in checks))
    for z in rows:
        print('r0=%3.0fM a=%+.1f: rmin %.4fM, peak %.4f g, I %.6f, exp(I) %.5f, %.4f days, spin shift %s min' %
              (z['r0_M'], z['a'], z['rmin_M'], z['peak_g'], z['impulse_over_c'], z['photon_mass_ratio'],
               z['t_days'], '%+.2f' % z['spin_shift_minutes'] if 'spin_shift_minutes' in z else '-'))
    print('ideal photon rocket, P = m a c: %.3e W per tonne at 1 g' % (1000*G0*C))
    assert all(c['status'] == 'PASS' for c in checks)


if __name__ == '__main__':
    main()
