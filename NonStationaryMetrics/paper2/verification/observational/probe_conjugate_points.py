"""Conjugate points of maintained-energy rails (section `sec:observational-context`).

A planar rail at maintained energy in a static potential is a geodesic of the optical metric
|dx|^2/f, f = v^2 = 2 Omega - C, and its travel time is the metric length. Coriolis forces are
normal to the velocity and are absorbed by the rail, so in a rotating frame the same holds with the
effective potential Omega and the Jacobi constant C. The Gaussian curvature is
    K = Lap(Omega) - 2 |grad Omega|^2 / f.
Along a ray, transverse Jacobi fields obey J'' + K J = 0 in the travel time t; a zero of J at t_c > 0
is a conjugate point, beyond which the rail is no longer a local minimum of the travel time.

Checks:
  1. symbolic: the curvature formula, and K = E mu/(r^2 (E r + mu)) for a Kepler field (no rotation);
  2. Earth-Moon circular restricted three-body problem, rays from L4 at 1.02 km/s: conjugate times,
     and the sign change of the daughter lag across t_c for the least favourable fixed-endpoint
     deformation (lowest Dirichlet eigenfunction; the initial direction is not held fixed);
  3. a hyperbolic rail winding around Jupiter (v_inf = 6 km/s): conjugate time and sign change;
  4. the Saturn rails of probe_saturn_comparison.py (Sun and four planets on prescribed circular
     orbits, time-dependent field): negative heliocentric energy, and no sign change of the in-plane
     and out-of-plane transverse endpoint fields of the rays from the Hill exit (centred differences
     in the launch angles at fixed epoch and energy, 4000 samples) within nine years. This is a
     sampled numerical check of the model, not a certificate.
Writes probe_conjugate_points.json.
"""
from pathlib import Path
import json
import math
import numpy as np
import sympy as sp
from scipy.integrate import solve_ivp
from scipy.linalg import eigh_tridiagonal
import probe_saturn_comparison as S

HERE = Path(__file__).resolve().parent
checks = []


def check(name, value, tolerance):
    checks.append(dict(name=name, residual=float(value), tolerance=tolerance,
                       status='PASS' if abs(value) <= tolerance else 'FAIL'))


# 1. Curvature of the optical metric -------------------------------------------------------------
x, y, C, E, mu = sp.symbols('x y C E mu', real=True)
Om = sp.Function('Omega')(x, y)
f = 2*Om - C
phi = -sp.log(f)/2                                     # metric e^{2 phi}|dx|^2
K = -sp.exp(-2*phi)*(sp.diff(phi, x, 2) + sp.diff(phi, y, 2))
claim = sp.diff(Om, x, 2) + sp.diff(Om, y, 2) - 2*(sp.diff(Om, x)**2 + sp.diff(Om, y)**2)/f
check('curvature K = Lap Omega - 2|grad Omega|^2/f (symbolic)', 0 if sp.simplify(K - claim) == 0 else 1, 0)
r = sp.sqrt(x**2 + y**2)
Kk = claim.subs(Om, mu/r).doit().subs(C, -2*E)
check('Kepler: K = E mu/(r^2 (E r + mu)) (symbolic)',
      0 if sp.simplify(Kk - E*mu/(r**2*(E*r + mu))) == 0 else 1, 0)

# 2. Earth-Moon CR3BP, rays from L4 --------------------------------------------------------------
MU = 0.012150585; LKM = 384400.; TU = 1/2.6617e-6; VU = LKM*1e3/TU; AU_ = VU/TU


def om(x, y):
    return .5*(x*x + y*y) + (1 - MU)/np.hypot(x + MU, y) + MU/np.hypot(x - 1 + MU, y)


def grad_lap(x, y):
    r1 = np.hypot(x + MU, y); r2 = np.hypot(x - 1 + MU, y)
    gx = x - (1 - MU)*(x + MU)/r1**3 - MU*(x - 1 + MU)/r2**3
    gy = y - (1 - MU)*y/r1**3 - MU*y/r2**3
    return gx, gy, 2 + (1 - MU)/r1**3 + MU/r2**3


def cr3bp_rhs(t, s, Cj):
    x, y, th, J, Jd = s; gx, gy, lap = grad_lap(x, y); f = 2*om(x, y) - Cj; v = math.sqrt(f)
    K = lap - 2*(gx*gx + gy*gy)/f
    return [v*math.cos(th), v*math.sin(th), -(-gx*math.sin(th) + gy*math.cos(th))/v, Jd, -K*J]


def first_zero(rhs, s0, tmax, args=()):
    ev = lambda t, s, *a: s[3] if t > 1e-3*tmax else 1.
    ev.terminal = True; ev.direction = -1
    sol = solve_ivp(rhs, (0, tmax), s0, args=args, events=ev, rtol=1e-11, atol=1e-13, dense_output=True)
    return (sol.t_events[0][0] if sol.t_events[0].size else None), sol


def lag_sign_test(path, speed, curvature, cost, t_end, scale, offsets):
    """Daughter lag along the lowest Dirichlet eigenfunction of -d^2/dt^2 - K (least favourable shape)."""
    t = np.linspace(0, t_end, 4001); P, th = path(t); v = speed(P); K = curvature(P, v)
    h = t[1] - t[0]
    w, V = eigh_tridiagonal(2/h**2 - K[1:-1], -np.ones(len(t) - 3)/h**2, select='i', select_range=(0, 0))
    eta = np.zeros(len(t)); eta[1:-1] = V[:, 0]; eta /= np.abs(eta).max()
    N = np.stack([-np.sin(th), np.cos(th)], 1); vm = v[np.argmax(np.abs(eta))]; T0 = cost(P)
    return float(w[0]), [float(cost(P + (d/scale/vm)*(v*eta)[:, None]*N) - T0) for d in offsets]


x4, y4 = .5 - MU, math.sqrt(3)/2; C4 = 2*om(x4, y4) - 1.0        # speed 1 unit = 1.023 km/s at L4


def cr3bp_cost(P):
    d = np.diff(P, axis=0); m = .5*(P[1:] + P[:-1])
    return np.sum(np.hypot(d[:, 0], d[:, 1])/np.sqrt(2*om(m[:, 0], m[:, 1]) - C4))


cr3bp = {}
for deg, ref in [(290, 4.2628), (326, 13.347)]:
    tc, sol = first_zero(cr3bp_rhs, [x4, y4, math.radians(deg), 0., 1.], 6., (C4,))
    xs = sol.sol(np.linspace(0, tc, 2001))
    gx, gy, _ = grad_lap(xs[0], xs[1]); v = np.sqrt(2*om(xs[0], xs[1]) - C4)
    thrust = np.abs(2*(gx*np.sin(xs[2]) - gy*np.cos(xs[2])) + 2*v)       # -2 grad.N + Coriolis 2v
    row = dict(direction_deg=deg, t_c_days=tc*TU/86400, peak_thrust_m_s2=float(thrust.max()*AU_),
               delta_v_to_tc_km_s=float(np.trapz(thrust, np.linspace(0, tc, 2001))*VU/1e3))
    check(f'CR3BP L4 ray {deg} deg: conjugate time (days) vs independent Mathematica value {ref}',
          row['t_c_days'] - ref, 2e-3)
    lags = {}
    for frac in (0.8, 1.2):
        sl = solve_ivp(lambda t, s: cr3bp_rhs(t, list(s) + [0, 0], C4)[:3], (0, frac*tc),
                       [x4, y4, math.radians(deg)], rtol=1e-12, atol=1e-13, dense_output=True)
        w, L = lag_sign_test(lambda t: (sl.sol(t)[:2].T, sl.sol(t)[2]),
                             lambda P: np.sqrt(2*om(P[:, 0], P[:, 1]) - C4),
                             lambda P, v: (lambda g: g[2] - 2*(g[0]**2 + g[1]**2)/v**2)(grad_lap(P[:, 0], P[:, 1])),
                             cr3bp_cost, frac*tc, LKM, (1000., 5000.))
        lags[str(frac)] = dict(lowest_eigenvalue=w, lag_s_d1000km=L[0]*TU, lag_s_d5000km=L[1]*TU)
    check(f'CR3BP L4 ray {deg} deg: daughter late before t_c (0.8 t_c)', min(0., lags['0.8']['lag_s_d1000km']), 0.)
    check(f'CR3BP L4 ray {deg} deg: daughter early beyond t_c (1.2 t_c)', max(0., lags['1.2']['lag_s_d1000km']), 0.)
    row['lags'] = lags; cr3bp[str(deg)] = row

# 3. Hyperbolic rail around Jupiter ----------------------------------------------------------------
GJ, RJ = 1.26686534e17, 7.1492e7; EJ = 6000.**2/2; r0 = 100*RJ


def kep_rhs(t, s):
    x, y, th, J, Jd = s; r = math.hypot(x, y); v = math.sqrt(2*EJ + 2*GJ/r)
    gx, gy = -GJ*x/r**3, -GJ*y/r**3
    return [v*math.cos(th), v*math.sin(th), -(-gx*math.sin(th) + gy*math.cos(th))/v, Jd,
            -EJ*GJ/(r**2*(EJ*r + GJ))*J]


th0 = math.pi + math.atan2(.5*RJ, r0)
tcJ, solJ = first_zero(kep_rhs, [r0, 0., th0, 0., 1.], 40*86400.)
xs = solJ.sol(np.linspace(0, tcJ, 4001)); rr = np.hypot(xs[0], xs[1])
ang = np.unwrap(np.arctan2(xs[1], xs[0]))
jup = dict(v_inf_km_s=6., start_R=100, impact_offset_R=.5, t_c_days=tcJ/86400, pericentre_R=float(rr.min()/RJ),
           swept_deg=float(abs(np.degrees(ang[-1] - ang[0]))), rule_of_2_thrust_at_pericentre_m_s2=2*GJ/(rr.min())**2)
check('Jupiter hyperbolic rail: conjugate time (days) vs independent Mathematica value 9.0144', jup['t_c_days'] - 9.0144, 2e-3)


def kep_cost(P):
    d = np.diff(P, axis=0); m = .5*(P[1:] + P[:-1]); rm = np.hypot(m[:, 0], m[:, 1])
    return np.sum(np.hypot(d[:, 0], d[:, 1])/np.sqrt(2*EJ + 2*GJ/rm))


jl = {}
for frac in (0.8, 1.2):
    sl = solve_ivp(lambda t, s: kep_rhs(t, list(s) + [0, 0])[:3], (0, frac*tcJ), [r0, 0., th0],
                   rtol=1e-12, atol=1e-6, dense_output=True)
    w, L = lag_sign_test(lambda t: (sl.sol(t)[:2].T, sl.sol(t)[2]),
                         lambda P: np.sqrt(2*EJ + 2*GJ/np.hypot(P[:, 0], P[:, 1])),
                         lambda P, v: EJ*GJ/(np.hypot(P[:, 0], P[:, 1])**2*(EJ*np.hypot(P[:, 0], P[:, 1]) + GJ)),
                         kep_cost, frac*tcJ, 1., (.1*RJ, .5*RJ))
    jl[str(frac)] = dict(lowest_eigenvalue=w, lag_s_d01R=L[0], lag_s_d05R=L[1])
check('Jupiter rail: daughter late before t_c (0.8 t_c)', min(0., jl['0.8']['lag_s_d01R']), 0.)
check('Jupiter rail: daughter early beyond t_c (1.2 t_c)', max(0., jl['1.2']['lag_s_d01R']), 0.)
jup['lags'] = jl

# 4. Saturn rails of the proposed comparison -------------------------------------------------------
options = []
for launch, peri in S.windows():
    fb, yy = S.flyby(peri); sol, quant, T, yf = S.ray(yy[:3], yy[3:], fb['exit_epoch'], 4.)
    options.append((S.rail_summary(sol, quant, T, fb['exit_epoch'])['min_Jupiter_distance_AU'], launch, fb, yy))
_, _, fb, yy = min([z for z in options if z[0] > 2.5 and min(z[2]['ring_plane_crossings_RS']) > 3.], key=lambda z: z[1])
A, vA, ep = yy[:3], yy[3:], fb['exit_epoch']
E_helio = float(np.dot(vA, vA)/2 - S.MU/np.linalg.norm(A))
check('Saturn rails: heliocentric energy negative, so K < 0 in the solar field', max(0., E_helio), 0.)
d0 = vA/np.linalg.norm(vA); e1 = np.cross([0, 0, 1.], d0); e1 /= np.linalg.norm(e1); e2 = np.cross(d0, e1)
horizon = 9.
base = S.ray(A, vA, ep, 4., endtime=horizon)[0]
ts = np.linspace(1e-3, horizon, 4000)


def transverse_fields(eps):
    # centred differences in the two launch angles; in- and out-of-plane blocks decouple by ecliptic symmetry
    rays = {k: S.ray(A, vA, ep, 4., initial_direction=(d0 + sg*eps*e)/np.linalg.norm(d0 + sg*eps*e), endtime=horizon)[0]
            for k, (sg, e) in {'a+': (1, e1), 'a-': (-1, e1), 'b+': (1, e2), 'b-': (-1, e2)}.items()}
    Jin = []; Jout = []
    for t in ts:
        s0 = base.sol(t); ev = s0[3:6]/np.linalg.norm(s0[3:6])
        nin = np.cross([0, 0, 1.], ev); nin /= np.linalg.norm(nin); nout = np.cross(ev, nin)
        Jin.append(np.dot((rays['a+'].sol(t)[:3] - rays['a-'].sol(t)[:3])/(2*eps), nin))
        Jout.append(np.dot((rays['b+'].sol(t)[:3] - rays['b-'].sol(t)[:3])/(2*eps), nout))
    return np.array(Jin), np.array(Jout)


Jin, Jout = transverse_fields(1e-6)
Jin5, Jout5 = transverse_fields(1e-5)
conv = max(np.max(np.abs(Jin5 - Jin))/np.max(np.abs(Jin)), np.max(np.abs(Jout5 - Jout))/np.max(np.abs(Jout)))
check('Saturn rails: centred differences at 1e-5 and 1e-6 agree (relative)', conv, 1e-5)
zeros = lambda J: int(np.sum(np.sign(J[1:]) != np.sign(J[:-1])))
rails = {str(R): S.ray(A, vA, ep, R)[2] for R in (4., 3.)}
check('Saturn rails: no zero of the in-plane Jacobi field within 9 years', zeros(Jin), 0)
check('Saturn rails: no zero of the out-of-plane Jacobi field within 9 years', zeros(Jout), 0)
saturn = dict(heliocentric_energy_AU2_yr2=E_helio, rail_duration_years=rails, horizon_years=horizon,
              difference_convergence=float(conv), min_in_plane_over_t=float(np.min(np.abs(Jin)/ts)),
              min_out_of_plane_over_t=float(np.min(np.abs(Jout)/ts)),
              in_plane_zeros=zeros(Jin), out_of_plane_zeros=zeros(Jout))

out = dict(caveat='planar optical-metric model; the Earth-Moon and Jupiter cases are illustrations, '
                  'not mission designs', cr3bp_units=dict(L_km=LKM, T_days=TU/86400, V_km_s=VU/1e3, A_m_s2=AU_),
           earth_moon_L4=cr3bp, jupiter_hyperbolic=jup, saturn_rails=saturn, checks=checks)
(HERE/'probe_conjugate_points.json').write_text(json.dumps(out, indent=2))
print('checks', len(checks), 'PASS', sum(c['status'] == 'PASS' for c in checks))
for k, row in cr3bp.items():
    print('L4 ray %s deg: t_c %.3f d, peak thrust %.3g m/s^2, lags 0.8/1.2 t_c (1000 km): %+.2f / %+.2f s' %
          (k, row['t_c_days'], row['peak_thrust_m_s2'], row['lags']['0.8']['lag_s_d1000km'], row['lags']['1.2']['lag_s_d1000km']))
print('Jupiter: t_c %.3f d, pericentre %.2f R, swept %.0f deg, lags 0.8/1.2 t_c (0.1 R): %+.2f / %+.2f s' %
      (jup['t_c_days'], jup['pericentre_R'], jup['swept_deg'], jl['0.8']['lag_s_d01R'], jl['1.2']['lag_s_d01R']))
print('Saturn rails: E = %.3f AU^2/yr^2, durations %s yr, Jacobi zeros in/out of plane within 9 yr: %d/%d' %
      (E_helio, {k: round(v, 3) for k, v in rails.items()}, zeros(Jin), zeros(Jout)))
assert all(c['status'] == 'PASS' for c in checks)
