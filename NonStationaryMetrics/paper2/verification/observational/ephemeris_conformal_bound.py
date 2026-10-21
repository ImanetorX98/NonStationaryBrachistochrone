"""Solar-System bound on a locally cosmological conformal factor (section `sec:observational-context`).

Slow-motion limit of the timelike action of g = A(eta)^2 g_Schw with constant seed mass m:
with cosmic time dt = A deta and physical radius R = A r, to first order in m and v^2,
    -A sqrt(f - rdot^2/f - r^2 phidot^2) deta = [-1 + m A/R + ((Rdot - H R)^2 + R^2 phidot_t^2)/2] dt,
so the Newtonian parameter is m A and grows at the rate H = Adot/A. Under the cosmological
identification H = H0, against the EPM2019 value of d(GM_sun)/dt / GM_sun
(Pitjeva et al. 2021, A&A 647, A141, eq. (1)). Writes ephemeris_conformal_bound.json.
"""
from pathlib import Path
import json
import sympy as sp

HERE = Path(__file__).resolve().parent
m, r, A, rd, pd, e, R, Rd, H, Pt = sp.symbols('m r A rdot phidot epsilon R Rdot H Phidot_t', real=True)
f = 1 - 2*m/r
lag = -A*sp.sqrt(f - rd**2/f - r**2*pd**2)
lag = lag.subs({m: e*m, rd: sp.sqrt(e)*rd, pd: sp.sqrt(e)*pd})
first = sp.series(lag, e, 0, 2).removeO().subs(e, 1)
# deta = dt/A, r = R/A, dr/deta = Rdot - H R, dphi/deta = A dphi/dt
per_dt = sp.expand(first/A).subs({r: R/A, rd: Rd - H*R, pd: A*Pt}, simultaneous=True)
target = -1 + m*A/R + ((Rd - H*R)**2 + R**2*Pt**2)/2
residual = sp.simplify(per_dt - target)

MPC = 3.0856775814913673e22; YEAR = 365.25*86400
H0 = 67.4e3/MPC*YEAR                     # per year, Planck 2018 value used in eq. (obs-drift)
ephemeris, sigma3 = -10.2e-14, 1.4e-14   # per year, 3 sigma
out = dict(slow_motion_residual=str(residual), newtonian_parameter='m A',
           H0_per_year=H0, EPM2019_dGM_over_GM_per_year=[ephemeris, sigma3],
           ratio_H0_to_ephemeris=H0/abs(ephemeris), H0_in_ephemeris_3sigma_units=(H0 - ephemeris)/sigma3,
           checks=[dict(name='slow-motion limit gives m A', status='PASS' if residual == 0 else 'FAIL'),
                   dict(name='H0 outside the EPM2019 interval', status='PASS' if H0 > ephemeris + sigma3 else 'FAIL')])
(HERE/'ephemeris_conformal_bound.json').write_text(json.dumps(out, indent=2))
print('slow-motion residual:', residual)
print('H0 = %.3e /yr; EPM2019 (%.1f +- %.1f)e-14 /yr; H0/|value| = %.0f' % (H0, ephemeris*1e14, sigma3*1e14, H0/abs(ephemeris)))
assert all(c['status'] == 'PASS' for c in out['checks'])
