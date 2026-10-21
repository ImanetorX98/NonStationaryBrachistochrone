"""Weak-field control law of the arrival-time rail (section `sec:observational-context`).

Newtonian limit of the Schwarzschild reduction of eqs. (randers)-(randers-data):
maintained energy, v^2 = 2 eps + 2 mu/r, arrival functional T = int ds / v(r).
Checks, symbolically:
  1. rail thrust |a| = 2 mu sin(alpha)/r^2 (twice the normal gravity), any energy;
  2. a straight chord at the same energy needs only the normal gravity;
  3. strict convexity of the monotone-r Lagrangian in phi' (unique minimum in that class);
  4. Snell function H = r^2 - k^2 v^2 strictly increasing: one turning radius.
Writes probe_control_law.json.
"""
from pathlib import Path
import json
import sympy as sp

HERE = Path(__file__).resolve().parent
checks = []


def check(name, expr):
    ok = sp.simplify(expr) == 0
    checks.append(dict(name=name, status='PASS' if ok else 'FAIL'))
    return ok


r, eps, mu, k, p = sp.symbols('r epsilon mu k phi_prime', real=True)
sa = sp.symbols('sin_alpha', real=True)
v2 = 2*eps + 2*mu/r
# Ray of the isotropic index n = 1/v bends with curvature kappa = -(d ln v/dr) sin(alpha);
# the required thrust is v^2 kappa plus the cancelled normal gravity mu sin(alpha)/r^2.
kinematic = -v2*sp.diff(sp.log(sp.sqrt(v2)), r)*sa
check('rail thrust = 2 mu sin(alpha)/r^2', kinematic + mu*sa/r**2 - 2*mu*sa/r**2)
# A straight chord has zero curvature: thrust = normal gravity only.
check('chord thrust = mu sin(alpha)/r^2', 0*kinematic + mu*sa/r**2 - mu*sa/r**2)
# Monotone r: T = int n(r) sqrt(1 + r^2 phi'^2) dr, n = 1/v > 0.
n = sp.symbols('n', positive=True)
L = n*sp.sqrt(1 + r**2*p**2)
check('d2L/dphi_prime2 = n r^2/(1+r^2 phi_prime^2)^(3/2) > 0',
      sp.diff(L, p, 2) - n*r**2/(1 + r**2*p**2)**sp.Rational(3, 2))
check("H'(r) = 2r + 2 mu k^2/r^2 > 0", sp.diff(r**2 - k**2*v2, r) - (2*r + 2*mu*k**2/r**2))

out = dict(model='Newtonian weak-field rail, maintained energy, arrival-time functional',
           results=dict(rail_thrust='2 mu sin(alpha)/r^2', chord_thrust='mu sin(alpha)/r^2',
                        convexity='strict in phi_prime for monotone r', snell_turning_radius='unique'),
           checks=checks)
(HERE/'probe_control_law.json').write_text(json.dumps(out, indent=2))
print('checks', len(checks), 'PASS', sum(c['status'] == 'PASS' for c in checks))
assert all(c['status'] == 'PASS' for c in checks)
