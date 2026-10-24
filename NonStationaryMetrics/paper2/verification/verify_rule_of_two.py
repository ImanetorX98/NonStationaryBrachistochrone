"""verify_rule_of_two.py

Relativistic rule of two (`prop:rule-of-two`, Solar-System paragraph of `sec:observational-context`).
In a static spacetime g = -f dt^2 + gamma, the arrival-time rail at maintained charge Ehat = f dt/dtau has proper
acceleration orthogonal to its velocity relative to the static observers, of magnitude 2 g_perp, with g_perp the
component normal to the spatial path of the static observers' acceleration grad(ln sqrt f).

1. Symbolic step: with v^2 = 1 - f/Ehat^2, n = 1/(v sqrt f), kappa = d_perp ln n and gamma_L^2 = Ehat^2/f,
   gamma_L^2 (v^2 kappa + g_perp) = 2 g_perp, with g_perp = d_perp ln sqrt f (f and its normal derivative generic).
2. Independent check: exact four-acceleration from the Christoffel symbols of equatorial Schwarzschild along the rail of
   eq:randers at a = 0 (Ehat = 1, 1.3; r_min = 10M, 50M, 1e5 M), compared with 2 g_s sin(alpha),
   g_s = M/(r^2 sqrt f), alpha the angle between path and radial direction in the static frame; and zero Killing power
   g(b, d_t) = 0.
3. Mutation: replacing the factor 2 by 1.9 makes check 2 fail.
Coverage: full for the identity in step 1; sampled for step 2 (it checks the kinematic decomposition used in step 1).

Run: python3 verify_rule_of_two.py   (exit status 1 if a check fails)
"""
import sys
import sympy as sp
import mpmath as mp

checks = []


def report(name, ok):
    checks.append((name, bool(ok)))
    print(('PASS ' if ok else 'FAIL ') + name)


# 1. symbolic identity
F, dF, E = sp.symbols('f df_perp Ehat', positive=True)
v2 = 1 - F/E**2
lnn = -sp.log(sp.sqrt(v2)) - sp.log(sp.sqrt(F))
kappa = sp.diff(lnn, F)*dF                 # d_perp ln n by the chain rule
gperp = sp.diff(sp.log(sp.sqrt(F)), F)*dF  # d_perp ln sqrt f
bperp = (E**2/F)*(v2*kappa + gperp)
report('gamma^2 (v^2 kappa + g_perp) = 2 g_perp', sp.simplify(bperp - 2*gperp) == 0)

# 2. exact four-acceleration along the Schwarzschild rail
r, j = sp.symbols('r j', positive=True)
f = 1 - 2/r
g = sp.diag(-f, 1/f, r**2)
gi = g.inv()
X = [sp.Symbol('t'), r, sp.Symbol('p')]
Gam = [[[sp.simplify(sum(gi[i, l]*(sp.diff(g[l, k], X[h]) + sp.diff(g[l, h], X[k]) - sp.diff(g[k, h], X[l]))/2
                         for l in range(3))) for h in range(3)] for k in range(3)] for i in range(3)]
Dl = r**2 - 2*r
P = E**2*r**2/(f*Dl*(E**2 - f)); Q = E**2*Dl/(f**2*(E**2 - f))      # eq:randers-data at a = 0, frozen A = 1
dpdr = sp.sqrt(P)/sp.sqrt(Q**2/j**2 - Q)                            # outgoing branch, j = sqrt(Q(r_min))
ut = E/f
ur = sp.sqrt((f*ut**2 - 1)/(1/f + r**2*dpdr**2))
up = dpdr*ur
U = [ut, ur, up]
acc = [ur*sp.diff(U[i], r) + sum(Gam[i][k][h]*U[k]*U[h] for k in range(3) for h in range(3)) for i in range(3)]
b2 = sum(g[i, i]*acc[i]**2 for i in range(3))
gs = 1/(r**2*sp.sqrt(f))
sina = r*dpdr/sp.sqrt(1/f + r**2*dpdr**2)
report('zero Killing power g(b, d_t) = 0 (symbolic)', sp.simplify(sum(g[0, k]*acc[k] for k in range(3))) == 0)
num = sp.lambdify((r, E, j), [b2, gs, sina], 'mpmath')
mp.mp.dps = 40


def worst(factor):
    w = mp.mpf(0)
    for rmin, Ev in [(10, 1), (10, mp.mpf('1.3')), (50, mp.mpf('1.05')), (10**5, 1)]:
        jv = mp.sqrt(Q.subs({r: rmin, E: Ev}))
        for x in ('1.0001', '1.3', '2'):
            rv = rmin*mp.mpf(x)
            bb, gg, ss = num(rv, mp.mpf(Ev), mp.sqrt(jv**2))
            w = max(w, abs(mp.sqrt(bb)/(factor*gg*ss) - 1))
    return w


w2 = worst(2)
print('largest relative deviation of |b| from 2 g_s sin(alpha):', mp.nstr(w2, 3))
report('|b| = 2 g_s sin(alpha) along the rails (40 digits)', w2 < mp.mpf('1e-25'))
report('mutation: factor 1.9 is rejected', worst(mp.mpf('1.9')) > mp.mpf('1e-3'))

npass = sum(ok for _, ok in checks)
print(f'{npass}/{len(checks)} checks passed')
sys.exit(0 if npass == len(checks) else 1)
