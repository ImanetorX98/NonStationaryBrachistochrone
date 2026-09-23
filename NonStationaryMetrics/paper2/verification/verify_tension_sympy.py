#!/usr/bin/env python3
"""Tension field of the optical submersion -- rebuilt from the connections.

Section 2.2, eq. (tension-residual).  Referee 1, major comment 10 and Protocol 4.

WHAT CHANGED, AND WHY IT MATTERS.  An earlier version of this file checked the
sum mu + grad_H log Lambda and called it the tension field.  That sum is correct
arithmetic, but it is the *Riemannian* assembly tau = -dpi(mu + grad_H log Lambda),
and our total space is Lorentzian with a timelike fibre: the vertical direction
enters the trace with the opposite sign.  The Mathematica route in
verify_submersion_link.wls shared the same premise, so the two computer-algebra
systems agreed on a wrongly identified object rather than cross-checking it.

This file therefore computes the tension field of the map directly from the
Christoffel symbols of the total metric g and of the base metric a,

    tau^i(pi) = g^{mu nu} [ -Gamma^{i}_{mu nu}(g)
                            + Gamma^i_{jk}(a) d_mu pi^j d_nu pi^k ],

for pi(t,r,theta,phi) = (r,theta,phi), and only then compares with the two
ingredients.  The discrepancy is exactly 2M/r^2.

Checks, for a static conformal factor over a Schwarzschild seed:

  1. fibre mean curvature of the orbits of W = d_t                  (unchanged)
  2. horizontal gradient of log of the Perlick dilation             (unchanged)
  3. the Lorentzian tension field, from the connections             (new)
  4. it differs from the Riemannian assembly by 2M/r^2              (new)
  5. it is not identically zero -- so the projection is NOT a harmonic
     morphism -- but it is not of one sign either: it tends to 1/(2M) at the
     stationary limit and vanishes on r = 10M when Ehat^2 = 6/5       (new)

Run:  python3 verify_tension_sympy.py
"""
import sys
import sympy as sp

t, r, th, ph, M, E = sp.symbols("t r theta varphi M Ehat", positive=True)
f = 1 - 2 * M / r

FAILURES = 0


def check(label, cond):
    global FAILURES
    ok = bool(cond)
    print(f"{'PASS  ' if ok else 'FAIL  '}{label}")
    if not ok:
        FAILURES += 1


def christoffel(gm, x):
    n = len(x)
    gi = gm.inv()
    return [[[sp.simplify(sum(
        gi[A, d]*(sp.diff(gm[d, b], x[c]) + sp.diff(gm[d, c], x[b])
                  - sp.diff(gm[b, c], x[d])) for d in range(n))/2)
        for c in range(n)] for b in range(n)] for A in range(n)]


print("=" * 62)
print("Lorentzian tension field of the optical submersion (SymPy)")
print("=" * 62)

# --- 1. fibre mean curvature -------------------------------------------
nabla_WW = sp.simplify(sp.Rational(1, 2) * f * sp.diff(f, r))
check("nabla_W W = (M/r^2) f d_r", sp.simplify(nabla_WW - M * f / r**2) == 0)
mu = sp.simplify(nabla_WW / f)
check("mu = (M/r^2) d_r", sp.simplify(mu - M / r**2) == 0)

# --- 2. dilation and its horizontal gradient ---------------------------
Lambda2 = E**2 / (f * (E**2 - f))
gradH = sp.simplify(f * sp.diff(sp.log(sp.sqrt(Lambda2)), r))
check("grad_H log Lambda closed form",
      sp.simplify(sp.together(gradH + (M/r**2)*(E**2 - 2*f)/(E**2 - f))) == 0)

# --- 3. the tension field, from the connections ------------------------
X4 = [t, r, th, ph]
g4 = sp.diag(-f, 1/f, r**2, r**2*sp.sin(th)**2)
X3 = [r, th, ph]
a3 = Lambda2*sp.diag(1/f, r**2, r**2*sp.sin(th)**2)
G4, G3 = christoffel(g4, X4), christoffel(a3, X3)
gi4 = g4.inv()

tau = []
for i in range(3):
    vert = -sum(gi4[m, n]*G4[i + 1][m][n] for m in range(4) for n in range(4))
    horiz = sum(gi4[j + 1, k + 1]*G3[i][j][k] for j in range(3) for k in range(3))
    tau.append(sp.simplify(vert + horiz))

tau_closed = M/r**2*(2*E**2 - 3*f)/(E**2 - f)
check("tau^r = (M/r^2)(2 Ehat^2 - 3f)/(Ehat^2 - f)",
      sp.simplify(tau[0] - tau_closed) == 0)
check("tau^theta = tau^phi = 0", tau[1] == 0 and tau[2] == 0)

# --- 4. the discrepancy with the Riemannian assembly -------------------
riemannian = sp.simplify(-(mu + gradH))
check("tau^r - [-(mu + grad_H log Lambda)] = 2M/r^2",
      sp.simplify(tau[0] - riemannian - 2*M/r**2) == 0)

# --- 5. non-vanishing, but not of one sign -----------------------------
check("tau is not identically zero", sp.simplify(tau[0]) != 0)
check("tau^r -> 1/(2M) as r -> 2M",
      sp.simplify(sp.limit(tau_closed, r, 2*M, "+") - 1/(2*M)) == 0)
root = sp.solve(sp.simplify(tau_closed.subs(E, sp.sqrt(sp.Rational(6, 5)))), r)
check(f"tau^r vanishes at r = 10M when Ehat^2 = 6/5  (got {root})",
      any(sp.simplify(x - 10*M) == 0 for x in root))
check("so the old claims 'strictly positive on the exterior' and 'vanishes only "
      "as f -> 0' are both false",
      sp.N(tau_closed.subs({M: 1, r: 3, E: sp.sqrt(sp.Rational(6, 5))})) > 0
      and sp.N(tau_closed.subs({M: 1, r: 20, E: sp.sqrt(sp.Rational(6, 5))})) < 0)
check("the sign change is confined to 1 < Ehat^2 < 3/2: for Ehat^2 >= 3/2 the "
      "numerator 2 Ehat^2 - 3f has no zero with f < 1",
      sp.solve(sp.Eq(2*sp.Rational(3, 2) - 3*f, 0), r) == []
      and sp.solve(sp.Eq(2*2 - 3*f, 0), r) == []
      and sp.solve(sp.Eq(2*sp.Rational(6, 5) - 3*f, 0), r) == [10*M])

# --- 6. a genuinely independent route: the composition law ---------------
# tau is what makes Box_g(u o pi) - Lambda^2 (Delta_a u) o pi fail to vanish, for
# ANY u on the base.  Testing that exercises tau through its defining role rather
# than recomputing the same Christoffel contraction, and it is the argument used
# in the text to settle harmonicity directly.  Both sides are written out for a
# radial u, where sqrt|g| = r^2 sin(theta) and sqrt|a| = Lambda^3 r^2/sqrt(f).
u = sp.Function("u")(r)
up = sp.diff(u, r)
L = sp.sqrt(Lambda2)

box_g = sp.together(sp.diff(r**2*f*up, r)/r**2)
lap_a = sp.together(sp.sqrt(f)/(L**3*r**2)*sp.diff(L*r**2*sp.sqrt(f)*up, r))
check("Box_g(u o pi) - Lambda^2 (Delta_a u) o pi = tau^r u' for arbitrary radial u",
      sp.simplify(sp.expand(box_g - Lambda2*lap_a - tau_closed*up)) == 0)

# the concrete base-harmonic function of the text, through its derivative: a
# radial u is a-harmonic exactly when sqrt|a| a^{rr} u' is constant in r.
uh = 1/(L*r**2*sp.sqrt(f))
flux = sp.simplify(L**3*r**2/sp.sqrt(f) * (f/Lambda2) * uh)
check("u' = 1/(Lambda r^2 sqrt f) gives a constant radial flux, so Delta_a u = 0",
      sp.simplify(sp.diff(flux, r)) == 0)
check("hence Box_g(u o pi) = tau^r u', not identically zero",
      sp.simplify(tau_closed*uh) != 0)

print("=" * 62)
print("ALL CHECKS PASSED" if FAILURES == 0 else f"{FAILURES} CHECK(S) FAILED")
print("=" * 62)
sys.exit(1 if FAILURES else 0)
