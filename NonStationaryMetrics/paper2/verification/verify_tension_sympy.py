#!/usr/bin/env python3
"""Tension field of the optical submersion -- independent SymPy route.

Second computer-algebra route for the identities of section 2.2, so that the
claims tagged evidence level (I) there are exact zeros in two independent
systems (SymPy here, Mathematica in verify_submersion_link.wls) rather than
one.  Referee 1, major comment 10 and Protocol 4.

Checks, for a static conformal factor over a Schwarzschild seed:

  1. fibre mean curvature of the orbits of W = d_t
  2. horizontal gradient of log of the Perlick dilation
  3. the tension residual mu + grad_H log Lambda, and that it is strictly
     positive on the exterior -- so the projection is NOT a harmonic morphism
  4. the residual vanishes only in the limit f -> 0

Run:  python3 verify_tension_sympy.py
"""
import sys
import sympy as sp

r, M, E, s = sp.symbols("r M Ehat s", positive=True)
f = 1 - 2 * M / r

FAILURES = 0


def check(label, cond):
    global FAILURES
    ok = bool(cond)
    print(f"{'PASS  ' if ok else 'FAIL  '}{label}")
    if not ok:
        FAILURES += 1


print("=" * 58)
print("tension field of the optical submersion (SymPy)")
print("=" * 58)

# --- 1. fibre mean curvature -------------------------------------------
# g = -f dt^2 + dr^2/f + r^2 dOmega^2,  W = d_t,  |W| = sqrt(f)
# nabla_t d_t = Gamma^r_{tt} d_r = (1/2) f f' d_r
nabla_WW = sp.simplify(sp.Rational(1, 2) * f * sp.diff(f, r))
check("nabla_W W = (M/r^2) f d_r",
      sp.simplify(nabla_WW - M * f / r**2) == 0)

# unit fibre tangent Wh = W/sqrt(f); mu = nabla_Wh Wh  (already horizontal)
mu = sp.simplify(nabla_WW / f)
check("mu = (M/r^2) d_r", sp.simplify(mu - M / r**2) == 0)

# --- 2. dilation and its horizontal gradient ---------------------------
Lambda2 = E**2 / (f * (E**2 - f))           # Perlick fixed-energy factor
gradH = sp.simplify(f * sp.diff(sp.log(sp.sqrt(Lambda2)), r))   # g^{rr} = f
gradH_closed = -(M / r**2) * (E**2 - 2 * f) / (E**2 - f)
check("grad_H log Lambda closed form",
      sp.simplify(sp.together(gradH - gradH_closed)) == 0)

# --- 3. the residual ---------------------------------------------------
resid = sp.simplify(sp.together(mu + gradH))
resid_closed = (M / r**2) * f / (E**2 - f)
check("mu + grad_H log Lambda = (M/r^2) f/(Ehat^2 - f)",
      sp.simplify(sp.together(resid - resid_closed)) == 0)

check("residual is not identically zero", sp.simplify(resid) != 0)

# strict positivity on the exterior.  Put r = 2M(1+s) with s>0 and
# Ehat^2 = 1+Y with Y>0; the residual must become a ratio of polynomials in
# (s,Y) with no negative coefficient, which settles the sign with no numerics.
Y = sp.symbols("Y", positive=True)
resid_sY = sp.simplify(sp.together(
    resid_closed.subs(r, 2 * M * (1 + s)).subs(E**2, 1 + Y).subs(E, sp.sqrt(1 + Y))))
num, den = sp.fraction(sp.cancel(resid_sY))
check("residual reduces to s / [4M (1+s)^2 (1 + Y + sY)]",
      sp.simplify(sp.cancel(resid_sY - s / (4 * M * (1 + s)**2 * (1 + Y + s * Y)))) == 0)
neg = [c for poly in (sp.Poly(sp.expand(num), s, Y, M), sp.Poly(sp.expand(den), s, Y, M))
       for c in poly.coeffs() if c.is_number and c < 0]
check("no negative coefficient in numerator or denominator", not neg)

for rv, Ev in ((3, sp.Rational(12, 10)), (6, sp.Rational(14, 10)),
               (20, sp.Rational(11, 10))):
    val = sp.N(resid_closed.subs({M: 1, r: rv, E: Ev}))
    check(f"residual > 0 at r={rv}M, Ehat={float(Ev)}  ({val:.6e})", val > 0)

# --- 4. behaviour at the boundary --------------------------------------
check("residual -> 0 as r -> 2M",
      sp.limit(resid_closed.subs(M, 1), r, 2, "+") == 0)

print("=" * 58)
print("ALL CHECKS PASSED" if FAILURES == 0 else f"{FAILURES} CHECK(S) FAILED")
print("=" * 58)
sys.exit(1 if FAILURES else 0)
