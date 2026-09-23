#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
verify_prograde_asymmetry.py -- the prograde/retrograde asymmetry of the
separation weight, proved rather than sampled.

WHAT THIS REPLACES.  The manuscript and the response reported, from one
integrated extremal at M=1, a=0.9, E=1.4, that W - K = B^2 + dB(N) is positive at
every sample, and separately that the two azimuthal senses give different minima
of W (-2.27e-03 prograde against -6.59e-03 retrograde).  Both were numerical
observations on a finite family, and were labelled as such.

They are consequences of a sign, and the sign is provable.

THE ARGUMENT, in three steps.

  (1) B has the closed form  B = 2 M a (E^2 - f)/(E^2 r^3 sqrt f),  f = 1 - 2M/r.
      Differentiating,

          B'(r) = -2 M a * [ 3(E^2-1)r^2 - (5E^2-13) M r - 14 M^2 ]
                  / ( E^2 r^{9/2} (r-2M)^{3/2} ) ,

      whose denominator is positive on the exterior.

  (2) The bracket is positive there.  Substituting r = 2M(1+x) with x >= 0 and
      E^2 = 1+u with u > 0 turns it into

          bracket/M^2 = 12 u x^2 + 14 u x + 2 u + 16 x + 2 ,

      every coefficient positive, hence strictly positive.  Therefore

          B'(r) < 0   on the whole exterior, for every E > 1 and every a > 0.

      This is a certificate, not a scan: the positivity is read off the
      coefficients, so it holds uniformly in (r, E) together, and at the
      boundary x = 0 the constant term 2 keeps it away from zero.

  (3) With N = J T and the conventions of verify_rotating_separation.py,
      N^r = -sqrt(Q/P) T^phi, so

          dB(N) = B'(r) N^r = -B'(r) sqrt(Q/P) T^phi ,

      and since B' < 0,   sign dB(N) = sign T^phi.

CONSEQUENCE.  For every PROGRADE extremal (T^phi > 0), dB(N) > 0, hence

          W - K = B^2 + dB(N) > 0   pointwise, with no sampling,

so the separation weight is strictly less negative than the Gauss curvature of
the same Riemannian metric.  For retrograde extremals dB(N) < 0 and the sign of
W - K is a genuine competition with B^2, which is why the two senses differ.

That is the mechanism behind the asymmetry, and it is the same parity that the
separatrix charges J_c^pm exhibit: b is odd in a, so B is odd in a, and reversing
the sense of travel has the same effect on dB(N) as reversing the spin.

WHAT IS STILL NOT A THEOREM.  The absence of conjugate points.  W < 0 on an arc
gives no focusing there, but W - K > 0 alone does not bound W away from zero, and
the retrograde case is not covered at all.  Conjugate points still have to be
located, and the runs in verify_rotating_separation.py remain an exploration.

Run:  python3 verify_prograde_asymmetry.py
"""
import sys
import sympy as sp

r, M, a, E, x, u = sp.symbols('r M a E x u', positive=True)

FAILURES = 0


def check(label, cond):
    global FAILURES
    ok = bool(cond)
    print(f"{'PASS  ' if ok else 'FAIL  '}{label}")
    if not ok:
        FAILURES += 1


print("=" * 72)
print("Prograde/retrograde asymmetry of W - K, from the sign of B'")
print("=" * 72)

f = 1 - 2*M/r
B = 2*M*a*(E**2 - f)/(E**2*r**3*sp.sqrt(f))

# --- (1) the derivative, and its denominator -----------------------------
Bp = sp.simplify(sp.diff(B, r))
bracket = 3*(E**2 - 1)*r**2 - (5*E**2 - 13)*M*r - 14*M**2
Bp_form = -2*M*a*bracket/(E**2*r**sp.Rational(9, 2)*(r - 2*M)**sp.Rational(3, 2))
check("B'(r) = -2 M a [3(E^2-1)r^2 - (5E^2-13)Mr - 14M^2] / "
      "(E^2 r^{9/2} (r-2M)^{3/2})",
      sp.simplify(sp.expand(sp.together(Bp - Bp_form))) == 0)

# --- (2) the positivity certificate for the bracket ----------------------
g = sp.expand(sp.simplify(bracket.subs({r: 2*M*(1 + x)}).subs(E**2, 1 + u)))
poly = sp.Poly(sp.expand(g/M**2), x, u)
coeffs = poly.coeffs()
check(f"bracket/M^2 = 12ux^2 + 14ux + 2u + 16x + 2, coefficients {sorted(set(coeffs))} "
      "all positive",
      sp.simplify(poly.as_expr() - (12*u*x**2 + 14*u*x + 2*u + 16*x + 2)) == 0
      and all(c > 0 for c in coeffs))
check("the constant term is 2, so the bracket does not degenerate at r = 2M",
      sp.simplify(poly.as_expr().subs({x: 0, u: 0})) == 2)
check("hence B' < 0 strictly on r > 2M, for every E > 1 and every a > 0 "
      "(denominator positive there)",
      all(sp.N(Bp.subs({M: 1, a: av, E: ev, r: rv})) < 0
          for av in (sp.Rational(1, 10), sp.Rational(9, 10))
          for ev in (sp.Rational(101, 100), sp.Rational(7, 5), 3)
          for rv in (sp.Rational(21, 10), 3, 10, 100)))

# --- (3) the consequence for dB(N) ---------------------------------------
# N^r = -sqrt(Q/P) T^phi  =>  dB(N) = B' N^r = -B' sqrt(Q/P) T^phi.
# B' < 0  =>  sign dB(N) = sign T^phi.
n2 = E**2/(E**2 - f)
P = n2*r**2/(f*(r**2 - 2*M*r + a**2))
Q = n2*(r**2 - 2*M*r + a**2)/f**2
Tphi = sp.symbols('T_phi', real=True)
dBN = Bp*(-sp.sqrt(Q/P)*Tphi)
sub = {M: 1, a: sp.Rational(9, 10), E: sp.Rational(7, 5)}
check("prograde (T^phi > 0): dB(N) > 0 at sampled radii",
      all(sp.N(dBN.subs(sub).subs({Tphi: 1, r: rv})) > 0
          for rv in (sp.Rational(21, 10), 3, 6, 20, 100)))
check("retrograde (T^phi < 0): dB(N) < 0 at the same radii",
      all(sp.N(dBN.subs(sub).subs({Tphi: -1, r: rv})) < 0
          for rv in (sp.Rational(21, 10), 3, 6, 20, 100)))
check("therefore W - K = B^2 + dB(N) > 0 for EVERY prograde extremal, "
      "with no sampling over curves",
      all(sp.N((B**2 + dBN).subs(sub).subs({Tphi: 1/sp.sqrt(Q.subs(sub).subs(r, rv)),
                                            r: rv})) > 0
          for rv in (3, 6, sp.Rational(9275, 1000), 20)))

# --- parity, tying it to the separatrix asymmetry ------------------------
check("B is odd in a, so reversing the sense of travel acts on dB(N) exactly "
      "as reversing the spin",
      sp.simplify(B.subs(a, -a) + B) == 0)

print()
print("  proved:   B' < 0 on the exterior  =>  sign dB(N) = sign T^phi")
print("            prograde  ->  W - K > 0 pointwise, for every such extremal")
print("            retrograde ->  dB(N) < 0; the sign of W - K is a competition")
print("  NOT proved: absence of conjugate points, in either sense.")
print("=" * 72)
print("ALL CHECKS PASSED" if FAILURES == 0 else f"{FAILURES} CHECK(S) FAILED")
print("=" * 72)
sys.exit(1 if FAILURES else 0)
