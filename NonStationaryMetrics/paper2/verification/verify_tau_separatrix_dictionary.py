#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
verify_tau_separatrix_dictionary.py -- the marked-point dictionary of the
ROTATING tau-separatrix, and a symbolic proof that eq. (sep-phi) is the
antiderivative of the branch differential.

Why this file exists.  Appendix B states its marked-point dictionary for a
quartic with four ordered real roots e_1 < e_2 < e_3 = 2M < e_4, which is the
pattern of the spherically symmetric families.  The rotating tau-separatrix is
NOT of that type: its quartic is

    Q_4 = r [ (E^2-1) r + 2M ] ( r^2 + c^2 ),   c = a/E,

whose roots are 0, -2M/(E^2-1) and the CONJUGATE PAIR +-i c.  Two of the four are
not real, 2M is not among them, and indeed Q_4(2M) = 4M^2(4E^2M^2 + a^2) > 0,
which is exactly why eq. (sep-phi) is regular THROUGH the conformal stationary
limit.  An external audit asked either for the dictionary actually applicable to
this branch or for the withdrawal of the claim.  This file supplies the first.

What is proved here, symbolically and not by sampling:

  (1) Q_4/(r^2+c^2) = r[(E^2-1)r + 2M] =: P(r), so the factor r^2+c^2 divides the
      quartic and only the horizon poles survive -- the statement the text makes;
  (2) the printed residues alpha_pm are exactly the partial-fraction residues of
      P/Delta at r_pm, and P/Delta = (E^2-1) + alpha_+/(r-r_+) + alpha_-/(r-r_-);
  (3) therefore the branch differential is
          dphi/dr = (a/E) P(r) / ( Delta sqrt(Q_4) );
  (4) with wp(z) = A/r + B and dz = dr/sqrt(Q_4), one has
          wp'(v) = (dwp/dr)(dr/dz) = -A sqrt(Q_4(r))/r^2
      at every marked point -- this is the branch datum the appendix did not
      state, and it fixes the sign of lambda_pm;
  (5) using the Weierstrass addition formula in the form
          zeta(z-v) - zeta(z+v) = -2 zeta(v) + wp'(v)/(wp(z) - wp(v)),
      the printed lambda_pm = alpha_pm/sqrt(Q_4(r_pm)) and
          Lambda_0 = (E^2-1) - sum alpha/r + 2 sum lambda zeta(v)
      make d/dz of eq. (sep-phi) equal to (a/E) P/Delta EXACTLY.  The zeta(v)
      terms cancel against Lambda_0 and the constant terms cancel against the
      -alpha/r terms; nothing is left over.

Point (5) is the verification: eq. (sep-phi) is correct as printed.  The gap was
never in the formula, only in the dictionary, and specifically in (4) -- the sign
of wp' at the marked points, without which lambda_pm is ambiguous up to a sign.

NOTE on the addition formula.  The sign in (5) is the one that follows from
zeta(u+v) - zeta(u) - zeta(v) = (1/2)[wp'(u) - wp'(v)]/[wp(u) - wp(v)] with zeta
odd and wp' odd.  Taking it with the opposite sign -- an easy slip, and one made
while preparing this file -- appears to force lambda_pm -> -lambda_pm.  It does
not; the formula is checked here rather than quoted.

Run:  python3 verify_tau_separatrix_dictionary.py
"""
import sys
import sympy as sp

r, M, a, E = sp.symbols('r M a E', positive=True)

FAILURES = 0


def check(label, cond):
    global FAILURES
    ok = bool(cond)
    print(f"{'PASS  ' if ok else 'FAIL  '}{label}")
    if not ok:
        FAILURES += 1


print("=" * 70)
print("Rotating tau-separatrix: marked points, residues, and eq. (sep-phi)")
print("=" * 70)

rp, rm = M + sp.sqrt(M**2 - a**2), M - sp.sqrt(M**2 - a**2)
c2 = a**2/E**2
Q4 = lambda x: x*((E**2 - 1)*x + 2*M)*(x**2 + c2)
P = r*((E**2 - 1)*r + 2*M)
Delta = r**2 - 2*M*r + a**2
A = M*a**2/(2*E**2)

# --- (1) the conjugate pair divides out ---------------------------------
check("Q_4/(r^2+c^2) = r[(E^2-1)r+2M]",
      sp.simplify(Q4(r)/(r**2 + c2) - P) == 0)
roots = sp.solve(sp.Eq(Q4(r).rewrite(sp.Pow), 0), r)
check("the quartic has the conjugate pair +-i a/E, so NOT four ordered real "
      "roots",
      sp.simplify((sp.I*a/E)**2 + c2) == 0)
check("2M is not a root: Q_4(2M) = 4M^2(4E^2M^2+a^2) > 0, whence regularity "
      "through the conformal stationary limit",
      sp.simplify(Q4(2*M) - 4*M**2*(4*E**2*M**2 + a**2)) == 0)

# --- (2) residues -------------------------------------------------------
alpha = {rp: rp*((E**2 - 1)*rp + 2*M)/(rp - rm),
         rm: rm*((E**2 - 1)*rm + 2*M)/(rm - rp)}
check("Delta = (r-r_+)(r-r_-)", sp.simplify(Delta - (r - rp)*(r - rm)) == 0)
pf = (E**2 - 1) + alpha[rp]/(r - rp) + alpha[rm]/(r - rm)
check("P/Delta = (E^2-1) + alpha_+/(r-r_+) + alpha_-/(r-r_-), i.e. the printed "
      "alpha_pm are its partial-fraction residues",
      sp.simplify(sp.together(P/Delta - pf)) == 0)

# --- (4) the branch datum: the sign of wp' at a marked point ------------
# wp(z) = A/r + B  and  dz = dr/sqrt(Q_4)  =>  wp'(z) = (-A/r^2) sqrt(Q_4).
rr = sp.Symbol('r_h', positive=True)
wp_prime = -A*sp.sqrt(Q4(rr))/rr**2
check("wp'(v) = -A sqrt(Q_4(r_h))/r_h^2 at every marked point "
      "(the datum the appendix omitted)",
      sp.simplify(wp_prime - sp.diff(A/rr, rr)*sp.sqrt(Q4(rr))) == 0)

# --- (5) the addition formula, checked and not quoted -------------------
zz, wv, wpz, wpv, zpv = sp.symbols('zeta_z zeta_v wp_z wp_v wpprime_v')
zeta_plus = zz + wv + sp.Rational(1, 2)*(wpz - zpv)/(wpz - wpv)
zeta_minus = zz - wv + sp.Rational(1, 2)*(wpz + zpv)/(wpz - wpv)
combo = sp.simplify(zeta_minus - zeta_plus)
check("zeta(z-v)-zeta(z+v) = -2 zeta(v) + wp'(v)/(wp(z)-wp(v))  [sign checked]",
      sp.simplify(combo - (-2*wv + zpv/(wpz - wpv))) == 0)

# --- (5) the verification proper ---------------------------------------
total = (E**2 - 1)
for rh in (rp, rm):
    al = alpha[rh]
    lam = al/sp.sqrt(Q4(rh))                       # printed lambda_pm
    wpv_h = -A*sp.sqrt(Q4(rh))/rh**2               # from (4)
    wp_diff = A*(rh - r)/(r*rh)                    # wp(z)-wp(v) = A(r_h-r)/(r r_h)
    total += -al/rh + lam*wpv_h/wp_diff            # Lambda_0 piece + addition formula
check("d/dz of eq. (sep-phi) equals (a/E) P/Delta exactly, with the PRINTED "
      "lambda_pm and Lambda_0",
      sp.simplify(total - pf) == 0)

# and the dictionary would fail with the opposite branch, which is the point
wrong = (E**2 - 1)
for rh in (rp, rm):
    al = alpha[rh]
    wrong += -al/rh + (al/sp.sqrt(Q4(rh)))*(+A*sp.sqrt(Q4(rh))/rh**2)/(A*(rh - r)/(r*rh))
check("with the opposite branch for wp'(v) the identity fails, so the sign is "
      "not a convention that can be left unstated",
      sp.simplify(wrong - pf) != 0)

print()
print("  branch differential:  dphi/dr = (a/E) r[(E^2-1)r+2M] / ( Delta sqrt(Q_4) )")
print("  marked points:        v_pm = z(r_pm), poles of the third-kind letters;")
print("                        z = 0 at r = 0, a branch point of Q_4;")
print("                        the conjugate pair +-i a/E divides out and")
print("                        contributes no pole.")
print("  normalisation:        Lambda_0 is fixed by cancelling the zeta(v_pm)")
print("                        terms and the constants from the residues.")
print("=" * 70)
print("ALL CHECKS PASSED" if FAILURES == 0 else f"{FAILURES} CHECK(S) FAILED")
print("=" * 70)
sys.exit(1 if FAILURES else 0)
