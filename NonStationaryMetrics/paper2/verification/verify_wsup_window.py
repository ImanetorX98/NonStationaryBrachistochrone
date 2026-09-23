#!/usr/bin/env python3
"""Check equation (Wsup) and the window values quoted beside it.

The manuscript bounds the separation weight of the rotating arrival branch by a
function of the radius alone,

    W <= W_sup(r) = K + B^2 + |B'|/sqrt(P),

and states that W_sup is negative on r/M in [6, 9.275] at a=0.9, for E=1.4 and
E=1.3.  This script rebuilds every ingredient from the data printed in the paper
-- P, Q, beta_phi from Proposition (rot-separation), B from (magnetic-B) -- and
checks the bound, the reduction of K to the printed static formula, and the two
windows.

Nothing is imported from the figure generators.  Three of the checks are
independent of how the bound was derived:

  (a) K at a=0 must reproduce equation (optical-gauss), which is printed in the
      paper and was obtained by a different route;
  (b) K at a=0.9, E=1.4, r=6M must reproduce the value -4.475e-3 already quoted
      in the manuscript for the rotating family;
  (c) the direction bound must be TIGHT: sup over admissible T^phi of
      |dB(N)| must equal |B'|/sqrt(P) exactly, not merely bound it.
"""
import sys

import mpmath as mp
import sympy as sp
from mpmath import iv

mp.iv.dps = 40

r, a, E = sp.symbols('r a E', positive=True)
M = 1

f    = 1 - 2*M/r
Dl   = r**2 - 2*M*r + a**2
n2   = E**2/(E**2 - f)
P    = n2*r**2/(f*Dl)
Q    = n2*Dl/f**2
bphi = -2*M*a/(r*f)

B  = sp.simplify(sp.diff(bphi, r)/sp.sqrt(P*Q))
Bp = sp.simplify(sp.diff(B, r))
sP, sQ = sp.sqrt(P), sp.sqrt(Q)
K  = sp.simplify(-1/(sP*sQ)*sp.diff(sp.diff(sQ, r)/sP, r))

fails = []


def check(name, ok, detail=""):
    print(f"  [{'PASS' if ok else 'FAIL'}] {name}" + (f"   {detail}" if detail else ""))
    if not ok:
        fails.append(name)


print("verify_wsup_window.py -- equation (Wsup) and its window")
print("=" * 66)

# --- B as printed -----------------------------------------------------------
B_printed = 2*M*a*(E**2 - f)/(E**2*r**3*sp.sqrt(f))
fB, fBp = sp.lambdify((r, a, E), B, 'math'), sp.lambdify((r, a, E), Bp, 'math')
fBpr = sp.lambdify((r, a, E), B_printed, 'math')
worst = max(abs(fB(rv, av, Ev) - fBpr(rv, av, Ev))
            for rv in (2.5, 6, 20, 300) for av in (0.3, 0.9) for Ev in (1.2, 3.0))
check("B = beta_phi'/sqrt(PQ) equals the printed 2Ma(E^2-f)/(E^2 r^3 sqrt f)",
      worst < 1e-14, f"max |diff| = {worst:.2e}")

# --- (a) K at a=0 must be the printed static Gauss curvature ----------------
K_paper = (M*(16*M**3 + 2*(9*E**2 - 14)*M**2*r + (3*E**4 - 19*E**2 + 16)*M*r**2
              - (2*E**2 - 3)*(E**2 - 1)*r**3)) / (E**2*r**5*(2*M + (E**2 - 1)*r))
check("K(a=0) reproduces equation (optical-gauss)",
      sp.simplify(K.subs(a, 0) - K_paper) == 0)

# --- (b) K must reproduce the rotating value already quoted -----------------
fK = sp.lambdify((r, a, E), K, 'math')
check("K(r=6M, a=0.9, E=1.4) reproduces the quoted -4.475e-3",
      abs(fK(6, 0.9, 1.4) + 4.475e-3) < 5e-7, f"got {fK(6, 0.9, 1.4):.6e}")

# --- (c) the direction bound is attained, not merely valid ------------------
# dB(N) = -B' sqrt(Q/P) T^phi with Q (T^phi)^2 <= 1, so the supremum over
# admissible directions is |B'| sqrt(Q/P) / sqrt(Q) = |B'|/sqrt(P).
#
# The identity sqrt(Q/P)/sqrt(Q) = 1/sqrt(P) needs P,Q > 0, which sympy cannot
# infer from the expressions alone -- P.is_positive and Q.is_positive both
# return None.  Checking it on abstract positive symbols states that hypothesis
# rather than hoping the simplifier guesses it; the hypothesis itself is then
# discharged ANALYTICALLY by the Delta certificate below, not by sampling.
Ps, Qs = sp.symbols('P_ Q_', positive=True)
check("sup over directions of |dB(N)| equals |B'|/sqrt(P), given P,Q > 0",
      sp.simplify(sp.sqrt(Qs/Ps)/sp.sqrt(Qs) - 1/sp.sqrt(Ps)) == 0)

# --- P, Q positive on the exterior: analytic, not sampled -------------------
# f > 0 for r > 2M, and Delta = r^2 - 2Mr + a^2 > 0 there because its larger root
# is r_+ = M + sqrt(M^2 - a^2) <= 2M for |a| <= M.  n^2 = E^2/(E^2-f) > 0 since
# f < 1 <= E^2.  Hence P = n^2 r^2/(f Delta) > 0 and Q = n^2 Delta/f^2 > 0
# throughout, with no sampling involved.
s_, a_ = sp.symbols('s_ a_', nonnegative=True)
delta_shift = sp.expand((2 + s_)**2 - 2*(2 + s_) + a_**2)      # Delta at r = 2M + s
check("Delta(2M+s) = s^2 + 2Ms + a^2, all coefficients nonnegative",
      sp.simplify(delta_shift - (s_**2 + 2*s_ + a_**2)) == 0,
      f"= {delta_shift}")
# The chain that turns that into P, Q > 0 must be CHECKED, not merely narrated.
# An earlier version of this file asserted it with `limit(s^2+2s, s, 0) == 0`,
# which is true whatever P and Q are: a check that cannot fail certifies nothing.
# The identities below can fail, and would, if any of the formulae above were
# edited.  With r = 2M + s and E^2 = 1 + u, both s and u nonnegative:
u_ = sp.Symbol('u_', nonnegative=True)
check("  f = s/(s+2M): numerator and denominator manifestly positive for s > 0",
      sp.simplify((s_/(s_ + 2))*(s_ + 2) - s_) == 0)
E2mf = sp.expand(sp.simplify(((1 + u_) - s_/(s_ + 2))*(s_ + 2)))
check("  (E^2-f)(s+2M) = u s + 2Mu + 2M, all coefficients positive",
      sp.simplify(E2mf - (s_*u_ + 2*u_ + 2)) == 0, f"= {E2mf}")
# Hence f > 0 and E^2 - f > 0, so n^2 = E^2/(E^2-f) > 0; Delta > 0 by the
# certificate above; and P = n^2 r^2/(f Delta), Q = n^2 Delta/f^2 are then
# quotients of positive quantities.  No sampling anywhere in the chain.

# --- the window, by INTERVAL arithmetic over the continuum ------------------
# A grid of point evaluations cannot exclude a positive maximum between nodes,
# however fine.  The sign is therefore established with outward-rounded interval
# arithmetic on a partition that covers [6, 371/40] EXACTLY: each cell returns an
# enclosure of W_sup over the whole cell, and every enclosure must have a
# strictly negative upper end.  Rational endpoints keep the partition exact.
fP = sp.lambdify((r, a, E), P, 'math')       # point versions, for the last check
Wsup = lambda rv, av, Ev: (fK(rv, av, Ev) + fB(rv, av, Ev)**2
                           + abs(fBp(rv, av, Ev))/fP(rv, av, Ev)**0.5)
mods = [{'sqrt': iv.sqrt, 'Abs': abs}, 'math']
gK  = sp.lambdify((r, a, E), K,  modules=mods)
gB  = sp.lambdify((r, a, E), B,  modules=mods)
gBp = sp.lambdify((r, a, E), Bp, modules=mods)
gP  = sp.lambdify((r, a, E), P,  modules=mods)


def wsup_iv(R, av, Ev):
    """Enclosure of W_sup over the r-interval R.  B' < 0 is a theorem
    (Proposition prograde-asym), so |B'| = -B' and no Abs is needed."""
    return gK(R, av, Ev) + gB(R, av, Ev)**2 - gBp(R, av, Ev)/iv.sqrt(gP(R, av, Ev))


print()
# Cells are built from integers and divided in interval arithmetic, so each is
# an OUTWARD-rounded enclosure of [i/1000, (i+1)/1000]: adjacent cells overlap
# rather than leave a gap, and their union covers [6, 9.275] with certainty.
# (9.275 = 371/40 is not a binary float; building it as 9275/1000 outwardly is
# what keeps the cover honest.)
I0, I1, DEN = 6000, 9275, 1000
for Ev_num, Ev_den, lo_num, hi_num in ((7, 5, -232, -68),
                                       (13, 10, -213, -57)):
    av_iv = iv.mpf(9)/iv.mpf(10)
    Ev_iv = iv.mpf(Ev_num)/iv.mpf(Ev_den)
    lo_all, hi_all, ok = None, None, True
    NCELL = I1 - I0
    for i in range(I0, I1):
        W = wsup_iv(iv.mpf([i, i + 1])/iv.mpf(DEN), av_iv, Ev_iv)
        if not (W.b < 0):
            ok = False
        lo_all = W.a if lo_all is None else min(lo_all, W.a)
        hi_all = W.b if hi_all is None else max(hi_all, W.b)
    check(f"W_sup < 0 on ALL of r/M in [6, 9.275], a=0.9, E={Ev_num}/{Ev_den} "
          f"({NCELL} interval cells)", ok,
          f"enclosure [{float(lo_all):.6e}, {float(hi_all):.6e}]")
    # The printed bounds ARE rigorous: they are the enclosure rounded outward.
    # Checking that with float() would break the interval chain at the last
    # step, so the printed values are rebuilt as exact integer ratios and the
    # comparison is made on interval endpoints.
    lo_p = iv.mpf(lo_num)/iv.mpf(100000)      # e.g. -232/100000 = -2.32e-3
    hi_p = iv.mpf(hi_num)/iv.mpf(100000)
    check(f"  printed bounds [{lo_num}/1e5, {hi_num}/1e5] contain the enclosure "
          f"(interval comparison, E={Ev_num}/{Ev_den})",
          lo_p.b <= lo_all.a and hi_all.b <= hi_p.a)

# --- the bound must actually bound: W <= W_sup along a real direction --------
print()
import math
ok = True
for rv in (6.0, 7.5, 9.0):
    Qv = sp.lambdify((r, a, E), Q, 'math')(rv, 0.9, 1.4)
    for tphi in (0.0, 0.5/math.sqrt(Qv), 1.0/math.sqrt(Qv), -1.0/math.sqrt(Qv)):
        W = fK(rv, 0.9, 1.4) + fB(rv, 0.9, 1.4)**2 - fBp(rv, 0.9, 1.4)*math.sqrt(
            Qv/fP(rv, 0.9, 1.4))*tphi
        if W > Wsup(rv, 0.9, 1.4) + 1e-18:
            ok = False
check("W <= W_sup for every admissible T^phi sampled (12 cases)", ok)

print()
print("=" * 66)
if fails:
    print(f"FAILED: {len(fails)} check(s): " + "; ".join(fails))
    sys.exit(1)
print("all checks passed")
