# -*- coding: utf-8 -*-
"""
Numerical validation of the sigma/zeta ingredient of the closed form
of the separatrix (kerr_separatrix_weierstrass.py).

The derivation chain has two links:
 (i)  reduction: dphi = (a/E)[(E^2-1) + a_+/(r-r_+) + a_-/(r-r_-)] dz
      -> already validated to 5.6e-17 (quadrature in r vs z);
 (ii) classical antiderivative (Whittaker-Watson 20.53):
      int dz/(P(z) - P(v)) = (1/P'(v)) [ 2 zeta(v) z
                                        + ln sigma(z-v)/sigma(z+v) ]
      -> validated HERE on a test lattice, with sigma, zeta, P built
      independently from the Jacobi theta functions (mpmath.jtheta).
      SCOPE OF THIS TEST.  The identity is a CLASSICAL THEOREM
      (Whittaker-Watson 20.53), not something this script proves.
      What the test establishes is that the IMPLEMENTATION of sigma, zeta and P
      built here from the Jacobi theta functions satisfies it, on the lattice tried.
      A rectangular lattice has real invariants, so it is also a special
      case: it does not sample the general complex lattice.  The previous
      version of this comment said that the check "establishes it in
      general", which exceeds both what the test does and what would be needed:
      for generality the theorem is cited.

Theta construction (WW ch. 20-21):
  v = pi z/(2 w1),  q = exp(i pi w2/w1)
  sigma(z) = (2w1/pi) exp(eta1 z^2/(2w1)) theta1(v)/theta1'(0)
  zeta(z)  = eta1 z/w1 + (pi/(2w1)) theta1'(v)/theta1(v)
  eta1     = -(pi^2/(12 w1)) theta1'''(0)/theta1'(0)
  P(z)     = -zeta'(z)
"""

import mpmath as mp

mp.mp.dps = 30
w1 = mp.mpf(1)
w2 = mp.mpc(0, '0.7')
q = mp.exp(1j * mp.pi * w2 / w1)

th1 = lambda v: mp.jtheta(1, v, q)
th1p = lambda v: mp.jtheta(1, v, q, 1)
th1ppp = lambda v: mp.jtheta(1, v, q, 3)
eta1 = -(mp.pi**2 / (12 * w1)) * th1ppp(0) / th1p(0)

def sigma(z):
    v = mp.pi * z / (2 * w1)
    return (2 * w1 / mp.pi) * mp.exp(eta1 * z**2 / (2 * w1)) \
        * th1(v) / th1p(0)

def zeta(z):
    v = mp.pi * z / (2 * w1)
    return eta1 * z / w1 + (mp.pi / (2 * w1)) * th1p(v) / th1(v)

def wp(z):
    return -mp.diff(zeta, z)

def wpp(z):
    return mp.diff(wp, z)

print("=" * 68)
print("[A] consistency of the implementation: (P')^2 = 4P^3 - g2 P - g3")
print("=" * 68)
# g2, g3 from the linear system at two points, checked at others
z1, z2 = mp.mpf('0.31'), mp.mpc('0.22', '0.17')
M_ = mp.matrix([[-wp(z1), -1], [-wp(z2), -1]])
b_ = mp.matrix([wpp(z1)**2 - 4 * wp(z1)**3, wpp(z2)**2 - 4 * wp(z2)**3])
sol = mp.lu_solve(M_, b_)
g2n, g3n = sol[0], sol[1]
worst = 0
for zt in (mp.mpf('0.41'), mp.mpc('0.13', '0.29'), mp.mpc('0.55', '0.1')):
    res = wpp(zt)**2 - (4 * wp(zt)**3 - g2n * wp(zt) - g3n)
    worst = max(worst, abs(res))
print(f"  g2 = {mp.nstr(g2n, 8)}   g3 = {mp.nstr(g3n, 8)}")
print(f"  largest residual of the Weierstrass identity: {mp.nstr(worst, 3)}")

print()
print("=" * 68)
print("[B] antiderivative WW 20.53 (the core of the closed form)")
print("=" * 68)
v0 = mp.mpc('0.28', '0.21')          # generic 'pole'

def antider(z):
    return (1 / wpp(v0)) * (2 * zeta(v0) * z
                            + mp.log(sigma(z - v0) / sigma(z + v0)))

# (1) pointwise derivative vs integrand
worst_d = 0
for zt in (mp.mpf('0.15'), mp.mpf('0.33'), mp.mpc('0.4', '0.12')):
    d_num = mp.diff(antider, zt)
    d_teo = 1 / (wp(zt) - wp(v0))
    worst_d = max(worst_d, abs(d_num - d_teo))
print(f"  max |d/dz antiderivative - 1/(P - P(v))|: {mp.nstr(worst_d, 3)}")

# (2) definite integral vs difference of the antiderivative
za, zb = mp.mpf('0.12'), mp.mpf('0.44')
I_quad = mp.quad(lambda t: 1 / (wp(t) - wp(v0)), [za, zb])
I_anti = antider(zb) - antider(za)
print(f"  quadrature        = {mp.nstr(I_quad, 15)}")
print(f"  antiderivative    = {mp.nstr(I_anti, 15)}")
print(f"  diff              = {mp.nstr(abs(I_quad - I_anti), 3)}")
print("\n  => (i) reduction to 5.6e-17 + (ii) implementation of the antiderivative checked:")
print("     both links of the sigma/zeta closed form phi(r) are validated (the identity itself is WW 20.53).")
print("DONE.")
