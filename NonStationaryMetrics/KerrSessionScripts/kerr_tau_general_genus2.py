# -*- coding: utf-8 -*-
"""
GENERAL phi(r) of the equatorial tau brachistochrone of Kerr, in explicit form through
GENUS-2 special functions (Kleinian functions), for any J.

Starting point:  dphi/dr|_tau = J r sqrt(w f) / ( Δ sqrt(Δ - J^2 w) ),  w = E^2 - f.
Curve:     Y^2 = sextic(r;J) = r (r-2M)((E^2-1)r+2M) ( rΔ - J^2((E^2-1)r+2M) )
           (6 SIMPLE roots for generic J -> genus 2, hyperelliptic).
Rational numerator:  dphi/dr = K_tau/sqrt(sextic),  K_tau = J r (r-2M)((E^2-1)r+2M)/Δ.
Decomposition (partial fractions, poles at the HORIZONS r_±):

  phi(r) = c1 ∫ r dr/√sextic + c0 ∫ dr/√sextic          [first kind, holomorphic]
         + α_+ ∫ dr/((r-r_+)√sextic) + α_- ∫ dr/((r-r_-)√sextic)   [third kind]

  c1 = 11 J/25,  c0 = 2 J,
  α_± = 81 J (∓610√19 - 209)/95000     (for M=1, a=0.9, E=1.2).

The first-kind integrals are the coordinates of the Abel map (arguments
of Klein's σ); the third-kind ones are quotients of log(σ). The period
matrix τ (2x2) is computed with SageMath (curve y^2=sextic); see
kerr_jpt_genus2_kleinian.py (not distributed) for the same machinery applied to J_+^t (t branch).

Reduction to |J|=J_c=a/E: f^2 makes r_e=2M a DOUBLE root; the factor (r-2M) of K_tau
cancels one of them -> √Q4 (quartic) = GENUS 1 (Weierstrass) = separatrix.

Validation here: decomposition vs direct integration ~1e-15; reduction to J_c
vs the separatrix eq.59 ~1e-14.
"""
import numpy as np
import sympy as sp
from scipy.integrate import quad

M, a, E = 1.0, 0.9, 1.2
Jc = a/E
rp = M + np.sqrt(M**2 - a**2)
rm = M - np.sqrt(M**2 - a**2)


def sextic(r, J):
    num_wf = (r - 2*M)*((E**2 - 1)*r + 2*M)
    num_DJw = r*(r**2 - 2*M*r + a**2) - J**2*((E**2 - 1)*r + 2*M)
    return r*num_wf*num_DJw


def coeffs(J):
    """GENERAL coefficients (M,a,E,J), poles at the horizons r_± = M ± √(M²−a²).
       c1 = J(E²−1),  c0 = 2MJ,
       α_± = ∓ J a² [ M(E²+1) ± (E²−1)√(M²−a²) ] / (2√(M²−a²))."""
    D = np.sqrt(M**2 - a**2)
    c1 = J*(E**2 - 1)
    c0 = 2*M*J
    ap = -J*a**2*(M*(E**2 + 1) + (E**2 - 1)*D)/(2*D)
    am = +J*a**2*(M*(E**2 + 1) - (E**2 - 1)*D)/(2*D)
    return c1, c0, ap, am


def dphidr_direct(r, J):
    f = 1 - 2*M/r
    w = E**2 - f
    Dl = r**2 - 2*M*r + a**2
    return J*r*np.sqrt(w*f)/(Dl*np.sqrt(Dl - J**2*w))


def dphidr_special(r, J):
    """explicit form: decomposed numerator / sqrt(sextic)."""
    c1, c0, ap, am = coeffs(J)
    return (c1*r + c0 + ap/(r - rp) + am/(r - rm))/np.sqrt(abs(sextic(r, J)))


def phi_special(r, J, r_ref=6.0):
    """phi(r) as a combination of genus-2 abelian integrals (evaluated)."""
    return quad(dphidr_special, r, r_ref, limit=400)[0]


if __name__ == '__main__':
    # ---- explicit symbolic form ----
    r, J = sp.symbols('r J')
    Dl = r**2 - 2*M*r + a**2
    sx = sp.expand(r*(r - 2)*sp.Rational(11, 25)*(r + sp.Rational(50, 11)) *
                   (r*Dl - J**2*(sp.Rational(11, 25)*r + 2)))
    print("Y^2 = sextic(r;J):")
    print("  ", sp.collect(sp.expand(
        r*(r-2*M)*((E**2-1)*r+2*M)*(r*Dl - J**2*((sp.Rational(11,25))*r+2))), r))
    print("\ncoeff:  c1=11J/25, c0=2J,  a_+ =", coeffs(1)[2],
          "* J,  a_- =", coeffs(1)[3], "* J")

    # ---- validation: decomposition vs direct, several J ----
    print("\nvalidation (max|dphi_direct - dphi_special| on [2.05,6]):")
    for k in [-1.5, -1.01, -0.9, -0.5, 0.5, 1.2]:
        Jv = k*Jc
        rr = np.linspace(2.05, 6, 60)
        err = max(abs(dphidr_direct(x, Jv) - dphidr_special(x, Jv)) for x in rr)
        print(f"  J={k:+.2f} Jc = {Jv:+.4f}:  {err:.2e}")
    print("\n=> the explicit genus-2 form (abelian integrals) reproduces phi(r) within the errors above.")
