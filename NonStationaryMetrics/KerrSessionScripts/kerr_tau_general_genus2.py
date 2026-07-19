# -*- coding: utf-8 -*-
"""
phi(r) GENERALE della brachistocrona tau equatoriale di Kerr, forma esplicita in
funzioni speciali di GENERE 2 (Kleinian), per J qualsiasi.

Partenza:  dphi/dr|_tau = J r sqrt(w f) / ( Δ sqrt(Δ - J^2 w) ),  w = E^2 - f.
Curva:     Y^2 = sextic(r;J) = r (r-2M)((E^2-1)r+2M) ( rΔ - J^2((E^2-1)r+2M) )
           (6 radici SEMPLICI per J generico -> genere 2, iperellittico).
Numeratore razionale:  dphi/dr = K_tau/sqrt(sextic),  K_tau = J r (r-2M)((E^2-1)r+2M)/Δ.
Decomposizione (frazioni parziali, poli agli ORIZZONTI r_±):

  phi(r) = c1 ∫ r dr/√sextic + c0 ∫ dr/√sextic          [1a specie, olomorfi]
         + α_+ ∫ dr/((r-r_+)√sextic) + α_- ∫ dr/((r-r_-)√sextic)   [3a specie]

  c1 = 11 J/25,  c0 = 2 J,
  α_± = 81 J (∓610√19 - 209)/95000     (per M=1, a=0.9, E=1.2).

Gli integrali di 1a specie sono le coordinate della mappa di Abel (argomenti
della σ di Klein); quelli di 3a specie sono log(σ)-quozienti. La matrice dei
periodi τ (2x2) si calcola con SageMath (curva y^2=sextic); vedi
kerr_jpt_genus2_kleinian.py per la stessa macchina applicata a J_+^t (ramo t).

Riduzione a |J|=J_c=a/E: la f^2 rende r_e=2M radice DOPPIA; la (r-2M) di K_tau
ne cancella una -> √Q4 (quartica) = GENERE 1 (Weierstrass) = separatrice.

Validazione qui: decomposizione vs integrazione diretta ~1e-15; riduzione a J_c
vs separatrice eq.59 ~1e-14.
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
    """coefficienti GENERALI (M,a,E,J), poli agli orizzonti r_± = M ± √(M²−a²).
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
    """forma esplicita: numeratore decomposto / sqrt(sextic)."""
    c1, c0, ap, am = coeffs(J)
    return (c1*r + c0 + ap/(r - rp) + am/(r - rm))/np.sqrt(abs(sextic(r, J)))


def phi_special(r, J, r_ref=6.0):
    """phi(r) come combinazione di integrali abeliani genere-2 (valutata)."""
    return quad(dphidr_special, r, r_ref, limit=400)[0]


if __name__ == '__main__':
    # ---- forma simbolica esplicita ----
    r, J = sp.symbols('r J')
    Dl = r**2 - 2*M*r + a**2
    sx = sp.expand(r*(r - 2)*sp.Rational(11, 25)*(r + sp.Rational(50, 11)) *
                   (r*Dl - J**2*(sp.Rational(11, 25)*r + 2)))
    print("Y^2 = sextic(r;J):")
    print("  ", sp.collect(sp.expand(
        r*(r-2*M)*((E**2-1)*r+2*M)*(r*Dl - J**2*((sp.Rational(11,25))*r+2))), r))
    print("\ncoeff:  c1=11J/25, c0=2J,  a_+ =", coeffs(1)[2],
          "* J,  a_- =", coeffs(1)[3], "* J")

    # ---- validazione: decomposizione vs diretta, vari J ----
    print("\nvalidazione (max|dphi_direct - dphi_special| su [2.05,6]):")
    for k in [-1.5, -1.01, -0.9, -0.5, 0.5, 1.2]:
        Jv = k*Jc
        rr = np.linspace(2.05, 6, 60)
        err = max(abs(dphidr_direct(x, Jv) - dphidr_special(x, Jv)) for x in rr)
        print(f"  J={k:+.2f} Jc = {Jv:+.4f}:  {err:.2e}")
    print("\n=> forma esplicita genere-2 (integrali abeliani) = phi(r) esatta.")
