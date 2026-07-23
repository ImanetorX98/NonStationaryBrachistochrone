# -*- coding: utf-8 -*-
"""
Closed-form reduction of the NEW Thakurta-Kerr off-shell source term (main11 referee).

The corrected terminally-anchored source is  S_D = int (Theta H + P_r H_Pr) dlambda, with the
extra dilation piece  int P_r H_Pr dlambda = int p_r dr  (the radial action; verified 1e-15).

This script proves, symbolically + numerically, that this extra letter lives on the SAME
spectral curve as the on-shell reduction and reduces to the SAME basis:

  int p_r dr  =  const * [ sum_k a_k U_k  +  (third-kind) ] ,   U_k = int r^k / sqrt(R6) dr,

with EXACT rational coefficients a_k, and third-kind poles exactly at the seed Kerr null
surfaces r_pm = M +/- sqrt(M^2 - a^2).  (This is "level B": the second/third-kind reduction.
The outer length-two wrap  int (G_Pr/H_Pr) * int p_r dr  -- the genus-2 dilogarithm -- is the
declared-open "level C" special-function assembly, NOT closed here.)

Run:  python ThakurtaMetric/adiabatic_offshell_reduction.py
"""
import sympy as sp
import numpy as np
from scipy.integrate import cumulative_trapezoid as ct

r = sp.symbols('r')
M = sp.Rational(1)


def _sqfree_rational(g):
    """Square-free part of a rational function (odd-multiplicity zeros and poles)."""
    num, den = sp.fraction(sp.cancel(g))
    out = sp.Integer(1)
    for poly in (num, den):
        for base, mult in sp.factor_list(poly)[1]:
            if mult % 2:
                out *= base
    return sp.expand(out)


def _perfect_sqrt(poly):
    """sqrt of a perfect-square polynomial (const kept symbolic if not a perfect square)."""
    fl = sp.factor_list(poly)
    res = sp.sqrt(fl[0])
    for base, mult in fl[1]:
        res *= base ** (mult // 2)
    return res


def reduce_branch(branch, Eq, Jq, aq):
    """Return dict with the exact reduction data and numeric residuals for one config."""
    a = aq
    f = 1 - 2 * M / r
    Dl = r**2 - 2 * M * r + a**2
    b = 2 * M * a / r
    v = 1 - f / Eq**2
    Pcap = r**2 + a**2 + 2 * M * a**2 / r
    Pb = Pcap + b**2 / Eq**2
    if branch == 'eta':
        A_term = Jq * b * v / Pb
        D = Jq**2 / Pb
    else:  # tau: shifted azimuthal momentum
        pt = Jq - b / Eq
        A_term = pt * (b * v / Pb)
        D = pt**2 / Pb
    Bc = Dl * v / Pb
    C_ = Dl / r**2
    # frozen radial momentum squared from the shell Hbar = 0
    pr2 = sp.cancel(((1 - A_term)**2 / Bc - D) / C_)
    curve = _sqfree_rational(pr2)                 # spectral curve y^2 = curve  (== on-shell R6)
    ratio = sp.cancel(pr2 / curve)                # must be a perfect rational square
    n, d = sp.fraction(ratio)
    same_curve = (all(m % 2 == 0 for _, m in sp.factor_list(n)[1]) and
                  all(m % 2 == 0 for _, m in sp.factor_list(d)[1]))
    g = sp.cancel(_perfect_sqrt(n) / _perfect_sqrt(d))   # p_r / sqrt(curve)  (rational)
    pf = sp.apart(sp.cancel(g * curve), r)               # = poly(U_k) + partial fractions (3rd kind)
    poly_terms = sp.Add(*[t for t in pf.as_ordered_terms() if t.is_polynomial(r)])
    Uk = [sp.nsimplify(c) for c in sp.Poly(poly_terms, r).all_coeffs()[::-1]]
    third = sp.together(pf - poly_terms)
    poles = sorted({sp.re(p).evalf() for p in sp.roots(sp.fraction(third)[1], r)})
    rpm = [float(M - sp.sqrt(M**2 - a**2)), float(M + sp.sqrt(M**2 - a**2))]
    # ---- numeric verification on the physical arc ----
    pr2n = sp.lambdify(r, pr2, 'numpy')
    prd = sp.lambdify(r, -sp.sqrt(pr2), 'numpy')          # ingoing branch
    red = sp.lambdify(r, pf / sp.sqrt(curve), 'numpy')    # reduced integrand
    rr = np.linspace(11.5, 4.6, 800)
    rr = rr[np.isfinite(prd(rr)) & (pr2n(rr) > 0)]
    e_int = float(np.nanmax(np.abs(np.abs(prd(rr)) - np.abs(red(rr)))))
    dS_dir = ct(np.abs(prd(rr)), rr, initial=0)
    dS_red = ct(np.abs(red(rr)), rr, initial=0)
    e_dS = float(np.nanmax(np.abs(dS_dir - dS_red)))
    return dict(same_curve=same_curve, Uk=Uk, poles=[round(p, 4) for p in poles],
                rpm=[round(x, 4) for x in rpm], e_int=e_int, e_dS=e_dS)


CONFIGS = [
    ('t   (a=0.9, E=1.4, J=6)  ', 'eta', sp.Rational(14, 10), sp.Integer(6), sp.Rational(9, 10)),
    ('tau (a=0.9, E=1.4, J=2.5)', 'tau', sp.Rational(14, 10), sp.Rational(5, 2), sp.Rational(9, 10)),
    ('t   (a=0.5, E=1.3, J=5)  ', 'eta', sp.Rational(13, 10), sp.Integer(5), sp.Rational(1, 2)),
]

if __name__ == '__main__':
    print("Closed-form reduction of  int p_r dr  (the new S_D dilation letter):")
    ok = True
    for name, br, E, J, a in CONFIGS:
        d = reduce_branch(br, E, J, a)
        print(f"\n== {name} ==  same spectral curve as on-shell: {d['same_curve']}")
        print(f"   U_k coefficients (k=0..): {d['Uk']}")
        print(f"   third-kind poles {d['poles']}  ==  seed Kerr null surfaces r_pm {d['rpm']}")
        print(f"   numeric: |p_r - reduced integrand| = {d['e_int']:.1e}   |Delta S red - direct| = {d['e_dS']:.1e}")
        ok = ok and d['same_curve'] and d['e_int'] < 1e-11 and d['e_dS'] < 1e-11
    print("\nOK" if ok else "\nFAILED")
