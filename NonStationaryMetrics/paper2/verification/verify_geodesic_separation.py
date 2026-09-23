#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
verify_geodesic_separation.py -- the geodesic separation on the fixed-energy
optical surface, i.e. the Jacobi (geodesic-deviation) equation the Round-2
referee asks to have calculated.

The referee's cited source (Lecian, `10.13140/RG.2.2.24140.58246`, section 7)
states the spacetime form, D^2 xi^mu / dtau^2 = R^mu_{nu alpha beta} u^nu u^alpha
xi^beta.  On a two-dimensional surface that reduces to the scalar Jacobi equation

    xi'' (s) + K(s) xi(s) = 0,        s = arclength of the optical metric,

with K the Gauss curvature.  Everything below is rebuilt from a = Lambda^2 h:
eq. (eq:optical-gauss) of the manuscript is checked, not assumed.

Proved here, for the frozen non-rotating branch with Ehat^2 >= 3/2:
  (i)  K < 0 everywhere outside r = 2M;
  (ii) hence a Jacobi field with xi(0)=0, xi'(0)=1 obeys xi'' = |K| xi > 0, so it
       is convex and strictly increasing: no conjugate point, for any arc length;
  (iii) the separation grows at least linearly, xi(s) >= s, and at least like
       sinh(k s)/k on any subarc where |K| >= k^2.
Below the threshold the far-field curvature changes sign and (i) fails, which is
why the threshold is where it is.
"""
import sympy as sp

r, ph, s, M, Eh, k = sp.symbols('r varphi s M Ehat k', positive=True)
f = 1 - 2*M/r
Lam2 = Eh**2/(f*(Eh**2 - f))

ok = True


def check(label, cond):
    global ok
    ok = ok and bool(cond)
    print(f'{"PASS" if cond else "FAIL"}  {label}')


# -- Gauss curvature of a = Lambda^2 (dr^2/f + r^2 dphi^2) -------------------
E = Lam2/f                      # a_rr
G = Lam2*r**2                   # a_phiphi   (no phi dependence)
W = sp.sqrt(E*G)
K = sp.simplify(-1/(2*W)*sp.diff(sp.diff(G, r)/W, r))
K = sp.simplify(sp.factor(K))

printed = (M*(16*M**3 + 2*(9*Eh**2 - 14)*M**2*r + (3*Eh**4 - 19*Eh**2 + 16)*M*r**2
              - (2*Eh**2 - 3)*(Eh**2 - 1)*r**3)
           / (Eh**2*r**5*(2*M + (Eh**2 - 1)*r)))
check('K reproduces eq. (eq:optical-gauss) of the manuscript',
      sp.simplify(K - printed) == 0)

check('Ehat -> infinity gives the Fermat value -(2M/r^3)(1 - 3M/2r)',
      sp.simplify(sp.limit(K, Eh, sp.oo) + 2*M/r**3*(1 - 3*M/(2*r))) == 0)
check('r^3 K -> -M(2 Ehat^2 - 3)/Ehat^2 as r -> infinity',
      sp.simplify(sp.limit(r**3*K, r, sp.oo) + M*(2*Eh**2 - 3)/Eh**2) == 0)

# -- (i) sign of K outside r = 2M at and above the threshold -----------------
# write r = 2M(1+x), x>0, and Ehat^2 = 3/2 + e, e>=0; show the numerator < 0.
x, e = sp.symbols('x epsilon', positive=True)
num = sp.numer(sp.together(K))
num_sub = sp.expand(sp.simplify(num.subs({r: 2*M*(1 + x), Eh**2: sp.Rational(3, 2) + e})
                                .subs(Eh, sp.sqrt(sp.Rational(3, 2) + e))))
poly = sp.Poly(sp.expand(-num_sub/M**4), x, e)
coeffs = poly.coeffs()
check('-numerator(K) is a polynomial in (x, epsilon) with all coefficients > 0',
      all(c > 0 for c in coeffs))
print('   coefficients:', sorted(set(coeffs)))

check('at threshold and horizon corner, K(2M, 3/2) < 0',
      sp.simplify(K.subs({r: 2*M, Eh: sp.sqrt(sp.Rational(3, 2))})) < 0)

# -- (ii)-(iii) Jacobi field ------------------------------------------------
# xi'' = -K xi with K < 0: convexity gives the two bounds, verified on the
# comparison solutions rather than asserted.
xi = sp.Function('xi')
comp = sp.sinh(k*s)/k
check('sinh(ks)/k solves the comparison equation xi\'\' - k^2 xi = 0 with xi(0)=0, xi\'(0)=1',
      sp.simplify(sp.diff(comp, s, 2) - k**2*comp) == 0
      and sp.limit(comp, s, 0) == 0 and sp.simplify(sp.diff(comp, s).subs(s, 0)) == 1)
check('the linear bound xi = s is the k -> 0 limit of the comparison solution',
      sp.simplify(sp.limit(comp, k, 0) - s) == 0)

# a concrete curvature scale for the manuscript's own configuration
for Ev in [sp.sqrt(sp.Rational(3, 2)), sp.Rational(7, 5), 2]:
    Kv = K.subs({M: 1, Eh: Ev})
    vals = [(rv, float(Kv.subs(r, rv))) for rv in (2.5, 3, 6, 10, 30)]
    print(f'   Ehat = {float(Ev):.4f}:  ' +
          '  '.join(f'K({rv})={v:+.3e}' for rv, v in vals))

print('\nALL CHECKS PASS' if ok else '\nSOME CHECKS FAILED')
raise SystemExit(0 if ok else 1)
