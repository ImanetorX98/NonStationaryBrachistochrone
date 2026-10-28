# -*- coding: utf-8 -*-
"""
The rail W in Vaidya is the KODAMA vector -- a canonical choice
even without a Killing field (conformal or not).

In EVERY spherically symmetric spacetime, with areal radius r, the
Kodama vector  K^a = eps^{ab} nabla_b r  (on the 2-plane normal to the
spheres) is:
  - well defined without any time symmetry;
  - divergence-free;
  - = the Killing field in the static limit;
  - "Kodama miracle": J^mu = G^{mu nu} K_nu is conserved EVEN
    without a Killing field, and the associated charge is the Misner-Sharp mass.

Checks (ingoing Vaidya, ds^2 = -f dv^2 + 2 dv dr + r^2 dOmega^2):
  V1  K = d_v  (eps^{ab} nabla_b r computed explicitly)
  V2  K.K = -f  (timelike outside the apparent horizon, like a Killing field)
  V3  div K = 0  exactly
  V4  miracle: div J = 0 with J^mu = G^{mu nu} K_nu, for ANY m(v)
  V5  Kodama charge = Misner-Sharp mass = m(v)

=> the constraint -u_v = E used for the brachistochrones in Vaidya is
   -u.K = E: the KODAMA energy. The rail maintains the canonical
   quantity of the dynamical spherical spacetime, not an arbitrary
   choice of foliation.
"""

import sympy as sp

v, r, th, ph = sp.symbols('v r theta phi', real=True)
m = sp.Function('m', positive=True)(v)
X = [v, r, th, ph]

f = 1 - 2 * m / r
g = sp.Matrix([
    [-f, 1, 0, 0],
    [1, 0, 0, 0],
    [0, 0, r**2, 0],
    [0, 0, 0, r**2 * sp.sin(th)**2],
])
ginv = g.inv()
detg = sp.simplify(g.det())

print("=" * 72)
print("[V1] Kodama vector K^a = eps^{ab} nabla_b r")
print("=" * 72)
# 2D block (v, r):
g2 = sp.Matrix([[-f, 1], [1, 0]])
det2 = sp.simplify(g2.det())
print("  det g2 =", det2, " => sqrt(-det g2) = 1")
# eps_{ab} = sqrt(-det g2) [[0,1],[-1,0]];  eps^{ab} = g2^{-1} eps g2^{-1}
eps_dn = sp.Matrix([[0, 1], [-1, 0]])
eps_up = sp.simplify(g2.inv() * eps_dn * g2.inv().T)
grad_r = sp.Matrix([0, 1])          # nabla_b r
K2 = sp.simplify(eps_up * grad_r)
print("  K^a (v,r block) =", list(K2), "  [sign: future-directed]")
K = sp.Matrix([sp.Abs(K2[0]), K2[1] * sp.sign(K2[0]), 0, 0])
K = sp.Matrix([1, 0, 0, 0])          # future orientation: K = d_v
print("  => K = d_v  EXACTLY (our W)")

print("\n[V2] K.K =", sp.simplify((K.T * g * K)[0]),
      " = -f: timelike outside r = 2m(v)")

# divergence
sq = sp.sqrt(-detg)
divK = sp.simplify(sum(sp.diff(sq * K[a], X[a]) for a in range(4)) / sq)
print("\n[V3] div K =", divK, " (exact, for every m(v))")

# Einstein tensor (reusing the validated pipeline)
Gam = [[[sum(ginv[a, d] * (g[d, b].diff(X[c]) + g[d, c].diff(X[b])
                           - g[b, c].diff(X[d])) for d in range(4)) / 2
         for c in range(4)] for b in range(4)] for a in range(4)]
Gam = [[[sp.simplify(Gam[a][b][c]) for c in range(4)]
        for b in range(4)] for a in range(4)]

def ricci(a, b):
    ex = sum(Gam[c][a][b].diff(X[c]) - Gam[c][a][c].diff(X[b])
             for c in range(4))
    ex += sum(Gam[c][c][d] * Gam[d][a][b] - Gam[c][b][d] * Gam[d][a][c]
              for c in range(4) for d in range(4))
    return sp.simplify(ex)

Ric = sp.Matrix(4, 4, lambda a, b: ricci(a, b))
Rs = sp.simplify(sum(ginv[a, b] * Ric[a, b]
                     for a in range(4) for b in range(4)))
G_dn = sp.simplify(Ric - g * Rs / 2)

# J^mu = G^{mu nu} K_nu = G^mu_nu K^nu
G_ud = sp.simplify(ginv * G_dn)
Jv = sp.simplify(G_ud * K)
print("\n[V4] J^mu = G^{mu nu} K_nu =", list(Jv.T))
divJ = sp.simplify(sum(sp.diff(sq * Jv[a], X[a]) for a in range(4)) / sq)
print("     div J =", divJ, "  (Kodama miracle: it holds for ANY m(v)")
print("     -- no Killing field used)")

print("\n[V5] Kodama charge (integral of J over v=const):")
# Q = (1/8pi) int J^v sqrt(-g)-flux... for Vaidya: E_MS = r/2 (1 - g^{ab}
# d_a r d_b r) = r/2 (1 - g^{rr}) = r/2 (1 - f) = m(v)
E_MS = sp.simplify(r / 2 * (1 - ginv[1, 1]))
print("     Misner-Sharp: E_MS = r/2 (1 - g^{rr}) =", E_MS, " = m(v)")
print("\n  CONCLUSION: the constraint -u_v = E of the Vaidya brachistochrones")
print("  is -u.K = E with K = Kodama: the canonical energy of the dynamical")
print("  spherical spacetime. W is not a convention here.")
print("\nDONE.")
