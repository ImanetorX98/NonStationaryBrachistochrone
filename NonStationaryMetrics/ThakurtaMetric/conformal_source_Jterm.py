# -*- coding: utf-8 -*-
"""
Missing term d_J F in the conformal adiabatic source of Thakurta-Kerr.

The conformal weights E_eff = Ehat/A and J_eff = J/A are those of the canonical derivation
of the manuscript (Euler source (E d_E + J d_J)F); they are not re-derived here.
Since the frozen shape is the Kerr one, F = F(r; E_eff, J_eff), the first-order
variation (A = A_0(1+eps eta), eps = A'/A) contains TWO terms:

    delta phi = -int eta [ Ehat d_E F + J d_J F ] dr        (E_0=Ehat, J_0=J at A_0=1)

This script checks that the term J d_J F reduces in the same basis as d_E F and
tabulates the coefficients c_k and b_i with and without it.
An earlier version used only -Ehat int eta d_E F dr, without the term J d_J F.
Structure: (E d_E + J d_J) F is the scaling operator; its numerator reduces
in the SAME basis of building blocks as d_E F (no new building block), so the
correction enters as coefficients c_k^J (generic) or b_i^J (separatrix), additive
to those in E.

NB: Vaidya (parameter m, not conformal) is NOT affected: E and J are both conserved,
the source is d_m F, and there is no d_J F term.

Check at the level of the complete ASSEMBLY (not only the coefficients): substituting the
Euler source (E d_E + J d_J)F in tk_t_sep_blockassembly.py, the block assembly
with ALL the special functions (Weierstrass zeta, wp through theta + genus-2 dilogarithm,
143 products) reproduces the direct integral to ~1e-8 on the TK-t separatrix --
the same precision as the E term alone. The closed form is therefore complete also with the
J term: no new special function, only the coefficients c_k^tot=c_k^E+c_k^J.

Output: c_k^J and the corrected b_i for the TK t and tau branches, with a closure check.
Branch numerators: t  K_t = r D_E (J(r-2M)+2Ma)/Delta;  tau  K_tau = Jc r (r-2M) D_E/Delta
(the spherical pattern J D_E holds only for a = 0).
The numerical parameters are substituted right after the derivatives in E and J, and the critical J is
obtained from the real roots of the factors of the resultant: the same mathematics, a few seconds
instead of hours of symbolic simplification.
"""
import sympy as sp
import numpy as np

r, M, a, E, J, Jc, rd = sp.symbols('r M a E J Jc r_d', real=True)


# ---------------------------------------------------------------- generic
def reduce_2nd(F, R6, sub, label):
    """Reduce d_E F, d_J F and (E d_E + J d_J)F in the second-kind basis r^k/sqrt(R6)."""
    out = {}
    for tag, src in (('E', E * sp.diff(F, E)), ('J', J * sp.diff(F, J)),
                     ('tot', E * sp.diff(F, E) + J * sp.diff(F, J))):
        N = sp.expand(sp.cancel((src * R6**sp.Rational(3, 2)).subs(sub)))
        R6s = sp.expand(R6.subs(sub))
        Ac = [sp.Symbol(f'A{i}') for i in range(6)]
        ck = [sp.Symbol(f'c{i}') for i in range(5)]
        A5 = sum(Ac[i] * r**i for i in range(6))
        Mp = sum(ck[i] * r**i for i in range(5))
        eq = sp.expand(2 * N - (2 * R6s * sp.diff(A5, r) - A5 * sp.diff(R6s, r) + 2 * R6s * Mp))
        sol = sp.solve(sp.Poly(eq, r).all_coeffs(), Ac + ck, dict=True)
        assert sol, f"{label}/{tag}: the second-kind reduction does NOT close"
        out[tag] = [float(sol[0][ck[i]]) for i in range(5)]
    return out


print("=" * 70)
print("GENERIC t branch (any J): c_k^E, c_k^J, c_k^tot  [they close in the same basis]")
print("=" * 70)
DE = (E**2 - 1) * r + 2 * M
Dl = r**2 - 2 * M * r + a**2
Q2 = (2*E**2*J**2*M*r - E**2*J**2*r**2 - 4*E**2*J*M*a*r + 2*E**2*M*a**2*r + E**2*a**2*r**2
      + E**2*r**4 + 4*J**2*M**2 - 4*J**2*M*r + J**2*r**2 - 8*J*M**2*a + 4*J*M*a*r + 4*M**2*a**2)
R6g = sp.expand(r * Q2 * DE)
Ktg = r * DE * (J * (r - 2 * M) + 2 * M * a) / Dl
subg = {M: 1, a: sp.Rational(9, 10), E: sp.Rational(6, 5), J: sp.Rational(5, 2)}
red = reduce_2nd(Ktg / sp.sqrt(R6g), R6g, subg, "t-gen")
print("  c_k^E   =", [round(x, 5) for x in red['E']], "  (reduction of E d_E F)")
print("  c_k^J   =", [round(x, 5) for x in red['J']], "  <-- NEW (reduction of J d_J F)")
print("  c_k^tot =", [round(x, 5) for x in red['tot']], "  (reduction of (E d_E + J d_J)F)")
chk = [red['E'][i] + red['J'][i] for i in range(5)]
print("  c^E + c^J =", [round(float(x), 5) for x in chk], "  (must equal c_k^tot, by linearity)")


# ---------------------------------------------------------------- separatrix
def bcoeffs(S, N):
    d = lambda ex, n: sp.diff(ex, r, n).subs(r, rd)
    S2, S3, S4 = d(S, 2), d(S, 3), d(S, 4)
    Q4rd, Q4prd, Q4pprd = S2 / 2, S3 / 6, S4 / 12
    s = sp.sqrt(Q4rd); a1 = Q4prd / (4 * s); a2 = Q4pprd / 12
    Nr, Np, Npp = N.subs(r, rd), sp.diff(N, r).subs(r, rd), sp.diff(N, r, 2).subs(r, rd)
    h0 = Nr / Q4rd
    Fp = (Np * Q4rd - Nr * Q4prd) / Q4rd**2
    Fpp = ((Npp * Q4rd - Nr * Q4pprd) * Q4rd - 2 * Q4prd * (Np * Q4rd - Nr * Q4prd)) / Q4rd**3
    h1 = Fp * s; h2 = sp.Rational(1, 2) * (Fpp * s**2 + Fp * (Q4prd / 2))
    return {'b3': h0 / s**3, 'b2': (h1 - 3 * a1 * h0) / s**3,
            'b1': (h2 - 3 * a1 * h1 + (6 * a1**2 - 3 * a2) * h0) / s**3}


def sep_bi(S, K, params, name, Jmin=5):
    Jv = sp.Symbol('Jvar')
    Sj, Kj = S.subs(Jc, Jv), K.subs(Jc, Jv)
    F = Kj / sp.sqrt(Sj)
    Sub0 = {M: params['M'], a: params['a'], E: params['E']}
    N_E = sp.expand(sp.cancel((sp.diff(F, E) * Sj**sp.Rational(3, 2)).subs(Sub0))).subs(Jv, Jc)
    N_tot = sp.expand(sp.cancel(((E * sp.diff(F, E) + Jv * sp.diff(F, Jv)) * Sj**sp.Rational(3, 2)).subs(Sub0))).subs(Jv, Jc)
    S = sp.expand(S.subs(Sub0))
    Sn = S.subs(Sub0)
    res = sp.Poly(sp.resultant(Sn, sp.diff(Sn, r), r), Jc)
    reals = []
    for fac, _ in res.factor_list()[1]:
        reals += [float(z) for z in fac.real_roots()]
    Jcv = min(z for z in reals if z > Jmin)
    rts = np.roots([complex(c) for c in sp.Poly(Sn.subs(Jc, Jcv), r).all_coeffs()])
    prs = [(i, j) for i in range(len(rts)) for j in range(i + 1, len(rts)) if abs(rts[i] - rts[j]) < 1e-5]
    rdv = float(np.real((rts[prs[0][0]] + rts[prs[0][1]]) / 2))
    sub = {M: params['M'], a: params['a'], E: params['E'], Jc: Jcv, rd: rdv}
    bE, bT = bcoeffs(S, N_E), bcoeffs(S, N_tot)
    print(f"\n  {name}:  Jc={Jcv:.5f}  r_d={rdv:.5f}   (NB r_d<0: cfr. separatrix_classification.py)")
    E0 = float(params['E'])
    for k in ['b1', 'b2', 'b3']:
        v0, v1 = float(bE[k].subs(sub)), float(bT[k].subs(sub))
        en = E0 * v0                    # energy-only source Ehat d_E F
        print(f"    {k}:  d_E F={v0:+.6f}   Ehat d_E F={en:+.6f}   total={v1:+.6f}   "
              f"axial J d_J F={v1 - en:+.6f}  ({100 * (v1 - en) / en:+.2f}% of Ehat d_E F)")


print()
print("=" * 70)
print("SEPARATRIX: paper b_i (d_E only) vs corrected (E d_E + J d_J)")
print("=" * 70)
S_tau = sp.expand(r * (r - 2 * M) * DE * (r * Dl - Jc**2 * DE))
sep_bi(S_tau, Jc * r * (r - 2 * M) * DE / Dl, {'M': 1, 'a': sp.Rational(9, 10), 'E': sp.Rational(6, 5)}, "TK tau")
Q2c = Q2.subs(J, Jc)
sep_bi(sp.expand(r * Q2c * DE), r * DE * (Jc * (r - 2 * M) + 2 * M * a) / Dl,
       {'M': 1, 'a': sp.Rational(9, 10), 'E': sp.Rational(6, 5)}, "TK t")
print()
print("Vaidya: NOT affected (parameter m, not conformal; E, J conserved; source d_m F).")
print("DONE.")
