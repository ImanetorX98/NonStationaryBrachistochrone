# -*- coding: utf-8 -*-
"""
Classification of the degenerations of the sextic: algebraic locus vs physical separatrix.

Critical correction: an earlier version called "separatrix" a locus of a
double root of the curve which, at the stated parameters, lies at NEGATIVE radius
(r_d<0): it is an algebraic degeneration (genus 2->1), NOT the physical separatrix
(unstable circular orbit, capture/scattering boundary), which requires a double
root in the CONTROL DOMAIN r_d>2M where two real turning points merge.
The criterion is r_d>2M and not r_d>r_+: between r_+ and 2M the selector W=d_eta is
spacelike and the compact control problem is not posed (Paper II, compactness
lemma); a double root there is an object of type (A), not a separatrix.

Distinction to keep in the paper:
  J_deg      : discriminant locus Res_r(curve, d_r curve)=0 (genus degeneration).
               r_d can be <0, complex or interior: the closed form (Weierstrass)
               stays valid, but it is NOT a physical separatrix.
  J_cap(v0)  : dynamical penetration threshold, a distinct PHYSICAL object.
  J_c^phys   : double root in the control domain (r_d>2M) = physical separatrix
               (merger of two real turning points).

Output: for each branch, ALL the real double roots, classified; for the t branch
of Thakurta-Kerr the outer physical separatrix and its b_i (corrected source
E d_E + J d_J, cf. conformal_source_Jterm.py).
"""
import sympy as sp
import numpy as np
import mpmath as mp

r, J, Es = sp.symbols('r J Es', real=True)
M = 1
mp.mp.dps = 25


def classify(curve, aval, label):
    """Res_r(curve, d_r curve)=0 -> all the J with a double root; classify r_d."""
    rp = 1 + np.sqrt(1 - aval**2) if aval < 1 else 2.0     # null surface of the seed r_+
    re = 2.0 * M                                           # stationary limit: boundary of the domain
    Res = sp.resultant(curve, sp.diff(curve, r), r)
    print(f"\n===== {label}  (r_+={rp:.4f}) =====")
    rows = []
    for js in sp.solve(sp.Eq(Res, 0), J):
        jc = complex(js)
        if abs(jc.imag) > 1e-8 or abs(jc.real) < 1e-6:
            continue
        jv = jc.real
        rts = np.roots([complex(c) for c in sp.Poly(curve.subs(J, jv), r).all_coeffs()])
        rd = None
        for i in range(len(rts)):
            for k in range(i + 1, len(rts)):
                if abs(rts[i] - rts[k]) < 1e-4 and abs(rts[i].imag) < 1e-6:
                    rd = (rts[i].real + rts[k].real) / 2
        if rd is None:
            continue
        # physical separatrix = merger of two real outer turning points
        real_ext = sorted([z.real for z in rts if abs(z.imag) < 1e-6 and z.real > re])
        phys = rd > re and any(abs(x - rd) < 1e-3 for x in real_ext)
        if rd > re:
            region = "control domain r>2M"
        elif rd > rp:
            region = "r_+<r<2M, W spacelike: type (A)"
        elif rd > 1e-9:
            region = "interior r<r_+"
        elif rd > -1e-9:
            region = "r_d=0, trivial"
        else:
            region = "NEGATIVE"
        rows.append((jv, rd, region, phys))
    for jv, rd, region, phys in sorted(rows):
        tag = "  <== PHYSICAL SEPARATRIX" if phys else "  (algebraic degeneration J_deg)"
        print(f"  J = {jv:+.5f}   r_d = {rd:+.5f}   [{region}]{tag}")
    return rows


DE = (Es**2 - 1) * r + 2 * M
Dl = lambda av: r**2 - 2 * M * r + av**2

# ---- #7: Schwarzschild tau (a=0), E=1.4 ----
S_schw = sp.expand(r * (r - 2 * M) * DE * (r * Dl(0) - J**2 * DE)).subs(Es, sp.Rational(7, 5))
classify(S_schw, 0.0, "#7 Schwarzschild tau, E=1.4")

# ---- #8: Thakurta-Kerr t-branch, a=0.9, E=1.2 ----
a = sp.Rational(9, 10)
Q2 = (2*Es**2*J**2*M*r - Es**2*J**2*r**2 - 4*Es**2*J*M*a*r + 2*Es**2*M*a**2*r + Es**2*a**2*r**2
      + Es**2*r**4 + 4*J**2*M**2 - 4*J**2*M*r + J**2*r**2 - 8*J*M**2*a + 4*J*M*a*r + 4*M**2*a**2)
R6 = sp.expand(r * Q2 * DE).subs(Es, sp.Rational(6, 5))
classify(R6, 0.9, "#8 Thakurta-Kerr t-branch, a=0.9, E=1.2")

# ---- PHYSICAL separatrix of the t branch: refine the outer (J_c, r_d) and compute the corrected b_i ----
print("\n=== t branch: outer PHYSICAL SEPARATRIX, b_i with source (E d_E + J d_J) ===")
R6p = sp.lambdify((r, J), R6, 'mpmath'); R6r = sp.lambdify((r, J), sp.diff(R6, r), 'mpmath')
sol = mp.findroot(lambda rd, jc: [R6p(rd, jc), R6r(rd, jc)], (mp.mpf('3.514'), mp.mpf('-8.0535')))
rd_v, Jc_v = float(sol[0]), float(sol[1])
print(f"  refined: J_c^phys = {Jc_v:.8f}   r_d = {rd_v:.8f}  (in the control domain, r_d>2M)")

Dl9 = r**2 - 2 * M * r + sp.Rational(81, 100)
Q2s = (2*Es**2*J**2*M*r - Es**2*J**2*r**2 - 4*Es**2*J*M*a*r + 2*Es**2*M*a**2*r + Es**2*a**2*r**2
       + Es**2*r**4 + 4*J**2*M**2 - 4*J**2*M*r + J**2*r**2 - 8*J*M**2*a + 4*J*M*a*r + 4*M**2*a**2)
R6s = sp.expand(r * Q2s * ((Es**2 - 1) * r + 2 * M))
Kts = r * ((Es**2 - 1) * r + 2 * M) * (J * (r - 2 * M) + 2 * M * a) / Dl9
Fs = Kts / sp.sqrt(R6s)
Ntot = sp.expand(sp.cancel((Es * sp.diff(Fs, Es) + J * sp.diff(Fs, J)) * R6s**sp.Rational(3, 2)).subs(Es, sp.Rational(6, 5)))
S6 = R6s.subs(Es, sp.Rational(6, 5))


def bcoeffs(S, N, rd):
    d = lambda ex, n: sp.diff(ex, r, n).subs(r, rd)
    Q4rd, Q4prd, Q4pprd = d(S, 2) / 2, d(S, 3) / 6, d(S, 4) / 12
    s = sp.sqrt(Q4rd); a1 = Q4prd / (4 * s); a2 = Q4pprd / 12
    Nr, Np, Npp = N.subs(r, rd), sp.diff(N, r).subs(r, rd), sp.diff(N, r, 2).subs(r, rd)
    h0 = Nr / Q4rd; Fp = (Np * Q4rd - Nr * Q4prd) / Q4rd**2
    Fpp = ((Npp * Q4rd - Nr * Q4pprd) * Q4rd - 2 * Q4prd * (Np * Q4rd - Nr * Q4prd)) / Q4rd**3
    h1 = Fp * s; h2 = sp.Rational(1, 2) * (Fpp * s**2 + Fp * (Q4prd / 2))
    return {'b3': h0 / s**3, 'b2': (h1 - 3 * a1 * h0) / s**3,
            'b1': (h2 - 3 * a1 * h1 + (6 * a1**2 - 3 * a2) * h0) / s**3}


bc = bcoeffs(S6.subs(J, Jc_v), Ntot.subs(J, Jc_v), rd_v)
for k in ['b1', 'b2', 'b3']:
    print(f"  {k} = {float(bc[k]):+.6f}")
print("\nOutcome: the PHYSICAL separatrix of the t branch is J_c^phys=-8.054, r_d=+3.514 (r_d>2M),")
print("       the only one in the control domain. J=+2.936 (r_+<r_d=1.512<2M) is an object")
print("       of type (A); J=+19.089/-18.671 (r_d<0), J=+1.267 (r_d<r_+) and J=a (r_d=0)")
print("       are algebraic degenerations. Schwarzschild tau (E=1.4)")
print("       has NO outer physical separatrix: its r_d=-3.36 is only J_deg.")
