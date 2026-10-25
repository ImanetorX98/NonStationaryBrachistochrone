# -*- coding: utf-8 -*-
"""
Analytic reduction of the fixed-endpoint no-inversion statement.

The paper states that, with FIXED ENDPOINTS, the brachistochrone of the t branch stays
SHALLOWER than the tau one: r_min^t > r_min^tau. The pointwise
bound n_t/n_tau=E/f>1 is NOT enough to order the minima of two different functionals.
Correct: at the same J the order is even the OPPOSITE (t deeper); the fixed-endpoint
result is an inversion driven by the angle map J->Phi.

This script shows the reduction to two lemmas (it proves Lemma A and only samples Lemma B):
  Parametrise each branch by the turning radius r_min (which fixes J through the
  turning condition: rDelta=J^2 DE for tau, Q2=0 for t). Let
     Phi_br(r_min) = int_{r_min}^{r0} F_br(r; J_br(r_min)) dr   (half-angle).
  (A) At equal r_min, the t integrand dominates POINTWISE: F_t(r;J_t)>F_tau(r;J_tau)
      for r in (r_min,r0). Since J_{t,tau}(r_min) are ALGEBRAIC, this is a
      polynomial (not transcendental) inequality => Phi_t(r_min)>Phi_tau(r_min).
  (B) Phi_tau strictly decreasing on [r_min^t, r0).
  Together: assume that a real tau curve of the compared family (one turning point)
  with the same half-angle exists. Then at fixed Phi,
           Phi_tau(r_min^tau)=Phi_t(r_min^t)>Phi_tau(r_min^t)
           => r_min^tau < r_min^t  (thesis).
  The argument excludes r_min^tau >= r_min^t when such a tau curve exists; it does
  not construct it. When Phi_t(r_min^t) exceeds the maximum of Phi_tau, the compared
  family has no tau curve with that half-angle: this is not an inversion.
  (A) is proved symbolically below, with M free. (B) is NOT proved by
  this script: it holds on the regions of Proposition prop:lemmaB of Paper II
  (proved in appendix app:lemmaB, checked by
  paper2/verification/verify_lemma_B.py) and not in general, because at low energy
  Phi_tau has several maxima. The numerical blocks on (B) below are diagnostics on
  samples, not a proof.
"""
import numpy as np
import sympy as sp
from scipy.integrate import quad


def branch_data(a, E, M=1.0, r0=10.0):
    r, J = sp.symbols('r J', real=True)
    DE = (E**2 - 1) * r + 2 * M
    Delta = r**2 - 2 * M * r + a**2
    Stau = sp.expand(r * (r - 2 * M) * DE * (r * Delta - J**2 * DE))
    Ktau = J * r * (r - 2 * M) * DE / Delta
    Ftau = sp.lambdify((r, J), sp.Abs(Ktau) / sp.sqrt(Stau), 'numpy')
    Q2 = (2*E**2*J**2*M*r - E**2*J**2*r**2 - 4*E**2*J*M*a*r + 2*E**2*M*a**2*r + E**2*a**2*r**2
          + E**2*r**4 + 4*J**2*M**2 - 4*J**2*M*r + J**2*r**2 - 8*J*M**2*a + 4*J*M*a*r + 4*M**2*a**2)
    R6 = sp.expand(r * Q2 * DE)
    Kt = r * DE * (J * (r - 2 * M) + 2 * M * a) / Delta
    Ft = sp.lambdify((r, J), sp.Abs(Kt) / sp.sqrt(R6), 'numpy')
    J2tau = sp.lambdify(r, r * Delta / DE, 'numpy')            # J_tau(r_min): rDelta=J^2 DE
    Ac, Bc, Cc = [sp.lambdify(r, c, 'numpy') for c in sp.Poly(Q2, J).all_coeffs()]

    def Jt(rm):                                                # J_t(r_min): prograde root of Q2=0
        A, B, C = Ac(rm), Bc(rm), Cc(rm); d = B * B - 4 * A * C
        if d < 0: return np.nan
        c = [x for x in ((-B + np.sqrt(d)) / (2 * A), (-B - np.sqrt(d)) / (2 * A)) if x > 0]
        return min(c) if c else np.nan
    return Ftau, Ft, J2tau, Jt, r0


def phis(a, E):
    Ftau, Ft, J2tau, Jt, r0 = branch_data(a, E)

    def Phi_tau(rm):
        return quad(lambda x: Ftau(x, np.sqrt(J2tau(rm))), rm + 1e-6, r0, limit=200)[0]

    def Phi_t(rm):
        Jv = Jt(rm)
        return np.nan if np.isnan(Jv) else quad(lambda x: Ft(x, Jv), rm + 1e-6, r0, limit=200)[0]
    return Phi_tau, Phi_t, Ftau, Ft, J2tau, Jt, r0


print("Check of the two lemmas for several (a,E):")
print("SAMPLED DIAGNOSTIC (not a proof): pointwise Lemma A and sampled monotonicity of Phi_tau")
print("  a    E     LemmaA (F_t>F_tau ptwise)  Phi_tau sampled, decreasing   sampled A and B checks")
for a, E in [(0.3, 1.2), (0.6, 1.2), (0.9, 1.2), (0.9, 1.05), (0.9, 1.5), (0.5, 1.3)]:
    Phi_tau, Phi_t, Ftau, Ft, J2tau, Jt, r0 = phis(a, E)
    A_ok = B_ok = True; prev = None
    for rm in np.linspace(3, 8, 12):
        Jta, Jtt = np.sqrt(J2tau(rm)), Jt(rm)
        if np.isnan(Jtt): continue
        # (A) pointwise integrand
        for rr in np.linspace(rm + 0.2, r0, 8):
            if Ft(rr, Jtt) <= Ftau(rr, Jta): A_ok = False
        # (B) monotonicity of Phi_tau, the function the deduction needs
        pt = Phi_tau(rm)
        if prev is not None and pt > prev: B_ok = False
        prev = pt
    sampled_AB = A_ok and B_ok
    print(f"  {a:.1f}  {E:.2f}   {str(A_ok):>5}                      {str(B_ok):>5}"
          f"                     {str(sampled_AB):>5}")
print("The last column only combines the two sampled checks. It does not solve the matching")
print("Phi_tau(r_min^tau) = Phi_t(r_min^t), so it neither tests that a tau curve with that half-angle")
print("exists nor proves that Phi_tau decreases between grid points (12 points on [3,8]).")
print("\nLemma A is ALGEBRAIC (J_{t,tau}(r_min) are roots of polynomials) => a symbolic proof")
print("is possible through polynomial positivity (done below). For (B) see verify_lemma_B.py.")


# ============================================================================
# SYMBOLIC PROOF of Lemma A at the TURNING POINT (dominant singular part).
# The sign of F_t^2 - F_tau^2 near r_min is given by the RESIDUES: F^2 ~ c/(r-r_min).
# c_t > c_tau  <=>  D := numerator of (c_t-c_tau) > 0. With the turning conditions
# (J_tau^2 = r_d Delta_d/DE_d, Q2(r_d,J_t)=0) it reduces to:
#     D = 4 r_d Delta(r_d) * B,   B = a DE(r_d) (J* - J_t),
#     J* = [E^2 r_d(a^2+r_d^2)+2Ma^2]/(a DE_d)  (>0).
# Chain of signs (r_d>2M, r_d>r_+, E>1, a>0), all factors explicit:
#   A   = coefficient of J_t^2 in Q2(r_d) = -(r_d-2M) DE_d     < 0   (downward parabola)
#   J_v = vertex = -B_q/2A           = -2Ma/(r_d-2M)           < 0
#   Q2(r_d,J*) = -E^2 r_d^3 (E^2 r_d^2+a^2) Delta_d/[a^2 DE_d] < 0
#   => J*>0>J_v and Q2(J*)<0 with A<0  => J* beyond the large root => J_t<J*
#   => B>0 => D>0 => c_t>c_tau => F_t^2>F_tau^2 at the turning point.  QED (singular part)
# The INTERIOR positivity (r away from r_min) reduces to N_G(r)>0 (degree 5),
# checked numerically above.
def symbolic_residue_proof():
    rd, a, E, M = sp.symbols('r_d a E M', positive=True)   # M kept symbolic
    Jt = sp.symbols('J_t', positive=True)
    Q2d = (2*E**2*Jt**2*M*rd - E**2*Jt**2*rd**2 - 4*E**2*Jt*M*a*rd + 2*E**2*M*a**2*rd
           + E**2*a**2*rd**2 + E**2*rd**4 + 4*Jt**2*M**2 - 4*Jt**2*M*rd + Jt**2*rd**2
           - 8*Jt*M**2*a + 4*Jt*M*a*rd + 4*M**2*a**2)
    P = sp.Poly(Q2d, Jt)
    A, Bq = P.coeff_monomial(Jt**2), P.coeff_monomial(Jt)
    DEd = (E**2-1)*rd + 2*M
    Jstar = (E**2*rd*(a**2+rd**2) + 2*M*a**2)/(a*DEd)
    assert sp.simplify(A - (-(rd-2*M)*DEd)) == 0
    assert sp.simplify(-Bq/(2*A) - (-2*M*a/(rd-2*M))) == 0
    assert sp.simplify(Q2d.subs(Jt, Jstar)
                       + E**2*rd**3*(E**2*rd**2+a**2)*(rd**2-2*M*rd+a**2)/(a**2*DEd)) == 0
    Q2star = sp.factor(Q2d.subs(Jt, Jstar))
    print("  A            =", sp.factor(A), "   (<0 for r_d>2M)")
    print("  vertex J_v   =", sp.factor(-Bq/(2*A)), "   (<0)")
    print("  J*           =", sp.factor(Jstar), "   (>0)")
    print("  Q2(r_d,J*)   =", Q2star, "   (<0 for r_d>r_+)")
    print("  => J_t < J* => B>0 => residue ordering proven.  QED")


if __name__ == '__main__':
    print("\n=== Symbolic proof of Lemma A at the turning point ===")
    symbolic_residue_proof()


# ============================================================================
# COMPLETE SYMBOLIC PROOF of Lemma A (F_t^2 > F_tau^2 for every r>r_d).
# Steps (M symbolic; valid for r_d>2M, E>1, a>0):
#   sign(F_t^2-F_tau^2) = sign(N_G),  N_G = (J_t(r-2M)+2Ma)^2 P_tau - J_tau^2(r-2M)Q2.
#   Reduced modulo the turning conditions [J_tau^2=r_d Delta_d/DE_d, Q2(r_d,J_t)=0]:
#     N_G = (r - r_d) * [4 M r Delta(r)/((r_d-2M) DE_d)] * W(r),   W(r) LINEAR in r.
#   W(r)>0 is needed for r>r_d. W(r_d) and the slope w1 are both linear in J_t:
#     W(r_d) = (r_d-2M) * B,   B = a DE_d (J*  - J_t),   J*  = [E^2 r_d(a^2+r_d^2)+2Ma^2]/(a DE_d)
#     w1     = a DE_d (J** - J_t),                        J** = [E^2 r_d^2(r_d-M)+E^2 a^2 r_d + Ma^2]/(a DE_d)
#   In both, J_t (prograde root of Q2(r_d,.)) is < J* and < J** because:
#     A = coefficient of J_t^2 in Q2(r_d) = -(r_d-2M)DE_d < 0   (downward parabola)
#     vertex J_v = -2Ma/(r_d-2M) < 0,  J*, J** > 0         (to the right of the vertex)
#     Q2(r_d,J*)  = -E^2 r_d^3 (E^2 r_d^2+a^2) Delta_d / [a^2 DE_d] < 0
#     Q2(r_d,J**) = -(r_d-2M)(E^2 r_d^2+a^2)(E^2 r_d^2[(r_d-M)^2+a^2]+M^2a^2)/[a^2 DE_d] < 0
#   => J* , J** beyond the large root => J_t < J*, J_t < J** => B>0, w1>0.
#   W linear with W(r_d)>0 and slope w1>0  =>  W(r)>0 for every r>r_d.  QED (Lemma A).
def symbolic_full_proof():
    rd, a, E, M = sp.symbols('r_d a E M', positive=True)   # M kept symbolic
    Jt = sp.symbols('J_t', positive=True)
    DEd = (E**2 - 1) * rd + 2 * M
    Q2d = (2*E**2*Jt**2*M*rd - E**2*Jt**2*rd**2 - 4*E**2*Jt*M*a*rd + 2*E**2*M*a**2*rd
           + E**2*a**2*rd**2 + E**2*rd**4 + 4*Jt**2*M**2 - 4*Jt**2*M*rd + Jt**2*rd**2
           - 8*Jt*M**2*a + 4*Jt*M*a*rd + 4*M**2*a**2)
    A = sp.Poly(Q2d, Jt).coeff_monomial(Jt**2)
    Jstar = (E**2*rd*(a**2+rd**2) + 2*M*a**2)/(a*DEd)                     # per W(r_d)>0
    Jss = (E**2*rd**2*(rd-M) + E**2*a**2*rd + M*a**2)/(a*DEd)             # per w1>0
    assert sp.simplify(A + (rd-2*M)*DEd) == 0                            # A = -(r_d-2M)DE_d
    Q2star = sp.factor(Q2d.subs(Jt, Jstar))
    Q2ss = sp.factor(Q2d.subs(Jt, Jss))
    # manifest sign of the two Q2 values (all factors signed for r_d>2M,E>1,a>0)
    print("  A = -(r_d-2M)DE_d :", sp.factor(A), " (<0)")
    print("  Q2(r_d,J*)  =", Q2star, " (<0)")
    print("  Q2(r_d,J**) =", Q2ss, " (<0)")
    assert sp.simplify(Q2ss + (rd-2*M)*(E**2*rd**2+a**2)
                       * (E**2*rd**2*((rd-M)**2+a**2) + M**2*a**2)/(a**2*DEd)) == 0
    print("  Q2(r_d,J**) = -(r_d-2M)(E^2 r_d^2+a^2)(E^2 r_d^2[(r_d-M)^2+a^2]+M^2 a^2)/(a^2 DE_d): True")
    # key factorization N_G = (r-r_d) 4 M r Delta(r)/[(r_d-2M) DE_d] W(r), W linear,
    # W(r_d) = (r_d-2M) a DE_d (J*-J_t), slope a DE_d (J**-J_t)
    r = sp.symbols('r', positive=True)
    Q2 = lambda rr, J: (E**2*rr**4 + (E**2*(a**2-J**2)+J**2)*rr**2
                        + 2*M*(J-a)*(E**2*(J-a)-2*J)*rr + 4*M**2*(J-a)**2)
    Dl = lambda rr: rr**2 - 2*M*rr + a**2
    DEr = (E**2-1)*r + 2*M
    Jtau2 = rd*Dl(rd)/DEd
    NG = (Jt*(r-2*M)+2*M*a)**2*(r*Dl(r) - Jtau2*DEr) - Jtau2*(r-2*M)*Q2(r, Jt)
    Wlin = (r-rd)*a*DEd*(Jss-Jt) + (rd-2*M)*a*DEd*(Jstar-Jt)
    diff = sp.together(NG - (r-rd)*4*M*r*Dl(r)/((rd-2*M)*DEd)*Wlin)
    num = sp.Poly(sp.expand(sp.numer(diff)), Jt)
    rem = sp.rem(num, sp.Poly(sp.expand(sp.numer(sp.together(Q2(rd, Jt)))), Jt))
    assert sp.expand(rem.as_expr()) == 0
    print("  N_G = (r-r_d) 4 M r Delta/[(r_d-2M) DE_d] W(r) mod Q2(r_d,J_t)=0, W linear: True")
    print("  => J_t<J* and J_t<J** => W(r_d)>0 and slope w1>0 => W(r)>0 for r>r_d.  Lemma A QED")


if __name__ == '__main__':
    print("\n=== Full symbolic proof of Lemma A (all r>r_d) ===")
    symbolic_full_proof()


# ============================================================================
# LEMMA B (monotonicity of Phi_tau) -- SAMPLED DIAGNOSTIC, not a proof.
# Phi_tau(r_min) is not globally monotonic. In the cases sampled here it has a single
# maximum r_pk ~ 3M: a fold of the fixed-endpoint half-angle map, NOT
# an unstable circular orbit of the tau branch (which does not exist outside: on the turning
# shell d/dr[r Delta - J^2 DE] = DE V' > 0). At low energy it can have several
# maxima (E^2 = 1.01, r0 = 1e4 M). The proved regions are those of
# Proposition prop:lemmaB of Paper II; this scan only illustrates them:
def lemma_B_regime():
    from scipy.integrate import quad
    M = 1.0
    def Phi_tau(rm, a, E, r0=10.0):
        Delta = lambda r: r*r - 2*M*r + a*a; DE = lambda r: (E*E-1)*r + 2*M
        J = np.sqrt(rm*Delta(rm)/DE(rm))
        def g(s):
            r = rm + s*s; P = r*Delta(r) - J*J*DE(r)
            return 2*J*np.sqrt(r*(r-2*M)*DE(r))/(Delta(r)*np.sqrt(P/(r-rm)))
        return quad(g, 0, np.sqrt(r0-rm), limit=200)[0]
    for a, E in [(0.2, 1.2), (0.5, 1.2), (0.9, 1.2), (0.9, 1.05), (0.9, 1.5)]:
        rms = np.linspace(2.2, 9, 120); ph = np.array([Phi_tau(rm, a, E) for rm in rms])
        ipk = int(np.argmax(ph)); rpk = rms[ipk]
        mono = all(ph[i+1] < ph[i] + 1e-9 for i in range(ipk, len(ph)-1))
        print(f"  a={a:.1f} E={E:.2f}:  r_pk={rpk:.2f}   Phi_tau decreasing for r_min>r_pk? {mono}")


if __name__ == '__main__':
    print("\n=== Lemma B: sampled diagnostic, not a proof (see verify_lemma_B.py) ===")
    lemma_B_regime()


# ============================================================================
# LEMMA B -- partial closed-form proof (grazing regime) + structure.
# Effective-potential form: Phi_tau = sqrt(Vmin) int_{rmin}^{r0} K/sqrt(V-Vmin) dr,
#   V=r Delta/DE,  K=sqrt(r(r-2M))/Delta,  Vmin=V(rmin),  W=K/V' (>0, W'<0).
# With V=Vmin+sigma^2 and an integration by parts:
#   dPhi/dVmin = (V0-2Vmin) W(Vmin)/(Sigma sqrt(Vmin))
#                + (2/sqrt(Vmin)) int_0^Sigma W'(Vmin+sigma^2) g(sigma) dsigma,
#   g(sigma) = (V0-2Vmin) sigma/Sigma + Vmin - sigma^2,  Sigma=sqrt(V0-Vmin).
# g(0)=Vmin>0, g(Sigma)=0, other root <0 => g>0 on [0,Sigma) => second term <0.
# If V(r0)<=2 V(rmin) (grazing): first term <=0 => dPhi/dVmin<0 PROVED.
# If V(r0)>2 V(rmin): first term >0, dominated by the integral (checked, but
# the inequality is tight at the transcendental threshold r_pk: symbolic proof open).
def lemma_B_partial():
    import sympy as sp
    r, a, E, s = sp.symbols('r a E sigma', positive=True); M = 1
    Vmin, V0, Sig = sp.symbols('V_min V_0 Sigma', positive=True)
    g = (V0 - 2*Vmin)*s/Sig + Vmin - s**2
    print("  g(sigma) =", g)
    print("  g(0) =", g.subs(s, 0), " (>0);  g(Sigma) =", sp.simplify(g.subs(s, Sig).subs(Sig**2, V0-Vmin)), " (=0)")
    roots = sp.solve(sp.Eq(g, 0), s)
    print("  roots of g:", [sp.simplify(rt.subs(Sig, sp.sqrt(V0-Vmin))) for rt in roots],
          " => one =Sigma, the other <0  =>  g>0 on [0,Sigma)")
    print("  => second term (2/sqrt Vmin) int W' g dsigma < 0 (W'<0).")
    print("  Grazing V0<=2Vmin: first term <=0 => dPhi/dVmin<0 PROVED in closed form.")


if __name__ == '__main__':
    print("\n=== Lemma B: closed-form partial proof (grazing regime) ===")
    lemma_B_partial()
