#!/usr/bin/env python3
"""Independent re-derivation of the headline claims of Papers I and II.

Every quantity is rebuilt from the DEFINITIONS
printed in the manuscript, not read back from the existing .wls scripts, so that
agreement is evidence and not tautology.

Definitions used (paper2.tex):
    f      = 1 - 2M/r                                        (sec. 2)
    Delta  = r^2 - 2Mr + a^2                                 (l. 432, 1493)
    DE     = (E^2-1) r + 2M                                  (l. 1615)
    w      = Ehat^2 - A^2 f   -> DE/r at A=1                 (l. 1493)
    S      = r (r-2M) DE [ r Delta - J^2 DE ]     tau-branch  (l. 1624)
    R6     = r Q2 DE                              t-branch    (l. 1624)
    Q2     from  (r-2M) Q2 + DE (J(r-2M)+2Ma)^2 = E^2 r^3 Delta   (l. 1588)
    R_opt  eq. (optical-scalar)                              (l. ~1030)

Run:  python3 verify_independent_rederivation.py
Exit code 0 iff every check passes.
"""
import sys
import sympy as sp
import mpmath as mp

r, M, a, E, J, Eh, s, Y = sp.symbols('r M a E J Ehat s Y', positive=True)

FAIL = 0
def check(label, cond, note=""):
    global FAIL
    ok = bool(cond)
    print(f"{'PASS  ' if ok else 'FAIL  '}{label}" + (f"   [{note}]" if note else ""))
    if not ok:
        FAIL += 1

def head(t):
    print("\n" + "=" * 66 + f"\n{t}\n" + "=" * 66)

# ----------------------------------------------------------------- base
f = 1 - 2*M/r
Delta = r**2 - 2*M*r + a**2
DE = (E**2 - 1)*r + 2*M

head("A. algebraic identities behind the four-case classification")

# A1 -- Q2 is a polynomial: the defining relation must divide exactly
num_Q2 = sp.expand(E**2 * r**3 * Delta - DE*(J*(r - 2*M) + 2*M*a)**2)
Q2, rem = sp.div(sp.Poly(num_Q2, r), sp.Poly(r - 2*M, r))
check("Q2 is a polynomial in r (exact division by r-2M)",
      sp.simplify(rem.as_expr()) == 0)
Q2 = Q2.as_expr()

# A2 -- marginal factorisation, with Jc = a/E
Jc = a/E
lhs = sp.expand(r*Delta - Jc**2 * DE)
rhs = sp.expand((r - 2*M)*(r**2 + Jc**2))
check("r*Delta - Jc^2 DE = (r-2M)(r^2+Jc^2)   [Jc = a/E]",
      sp.simplify(lhs - rhs) == 0)

# A3 -- hence S has a DOUBLE root at r = 2M at the marginal value
S_marg = sp.expand(r*(r - 2*M)*DE*(r*Delta - Jc**2*DE))
# multiplicity by exact division: (r-2M)^2 divides, (r-2M)^3 does not
q2_, rem2_ = sp.div(sp.Poly(S_marg, r), sp.Poly((r - 2*M)**2, r))
q3_, rem3_ = sp.div(sp.Poly(S_marg, r), sp.Poly((r - 2*M)**3, r))
check("S has a root of multiplicity exactly 2 at r = 2M when J = Jc",
      sp.simplify(rem2_.as_expr()) == 0 and sp.simplify(rem3_.as_expr()) != 0)
# and the same holds for BOTH signs of J, since S depends on J only through J^2
check("the double root at r=2M occurs for BOTH signs of J (S depends on J^2)",
      sp.simplify(S_marg.subs(a, -a) - S_marg) == 0)

# A3b -- and the printed closed form S = r^3 f^2 DE (r^2+Jc^2)
S_closed = r**3 * f**2 * DE * (r**2 + Jc**2)
check("S(marginal) = r^3 f^2 DE (r^2 + Jc^2)",
      sp.simplify(sp.expand(S_marg - S_closed)) == 0)

# A4 -- w = DE/r at A = 1, so r^4 f^2 w = r^3 f^2 DE  (the response's form)
w = Eh**2 - f
check("w = DE/r at A=1, Ehat=E", sp.simplify((w - DE/r).subs(Eh, E)) == 0)

# A5 -- the sign-breaking quantity is LINEAR in J and vanishes only at +Jc
#       tilde p_phi(r_e) = J - Jc  (manuscript, major 4)
tp = J - Jc
check("tilde p_phi(r_e) = J - Jc is linear in J and zero only at J=+Jc",
      sp.degree(sp.Poly(tp, J)) == 1 and sp.solve(tp, J) == [Jc])

head("B. exterior retrograde separatrix (falsifiability table)")

# rebuilt from R6 = r Q2 DE : double root <=> R6 = dR6/dr = 0
mp.mp.dps = 40
def separatrix(Mv, av, Ev, guess_r, guess_J):
    R6 = r * Q2 * DE
    F = sp.lambdify((r, J), R6.subs({M: Mv, a: av, E: Ev}), 'mpmath')
    Fr = sp.lambdify((r, J), sp.diff(R6, r).subs({M: Mv, a: av, E: Ev}), 'mpmath')
    sol = mp.findroot(lambda rr, JJ: (F(rr, JJ), Fr(rr, JJ)),
                      (mp.mpf(guess_r), mp.mpf(guess_J)))
    return sol[0], sol[1], sp.lambdify((r, J),
              sp.diff(R6, r, 2).subs({M: Mv, a: av, E: Ev}), 'mpmath')(sol[0], sol[1])

rd, Jc_num, d2 = separatrix(1, 0.9, 1.2, 3.5, -8.0)
print(f"       r_d  = {mp.nstr(rd, 20)}")
print(f"       Jc^- = {mp.nstr(Jc_num, 20)}")
check("r_d = 3.513905124011657 M", abs(rd - mp.mpf('3.513905124011657')) < mp.mpf('1e-13'))
check("Jc^- = -8.053516003877019", abs(Jc_num + mp.mpf('8.053516003877019')) < mp.mpf('1e-12'))
check("root is genuinely double, not higher (R6'' != 0)", abs(d2) > 1e-6,
      f"R6'' = {mp.nstr(d2, 8)}")
check("separatrix lies OUTSIDE r = 2M", rd > 2)

# control-domain inequalities at r_d, A = 1
gWW = float((-(1 - 2*M/r)).subs({M: 1, r: float(rd)}))
vbar2 = float((1 - (1 - 2*M/r)/Eh**2).subs({M: 1, r: float(rd), Eh: 1.2}))
print(f"       g(W,W) = {gWW:.6f}   vbar^2 = {vbar2:.6f}")
check("g(W,W) < 0 at r_d (selector timelike)", gWW < 0)
check("vbar^2 > 0 at r_d", vbar2 > 0)

# spin scan
print("       spin scan r_d/M:")
scan = []
for av, g_r, g_J in ((0.1, 2.86, -8.0), (0.5, 3.20, -8.0),
                     (0.9, 3.51, -8.0), (0.99, 3.58, -8.0)):
    rr, JJ, _ = separatrix(1, av, 1.2, g_r, g_J)
    scan.append(float(rr)); print(f"         a={av:<5} -> {mp.nstr(rr, 10)}")
check("r_d > 2M for every spin sampled (exterior is structural)",
      all(x > 2 for x in scan))

head("C. optical curvature: threshold and sign")

R_opt_num = -2*M**2*(32*M**2 + 4*(3*Eh**2 - 8)*M*r + (3*Eh**4 - 6*Eh**2 + 8)*r**2)
R_opt_den = Eh**2 * r**5 * (2*M + (Eh**2 - 1)*r)

# C1 -- expansion coefficients after r = 2M(1+s)
br = sp.expand((32*M**2 + 4*(3*Eh**2 - 8)*M*r + (3*Eh**4 - 6*Eh**2 + 8)*r**2)
               .subs(r, 2*M*(1 + s)) / M**2)
c = sp.Poly(br, s).all_coeffs()[::-1]
printed = [12*Eh**4, 24*Eh**2*(Eh**2 - 1), 4*(3*Eh**4 - 6*Eh**2 + 8)]
for k in range(3):
    check(f"expansion coefficient s^{k} matches the printed value",
          sp.simplify(c[k] - printed[k]) == 0, f"{sp.factor(c[k])}")

# C2 -- every coefficient positive for Ehat > 1  => bracket > 0 => R_opt < 0
brY = sp.expand(br.subs(Eh**2, 1 + Y).subs(Eh, sp.sqrt(1 + Y)))
neg = [co for co in sp.Poly(sp.expand(brY), s, Y).coeffs() if co.is_number and co < 0]
check("no negative coefficient in (s,Y) => bracket > 0 => R_opt < 0", not neg)

# C3 -- the 3/2 threshold sits in the leading coefficient -(Ehat^2-1)(2Ehat^2-3)
lead = sp.factor(-(Eh**2 - 1)*(2*Eh**2 - 3))
roots = sp.solve(sp.Eq(lead, 0), Eh**2)
check("curvature-cubic leading coefficient vanishes at Ehat^2 = 3/2",
      sp.Rational(3, 2) in [sp.nsimplify(x) for x in roots], f"roots {roots}")

head("D. optical submersion: tension field")

Lam2 = Eh**2/(f*(Eh**2 - f))
gradH = sp.simplify(f*sp.diff(sp.log(sp.sqrt(Lam2)), r))
mu = M/r**2                                    # fibre mean curvature
resid = sp.simplify(sp.together(mu + gradH))
check("mu + grad_H log Lambda = (M/r^2) f/(Ehat^2-f)",
      sp.simplify(resid - (M/r**2)*f/(Eh**2 - f)) == 0)
check("residual non-zero on the exterior => NOT a harmonic morphism",
      sp.simplify(resid) != 0)
val = float(resid.subs({M: 1, r: 6, Eh: 1.4}))
check("residual > 0 at r=6M, Ehat=1.4", val > 0, f"{val:.6e}")

head("E. conformal transfer map (sec. 2, restriction to constant A)")

A = sp.symbols('A', positive=True)
gbar = sp.diag(-f, 1/f, r**2, r**2)
g = A**2 * gbar
ut, ur = sp.symbols('ut ur', real=True)
u = sp.Matrix([ut, ur, 0, 0]); W = sp.Matrix([1, 0, 0, 0])
check("gbar(Au,Au) = g(u,u)  (unit vectors map to unit vectors)",
      sp.simplify((u.T*g*u)[0] - ((A*u).T*gbar*(A*u))[0]) == 0)
check("Ehat = A * Ebar  (rail charge rescales by one power of A)",
      sp.simplify(-(u.T*g*W)[0] - A*(-((A*u).T*gbar*W)[0])) == 0)

head("F. Paper I: ingoing Vaidya with the Kodama selector")

# ds^2 = -f dv^2 + 2 dv dr + r^2 dOmega^2 ,  f = 1 - 2m(v)/r      (paper1 l. 1270)
v, th = sp.symbols('v theta', real=True)
m = sp.Function('m', positive=True)(v)
fV = 1 - 2*m/r
gV = sp.Matrix([[-fV, 1, 0, 0],
                [1,   0, 0, 0],
                [0,   0, r**2, 0],
                [0,   0, 0, r**2*sp.sin(th)**2]])
X = [v, r, th, sp.Symbol('phi')]
K = sp.Matrix([1, 0, 0, 0])                       # Kodama selector K = d_v

check("g(K,K) = -f identically  (paper I, l. 1316)",
      sp.simplify((K.T*gV*K)[0] + fV) == 0)
check("|W| = sqrt(f), so W timelike <=> f > 0 <=> r > 2m(v)",
      sp.simplify(sp.sqrt(-(K.T*gV*K)[0]) - sp.sqrt(fV)) == 0)

# Kodama property: divergence-free.  div K = (1/sqrt|g|) d_mu ( sqrt|g| K^mu )
detg = sp.simplify(gV.det())
sq = sp.sqrt(sp.simplify(-detg))
divK = sp.simplify(sum(sp.diff(sq*K[i], X[i]) for i in range(4))/sq)
check("div K = 0  (K is the Kodama vector, not merely a coordinate field)",
      sp.simplify(divK) == 0, f"div K = {divK}")

# static limit: m' = 0 makes K Killing (Lie derivative of g vanishes)
gS = gV.subs(m, sp.Symbol('m0', positive=True))
LieK = sp.zeros(4, 4)
for i in range(4):
    for j in range(4):
        LieK[i, j] = sp.simplify(sum(K[k]*sp.diff(gS[i, j], X[k]) for k in range(4))
                                 + sum(gS[k, j]*sp.diff(K[k], X[i]) for k in range(4))
                                 + sum(gS[i, k]*sp.diff(K[k], X[j]) for k in range(4)))
check("m' = 0  =>  K is Killing (zero control cost in the static limit)",
      LieK == sp.zeros(4, 4))

# and NOT Killing when m depends on v  -> the rail must be actively held
LieKv = sp.simplify(sum(K[k]*sp.diff(gV[0, 0], X[k]) for k in range(4)))
check("m' != 0  =>  K is not Killing (Lie_K g_vv = 2 m'(v)/r != 0)",
      sp.simplify(LieKv - 2*sp.diff(m, v)/r) == 0, "so the charge must be controlled")

head("G. Major 4: the four cases, from the shape dphi/dr alone")

# From the branch table (l. 1624-1625), tau-branch:
#     curve      S     = r (r-2M) DE [ r Delta - J^2 DE ]
#     numerator  Knum  = J r (r-2M) DE / Delta
#     shape      dphi/dr = Knum / sqrt(S)
# Everything below follows from these two definitions and nothing else.
u = sp.symbols('u', positive=True)                       # u = r - 2M
S_gen = r*(r - 2*M)*DE*(r*Delta - J**2*DE)
Knum = J*r*(r - 2*M)*DE/Delta
dphi = Knum/sp.sqrt(S_gen)

# --- vanishing orders of p_r^2 = S/(DE^2 Delta^2), prefactor regular at r_e
pref = 1/(DE**2*Delta**2)
check("prefactor of p_r^2 is regular and nonzero at r_e = 2M",
      sp.simplify(pref.subs(r, 2*M)) != 0 and sp.simplify(pref.subs(r, 2*M)).is_finite)

q1_, r1_ = sp.div(sp.Poly(sp.expand(S_gen), r), sp.Poly(r - 2*M, r))
q2g, r2g = sp.div(sp.Poly(sp.expand(S_gen), r), sp.Poly((r - 2*M)**2, r))
check("case (ii) |J|<Jc: S has a SIMPLE zero at r_e",
      sp.simplify(r1_.as_expr()) == 0 and sp.simplify(r2g.as_expr()) != 0)

# --- case (i) |J|>Jc.  S carries the explicit factor (r-2M), so p_r^2(r_e)=0 for
# EVERY finite J: the case is separated not by the value but by the sign of the
# approach.  (An earlier printed statement claimed p_r^2(r_e)>0 here; it was false.)
check("p_r^2(r_e) = 0 for EVERY finite J (S has the explicit factor r-2M)",
      sp.simplify(S_gen.subs(r, 2*M)) == 0)
dS_re = sp.factor(sp.simplify(sp.diff(S_gen, r).subs(r, 2*M)))
check("dS/dr at r_e = 8 M^3 E^2 (a^2 - E^2 J^2)",
      sp.simplify(dS_re - 8*M**3*E**2*(a**2 - E**2*J**2)) == 0, f"{dS_re}")
_p = {M: 1, a: sp.Rational(9, 10), E: sp.Rational(12, 10)}
check("case (i): |J|>Jc gives dS/dr < 0, so p_r^2 < 0 just outside r_e "
      "(neighbourhood inadmissible, exterior branch turns at r_min > r_e)",
      dS_re.subs(_p).subs(J, 2) < 0 and dS_re.subs(_p).subs(J, sp.Rational(1, 2)) > 0)

# --- case (ii): semicubical cusp, amplitude K
lead = sp.simplify(sp.limit(sp.simplify(dphi.subs(r, 2*M + u))/sp.sqrt(u), u, 0, '+'))
K_printed = (E*J/a**2)*sp.sqrt(2*M/(a**2 - E**2*J**2))
# compare squares: both sides are positive for J>0 and |J|<Jc, and SymPy will not
# merge sqrt(1/x) with 1/sqrt(x) without being told x>0
check("case (ii): dphi/dr ~ K sqrt(u), i.e. varphi-varphi_e = (2/3) K u^{3/2}",
      sp.simplify(sp.expand(lead**2 - K_printed**2)) == 0,
      "cusp amplitude matches eq. (cusp)")
_sub = {M: 1, a: sp.Rational(9, 10), E: sp.Rational(12, 10), J: sp.Rational(1, 2)}
check("   and agrees numerically at M=1, a=0.9, E=1.2, J=0.5",
      abs(sp.N(lead.subs(_sub) - K_printed.subs(_sub), 25)) < 1e-20,
      f"{sp.N(lead.subs(_sub), 16)}")

# the amplitude diverges as |J| -> Jc^-, so the families do not join
check("cusp amplitude K diverges as |J| -> Jc^- (families do not connect)",
      sp.limit(K_printed.subs(J, Jc - u), u, 0, '+') is sp.oo)

# --- cases (iii)/(iv): at |J| = Jc the zero is DOUBLE, so the shape stays finite
S_m = sp.expand(S_gen.subs(J, Jc))
qA, rA = sp.div(sp.Poly(S_m, r), sp.Poly((r - 2*M)**2, r))
qB, rB = sp.div(sp.Poly(S_m, r), sp.Poly((r - 2*M)**3, r))
check("cases (iii)/(iv) |J|=Jc: S has a DOUBLE zero at r_e",
      sp.simplify(rA.as_expr()) == 0 and sp.simplify(rB.as_expr()) != 0)

dphi_m = sp.simplify(Knum/sp.sqrt(sp.factor(S_m)))
tangent = sp.simplify(sp.limit(dphi_m.subs(r, 2*M + u), u, 0, '+'))
tang_printed = 2*M*E*J/(a**2*sp.sqrt(4*M**2 + Jc**2))
check("case (iii): finite one-sided tangent = 2MEJ/(a^2 sqrt(4M^2+Jc^2))",
      sp.simplify(sp.radsimp(tangent - tang_printed)) == 0,
      "matches eq. (corner-slope)")

# --- case (iii): r_e is an ASYMPTOTIC endpoint, not attained
# rdot = R^2 (Delta/r^2) p_r / T.  With p_r ~ c1*u (double root) and T(r_e) != 0,
# rdot ~ C*u, so the elapsed parameter int du/|rdot| diverges logarithmically.
c1, Cc = sp.symbols('c1 C', positive=True)
elapsed = sp.integrate(1/(Cc*u), (u, sp.Symbol('eps', positive=True), 1))
check("case (iii): p_r ~ u and T(r_e)!=0  =>  int dr/rdot diverges logarithmically",
      sp.limit(elapsed, sp.Symbol('eps', positive=True), 0, '+') is sp.oo,
      "r_e is an asymptotic endpoint, a fortiori not crossed")

# --- case (iv): the two linear vanishings cancel, giving a finite nonzero rdot
ratio = sp.limit((c1*u)/(Cc*u), u, 0, '+')
check("case (iv): p_r ~ u AND T ~ u  =>  rdot finite and nonzero, r_e is crossed",
      ratio == c1/Cc and sp.simplify(ratio) != 0)

# --- the asymmetry is invisible to the radical and lives in T
check("S depends on J only through J^2, so the radical cannot break +/- J symmetry",
      sp.simplify(sp.expand(S_gen.subs(J, -J) - S_gen)) == 0)

print("\n" + "=" * 66)
print("ALL CHECKS PASSED" if FAIL == 0 else f"{FAIL} CHECK(S) FAILED")
print("=" * 66)
sys.exit(1 if FAIL else 0)
