#!/usr/bin/env python3
"""Independent SymPy cross-check of the formal core of Paper I (JMP variant).

Mirrors verify_paper1_core.wls in a second computer-algebra system.  Every
check prints [OK] only if the residual is an exact zero (or, where a
determinant is evaluated at sample rational parameter points, exactly the
predicted value).

    python3 verify_paper1_core.py

Claims covered, with their location in paper1_JMP.tex:

    Lemma I.A     rail indicatrix; Vaidya ellipse; support function = eq:Hv
    Theorem I.1   (ii) closed-form branch Hamiltonian
    Lemma I.4     det M = -32 (Ehat^2 - 1) disc_r S; one-dimensional kernel
                  at the genus-degeneration locus
    Lemma I.6     Euler homogeneity of S; N_tot(r_d) = 0
    Theorem I.5   Randers-Euler, self-similarity, Theta H_v, boundary form
"""
import sympy as sp

FAILS = []


def ok(label, value):
    """Print [OK] iff the residual is an exact zero.

    The second pass with powsimp(force=True) is needed only to let SymPy
    combine sqrt(A)*sqrt(B) into sqrt(A*B); on the regular domain
    (r > 2m > 0, Ehat > 1) both radicands are positive, so the rewriting is
    valid there.  Mathematica performs it directly from the assumptions.
    """
    z = sp.simplify(value)
    if z != 0:
        z = sp.simplify(sp.powsimp(sp.radsimp(z), force=True))
    good = (z == 0)
    print(f"  [{'OK' if good else 'FAIL'}]   {label}  ->  {z}")
    if not good:
        FAILS.append(label)


# plain symbols: r must be allowed to be negative (the double root is)
r, m, E, J, pr, uv, Vr, Vphi, al = sp.symbols('r m E J pr uv Vr Vphi alpha')

print("=" * 66)
print(" LEMMA I.A -- rail indicatrix (ingoing Vaidya)")
print("=" * 66)

f = 1 - 2 * m / r
w = E**2 - f

# eliminate u^v between g(u,u) = -1 and -u_v = Ehat
shell = -f * uv**2 + 2 * uv * (uv * Vr) + r**2 * (uv * Vphi / r)**2 + 1
rail = f * uv - uv * Vr - E
uvsol = sp.solve(rail, uv)[0]
ind = sp.simplify(sp.numer(sp.together(shell.subs(uv, uvsol))))

ok("indicatrix == r^2 [ (Vr+(E^2-f))^2 + E^2 Vphi^2 - E^2 w ]",
   sp.expand(ind - r**2 * ((Vr + (E**2 - f))**2 + E**2 * Vphi**2 - E**2 * w)))

th = sp.Symbol('theta')
ok("paper's eq:vaidya-ind parametrization satisfies it",
   sp.simplify(ind.subs({Vr: (f - E**2) + E * sp.sqrt(w) * sp.cos(th),
                         Vphi: sp.sqrt(w) * sp.sin(th)})))

Hv = pr * (f - E**2) - 1 + sp.sqrt(w) * sp.sqrt(E**2 * pr**2 + J**2 / r**2)
sup = pr * (f - E**2) + sp.sqrt((pr * E * sp.sqrt(w))**2 + (J * sp.sqrt(w) / r)**2)
ok("support function of the ellipse - 1 - H_v  [Theorem I.1(ii)]",
   sp.simplify(sup - 1 - Hv))

print("\n" + "=" * 66)
print(" LEMMA I.4 -- the 11x11 reduction system")
print("=" * 66)


def reduction_matrix(mv, Ev, Jv):
    D_E = (Ev**2 - 1) * r + 2 * mv
    S = sp.expand(r * (r - 2 * mv) * D_E * (r**2 * (r - 2 * mv) - Jv**2 * D_E))
    a = sp.symbols('a0:6')
    c = sp.symbols('c0:5')
    A = sum(a[i] * r**i for i in range(6))
    P = sp.expand(2 * sp.diff(A, r) * S - A * sp.diff(S, r)
                  + 2 * S * sum(c[k] * r**k for k in range(5)))
    unk = list(a) + list(c)
    poly = sp.Poly(P, r)
    rows = [[sp.expand(poly.coeff_monomial(r**n)).coeff(u) for u in unk]
            for n in range(11)]
    return sp.Matrix(rows), S


print("\ndet M against -32 (E^2-1) disc_r S at sample rational points:")
for mv, Ev, Jv in [(sp.Integer(1), sp.Rational(7, 5), sp.Rational(5, 2)),
                   (sp.Rational(3, 2), sp.Rational(6, 5), sp.Rational(4, 3)),
                   (sp.Integer(2), sp.Rational(9, 5), sp.Rational(7, 2))]:
    M, S = reduction_matrix(mv, Ev, Jv)
    pred = -32 * (Ev**2 - 1) * sp.discriminant(sp.Poly(S, r))
    ok(f"m={mv}, E={Ev}, J={Jv}:  det M - (-32 (E^2-1) disc S)",
       sp.nsimplify(M.det() - pred))
    g = sp.gcd(sp.Poly(S, r), sp.Poly(sp.diff(S, r), r)).as_expr()
    print(f"          gcd(S,S') = {g}   (1 => S squarefree => nonsingular)")

# degeneration locus at m = 1, E = 7/5
mv, Ev = sp.Integer(1), sp.Rational(7, 5)
rds = sp.solve(sp.Eq((Ev**2 - 1) * r**2 - mv * (Ev**2 - 4) * r - 4 * mv**2, 0), r)
J2 = lambda x: sp.simplify((3 * x**2 - 4 * mv * x) / (Ev**2 - 1))
rd = [x for x in rds if sp.N(J2(x)) > 0][0]
Jd = sp.sqrt(J2(rd))
print(f"\nr_d = {sp.N(rd, 16)}   J_deg = {sp.N(Jd, 16)}")
ok("J_deg - paper's closed form 5*sqrt(3011/3072 + 581*sqrt(249)/9216)",
   sp.radsimp(Jd - 5 * sp.sqrt(sp.Rational(3011, 3072)
                               + 581 * sp.sqrt(249) / sp.Integer(9216))))

M, S = reduction_matrix(mv, Ev, Jd)
Sn = sp.expand(S)
ok("S(r_d)", sp.simplify(Sn.subs(r, rd)))
ok("S'(r_d)", sp.simplify(sp.diff(Sn, r).subs(r, rd)))
ok("det M at J = J_deg", sp.simplify(M.det()))
print(f"  rank M at J = J_deg = {M.rank()}  (of 11: kernel exactly 1-dimensional)")

Q4 = sp.Poly(sp.div(sp.Poly(Sn, r), sp.Poly((r - rd)**2, r))[0], r)
Aker = sp.expand((r - rd) * Q4.as_expr())
ok("S - (r-r_d)^2 Q4", sp.expand(Sn - (r - rd)**2 * Q4.as_expr()))
print(f"  deg Q4 = {Q4.degree()}, deg A_kernel = {sp.degree(Aker, r)}"
      "  (must be 4 and 5)")
quot, rem = sp.div(sp.Poly(sp.expand(2 * sp.diff(Aker, r) * Sn
                                     - Aker * sp.diff(Sn, r)), r), sp.Poly(Sn, r))
ok("remainder of 2A'S - A S' modulo S", sp.simplify(rem.as_expr()))
ok("C - (A' - Q4)", sp.simplify(quot.as_expr() - (sp.diff(Aker, r) - Q4.as_expr())))
print(f"  deg C = {quot.degree()}  (must be <= 4)")

print("\n" + "=" * 66)
print(" LEMMA I.6 -- degeneration-family derivative")
print("=" * 66)
D_E = (E**2 - 1) * r + 2 * m
Ssym = sp.expand(r * (r - 2 * m) * D_E * (r**2 * (r - 2 * m) - J**2 * D_E))
ok("Euler:  r S_r + m S_m + J S_J - 6 S",
   sp.expand(r * sp.diff(Ssym, r) + m * sp.diff(Ssym, m)
             + J * sp.diff(Ssym, J) - 6 * Ssym))
ok("explicit scaling S(a r, a m, a J) - a^6 S",
   sp.expand(Ssym.subs({r: al * r, m: al * m, J: al * J}) - al**6 * Ssym))

K = J * D_E
N = Ssym * sp.diff(K, m) - K * sp.diff(Ssym, m) / 2
NJ = Ssym * sp.diff(K, J) - K * sp.diff(Ssym, J) / 2
Ntot = sp.expand(N + (J / m) * NJ)
sub = {m: mv, E: Ev, J: Jd}
ok("dS/dm + (J/m) dS/dJ  at r_d",
   sp.simplify((sp.diff(Ssym, m) + (J / m) * sp.diff(Ssym, J)).subs(sub).subs(r, rd)))
print(f"  N(r_d)   = {sp.N(sp.simplify(N.subs(sub).subs(r, rd)), 12)}   (nonzero)")
print(f"  N_J(r_d) = {sp.N(sp.simplify(NJ.subs(sub).subs(r, rd)), 12)}   (nonzero)")
ok("N_tot(r_d)", sp.simplify(Ntot.subs(sub).subs(r, rd)))

print("\n" + "=" * 66)
print(" THEOREM I.5 -- first-order response and the boundary form")
print("=" * 66)
ok("Randers-Euler  (J d_J + pr d_pr) H - (H+1)",
   sp.simplify((J * sp.diff(Hv, J) + pr * sp.diff(Hv, pr)) - (Hv + 1)))
ok("self-similarity  r H_r + m H_m + J H_J  (degree zero, off shell too)",
   sp.simplify(m * sp.diff(Hv, m) + r * sp.diff(Hv, r) + J * sp.diff(Hv, J)))
ok("explicit scaling H(a r, a m, a J) - H",
   sp.simplify(Hv.subs({r: al * r, m: al * m, J: al * J}) - Hv))
ok("Theta H_v - eq:ThetaHv",
   sp.simplify(m * sp.diff(Hv, m)
               - m * (-2 * pr / r
                      + sp.sqrt(E**2 * pr**2 + J**2 / r**2) / (r * sp.sqrt(w)))))
ok("d(r pr)/dlambda - (1 + Theta H) - H   [identically 0]",
   sp.simplify((pr * sp.diff(Hv, pr) - r * sp.diff(Hv, r))
               - (1 + m * sp.diff(Hv, m)) - Hv))
print("  => on the shell H = 0:  Theta H = d(r pr)/dlambda - 1,")
print("     hence  S_D(lambda) = [r pr]_0^lambda - lambda.   (eq:SD-vaidya)")

print("\n" + "=" * 66)
if FAILS:
    print(f" {len(FAILS)} CHECK(S) FAILED: " + "; ".join(FAILS))
    raise SystemExit(1)
print(" ALL CHECKS PASSED -- every residual is an exact zero.")
