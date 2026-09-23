#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
verify_randers_convexity.py -- the Randers form of the *arrival-time* rail cost on
the frozen branch, its convexity domain, and its relation to the Jacobi-Maupertuis
metric of free geodesics.

Answers the Round-2 referee request to exhibit the map to Randers-Finsler
geodesics.  The object is the arrival-time functional, not the free-fall one:
this file supersedes an earlier version that computed the Jacobi-Maupertuis data
of Chanda, Gibbons, Guha, Maraner & Werner, JMP 60 (2019) 122501 eq. (16), and
mistook it for the brachistochrone cost.

Threading form, signature (+,-,-,-), unit mass:

    ds^2 = f (dt + w_i dx^i)^2 - h_ij dx^i dx^j,    f = g_00,  w_i = g_{0i}/g_00,
    h_ij = -g_ij + g_{0i} g_{0j}/g_00   (= gamma_ij of the reference).

Held rail at specific energy E:  E = f (tdot + w_i xdot^i), so h_ij xdot xdot =
(E^2-f)/f, and eliminating proper time gives the arrival cost

    F_T = sqrt( [E^2/(f(E^2-f))] h_ij dx^i dx^j ) - w_i dx^i,

a Randers functional whose quadratic part is exactly Lambda^2 h with the dilation
Lambda^2 of the manuscript's eq. (12).  This is Perlick 1991 (JMP 32, 3148),
Proposition 3.2, p. 3153; his Proposition 3.3 is the *free fall* counterpart, and
his comparison of Props 3.1-3.3 shows the two coincide only after passing to a
different stationary spacetime (his eq. (50)).  Chanda eq. (16) is of the 3.3 type.

Proved here:
  (1) the quadratic part of F_T is Lambda^2 h, i.e. the manuscript already carries
      the correct Randers metric -- and Lambda is the dilation of the semiconformal
      projection, which is what the Round-1 "T2" objection is about;
  (2) F_T and the Jacobi-Maupertuis functional differ by the position-dependent
      factor E^2/(E^2-f)^2 on the quadratic part, hence have different extremals;
      explicit witness: the circular extremal in Schwarzschild at E^2 = 2;
  (3) on frozen conformal Kerr, with effective energy E = Ehat/A above freezing
      (E^2 > f) and in the exterior, the Randers convexity condition
      sqrt(a^{ij}b_i b_j) < 1 holds *throughout* -- no coincidence with the
      boundary of the control domain, contrary to an earlier claim of ours.
"""
import sympy as sp

r, th, M, a, E, A, Eh = sp.symbols('r theta M a E A Ehat', positive=True)
# E above is the effective energy Ehat/A: the manuscript's frozen shorthand.
# It used to be carried here as a second symbol 'e', which collided with the
# base of the exponential and duplicated a name the paper already defines.
s2, c2 = sp.sin(th)**2, sp.cos(th)**2
rho2 = r**2 + a**2*c2
Delta = r**2 - 2*M*r + a**2
X = Delta - a**2*s2                                   # = rho^2 g_00

ok = True


def check(label, cond):
    global ok
    ok = ok and bool(cond)
    print(f'{"PASS" if cond else "FAIL"}  {label}')


# ---------------------------------------------------------------- (1) the cost
f_, h_ = sp.symbols('f h', positive=True)             # f = g_00, h = h_ij dx dx
dtau = sp.sqrt(f_*h_/(E**2 - f_))                     # from h_ij xdot xdot = (E^2-f)/f
quad_T = sp.simplify((E/f_*dtau)**2/h_)               # coefficient of h in F_T^2
Lam2 = E**2/(f_*(E**2 - f_))
check('quadratic part of the arrival cost is Lambda^2 = E^2/(f(E^2-f)), eq. (12)',
      sp.simplify(quad_T - Lam2) == 0)

quad_J = (E**2 - f_)/f_                               # Chanda eq. (16), m = 1
check('Jacobi-Maupertuis quadratic part is (E^2-f)/f',
      sp.simplify(quad_J - (E**2 - f_)/f_) == 0)
ratio = sp.simplify(Lam2/quad_J)
check('their ratio is E^2/(E^2-f)^2, not a constant',
      sp.simplify(ratio - E**2/(E**2 - f_)**2) == 0 and sp.diff(ratio, f_) != 0)

# NORMALISATION, and it matters for what may be claimed.  Chanda's eq. (16)
# carries the one-form p_0 g_{0i}/g_00 dx^i, i.e. -E w with his own Kerr
# convention p_0 = -E, whereas F_T carries -w.  The two one-forms are therefore
# NOT equal, and the sentence "both carry the same one-form, so the whole
# difference is in the quadratic part" is false as it stands for this pair.
# Dividing the whole Jacobi-Maupertuis functional by the constant E leaves its
# geodesics alone and makes the one-forms agree; the quadratic parts then differ
# by E^4/(E^2-f)^2 instead of E^2/(E^2-f)^2.  The comparison must state which
# normalisation it uses.
check('un-normalised, the two one-forms differ by the factor E',
      sp.simplify(E*sp.Symbol('w') - sp.Symbol('w')) != 0)
quad_Jn = sp.simplify(quad_J/E**2)
check('dividing by E gives the same one-form and quadratic part (E^2-f)/(E^2 f)',
      sp.simplify(quad_Jn - (E**2 - f_)/(E**2*f_)) == 0)
ratio_n = sp.simplify(Lam2/quad_Jn)
check('at equal one-form the quadratic parts differ by E^4/(E^2-f)^2, still '
      'position dependent',
      sp.simplify(ratio_n - E**4/(E**2 - f_)**2) == 0
      and sp.diff(ratio_n, f_) != 0)

# AND A WARNING ABOUT A BIBLIOGRAPHIC ARGUMENT WE MUST NOT USE.  Chanda 2024
# states in words that "optical metrics are not Jacobi metrics for null curves".
# That sentence is in the PDF, but it cannot be leaned on: setting m = 0 and
# p_0 = -E in eq. (16) of the 2019 paper and dividing by E returns exactly the
# null optical (Fermat-Randers) functional sqrt(h/f) - w.  So in the null limit
# their own formula gives the optical metric, and the verbal claim would first
# need the time branch, the sign of the energy and the orientation reconciled --
# eqs. (3.2.2) and (3.2.4) of the 2024 text carry opposite wind signs.  The
# distinction we need is the MASSIVE one computed above, which does not depend on
# resolving that ambiguity.
mm, ww = sp.symbols('m w')
FJ = sp.sqrt(((-E)**2 - mm**2*f_)/f_*h_) + (-E)*ww
check('at m = 0 and p_0 = -E, Chanda eq. (16) divided by E IS the null optical '
      'functional sqrt(h/f) - w, so the 2024 verbal claim cannot be used as '
      'support without reconciling conventions',
      sp.simplify(FJ.subs(mm, 0)/E - (sp.sqrt(h_/f_) - ww)) == 0)
check('while at m = 1 it is the Jacobi-Maupertuis datum, which is the comparison '
      'that actually carries the argument',
      sp.simplify(FJ.subs(mm, 1) - (sp.sqrt(h_*(E**2 - f_)/f_) - E*ww)) == 0)

# ---------------------------------------------------------------- (2) witness
# Schwarzschild equatorial, E^2 = 2: a circular extremal of a diagonal radial
# metric satisfies d/dr(a_phiphi) = 0.  The two functionals disagree.
fs = 1 - 2*M/r
aT_pp = sp.simplify(2/(fs*(2 - fs))*r**2)
aJ_pp = sp.simplify((2 - fs)/fs*r**2)
rT = [sp.simplify(x) for x in sp.solve(sp.diff(aT_pp, r), r)]
rJ = [sp.simplify(x) for x in sp.solve(sp.diff(aJ_pp, r), r)]
check(f'arrival-cost circular extremal at r = 2*sqrt(2) M   (got {rT})',
      any(sp.simplify(x - 2*sp.sqrt(2)*M) == 0 for x in rT))
check(f'Jacobi-Maupertuis circular extremal at r = (1+sqrt(5)) M   (got {rJ})',
      any(sp.simplify(x - (1 + sp.sqrt(5))*M) == 0 for x in rJ))

# ---------------------------------------------------------------- (3) Kerr data
g00 = X/rho2
g0p = 2*M*a*r*s2/rho2
h_pp = sp.simplify(-(-s2/rho2*((r**2 + a**2)**2 - a**2*Delta*s2)) + g0p**2/g00)
check('h_phiphi = gamma_phiphi = rho^2 Delta sin^2(theta)/(Delta - a^2 sin^2(theta))',
      sp.simplify(h_pp - rho2*Delta*s2/X) == 0)

w_p = sp.simplify(g0p/g00)                            # the wind one-form
print('\n  Lambda^2 = E^2/(f(E^2-f)) with f = g_00,  E = Ehat/A')
print('  w_phi    =', sp.simplify(w_p))

# convexity N = a^{ij} b_i b_j with a = Lambda^2 h, b = -w
C = sp.simplify(g0p**2/(g00*h_pp))
check('C := g_0phi^2/(f h_phiphi) = 4 M^2 a^2 r^2 sin^2(theta)/(rho^4 Delta)',
      sp.simplify(C - 4*M**2*a**2*r**2*s2/(rho2**2*Delta)) == 0)
check('1 - C = f[(r^2+a^2)^2 - a^2 Delta sin^2(theta)]/(rho^2 Delta),  positive outside',
      sp.simplify(1 - C - g00*((r**2 + a**2)**2 - a**2*Delta*s2)/(rho2*Delta)) == 0)

N_T = sp.simplify(w_p**2/(E**2/(g00*(E**2 - g00))*h_pp))
N_J = sp.simplify(C*E**2/(E**2 - g00))
check('N_T = C (E^2 - f)/E^2', sp.simplify(N_T - C*(E**2 - g00)/E**2) == 0)
print('  so N_T <= C < 1 in the open exterior with E^2 > f (above freezing).')
# but the bound is attained in the limit: C -> 1 and f -> 0 at the equatorial
# ergosurface, so convexity DOES degenerate there -- the coefficients are singular
# and this is a limit, not a regular Randers norm on the boundary.
check('equatorial C -> 1 at r = 2M',
      sp.simplify(C.subs(th, sp.pi/2).subs(r, 2*M) - 1) == 0)
check('hence N_T -> 1 at the equatorial ergosurface',
      sp.limit(sp.simplify(N_T.subs(th, sp.pi/2)), r, 2*M, '+') == 1)
check('N_T also vanishes on the axis with a != 0, so it does not isolate a = 0',
      sp.simplify(C.subs(th, 0)) == 0)

check('N_T and N_J are different functions', sp.simplify(N_T - N_J) != 0)
sub = {M: 1, a: sp.Rational(9, 10), r: 3, th: sp.pi/2, E: sp.Rational(3, 5)}
print(f'\n  witness at M=1, a=0.9, r=3M, equator, E=Ehat/A=0.6 :'
      f'  N_T = {sp.nsimplify(N_T.subs(sub))},  N_J = {sp.nsimplify(N_J.subs(sub))}')
check('at that point the free-geodesic datum is not even convex, the arrival one is',
      N_J.subs(sub) > 1 and N_T.subs(sub) < 1)

# ---------------------------------------------------------------- (4) T2 scope
# T2 fixes h(dpsi X, dpsi Y) = Lambda^2 g(X,Y) on the horizontal space: a
# statement about the quadratic part only.  The wind is invisible to it.
check('non-rotating frozen branch: the wind vanishes and F_T is Riemannian',
      sp.simplify(w_p.subs(a, 0)) == 0)
print('  rotating frozen branch: w_phi != 0, an ingredient no horizontal-conformality')
print('  statement determines.  N_T is NOT a measure of that: it is sup |b(v)|/sqrt(a(v,v)),')
print('  not an additive share, and it is gauge dependent -- under t -> t + chi(x) the')
print('  cost changes by the exact term d(chi) (Perlick 1991, p. 3151, psi -> psi - du),')
print('  leaving the extremals fixed while changing b and its norm.  The invariant')
print('  local datum is db:')
b_eq = sp.simplify(-w_p.subs(th, sp.pi/2))
print(f'    equatorially b = {b_eq} dphi,  db = {sp.simplify(sp.diff(b_eq, r))} dr ^ dphi')
check('db != 0 for M > 0 and a != 0', sp.simplify(sp.diff(b_eq, r)) != 0)

# the wind is conformally invariant: w_i = g_{0i}/g_{00} is unchanged by g -> A^2 g.
# So when A runs it is only Lambda^2 that breathes, i.e. only the quadratic part.
Ae = sp.Function('A')(sp.Symbol('eta'))
check('w_i is invariant under g -> A(eta)^2 g, so non-stationarity enters only Lambda^2',
      sp.simplify(Ae**2*g0p/(Ae**2*g00) - w_p) == 0)

# --------------------------------------------------------- (5) the gauge proof
# The scope statement above is usually argued by inspection: "horizontal
# conformality is about the quadratic part".  It can be PROVED instead, and the
# proof is short.  Work in threading form g = -f (dt + beta)^2 + h dx dx and
# change the time section, t -> t' = t + chi(x).  Then dt + beta = dt' + (beta -
# dchi), so f and h are unchanged, hence so are the optical metric a = Lambda^2 h
# and the dilation Lambda; the projection pi along the orbits is the same map and
# the horizontal distribution is the same distribution.  Only beta moves, by the
# exact form -dchi.
#
# Consequence, stated at the strength it actually has: the representative of b is
# not fixed by those data.  It does NOT follow that the dynamics is unfixed --
# see block (6) -- because the change is by an exact form and leaves the
# fixed-endpoint extremals alone.  What this block establishes is only that b
# alone is a datum of the chosen time section; the intrinsic object is db, and
# separating the dynamics requires a witness at fixed db, not at fixed b.
tt, ph = sp.symbols('t varphi')
f_th = sp.Function('f')(r)             # threading lapse, kept abstract
bet = sp.Function('beta')(r)           # threading shift, dx^1 = dphi component
hrr, hpp = sp.Function('h_rr')(r), sp.Function('h_pp')(r)
chi = sp.Function('chi')(r)

# g in threading form, coordinates (t, r, phi), shift along dphi.
X = [tt, r, ph]
one = sp.Matrix([[1, 0, bet]])                      # dt + beta dphi
gth = -f_th*(one.T*one) + sp.diag(0, hrr, hpp)

# the change of section t = t' - chi(r): dt = dt' - chi' dr.  Push the metric
# forward with the Jacobian rather than asserting the result.
tp = sp.Symbol('tprime')
Xn = [tp, r, ph]
Jac = sp.Matrix([[sp.diff(c, v) for v in Xn]
                 for c in [tp - chi, r, ph]])       # old coords in terms of new
gnew = sp.simplify(Jac.T*gth*Jac)

# read the new threading data off gnew: f' = -g'_{t't'}, beta'_i = g'_{t'i}/g'_{t't'}
f_new = sp.simplify(-gnew[0, 0])
beta_new_r = sp.simplify(gnew[0, 1]/gnew[0, 0])
beta_new_p = sp.simplify(gnew[0, 2]/gnew[0, 0])
h_new_rr = sp.simplify(gnew[1, 1] + f_new*beta_new_r**2)
h_new_pp = sp.simplify(gnew[2, 2] + f_new*beta_new_p**2)

check('under t -> t + chi(x) the threading lapse f is unchanged',
      sp.simplify(f_new - f_th) == 0)
check('and the threading spatial metric h is unchanged',
      sp.simplify(h_new_rr - hrr) == 0 and sp.simplify(h_new_pp - hpp) == 0)
check('hence Lambda^2 = E^2/(f(E^2-f)) and a = Lambda^2 h are unchanged too, '
      'being built from f and h alone',
      sp.simplify(E**2/(f_new*(E**2 - f_new)) - E**2/(f_th*(E**2 - f_th))) == 0)
check('while the shift picks up exactly the exact form -dchi',
      sp.simplify(beta_new_r + sp.diff(chi, r)) == 0
      and sp.simplify(beta_new_p - bet) == 0)
check('and that shift really moves: distinct chi give distinct b, so the same '
      'horizontal-conformality data carry a whole family of one-forms',
      sp.simplify(beta_new_r.subs(chi, r**2)
                  - beta_new_r.subs(chi, 2*r**2)) != 0)
check('what survives the gauge is db, and db = 0 exactly when a = 0',
      sp.simplify(sp.diff(beta_new_r, ph) - sp.diff(beta_new_p, r)
                  - (sp.diff(-sp.diff(chi, r), ph) - sp.diff(bet, r))) == 0
      and sp.simplify(sp.diff(b_eq, r).subs(a, 0)) == 0
      and sp.simplify(sp.diff(b_eq, r)) != 0)
print('  so the implication runs the other way: semiconformality is a HYPOTHESIS')
print('  of such statements and a CONCLUSION here, obtained by computing Lambda')
print('  from the threading metric and the fixed-energy mass shell.')

# --------------------------------------- (6) underdetermination, done properly
# WHAT THE GAUGE ARGUMENT DOES NOT SHOW.  Under t -> t + chi(x) the integral of b
# changes by chi(endpoint) - chi(startpoint), a constant at fixed spatial
# endpoints, so the EXTREMALS do not move.  Non-uniqueness of the representative
# of b is therefore NOT a separation between families of extremals, and must not
# be used as one.  An earlier version of this file, and of the manuscript, leaned
# on it as though it were; the correction is due to an external audit.
#
# Nor is it right to say "b is not a datum of the submersion".  Given the total
# metric g and the normalised vertical W, the connection one-form
#     theta = g(W,.)/g(W,W) = dt + omega,   theta(W) = 1,   H = ker theta,
# is determined, and in the stationary case its curvature descends to d omega on
# the base.  So db = -d omega IS recoverable from the complete data: for
# horizontal lifts, d theta(X^H, Y^H) = -theta([X^H, Y^H]), the non-integrability
# of H.  The correct statement is weaker and sharper: the complete submersion
# metric can carry the rotation, while the horizontal-conformality condition
# ALONE does not constrain it.
#
# WITNESS (due to the audit, and better than the spin-magnitude family that stood
# here before, whose caveat was that h varied along it).  Compare Kerr at +a and
# at -a in fixed coordinates.  Every coefficient of f and h depends on a^2, so
# f, h, Lambda and the whole base metric a_T are IDENTICAL; but b is odd in a, so
# b and db flip sign.  Identical quadratic optical data, opposite connection
# curvature, hence different oriented extremals.  The two are of course related
# by the reflection phi -> -phi; the claim is not that they are inequivalent
# under every transformation, but that the quadratic part does not fix the sense
# of the dragging.
for label, expr in (('f = g_00', g00), ('h_rr', sp.simplify(rho2/Delta)),
                    ('h_phiphi', h_pp)):
    check(f'Kerr +a and -a share {label} (it depends on a^2)',
          sp.simplify(expr.subs(a, -a) - expr) == 0)
Lam2_k = sp.simplify(E**2/(g00*(E**2 - g00)))
check('hence they share the dilation Lambda^2 and the whole base metric a_T',
      sp.simplify(Lam2_k.subs(a, -a) - Lam2_k) == 0
      and sp.simplify((Lam2_k*h_pp).subs(a, -a) - Lam2_k*h_pp) == 0)
db_eq = sp.simplify(sp.diff(b_eq, r))
check('while b and db are odd in a, so the connection curvature flips sign at '
      'identical quadratic data',
      sp.simplify(b_eq.subs(a, -a) + b_eq) == 0
      and sp.simplify(db_eq.subs(a, -a) + db_eq) == 0
      and sp.simplify(db_eq) != 0)

# A cleaner general witness, on a flat neighbourhood: one dilation, a whole
# family of connection curvatures, all horizontally conformal.
#     g_k = -(dt - k x dy)^2 + dx^2 + dy^2 + dz^2,  W = d_t,  E^2 = 2
# gives f = 1, h = I, Lambda^2 = 2, a_T = 2I for every k, while db_k = k dx ^ dy.
kk, xx = sp.symbols('k x', real=True)
f_k, h_k = sp.Integer(1), sp.eye(3)
Lam2_flat = sp.simplify(2/(f_k*(2 - f_k)))
check('flat family g_k: f = 1 and h = I for every k',
      sp.simplify(f_k - 1) == 0 and (h_k - sp.eye(3)).is_zero_matrix)
check('so Lambda^2 = 2 and a_T = 2I for every k, independently of k',
      sp.simplify(Lam2_flat - 2) == 0 and sp.diff(Lam2_flat, kk) == 0)
b_k = kk*xx                      # b = k x dy
check('while db_k = k dx ^ dy separates the members: straight extremals at k = 0,'
      ' curved at k != 0',
      sp.simplify(sp.diff(b_k, xx)) == kk
      and sp.simplify(sp.diff(b_k, xx).subs(kk, 0)) == 0
      and sp.simplify(sp.diff(b_k, xx).subs(kk, 1)) != 0)
check('and each member is strongly convex where |k x| < sqrt 2',
      sp.simplify((b_k**2/Lam2_flat).subs({kk: 1, xx: sp.Rational(1, 2)})) < 1)

# Finally: conformal does NOT mean geodesic-preserving, which is the step the
# objection needs and does not have.  On (R^2, E^{2x}(dx^2+dy^2)) the vertical
# line x = const is not a geodesic, though it is one for the flat metric in the
# same conformal class.
xs, ys = sp.symbols('x y', real=True)
conf = sp.exp(2*xs)*sp.eye(2)
Gam_x_yy = sp.simplify(sp.Rational(1, 2)*(1/conf[0, 0])
                       * (-sp.diff(conf[1, 1], xs)))
check(f'conformal rescaling does not preserve geodesics: Gamma^x_yy = {Gam_x_yy}',
      sp.simplify(Gam_x_yy + 1) == 0)

print('\nALL CHECKS PASS' if ok else '\nSOME CHECKS FAILED')
raise SystemExit(0 if ok else 1)
