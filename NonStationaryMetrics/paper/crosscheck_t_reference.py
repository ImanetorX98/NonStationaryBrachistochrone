"""Python reference values for the Mathematica checks of the Thakurta-Kerr t branch (and of the Vaidya v branch).

Frozen Kerr (A = 1), M = 1, a = 9/10, E = 6/5, endpoint r0 = 20.  At each of the two double roots of the
t-branch quartic Q2 (prograde J ~ +19.089, r_d ~ -6.62; retrograde J ~ -18.671, r_d ~ -6.59) it computes the
first-order energy response

    delta phi_t(x) = int_{r0}^{x} (d_E F_t)(r) eta(r) dr,     eta(r) = int_{r0}^{r} rho_t/sqrt(R6) dr',

with R6 = r Q2 DE, F_t = K_t/sqrt(R6), K_t = r DE (J(r-2M)+2Ma)/Delta, and the t-branch drift clock
rho_t = (E^2 r^3 - (2Ma/r) K_t)/f of Paper II, eq. clock-t.  The derivative d_E F_t = N_t/R6^(3/2) uses the
t-branch numerator N_t = (d_E K_t) R6 - K_t (d_E R6)/2, derived symbolically at fixed J.

Two independent routes, which must agree:
  (A) the initial-value problem eta' = rho_t/sqrt(R6), (delta phi)' = (d_E F_t) eta, integrated with DOP853;
  (B) nested adaptive quadrature in mpmath at 30 digits.
It also checks that the polynomial part of rho_t equals the P_3 printed in Paper II after eq. clock-t.
r_d < 0 at both roots: they are algebraic genus degenerations, not physical separatrices; the integrals are
nevertheless real on the exterior arc [18, 20] used here (R6 > 0 there, checked).

Vaidya v branch (M = 1, E = 7/5, r0 = 12), at the double root J_c of the spectral curve S (r_d ~ -3.364 < 0):
    delta phi_v(11) = int_{r0}^{11} (d_m F) v dr,   F = J D_E/sqrt(S),   v = E U_3 + (r - r0) + 2m ln((r-2m)/(r0-2m)),
with U_3 = int_{r0}^{r} t^3/sqrt(S) dt and d_m at fixed J, by the same two routes.

The values printed here are the references quoted in crosscheck_tkt_flat.wl and crosscheck_vt_numeric.wl.

Run: python3 crosscheck_t_reference.py   (exit status 1 if a check fails)
"""
import sys
import mpmath as mp
import numpy as np
import sympy as sp
from scipy.integrate import solve_ivp

mp.mp.dps = 30
r, E, J = sp.symbols('r E J')
M, a, Ev, r0 = 1, sp.Rational(9, 10), sp.Rational(6, 5), 20
DE = (E**2 - 1)*r + 2*M
Dl = r**2 - 2*M*r + a**2
f = 1 - 2*M/r
Q2 = (2*E**2*J**2*M*r - E**2*J**2*r**2 - 4*E**2*J*M*a*r + 2*E**2*M*a**2*r + E**2*a**2*r**2 + E**2*r**4
      + 4*J**2*M**2 - 4*J**2*M*r + J**2*r**2 - 8*J*M**2*a + 4*J*M*a*r + 4*M**2*a**2)
R6 = r*Q2*DE
K = r*DE*(J*(r - 2*M) + 2*M*a)/Dl
Nt = sp.together(sp.diff(K, E)*R6 - K*sp.diff(R6, E)/2)
rho = sp.together((E**2*r**3 - (2*M*a/r)*K)/f)

fails = 0
def check(label, ok, detail=''):
    global fails
    fails += not ok
    print(f"  {'PASS' if ok else 'FAIL'}  {label}{('  ' + detail) if detail else ''}")

# polynomial part of rho_t against the P_3 printed in Paper II (symbolic in E, J)
num, den = sp.fraction(sp.cancel(rho))
quo, rem = sp.div(sp.Poly(num, r), sp.Poly(den, r))
P3_paper = E**2*r**3 + 2*E**2*M*r**2 + 4*E**2*M**2*r + 8*E**2*M**3 - 2*(E**2 - 1)*J*M*a
check("denominator of rho_t is a constant multiple of Delta (the ergosphere pole r = 2M is removable)",
      sp.simplify(den/Dl).is_constant())
check("polynomial part of rho_t equals P_3 of Paper II (after eq. clock-t)",
      sp.simplify(quo.as_expr() - P3_paper) == 0)
check("the horizon numerator R_Delta has degree one", rem.degree() == 1)

print("Thakurta-Kerr t branch, M=1, a=9/10, E=6/5, r0=20: direct delta phi_t at the two double roots of Q2")
Q2v = Q2.subs(E, Ev)
Q2f = sp.lambdify((r, J), Q2v, 'mpmath'); dQ2f = sp.lambdify((r, J), sp.diff(Q2v, r), 'mpmath')
for tag, guess in (('prograde', (-6.62, 19.089)), ('retrograde', (-6.588, -18.671))):
    rd, Jc = mp.findroot(lambda x, y: (Q2f(x, y), dQ2f(x, y)), guess)
    sub = {E: Ev, J: sp.Float(str(Jc), 40)}
    R6n = sp.lambdify(r, R6.subs(sub), 'mpmath')
    rhon = sp.lambdify(r, rho.subs(sub), 'mpmath')
    dEF = sp.lambdify(r, (Nt.subs(sub)), 'mpmath')
    dEFn = lambda x: dEF(x)/R6n(x)**mp.mpf(1.5)
    check(f"{tag}: R6 > 0 on [18, 20]", all(R6n(mp.mpf(x)) > 0 for x in np.linspace(18, 20, 41)))
    # (B) nested quadrature, 30 digits
    eta = lambda x: mp.quad(lambda t: rhon(t)/mp.sqrt(R6n(t)), [r0, x])
    dphiB = {x: mp.quad(lambda t: dEFn(t)*eta(t), [r0, x]) for x in (19, 18)}
    # (A) initial-value problem, double precision
    R6d = sp.lambdify(r, R6.subs(sub), 'numpy'); rhod = sp.lambdify(r, rho.subs(sub), 'numpy')
    Ntd = sp.lambdify(r, Nt.subs(sub), 'numpy')
    rhs = lambda x, y: [rhod(x)/np.sqrt(R6d(x)), Ntd(x)/R6d(x)**1.5*y[0]]
    sol = solve_ivp(rhs, [r0, 18], [0.0, 0.0], method='DOP853', rtol=1e-13, atol=1e-16, dense_output=True)
    dphiA = {x: sol.sol(x)[1] for x in (19, 18)}
    print(f"  {tag}: J = {mp.nstr(Jc, 15)}, r_d = {mp.nstr(rd, 12)}")
    for x in (19, 18):
        rel = abs(float(dphiB[x]) - dphiA[x])/abs(float(dphiB[x]))
        print(f"    delta phi_t({x}) = {mp.nstr(dphiB[x], 16)}   (quadrature)   {dphiA[x]:.15e}   (ODE)   rel. diff {rel:.1e}")
        check(f"{tag}, x = {x}: quadrature and ODE agree to 1e-10", rel < 1e-10)
        ok = mp.isfinite(dphiB[x]) and np.isfinite(dphiA[x])
        check(f"{tag}, x = {x}: finite", ok)

print("Vaidya v branch, M=1, E=7/5, r0=12: direct delta phi_v at the double root of S")
m, x = sp.symbols('m x')
Evv, r0v = sp.Rational(7, 5), 12
DEv = (Evv**2 - 1)*r + 2*m
Sv = r*(r - 2*m)*DEv*(r**2*(r - 2*m) - J**2*DEv)
Sv1 = sp.expand(Sv.subs(m, 1))
Cv = (r**2*(r - 2*m) - J**2*DEv).subs(m, 1)          # the factor of S that carries the double root
Cf = sp.lambdify((r, J), Cv, 'mpmath'); dCf = sp.lambdify((r, J), sp.diff(Cv, r), 'mpmath')
rdv, Jcv_mp = mp.findroot(lambda y, j: (Cf(y, j), dCf(y, j)), (-3.364, 7.0266))
Jcv = sp.Float(str(Jcv_mp), 40)
Fv = J*DEv/sp.sqrt(Sv)
dmF = sp.lambdify(r, sp.diff(Fv, m).subs({m: 1, J: Jcv}), 'mpmath')
Svn = sp.lambdify(r, Sv1.subs(J, Jcv), 'mpmath')
U3 = lambda y: mp.quad(lambda t: t**3/mp.sqrt(Svn(t)), [r0v, y])
vclock = lambda y: Evv*U3(y) + (y - r0v) + 2*mp.log((y - 2)/(r0v - 2))
dphivB = mp.quad(lambda t: dmF(t)*vclock(t), [r0v, 11])
dmFd = sp.lambdify(r, sp.diff(Fv, m).subs({m: 1, J: Jcv}), 'numpy'); Svd = sp.lambdify(r, Sv1.subs(J, Jcv), 'numpy')
rhsv = lambda y, z: [y**3/np.sqrt(Svd(y)), dmFd(y)*(float(Evv)*z[0] + (y - r0v) + 2*np.log((y - 2)/(r0v - 2)))]
solv = solve_ivp(rhsv, [r0v, 11], [0.0, 0.0], method='DOP853', rtol=1e-13, atol=1e-16)
dphivA = solv.y[1, -1]
relv = abs(float(dphivB) - dphivA)/abs(float(dphivB))
print(f"  J_c = {mp.nstr(Jcv_mp, 15)}, r_d = {mp.nstr(rdv, 12)}")
print(f"    delta phi_v(11) = {mp.nstr(dphivB, 16)}   (quadrature)   {dphivA:.15e}   (ODE)   rel. diff {relv:.1e}")
check("Vaidya, x = 11: quadrature and ODE agree to 1e-10", relv < 1e-10)
check("Vaidya, x = 11: finite", mp.isfinite(dphivB) and np.isfinite(dphivA))

print("ALL CHECKS PASSED" if fails == 0 else f"{fails} CHECK(S) FAILED")
sys.exit(1 if fails else 0)
