# -*- coding: utf-8 -*-
"""
The drift clock on the Thakurta--Kerr tau-branch.

The conformal factor runs with the conformal time, A = A(eta), so the first-order
drift of the frozen charges along an orbit is weighted by the conformal time
ELAPSED ALONG THAT ORBIT, eta(r), whichever clock the cost functional measures.
On the t-branch eta and t coincide at A = 1.  On the tau-branch they do not: the
proper time tau(r) is the clock the functional measures (it sets the timing), but
it is not the slow variable.  An earlier version weighted the tau-branch drift by
tau(r); this script establishes the correct weight and its closed form.

  (1) symbolic: the rail relation d tau/d eta = (f + (2Ma/r) dphi/deta)/E at A = 1
      gives   d eta/dr = [E r^3 - 2MaJ r D_E/Delta] / sqrt(S),
      and it reduces to E r^3/sqrt(S) = E U_3' at a = 0 (Schwarzschild coordinate
      time, the bulk part of the Vaidya advanced-time clock of Paper I);
  (2) independent: d eta/dr = 1/H_{p_r} along the Hamiltonian flow of H_tau
      (eta-parametrised: velocities d x/d eta, running cost d tau/d eta);
  (3) partial fractions: second-kind vector b = (-2MaJ(E^2-1), 0, 0, E, 0) and a
      horizon third-kind part -2MaJ [2ME^2 r - (E^2-1)a^2]/(Delta sqrt S);
      coefficient of 1/r at infinity E M (2E^2-3)/(E^2-1)^{3/2};
  (4) assembly: int F_lambda * eta dr = [W_kj reorganisation of the second-kind
      part] + [horizon part], for the energy source and for the full Euler source
      (E d_E + J d_J), against the direct double integral;
  (5) size of the error the proper-time weight made.

Exit status 0 iff every check passes.
"""
import sys
import numpy as np
import sympy as sp
from scipy.integrate import quad, solve_ivp
from scipy.optimize import brentq

FAIL = []


def check(label, ok, detail=""):
    print(f"  [{'PASS' if ok else 'FAIL'}] {label}" + (f"   {detail}" if detail else ""))
    if not ok:
        FAIL.append(label)


r, M, a, E, J = sp.symbols('r M a E J', positive=True)
f = 1 - 2*M/r
Dl = r**2 - 2*M*r + a**2
DE = (E**2 - 1)*r + 2*M
S = r*(r - 2*M)*DE*(r*Dl - J**2*DE)
Ktau = J*r*(r - 2*M)*DE/Dl                     # dphi/dr = Ktau/sqrt(S)
dtau_num = r**2*(r - 2*M)                      # dtau/dr = dtau_num/sqrt(S)

# ---------------------------------------------------------------- (1) symbolic
print("(1) symbolic form of the drift clock")
deta_num = sp.simplify((E/f)*dtau_num - (2*M*a/(r*f))*Ktau)
claim = E*r**3 - 2*M*a*J*r*DE/Dl
check("d eta/dr * sqrt(S) = E r^3 - 2MaJ r D_E/Delta",
      sp.simplify(deta_num - claim) == 0)
check("a = 0: d eta/dr = E r^3/sqrt(S)", sp.simplify(claim.subs(a, 0) - E*r**3) == 0)
check("it is NOT the proper-time density (differs unless E = f and a = 0)",
      sp.simplify(claim - dtau_num) != 0)
q, rem = sp.div(sp.Poly(sp.expand(r*DE), r), sp.Poly(sp.expand(Dl), r))
check("r D_E/Delta = (E^2-1) + [2ME^2 r - (E^2-1)a^2]/Delta",
      sp.simplify(q.as_expr() - (E**2 - 1)) == 0 and
      sp.simplify(rem.as_expr() - (2*M*E**2*r - (E**2 - 1)*a**2)) == 0)
x = sp.symbols('x', positive=True)
lead = sp.series((claim/sp.sqrt(S)).subs(r, 1/x), x, 0, 2).removeO()
c1 = sp.simplify(lead.coeff(x, 1))
check("coefficient of 1/r at infinity = E M (2E^2-3)/(E^2-1)^(3/2)",
      sp.simplify(c1 - E*M*(2*E**2 - 3)/(E**2 - 1)**sp.Rational(3, 2)) == 0)

# ------------------------------------------------ (2) against the Hamiltonian flow
print("\n(2) against the eta-parametrised Hamiltonian flow of H_tau")
rr, pr, Es, Js = sp.symbols('rr pr Es Js')
Mn, an = 1.0, 0.9


def htau(Mv, av):
    ff = 1 - 2*Mv/rr; Dd = rr**2 - 2*Mv*rr + av**2; bb = 2*Mv*av/rr; vv = 1 - ff/Es**2
    P = rr**2 + av**2 + 2*Mv*av**2/rr; pt = Js - bb/Es; Pb = P + bb**2/Es**2
    return pt*(bb*vv/Pb) + sp.sqrt(Dd*vv/Pb)*sp.sqrt((Dd/rr**2)*pr**2 + pt**2/Pb) - ff/Es


H = htau(Mn, an)
L = lambda X: sp.lambdify((rr, pr, Es, Js), X, 'numpy')
Hn, Hp, Hr, HJ = L(H), L(sp.diff(H, pr)), L(sp.diff(H, rr)), L(sp.diff(H, Js))
for (Ev, Jv, r0) in [(1.4, 2.5, 12.0), (1.2, 2.5, 10.0), (1.6, 4.0, 14.0)]:
    pg = np.linspace(-80, 80, 4001); hv = Hn(r0, pg, Ev, Jv)
    rts = [brentq(lambda p: Hn(r0, p, Ev, Jv), pg[i], pg[i+1]) for i in range(len(pg)-1)
           if np.isfinite(hv[i]*hv[i+1]) and hv[i]*hv[i+1] < 0]
    p0 = min(p for p in rts if Hp(r0, p, Ev, Jv) < 0)
    ev = lambda l, y: y[1]; ev.terminal = True; ev.direction = 1
    sol = solve_ivp(lambda l, y: [Hp(y[0], y[1], Ev, Jv), -Hr(y[0], y[1], Ev, Jv),
                                  HJ(y[0], y[1], Ev, Jv)],
                    [0, 400], [r0, p0, 0.0], rtol=1e-12, atol=1e-14, events=ev,
                    max_step=0.01, dense_output=True)
    num = sp.lambdify(r, (claim/sp.sqrt(S)).subs({M: Mn, a: an, E: Ev, J: Jv}), 'numpy')
    phinum = sp.lambdify(r, (Ktau/sp.sqrt(S)).subs({M: Mn, a: an, E: Ev, J: Jv}), 'numpy')
    taunum = sp.lambdify(r, (dtau_num/sp.sqrt(S)).subs({M: Mn, a: an, E: Ev, J: Jv}), 'numpy')
    rmin = sol.y[0].min(); worst = 0.0; worstphi = 0.0; besttau = np.inf
    for rv in np.linspace(rmin + 0.3, r0 - 0.3, 9):
        i = np.argmin(abs(sol.y[0] - rv)); rv = sol.y[0][i]; pv = sol.y[1][i]
        worst = max(worst, abs(abs(1/Hp(rv, pv, Ev, Jv)) - num(rv))/num(rv))
        besttau = min(besttau, abs(abs(1/Hp(rv, pv, Ev, Jv)) - taunum(rv))/num(rv))
        worstphi = max(worstphi, abs(abs(HJ(rv, pv, Ev, Jv)/Hp(rv, pv, Ev, Jv)) - phinum(rv))
                       / phinum(rv))
    check(f"E={Ev}, J={Jv}: |1/H_pr| = closed form", worst < 1e-8, f"rel {worst:.1e}")
    check(f"E={Ev}, J={Jv}: dphi/dr = Ktau/sqrt(S) (sanity)", worstphi < 1e-8,
          f"rel {worstphi:.1e}")
    # falsifiability: the proper-time density must FAIL the same comparison
    check(f"E={Ev}, J={Jv}: proper-time density is rejected by the same test",
          besttau > 1e-2, f"smallest rel gap {besttau:.2f}")

# --------------------------------------------------------- (3)+(4) the assembly
print("\n(3)-(4) weight-two assembly with the drift clock")


def assemble(sub, source, r0=12.0):
    """source: 'E' (E d_E) or 'Euler' (E d_E + J d_J).  Returns direct and assembled."""
    Ss = sp.expand(S.subs(sub)); Kf = Ktau
    if source == 'E':
        N = sp.expand(sp.simplify(E*sp.diff(Kf/sp.sqrt(S), E)*S**sp.Rational(3, 2)))
    else:
        N = sp.expand(sp.simplify((E*sp.diff(Kf/sp.sqrt(S), E) + J*sp.diff(Kf/sp.sqrt(S), J))
                                  * S**sp.Rational(3, 2)))
    Ns = sp.expand(N.subs(sub))
    ai = [sp.Symbol(f'a{i}') for i in range(6)]; ck = [sp.Symbol(f'c{i}') for i in range(5)]
    A = sum(ai[i]*r**i for i in range(6)); Mp = sum(ck[i]*r**i for i in range(5))
    ident = sp.expand(2*sp.diff(A, r)*Ss - A*sp.diff(Ss, r) + 2*Ss*Mp - 2*Ns)
    sol = sp.solve(sp.Poly(ident, r).all_coeffs(), ai + ck, dict=True)[0]
    exact = sp.Poly(sp.expand(ident.subs(sol)), r).is_zero
    c = [float(sol.get(ck[i], 0)) for i in range(5)]
    Acal = sp.lambdify(r, sum(sol.get(ai[i], 0)*r**i for i in range(6)), 'numpy')
    Mv, av, Ev, Jv = [float(sub[s]) for s in (M, a, E, J)]
    b = [-2*Mv*av*Jv*(Ev**2 - 1), 0.0, 0.0, Ev, 0.0]
    RD = lambda t: -2*Mv*av*Jv*(2*Mv*Ev**2*t - (Ev**2 - 1)*av**2)
    Dn = lambda t: t**2 - 2*Mv*t + av**2
    Sn = sp.lambdify(r, Ss, 'numpy'); Nn = sp.lambdify(r, Ns, 'numpy')
    sq = lambda t: np.sqrt(Sn(t))
    rs = np.linspace(2.2, r0, 20000); vs = Sn(rs); idx = np.where(np.diff(np.sign(vs)))[0]
    xf = max(brentq(Sn, rs[i], rs[i+1]) for i in idx) + 0.4
    Fl = lambda t: Nn(t)/Sn(t)**1.5
    U = lambda t, k: quad(lambda s: s**k/sq(s), r0, t, limit=200)[0]
    eta2 = lambda t: sum(b[j]*U(t, j) for j in range(5) if b[j] != 0.0)
    eta3 = lambda t: quad(lambda s: RD(s)/(Dn(s)*sq(s)), r0, t, limit=200)[0]
    W = lambda t, k, j: quad(lambda s: (U(s, k)*s**j - U(s, j)*s**k)/sq(s), r0, t, limit=150)[0]
    Ij = lambda t, j: quad(lambda s: Acal(s)*s**j/Sn(s), r0, t, limit=200)[0]
    # the clock itself, closed form against its direct quadrature
    direct_clock = quad(lambda s: (Ev*s**3 - 2*Mv*av*Jv*s*((Ev**2-1)*s + 2*Mv)/Dn(s))/sq(s),
                        r0, xf, limit=200)[0]
    clock_ok = abs(direct_clock - (eta2(xf) + eta3(xf)))
    Uv = [U(xf, k) for k in range(5)]
    Qsum = sum((c[k]*b[j] - c[j]*b[k])*W(xf, k, j) for k in range(5) for j in range(k+1, 5)
               if abs(c[k]*b[j] - c[j]*b[k]) > 1e-15)
    I2 = (-0.5*Qsum + eta2(xf)*(Acal(xf)/sq(xf) + 0.5*sum(c[k]*Uv[k] for k in range(5)))
          - sum(b[j]*Ij(xf, j) for j in range(5) if b[j] != 0.0))
    I3 = quad(lambda s: Fl(s)*eta3(s), r0, xf, limit=200)[0]
    Idir = quad(lambda s: Fl(s)*(eta2(s) + eta3(s)), r0, xf, limit=200)[0]
    Ptau = quad(lambda s: Fl(s)*quad(lambda t: t**2*(t - 2*Mv)/sq(t), r0, s, limit=200)[0],
                r0, xf, limit=200)[0]
    return exact, clock_ok, Idir, I2 + I3, I3, Ptau


for sub in [{M: 1, a: sp.Rational(9, 10), E: sp.Rational(7, 5), J: sp.Rational(5, 2)},
            {M: 1, a: sp.Rational(9, 10), E: sp.Rational(6, 5), J: sp.Rational(5, 2)}]:
    for src in ('E', 'Euler'):
        exact, cl, Idir, Iasm, I3, Ptau = assemble(sub, src)
        tag = f"E={sub[E]}, source {src}"
        check(f"{tag}: reduction identity exact", exact)
        check(f"{tag}: clock = second-kind + horizon part", cl < 1e-9, f"diff = {cl:.2e}")
        check(f"{tag}: assembled = direct", abs(Idir - Iasm) < 1e-8*max(1, abs(Idir)),
              f"diff = {abs(Idir - Iasm):.2e}")
        print(f"      direct {Idir:+.10f}   horizon part {I3:+.6f}   "
              f"proper-time weight would give {Ptau:+.6f} "
              f"({100*abs(Ptau - Idir)/abs(Idir):.0f}% off)")

print()
if FAIL:
    print(f"{len(FAIL)} check(s) FAILED:", *FAIL, sep="\n  ")
    sys.exit(1)
print("all checks passed")
