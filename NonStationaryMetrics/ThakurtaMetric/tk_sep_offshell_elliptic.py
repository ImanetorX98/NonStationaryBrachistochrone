# -*- coding: utf-8 -*-
"""
Elliptic collapse of the OFF-SHELL term at a degeneration of the t-branch curve
(Paper II, "Separatrix limit" paragraph and Remark rmk:collapse-scope).

At a double root r_d of the on-shell sextic S = r[(E^2-1)r+2M] Q2 one has
S = (r-r_d)^2 Q4, and the paper states three things, all checked here:

  (1) the curve degenerates:  S - (r-r_d)^2 Q4 == 0  (exact polynomial division);
  (2) the spurious poles cancel:  sum_k a_k(J_c) r_d^k + R/deng |_{r_d} = [S/deng](r_d) = 0,
      with a_k, R the radial-action reduction of eq. (ak) (adiabatic_offshell_coefficients.py);
  (3) each hyperelliptic letter collapses to the elliptic one,
          U_k = int r^k/sqrt(S) dr = V_k + r_d^k Pi_rd,
          V_k = int q_k(r)/sqrt(Q4) dr,  q_k = (r^k - r_d^k)/(r - r_d),
          Pi_rd = int dr/[(r-r_d) sqrt(Q4)],
      and so does the weight-two block W_jk = int (r^j U_k - r^k U_j)/sqrt(S) dr,
      which is recomputed from the collapsed letters and compared with the direct one.

The identities are algebraic and hold at EVERY double root; the numbers are spot checks.
They are run at two degenerations of the t-branch at M=1, a=9/10, E=6/5:

  J_c^- = -8.0535  r_d = +3.514   exterior retrograde separatrix, the physical one (r_d > 2M);
  J     = +2.9364  r_d = +1.512   prograde root inside r_e = 2M: a type-(A) object, NOT on the
                                  physical arc (an earlier comment here called it physical on the
                                  seed criterion r > r_+; the selector criterion is r > 2M).

The separatrix TRACKING of the dilation pole (the non-uniformity at r_d) is a separate, open
question and is not addressed here; see tk_sep_offshell_divergence.py.
"""
import sys
import mpmath as mp
import sympy as sp

mp.mp.dps = 30
FAIL = []


def check(name, ok, detail=""):
    print(f"  [{'PASS' if ok else 'FAIL'}] {name}" + (f"   {detail}" if detail else ""))
    if not ok:
        FAIL.append(name)


r, J = sp.symbols('r J')
M, a, E = sp.Integer(1), sp.Rational(9, 10), sp.Rational(6, 5)
Dl = r**2 - 2*M*r + a**2
DE = (E**2 - 1)*r + 2*M
Q2 = (2*E**2*J**2*M*r - E**2*J**2*r**2 - 4*E**2*J*M*a*r + 2*E**2*M*a**2*r + E**2*a**2*r**2
      + E**2*r**4 + 4*J**2*M**2 - 4*J**2*M*r + J**2*r**2 - 8*J*M**2*a + 4*J*M*a*r + 4*M**2*a**2)
S = sp.expand(r*DE*Q2)
deng = sp.expand(Dl*DE)
Qd, R = sp.div(sp.Poly(S, r), sp.Poly(deng, r))
a_k = Qd.all_coeffs()[::-1]                       # a_0..a_3, polynomials in J

# eq. (ak) of the paper, printed form
ak_printed = [2*M*(4*E**2*M**2 - J**2 + 2*(1 - E**2)*J*a),
              4*E**2*M**2 - (E**2 - 1)*J**2, 2*E**2*M, E**2]
print("eq. (ak) against the polynomial division")
check("a_0..a_3 as printed in eq. (ak)",
      all(sp.expand(x - y) == 0 for x, y in zip(a_k, ak_printed)))

# double roots: solve Q2 = dQ2/dr = 0 by Newton at 30 digits from the two known roots
fQ2 = sp.lambdify((J, r), Q2, 'mpmath')
fQ2r = sp.lambdify((J, r), sp.diff(Q2, r), 'mpmath')
targets = {"exterior retrograde, r_d>2M (physical)": (-8.0535, 3.514),
           "prograde, r_+<r_d<2M (type A)": (2.9364, 1.512)}

for label, (jg, rg) in targets.items():
    sol = mp.findroot(lambda x, y: [fQ2(x, y), fQ2r(x, y)], (mp.mpf(jg), mp.mpf(rg)))
    jv, rd = sp.Float(str(sol[0]), 30), sp.Float(str(sol[1]), 30)
    jf = float(jv)
    print(f"\n=== {label}:  J_c = {jf:+.10f}   r_d = {float(rd):+.10f}")
    Sj = sp.Poly(sp.expand(S.subs(J, jv)), r)
    Q4, rem = sp.div(Sj, sp.Poly((r - rd)**2, r))
    check("(1) S = (r-r_d)^2 Q4", max(abs(c) for c in rem.all_coeffs()) < 1e-25,
          f"remainder {float(max(abs(c) for c in rem.all_coeffs())):.1e}")
    # (2) residue cancellation at r_d
    val = sum(a_k[k].subs(J, jv)*rd**k for k in range(4)) + (R.as_expr()/deng).subs({J: jv, r: rd})
    Sd = (S/deng).subs({J: jv, r: rd})
    check("(2) sum a_k r_d^k + R/deng |_{r_d} = [S/deng](r_d) = 0",
          abs(val) < 1e-25 and abs(val - Sd) < 1e-25, f"|value| {float(abs(val)):.1e}")
    # (3) letter collapse on a compact sub-arc away from r_d
    fS = sp.lambdify(r, Sj.as_expr(), 'mpmath')
    fQ = sp.lambdify(r, Q4.as_expr(), 'mpmath')
    rdm = mp.mpf(str(sp.N(rd, 40)))
    r0, r1 = mp.mpf(6), mp.mpf(11)
    assert fS(r0) > 0 and fQ(r0) > 0 and fS(r1) > 0
    U = lambda k, x: mp.quad(lambda t: t**k/mp.sqrt(fS(t)), [r0, x])
    q_k = lambda k, t: sum(t**i*rdm**(k - 1 - i) for i in range(k))
    V = lambda k, x: mp.quad(lambda t: q_k(k, t)/mp.sqrt(fQ(t)), [r0, x])
    Pi = lambda x: mp.quad(lambda t: 1/((t - rdm)*mp.sqrt(fQ(t))), [r0, x])
    worst = max(abs(U(k, r1) - (V(k, r1) + rdm**k*Pi(r1)))/abs(U(k, r1)) for k in range(4))
    check("(3) U_k = V_k + r_d^k Pi_rd, k=0..3", worst < mp.mpf(10)**-20, f"max rel {float(worst):.1e}")
    Uc = lambda k, x: V(k, x) + rdm**k*Pi(x)
    for (j, k) in [(0, 1), (0, 2), (1, 3)]:
        Wd = mp.quad(lambda x: (x**j*U(k, x) - x**k*U(j, x))/mp.sqrt(fS(x)), [r0, r1])
        Wc = mp.quad(lambda x: (x**j*Uc(k, x) - x**k*Uc(j, x))/((x - rdm)*mp.sqrt(fQ(x))), [r0, r1])
        check(f"(3) W_{j}{k} from the collapsed letters", abs(Wd - Wc) <= mp.mpf(10)**-20*abs(Wd),
              f"rel {float(abs(Wd - Wc)/abs(Wd)):.1e}")
    # falsifiability: the collapse must fail if r_d is displaced
    bad = mp.quad(lambda t: 1/((t - rdm - mp.mpf('0.01'))*mp.sqrt(fQ(t))), [r0, r1])
    check("falsifiability: displacing r_d by 0.01 breaks (3)",
          abs(U(0, r1) - (V(0, r1) + bad)) > mp.mpf(10)**-6)

print()
if FAIL:
    print(f"{len(FAIL)} FAILED:", *FAIL, sep="\n  ")
    sys.exit(1)
print("ALL CHECKS PASSED")
