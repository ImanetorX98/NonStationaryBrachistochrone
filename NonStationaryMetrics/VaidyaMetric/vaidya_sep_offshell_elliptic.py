# -*- coding: utf-8 -*-
"""
Vaidya: elliptic collapse of the weight-two letters at the genus degeneration
(Paper I, "Algebraic genus-degeneration limit" and the off-shell paragraph citing this file).

The frozen sextic S = r(r-2m) D_E [r^2(r-2m) - J^2 D_E], D_E = (E^2-1)r + 2m, acquires a double
root r_d at J = J_c, S = (r-r_d)^2 Q4, and the paper states that each hyperelliptic letter
collapses to the elliptic one,

    U_k = int r^k/sqrt(S) dr = V_k + r_d^k Pi_rd,
    V_k = int q_k(r)/sqrt(Q4) dr,  q_k = (r^k - r_d^k)/(r - r_d),
    Pi_rd = int dr/[(r-r_d) sqrt(Q4)],

and that the weight-two W_jk = int (r^j U_k - r^k U_j)/sqrt(S) dr follow, so the off-shell
degeneration form is elliptic. All of it is checked here at m=1, E=7/5, where
J_c = 5 sqrt(3011/3072 + 581 sqrt(249)/9216) = 7.0266... and r_d = -3.3637... < 0: the double
root is a branch point OFF the physical arc, an algebraic degeneration J_deg, not a
separatrix. For Vaidya D = Theta carries no dilation letter, so there is no off-shell
tracking counterterm at r_d (contrast ThakurtaMetric/tk_sep_offshell_elliptic.py).
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
m, E = sp.Integer(1), sp.Rational(7, 5)
DE = (E**2 - 1)*r + 2*m
g = r**2*(r - 2*m) - J**2*DE
S = sp.expand(r*(r - 2*m)*DE*g)

Jc_closed = 5*sp.sqrt(sp.Rational(3011, 3072) + 581*sp.sqrt(249)/9216)
fg = sp.lambdify((J, r), g, 'mpmath')
fgr = sp.lambdify((J, r), sp.diff(g, r), 'mpmath')
sol = mp.findroot(lambda x, y: [fg(x, y), fgr(x, y)], (mp.mpf('7.03'), mp.mpf('-3.36')))
jv, rd = sp.Float(str(sol[0]), 30), sp.Float(str(sol[1]), 30)
print(f"J_c = {float(jv):.10f}   r_d = {float(rd):+.10f}")
check("J_c agrees with the closed form printed in Paper I",
      abs(jv - sp.N(Jc_closed, 30)) < 1e-25, f"diff {float(abs(jv - sp.N(Jc_closed, 30))):.1e}")
check("r_d < 0: a branch point off the physical arc (J_deg, not a separatrix)", rd < 0)

Sj = sp.Poly(sp.expand(S.subs(J, jv)), r)
Q4, rem = sp.div(Sj, sp.Poly((r - rd)**2, r))
check("S = (r-r_d)^2 Q4", max(abs(c) for c in rem.all_coeffs()) < 1e-24,
      f"remainder {float(max(abs(c) for c in rem.all_coeffs())):.1e}")

fS = sp.lambdify(r, Sj.as_expr(), 'mpmath')
fQ = sp.lambdify(r, Q4.as_expr(), 'mpmath')
rdm = mp.mpf(str(sp.N(rd, 40)))
r0, r1 = mp.mpf('9.2'), mp.mpf(14)   # S>0 beyond the simple root r=8.727 of the bracket
assert fS(r0) > 0 and fQ(r0) > 0 and fS(r1) > 0
U = lambda k, x: mp.quad(lambda t: t**k/mp.sqrt(fS(t)), [r0, x])
q_k = lambda k, t: sum(t**i*rdm**(k - 1 - i) for i in range(k))
V = lambda k, x: mp.quad(lambda t: q_k(k, t)/mp.sqrt(fQ(t)), [r0, x])
Pi = lambda x: mp.quad(lambda t: 1/((t - rdm)*mp.sqrt(fQ(t))), [r0, x])
worst = max(abs(U(k, r1) - (V(k, r1) + rdm**k*Pi(r1)))/abs(U(k, r1)) for k in range(5))
check("U_k = V_k + r_d^k Pi_rd, k=0..4", worst < mp.mpf(10)**-20, f"max rel {float(worst):.1e}")
Uc = lambda k, x: V(k, x) + rdm**k*Pi(x)
for (j, k) in [(0, 1), (0, 2), (1, 3)]:
    Wd = mp.quad(lambda x: (x**j*U(k, x) - x**k*U(j, x))/mp.sqrt(fS(x)), [r0, r1])
    Wc = mp.quad(lambda x: (x**j*Uc(k, x) - x**k*Uc(j, x))/((x - rdm)*mp.sqrt(fQ(x))), [r0, r1])
    check(f"W_{j}{k} from the collapsed letters", abs(Wd - Wc) <= mp.mpf(10)**-20*abs(Wd),
          f"rel {float(abs(Wd - Wc)/abs(Wd)):.1e}")
bad = mp.quad(lambda t: 1/((t - rdm - mp.mpf('0.01'))*mp.sqrt(fQ(t))), [r0, r1])
check("falsifiability: displacing r_d by 0.01 breaks the collapse",
      abs(U(0, r1) - (V(0, r1) + bad)) > mp.mpf(10)**-6)

print()
if FAIL:
    print(f"{len(FAIL)} FAILED:", *FAIL, sep="\n  ")
    sys.exit(1)
print("ALL CHECKS PASSED")
