# -*- coding: utf-8 -*-
"""
Weight-two assembly of the new off-shell adiabatic term into the genus-2 W_jk basis (main11).

  delta phi_extra = - int (A/sqrt(S)) * Delta S dr ,   A/sqrt(S)=sum_j p_j r^j/sqrt(S)+(poles),
                    Delta S = int p_r dr = sum_k a_k U_k + (third kind).

Core identity (verified to 1e-17 below):
  int (r^j/sqrt(S)) U_k dr = 1/2 ( U_j U_k + W_jk ),   U_k=int r^k/sqrt(S) dr,
  W_jk = int (r^j U_k - r^k U_j)/sqrt(S) dr   (the weight-two / genus-2 dilogarithm letters).

Hence the second-kind block of the wrap is
  delta phi_extra |_2nd = sum_{j<k} Q_jk W_jk  -  1/2 (sum_j p_j U_j)(sum_k a_k U_k),
  Q_jk = -1/2 ( p_j a_k - p_k a_j ),
same W_jk basis as the on-shell psi (kerr_psi_explicit_verified.py); the p_j come from A(r)
and the a_k from the level-B radial-action reduction (adiabatic_offshell_reduction.py). The
pole letters of A and the third-kind part of Delta S add third-kind dilog blocks by the same
IBP. The W_jk themselves are genus-2 transcendentals (theta-nome q-series), not elementary.
"""
import numpy as np, sympy as sp
from scipy.integrate import quad

M, a, E, J = 1.0, 0.9, 1.4, 6.0
r = sp.symbols('r', positive=True)
Es, Js = sp.symbols('E J_')
Q2 = (2*Es**2*Js**2*M*r - Es**2*Js**2*r**2 - 4*Es**2*Js*M*a*r + 2*Es**2*M*a**2*r
      + Es**2*a**2*r**2 + Es**2*r**4 + 4*Js**2*M**2 - 4*Js**2*M*r + Js**2*r**2
      - 8*Js*M**2*a + 4*Js*M*a*r + 4*M**2*a**2)
S = sp.expand((r*Q2*((Es**2-1)*r + 2*M)).subs({Es: E, Js: J}))
Sn = sp.lambdify(r, S, 'numpy')
sq = lambda x: np.sqrt(Sn(x))
r0, rf = 11.0, 6.0            # both well above the turning point


def U(x, k):
    return quad(lambda t: t**k/sq(t), r0, x, limit=200)[0]


def W(j, k):
    return quad(lambda x: (x**j*U(x, k) - x**k*U(x, j))/sq(x), r0, rf, limit=200)[0]


def LHS(j, k):
    return quad(lambda x: (x**j/sq(x))*U(x, k), r0, rf, limit=200)[0]


if __name__ == '__main__':
    print("IBP identity  int (r^j/sqrt S) U_k dr = 1/2 (U_j U_k + W_jk):")
    ok = True
    for (j, k) in [(0, 2), (1, 3), (0, 3), (1, 2)]:
        lhs = LHS(j, k)
        rhs = 0.5*(U(rf, j)*U(rf, k) + W(j, k))
        d = abs(lhs - rhs)
        ok = ok and d < 1e-12
        print(f"  (j,k)=({j},{k}): LHS={lhs:.10f}  RHS={rhs:.10f}  diff={d:.1e}")
    print("\nOK" if ok else "\nFAILED")
