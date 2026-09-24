# -*- coding: utf-8 -*-
"""Riemann period matrix of the generic t-branch curve y^2 = R6 (Paper II, eq. t-genus2).

M = 1, a = 9/10, E = 6/5, J = 3: the genus-2 orbit of kerr_jpt_genus2_figure.py, which
penetrates the ergosphere and reproduces the Hamiltonian flow to ~2e-13. Checks the two
Riemann relations: tau = tau^T and Im(tau) positive definite. The eigenvalues of Im(tau)
depend on the homology basis RiemannSurface chooses and are printed for the record only.

Run with:  sage -python kerr_tbranch_genus2_periods.sage.py
"""
import sys
from sage.all import QQ, PolynomialRing, lcm, expand, RealField, matrix
from sage.schemes.riemann_surfaces.riemann_surface import RiemannSurface

Rx = PolynomialRing(QQ, 'x'); x = Rx.gen()
M, a, E, J = QQ(1), QQ(9)/10, QQ(6)/5, QQ(3)
Q2 = (2*E**2*J**2*M*x - E**2*J**2*x**2 - 4*E**2*J*M*a*x + 2*E**2*M*a**2*x + E**2*a**2*x**2
      + E**2*x**4 + 4*J**2*M**2 - 4*J**2*M*x + J**2*x**2 - 8*J*M**2*a + 4*J*M*a*x + 4*M**2*a**2)
R6 = expand(x*((E**2 - 1)*x + 2*M)*Q2)
L = lcm([c.denominator() for c in R6.coefficients()])
Rxy = PolynomialRing(QQ, ['X', 'Y']); X, Y = Rxy.gens()
S = RiemannSurface(Y**2 - expand(L*R6)(X), prec=40)
tau = S.riemann_matrix()
sym = float((tau - tau.transpose()).norm())
ImT = matrix(RealField(40), 2, 2, [complex(tau[i][j]).imag for i in range(2) for j in range(2)])
eig = [float(e) for e in ImT.eigenvalues()]
print(f"genus = {S.genus}")
print(f"|tau - tau^T| = {sym:.2e}")
print(f"eigenvalues of Im(tau) = {[round(e, 3) for e in eig]}  (basis-dependent)")
ok = S.genus == 2 and sym < 1e-10 and min(eig) > 0
print("ALL CHECKS PASSED" if ok else "FAILED")
sys.exit(0 if ok else 1)
