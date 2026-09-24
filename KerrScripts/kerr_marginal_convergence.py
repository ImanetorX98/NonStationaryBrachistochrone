# -*- coding: utf-8 -*-
"""Uniform convergence of the tau-branch shapes to the marginal one as |J| -> J_c^-
(Paper II, Proposition prop:classification, case (iii) discussion).

phi_J(r) = int_{2M}^{r} K/sqrt(S) dr, anchored at r_e = 2M, with
K = J r (r-2M) D_E/Delta and S = r(r-2M) D_E (r Delta - J^2 D_E), at M=1, a=0.9, E=1.2
(J_c = a/E = 0.75). Prints sup_{[2M,6M]} |phi_J - phi_Jc| and the cusp amplitude
K = (E J/a^2) sqrt(2M/(a^2 - E^2 J^2)) of eq. (cusp) for J_c - J = 1e-2 ... 1e-6.
The sup falls while K diverges: convergence in C^0, not in C^1.
"""
import sys
import mpmath as mp

mp.mp.dps = 30
M, a, E = mp.mpf(1), mp.mpf('0.9'), mp.mpf('1.2')
Jc = a/E
DE = lambda r: (E**2 - 1)*r + 2*M
Dl = lambda r: r**2 - 2*M*r + a**2


def dphi(r, J):
    return J*r*(r - 2*M)*DE(r)/(Dl(r)*mp.sqrt(r*(r - 2*M)*DE(r)*(r*Dl(r) - J**2*DE(r))))


e0 = 2*M + mp.mpf('1e-25')
brk = [e0, 2*M + mp.mpf('1e-8'), 2*M + mp.mpf('1e-4'), 2*M + mp.mpf('1e-2'), mp.mpf('2.1')]
phi = lambda r, J: mp.quad(lambda t: dphi(t, J), brk + [r])
grid = [2 + 4*mp.mpf(i)/40 for i in range(1, 41)]
ref = [phi(r, Jc) for r in grid]
sups, Ks = [], []
for k in range(2, 7):
    J = Jc - mp.mpf(10)**-k
    sups.append(max(abs(phi(r, J) - f) for r, f in zip(grid, ref)))
    Ks.append(E*J/a**2*mp.sqrt(2*M/(a**2 - E**2*J**2)))
    print(f"  Jc - J = 1e-{k}:  sup|phi_J - phi_Jc| on [2M,6M] = {float(sups[-1]):.2e}   K = {float(Ks[-1]):.1f}")
ok = all(sups[i + 1] < sups[i] for i in range(4)) and all(Ks[i + 1] > Ks[i] for i in range(4))
print("ALL CHECKS PASSED" if ok else "FAILED")
sys.exit(0 if ok else 1)
