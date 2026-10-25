# -*- coding: utf-8 -*-
# Adiabatic Vaidya: reduction of dM F (the analogue of dE F for Thakurta-Kerr).
# Frozen = Schwarzschild (a=0), genus-2 curve S. dM F = N_M/S^(3/2) -> reduction
# dM F = d(A/sqrtS) + sum c_k^M r^k/sqrtS. Checked numerically to ~1e-15.
import sympy as sp

r, M, E, J = sp.symbols('r M E J', positive=True)
Emu = (E**2-1)*r + 2*M
S = sp.expand(r*(r-2*M)*Emu*(r**2*(r-2*M) - J**2*Emu))    # tau sextic, a=0
K = J*Emu                                                 # K_tau a=0

# N_M polynomial: dM F = N_M/S^(3/2), N_M = S dM K - K dM S / 2
N_M = sp.expand(S*sp.diff(K, M) - K*sp.diff(S, M)/2)
print('N_M degree', sp.degree(N_M, r), '(polynomial)')

# reduction: 2 N_M = 2S A' - A S' + 2S Mpoly, A deg5, Mpoly deg4
Ac = [sp.Symbol(f'A{i}') for i in range(6)]
ck = [sp.Symbol(f'c{i}') for i in range(5)]
Acal = sum(Ac[i]*r**i for i in range(6)); Mp = sum(ck[i]*r**i for i in range(5))
eq = sp.expand(2*N_M - (2*S*sp.diff(Acal, r) - Acal*sp.diff(S, r) + 2*S*Mp))
sol = sp.solve(sp.Poly(eq, r).all_coeffs(), Ac + ck, dict=True)[0]
print('reduction solved')

# numerical check
sub = {M: 1, E: sp.Rational(7, 5), J: sp.Rational(5, 2)}
dMF = sp.diff(K, M)/sp.sqrt(S) - K*sp.diff(S, M)/(2*S**sp.Rational(3, 2))
dMFn = sp.lambdify(r, dMF.subs(sub), 'numpy')
rhs = sp.diff(Acal.subs(sol)/sp.sqrt(S), r) + Mp.subs(sol)/sp.sqrt(S)
rhsn = sp.lambdify(r, rhs.subs(sub), 'numpy')
print('CHECK dM F == d(A/sqrtS) + M_poly/sqrtS:')
for rv in [11., 9., 7., 5.]:
    print(f'  r={rv}: diff = {abs(dMFn(rv)-rhsn(rv)):.2e}')
print('c_k^M @E=7/5 =', [round(float(sp.simplify(sol[ck[i]]).subs(sub)), 4)
                         for i in range(5)])
print('=> the Thakurta-Kerr machinery transferred to Vaidya (M(v) as the slow parameter).')
