# -*- coding: utf-8 -*-
# Adiabatico Vaidya: riduzione di dM F (analogo di dE F di Thakurta-Kerr).
# Frozen = Schwarzschild (a=0), curva genus-2 S. dM F = N_M/S^(3/2) -> riduzione
# dM F = d(A/sqrtS) + sum c_k^M r^k/sqrtS. Verificata numericamente a ~1e-15.
import sympy as sp

r, M, E, J = sp.symbols('r M E J', positive=True)
Emu = (E**2-1)*r + 2*M
S = sp.expand(r*(r-2*M)*Emu*(r**2*(r-2*M) - J**2*Emu))    # sestica tau, a=0
K = J*Emu                                                 # K_tau a=0

# N_M polinomiale: dM F = N_M/S^(3/2), N_M = S dM K - K dM S / 2
N_M = sp.expand(S*sp.diff(K, M) - K*sp.diff(S, M)/2)
print('N_M grado', sp.degree(N_M, r), '(polinomiale)')

# riduzione: 2 N_M = 2S A' - A S' + 2S Mpoly, A deg5, Mpoly deg4
Ac = [sp.Symbol(f'A{i}') for i in range(6)]
ck = [sp.Symbol(f'c{i}') for i in range(5)]
Acal = sum(Ac[i]*r**i for i in range(6)); Mp = sum(ck[i]*r**i for i in range(5))
eq = sp.expand(2*N_M - (2*S*sp.diff(Acal, r) - Acal*sp.diff(S, r) + 2*S*Mp))
sol = sp.solve(sp.Poly(eq, r).all_coeffs(), Ac + ck, dict=True)[0]
print('riduzione risolta')

# verifica numerica
sub = {M: 1, E: sp.Rational(7, 5), J: sp.Rational(5, 2)}
dMF = sp.diff(K, M)/sp.sqrt(S) - K*sp.diff(S, M)/(2*S**sp.Rational(3, 2))
dMFn = sp.lambdify(r, dMF.subs(sub), 'numpy')
rhs = sp.diff(Acal.subs(sol)/sp.sqrt(S), r) + Mp.subs(sol)/sp.sqrt(S)
rhsn = sp.lambdify(r, rhs.subs(sub), 'numpy')
print('VERIFICA dM F == d(A/sqrtS) + M_poly/sqrtS:')
for rv in [11., 9., 7., 5.]:
    print(f'  r={rv}: diff = {abs(dMFn(rv)-rhsn(rv)):.2e}')
print('c_k^M @E=7/5 =', [round(float(sp.simplify(sol[ck[i]]).subs(sub)), 4)
                         for i in range(5)])
print('=> macchina di Thakurta-Kerr trasferita a Vaidya (M(v) parametro lento).')
