# -*- coding: utf-8 -*-
# Adiabatic Vaidya: reduction of dM F (the analogue of dE F for Thakurta-Kerr).
# Frozen = Schwarzschild (a=0), genus-2 curve S. dM F = N_M/S^(3/2) -> reduction
# dM F = d(A/sqrtS) + sum c_k^M r^k/sqrtS. Checked numerically to ~1e-15.
# Stable verification fields are independent of the human-readable messages.
# Without verification_results.py (e.g. a script copied out of the repository) the
# program runs as before and exports nothing, unless NSB_RESULT_JSON asks for a result.
import sys as _result_sys
from pathlib import Path as _ResultPath
_result_source = _ResultPath(__file__).resolve()
_result_module = next((p / "NonStationaryMetrics" for p in _result_source.parents
                       if (p / "NonStationaryMetrics" / "verification_results.py").is_file()), None)
if _result_module is not None:
    _result_sys.path.insert(0, str(_result_module))
    from verification_results import ResultWriter as _ResultWriter
    _results = _ResultWriter(_result_source)
else:
    import os as _result_os
    if _result_os.environ.get("NSB_RESULT_JSON"):
        raise SystemExit("cannot export NSB_RESULT_JSON: verification_results.py not found above "
                         + str(_result_source))
    class _NoResults:
        def add(self, metric, sample, value): pass
        def finish(self): pass
    _results = _NoResults()

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
    _results.add('residual', f'reduction_r{rv:g}', abs(dMFn(rv)-rhsn(rv)))
    print(f'  r={rv}: diff = {abs(dMFn(rv)-rhsn(rv)):.2e}')
print('c_k^M @E=7/5 =', [round(float(sp.simplify(sol[ck[i]]).subs(sub)), 4)
                         for i in range(5)])
print('=> the Thakurta-Kerr machinery transferred to Vaidya (M(v) as the slow parameter).')

_results.finish()
