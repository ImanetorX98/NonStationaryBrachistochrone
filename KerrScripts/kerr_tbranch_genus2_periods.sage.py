# -*- coding: utf-8 -*-
"""Riemann period matrix of the generic t-branch curve y^2 = R6 (Paper II, eq. t-genus2).

M = 1, a = 9/10, E = 6/5, J = 3: the genus-2 orbit of kerr_jpt_genus2_figure.py, which
penetrates the ergosphere and reproduces the Hamiltonian flow to ~2e-13. Checks the two
Riemann relations: tau = tau^T and Im(tau) positive definite. The eigenvalues of Im(tau)
depend on the homology basis RiemannSurface chooses and are printed for the record only.

Run with the Python of the Sage environment (conda's `sage` has no `-python` option; see
NonStationaryMetrics/sage_env/SAGE_SETUP.md), e.g. via NonStationaryMetrics/sage_env/run_sage_checks.py, or
    <env>/bin/python kerr_tbranch_genus2_periods.sage.py
"""
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
_results.add('symmetry', 'norm', sym)
for _i, _e in enumerate(eig):
    _results.add('positive_eigenvalues', f'eigenvalue_{_i}', _e)
_results.finish()
sys.exit(0 if ok else 1)
