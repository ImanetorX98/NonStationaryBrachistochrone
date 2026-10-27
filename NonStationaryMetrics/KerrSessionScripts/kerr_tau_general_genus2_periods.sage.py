"""Riemann period matrix of the generic tau-branch curve y^2 = sextic (Paper II, genus-2 tau branch).

M = 1, a = 9/10, E = 6/5, J = -0.9 J_c = -27/40 (J_c = a/E); sextic = r (r-2M)((E^2-1)r+2M)(r Delta - J^2((E^2-1)r+2M)),
the curve of kerr_tau_general_genus2.py.  Prints the six branch points, the integral model, the genus and the Riemann
matrix tau computed by SageMath's RiemannSurface (40 bits), and checks the two Riemann relations: tau = tau^T and
Im(tau) positive definite.  The eigenvalues of Im(tau) depend on the homology basis RiemannSurface chooses and are
printed for the record only.  The t-branch companion is KerrScripts/kerr_tbranch_genus2_periods.sage.py.

Run with the Python of the Sage environment (conda's `sage` has no `-python` option; see
NonStationaryMetrics/sage_env/SAGE_SETUP.md), e.g. via NonStationaryMetrics/sage_env/run_sage_checks.py, or
    <env>/bin/python kerr_tau_general_genus2_periods.sage.py
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

import time
t0=time.time()
from sage.all import QQ, ZZ, PolynomialRing, lcm, expand
from sage.schemes.riemann_surfaces.riemann_surface import RiemannSurface
def log(m): print(f"[{time.time()-t0:6.1f}s] {m}", flush=True)

Rx=PolynomialRing(QQ,'x'); x=Rx.gen()
M=QQ(1); a=QQ(9)/10; E=QQ(6)/5; Jc=a/E; J=-QQ(9)/10*Jc   # -0.9 Jc = -27/40
log(f"J = -0.9 Jc = {J}")
num_wf  = (x-2*M)*((E**2-1)*x+2*M)
num_DJw = x*(x**2-2*M*x+a**2) - J**2*((E**2-1)*x+2*M)
sextic  = expand(x*num_wf*num_DJw)
log(f"sextic = {sextic}")
# roots
from sage.all import CDF
rts=sextic.roots(CDF, multiplicities=False)
log("branch points:")
for z in sorted(rts, key=lambda z: z.real()):
    print(f"     {complex(z).real:+.4f} {complex(z).imag:+.4f}i")
# integral model
L=lcm([c.denominator() for c in sextic.coefficients()])
R6i=expand(L*sextic)
log(f"L={L}; integral? {all(c in ZZ for c in R6i.coefficients())}")
Rxy=PolynomialRing(QQ,['X','Y']); X,Y=Rxy.gens()
f=Y**2-R6i(X)
log("build RiemannSurface prec=40 ...")
S=RiemannSurface(f, prec=40)
log(f"  genus = {S.genus}")
tau=S.riemann_matrix()
log("Riemann matrix tau =")
for row in tau:
    print("   " + "  ".join(f"{complex(z).real:+.5f}{complex(z).imag:+.5f}i" for z in row))
sym=(tau-tau.transpose()).norm()
from sage.all import RealField, matrix
ImT=matrix(RealField(40),2,2,[complex(tau[i][j]).imag for i in range(2) for j in range(2)])
log(f"|tau-tau^T| = {float(sym):.2e} (symmetric); eigenvalues of Im(tau) = {[float(e) for e in ImT.eigenvalues()]} (>0)")
_results.add('symmetry', 'norm', sym)
for _i, _e in enumerate(ImT.eigenvalues()):
    _results.add('positive_eigenvalues', f'eigenvalue_{_i}', _e)
log("DONE.")

_results.finish()
