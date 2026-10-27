# Closure of the HORIZON DILOGARITHM D_k = int U_k/(r-2m) dr (Vaidya, frozen Schwarzschild).
# ---------------------------------------------------------------------------
# STRUCTURE (exact integration by parts):
#   eta_h = dr/(r-2m).  As a 1-form on the curve y^2=S, in t=sqrt(r-2m):
#     r-2m=t^2, dr=2t dt  =>  eta_h = 2 dt/t  = SIMPLE pole, residue 2 at the
#     Weierstrass point above r=2m  => THIRD KIND (dlog letter at the HORIZON).
#   Elementary primitive: int eta_h = ln(r-2m).
#   => D_k = U_k ln(r-2m) - G_k,   G_k = int ln(r-2m) r^k/sqrt(S) dr   (STEP1)
#   G_k = genus-2 hyperelliptic dilogarithm with a dlog letter at the HORIZON
#         (companion of W_kj, which has the dipole letter at INFINITY).
# SEPARATRIX |J|=Jc: S acquires a double root -> the genus drops to 1 ->
#   G_k, D_k CLOSE in the ELLIPTIC DILOGARITHM (Bloch-Wigner), tabulated (STEP2).
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

import numpy as np
from scipy.integrate import quad
from scipy.optimize import brentq
import sympy as sp

M,E,J=1.0,1.4,2.5; r0=12.0
def DE(r): return (E**2-1)*r+2*M
def Sn(r):  return r*(r-2*M)*DE(r)*(r**2*(r-2*M)-J**2*DE(r))
sq=lambda x:np.sqrt(Sn(x))
rs=np.linspace(2.01,r0,20000); vv=Sn(rs); idx=np.where(np.diff(np.sign(vv)))[0]
rmin=max(brentq(Sn,rs[i],rs[i+1]) for i in idx if rs[i]>2.0); xf=rmin+0.4

def Uk(x,k): return quad(lambda t:t**k/sq(t),r0,x,limit=200)[0]
def Dk_dir(x,k): return quad(lambda t:Uk(t,k)/(t-2*M),r0,x,limit=150)[0]     # DIRECT
def Gk(x,k): return quad(lambda t:np.log(abs(t-2*M))*t**k/sq(t),r0,x,limit=150)[0]

# ---- STEP1: D_k = U_k ln(r-2m) - G_k  (integration-by-parts identity) ----
print('STEP1  D_k = U_k ln(r-2m) - G_k   (by parts, dlog letter at the horizon)')
print(' k :  D_k direct       U_k ln - G_k       diff')
c=np.log(abs(r0-2*M))   # note: U_k(r0)=0, so the boundary term at r0 vanishes
for k in range(5):
    d_dir=Dk_dir(xf,k)
    d_ibp=Uk(xf,k)*np.log(abs(xf-2*M)) - Gk(xf,k)
    _results.add('residual', f'horizon_k{k}', abs(d_dir-d_ibp))
    print(f' {k} : {d_dir:+.9f}   {d_ibp:+.9f}   {abs(d_dir-d_ibp):.2e}')

# ---- separatrix Jc: the bracket B(r)=r^2(r-2m)-J^2 DE has a double root ----
r=sp.symbols('r'); Js=sp.symbols('J',positive=True)
B=r**2*(r-2*M)-Js**2*((E**2-1)*r+2*M)
disc=sp.discriminant(sp.Poly(B,r),r)
Jc_sols=[s for s in sp.solve(sp.Eq(disc,0),Js) if s.is_real and s>0]
Jc=float(min(Jc_sols,key=lambda s:abs(float(s))))
print(f'\nSTEP2  separatrix Jc = {Jc:.6f}  (double root of the bracket)')
Bc=sp.Poly(B.subs(Js,Jc),r)
roots=sorted([complex(x).real for x in np.roots([float(c) for c in Bc.all_coeffs()]) if abs(complex(x).imag)<1e-6])
print(f'  roots of the bracket @Jc: {[round(x,5) for x in roots]}  (double = genus 1)')
print('  => on Jc: sqrt(S)=(r-r_d) sqrt(quartic), U_k elliptic,')
print('     D_k = ELLIPTIC DILOGARITHM (Bloch-Wigner), tabulated.')

_results.finish()
