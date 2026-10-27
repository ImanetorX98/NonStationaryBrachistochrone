# FULLY EXPLICIT FORM of delta phi_V/Mdot (adiabatic Vaidya).
# Every integral expanded in NAMED building blocks:
#   U_k      = int r^k/sqrt(S) dr                     (second kind, genus-2 polylogarithm)
#   W_3k     = int (U_3 dU_k - U_k dU_3)              (weight-2 polylogarithm, genus 2)
#   D_k      = int U_k/(r-2M) dr                      (HORIZON dilogarithm, weight 2)
#   L_2M     = int Acal^M r/((r-2M) sqrt(S)) dr       (third kind, horizon pole)   [note below]
#   I_poly   = int Acal^M r^3/S dr                    (ELEMENTARY: logs/rational fractions)
# ALGEBRAIC coefficients c_k^M, a_j^M (from the reduction dM F = d(Acal^M/sqrtS)+sum c_k^M r^k/sqrtS).
#
# ASSEMBLY (by parts + reduction, exact):
#   delta phi_V/Mdot = A_M v
#       - E*I_poly - (E/2) sum_k c_k^M (U_k U_3 - W_3k)
#       - sum_j a_j^M U_j - sum_k c_k^M (r U_k - U_{k+1}) - 2M*Lh - 2M sum_k c_k^M D_k
# with A_M = Acal^M/sqrtS + sum_k c_k^M U_k,  v = E U_3 + (r-r0) + 2M ln((r-2M)/(r0-2M)).
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

import sympy as sp, numpy as np
from scipy.integrate import quad
from scipy.optimize import brentq

M,E,J=1.0,1.4,2.5; r0=12.0
r=sp.symbols('r',positive=True); Ms,Es,Js=sp.symbols('M E J',positive=True)
Emu=(Es**2-1)*r+2*Ms
Ssym=sp.expand(r*(r-2*Ms)*Emu*(r**2*(r-2*Ms)-Js**2*Emu))
sub={Ms:1,Es:sp.Rational(7,5),Js:sp.Rational(5,2)}
Sn=sp.lambdify(r,Ssym.subs(sub),'numpy'); sq=lambda x:np.sqrt(Sn(x))

# ---- reduction: dM F = d(Acal/sqrtS) + sum c_k r^k/sqrtS ----
K=Js*Emu; NM=sp.expand(Ssym*sp.diff(K,Ms)-K*sp.diff(Ssym,Ms)/2)
Ac=[sp.Symbol(f'A{i}') for i in range(6)]; ck=[sp.Symbol(f'c{i}') for i in range(5)]
Acal=sum(Ac[i]*r**i for i in range(6)); Mp=sum(ck[i]*r**i for i in range(5))
sol=sp.solve(sp.Poly(sp.expand(2*NM-(2*Ssym*sp.diff(Acal,r)-Acal*sp.diff(Ssym,r)+2*Ssym*Mp)),r).all_coeffs(),Ac+ck,dict=True)[0]
cM=[float(sp.simplify(sol[ck[i]]).subs(sub)) for i in range(5)]
aM=[float(sp.simplify(sol.get(Ac[i],0)).subs(sub)) for i in range(6)]  # coefficients of Acal^M, deg5
def Acaln(x): return sum(aM[i]*x**i for i in range(6))
dMF=(sp.diff(K,Ms)/sp.sqrt(Ssym)-K*sp.diff(Ssym,Ms)/(2*Ssym**sp.Rational(3,2)))
dMFn=sp.lambdify(r,dMF.subs(sub),'numpy')
print('c_k^M =',[round(c,5) for c in cM])
print('a_j^M =',[round(a,5) for a in aM])

# ---- integration endpoint (inside, beyond the turning point) ----
rs=np.linspace(2.01,r0,20000); vv=Sn(rs); idx=np.where(np.diff(np.sign(vv)))[0]
rmin=max(brentq(Sn,rs[i],rs[i+1]) for i in idx if rs[i]>2.0); xf=rmin+0.4

# ---- building block ----
def Uk(x,k): return quad(lambda t:t**k/sq(t),r0,x,limit=200)[0]
def W3k(x,k): # int (U_3 dU_k - U_k dU_3) = int (U_3 r^k - U_k r^3)/sqrtS dr
    return quad(lambda t:(Uk(t,3)*t**k-Uk(t,k)*t**3)/sq(t),r0,x,limit=150)[0]
def Dk(x,k): return quad(lambda t:Uk(t,k)/(t-2*M),r0,x,limit=150)[0]
def Lh(x): return quad(lambda t:Acaln(t)*t/((t-2*M)*sq(t)),r0,x,limit=150)[0]
def Ipoly(x): return quad(lambda t:Acaln(t)*t**3/Sn(t),r0,x,limit=150)[0]
def AM(x): return Acaln(x)/sq(x)-Acaln(r0)/sq(r0)+sum(cM[k]*Uk(x,k) for k in range(5))
def vclock(x): return E*Uk(x,3)+(x-r0)+2*M*(np.log(x-2*M)-np.log(r0-2*M))

# ---- DIRECT ----
dphi_dir=quad(lambda x:dMFn(x)*vclock(x),r0,xf,limit=150)[0]

# ---- FULLY EXPLICIT ASSEMBLY ----
# A_M = Acal/sqrtS - alpha0 + sum c_k U_k,  alpha0 = Acal(r0)/sqrtS(r0) constant.
x=xf
a0 = Acaln(r0)/np.sqrt(Sn(r0))
bdy = AM(x)*vclock(x)
# E * int A_M r^3/sqrtS = E[ I_poly - a0 U_3 + sum c_k (1/2)(U_k U_3 - W_3k) ]
p2  = E*( Ipoly(x) - a0*Uk(x,3) + 0.5*sum(cM[k]*(Uk(x,k)*Uk(x,3)-W3k(x,k)) for k in range(5)) )
# int A_M r/(r-2M) = Lh - a0 (v - E U_3) + sum c_k[(r U_k - U_{k+1}) + 2M D_k]
p3  = ( Lh(x) - a0*(vclock(x)-E*Uk(x,3))
      + sum(cM[k]*(x*Uk(x,k)-Uk(x,k+1)) for k in range(5))
      + 2*M*sum(cM[k]*Dk(x,k) for k in range(5)) )
dphi_expl = bdy - p2 - p3

print(f'\ndelta phi_V/Mdot  DIRECT         = {dphi_dir:.10f}')
print(f'                  FULL EXPLICIT  = {dphi_expl:.10f}')
_results.add('residual', 'assembly', abs(dphi_dir-dphi_expl))
print(f'  diff = {abs(dphi_dir-dphi_expl):.2e}   <-- identity (algebraic coefficients given)')

# ---- values of the building blocks (for the paper) ----
print('\nBuilding block @ xf=%.4f:'%xf)
print('  U_k   =',[round(Uk(x,k),4) for k in range(6)])
print('  W_3k  =',[round(W3k(x,k),4) for k in range(5)])
print('  D_k   =',[round(Dk(x,k),4) for k in range(5)])
print('  L_2M  =',round(Lh(x),4),'   I_poly =',round(Ipoly(x),4))
print('  A_M   =',round(AM(x),4),'   v =',round(vclock(x),4))

_results.finish()
