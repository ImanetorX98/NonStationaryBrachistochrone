# Vaidya TAU BRANCH (proper time): first-order adiabatic correction, ANALYTIC.
# Frozen = Schwarzschild (a=0). MASS reduction dm F = N_m/S^{3/2}.
# Proper clock: dtau/dr = r^2(r-2m)/sqrt(S) => eta = U_3 - 2m U_2,  b=(0,0,-2m,1,0).
# Assembly (eq:psi-split), Q_kj = c_k^m b_j - c_j^m b_k, ALL algebraic (no fit):
#   delta phi_tau / mdot = int dm F * eta dr
#      = -1/2 sum_{k<j} Q_kj W_kj + eta(Acal^m/sqrtS + 1/2 sum c_k^m U_k) - sum_j b_j I_j
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

r=sp.symbols('r',positive=True); m,E,J=sp.symbols('m E J',positive=True)
sub={m:1,E:sp.Rational(7,5),J:sp.Rational(5,2)}; r0=12.0
DE=(E**2-1)*r+2*m
S=sp.expand(r*(r-2*m)*DE*(r**2*(r-2*m)-J**2*DE))
K=J*DE                                  # tau-branch numerator (a=0 limit of TK)
# --- mass reduction: dm F = N_m/S^{3/2},  N_m = S dm K - 1/2 K dm S ---
Nm=sp.expand(S*sp.diff(K,m)-sp.Rational(1,2)*K*sp.diff(S,m))
Ss=sp.expand(S.subs(sub)); Ns=sp.expand(Nm.subs(sub))
assert Ns.is_polynomial(r), "N_m is not a polynomial"
ai=[sp.Symbol(f'a{i}') for i in range(6)]; ck=[sp.Symbol(f'c{i}') for i in range(5)]
A=sum(ai[i]*r**i for i in range(6)); Mp=sum(ck[i]*r**i for i in range(5))
sol=sp.solve(sp.Poly(sp.expand(2*sp.diff(A,r)*Ss-A*sp.diff(Ss,r)+2*Ss*Mp-2*Ns),r).all_coeffs(),ai+ck,dict=True)[0]
cM=[float(sol[ck[i]]) for i in range(5)]
Acal=sum(sol.get(ai[i],0)*r**i for i in range(6))
resid=sp.Poly(sp.expand(2*sp.diff(A,r)*Ss-A*sp.diff(Ss,r)+2*Ss*Mp-2*Ns).subs(sol),r).all_coeffs()
print("c_k^m =",[round(x,5) for x in cM]," reduction residual =",float(max(abs(x) for x in resid)) if resid else 0.0)

# --- clock tau: b=(0,0,-2m,1,0) ---
b=[0,0,-2*sub[m],1,0]
print("b (clock tau, eta=U_3-2m U_2) =",[float(x) if not isinstance(x,int) else x for x in b])

# --- Q_kj = c_k^m b_j - c_j^m b_k (algebraic) ---
Qlist=[]
for k in range(5):
    for j in range(k+1,5):
        Q=cM[k]*b[j]-cM[j]*b[k]
        if abs(Q)>1e-14: Qlist.append((k,j,Q))
print("Q_kj:",[f'Q_{k}{j}={Q:+.5f}' for k,j,Q in Qlist])

# --- numerical building blocks ---
Sn=sp.lambdify(r,Ss,'numpy'); sq=lambda x:np.sqrt(Sn(x))
Acaln=sp.lambdify(r,Acal.subs(sub) if Acal.free_symbols-{r} else Acal,'numpy')
dMF=(sp.diff(K,m)/sp.sqrt(S)-K*sp.diff(S,m)/(2*S**sp.Rational(3,2)))
dMFn=sp.lambdify(r,dMF.subs(sub),'numpy')
rr=np.linspace(2.01,r0,20000); vv=Sn(rr); idx=np.where(np.diff(np.sign(vv)))[0]
rt=max(brentq(Sn,rr[i],rr[i+1]) for i in idx if rr[i]>2.0); xf=rt+0.4
def U(x,k): return quad(lambda t:t**k/sq(t),r0,x,limit=200)[0]
def W(x,k,j): return quad(lambda t:(U(t,k)*t**j-U(t,j)*t**k)/sq(t),r0,x,limit=150)[0]
def Ical(x,j): return quad(lambda t:Acaln(t)*t**j/Sn(t),r0,x,limit=200)[0]
def eta(x): return sum(b[j]*U(x,j) for j in range(5))

x=xf; Uv=[U(x,k) for k in range(5)]
print(f"\nturning point S=0 at r={rt:.4f}, evaluating at x={x:.4f}")
# --- DIRECT ---
I_dir=quad(lambda t:dMFn(t)*eta(t),r0,xf,limit=200)[0]
# --- ASSEMBLED ---
I_asm=(-0.5*sum(Q*W(x,k,j) for k,j,Q in Qlist)
       + eta(x)*(Acaln(x)/sq(x)+0.5*sum(cM[k]*Uv[k] for k in range(5)))
       - sum(b[j]*Ical(x,j) for j in range(5)))
print(f"\ndelta phi_tau/mdot  DIRECT    = {I_dir:.10f}")
print(f"                    ASSEMBLED = {I_asm:.10f}")
_results.add('residual', 'assembly', abs(I_dir-I_asm))
print(f"  difference = {abs(I_dir-I_asm):.2e}   <-- Vaidya tau branch, all algebraic")

_results.finish()
