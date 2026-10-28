# CLOSING THE LOOP: assembly of psi for the t-branch (on-curve), verified end to end.
# Clock t(r) = eta2 + eta3:  eta2 = int P3/sqrt(R6) = sum_j b_j U_j  (SECOND KIND),
#                            eta3 = int R_Delta/(Delta sqrt(R6))     (THIRD KIND horizon).
# By parts (identity eq:psi-split):  with I := int dE F * eta2 dr,
#   I = -1/2 sum_{k<j} Q_kj W_kj + eta2*(Acal/sqrt(R6)+1/2 sum c_k U_k) - sum_j b_j Ical_j
#   Q_kj=c_k b_j-c_j b_k,  W_kj=int(U_k dU_j-U_j dU_k),  Ical_j=int Acal r^j/R6 dr.
# Full: int dE F * t dr = I + int dE F * eta3 dr  (the horizon piece, log-sigma at r_pm).
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

r0=12.0
r=sp.symbols('r',positive=True); M,a,E,J=sp.symbols('M a E J',positive=True)
sub={M:1,a:sp.Rational(9,10),E:sp.Rational(6,5),J:sp.Rational(5,2)}
DE=(E**2-1)*r+2*M; Delta=r**2-2*M*r+a**2
Q2=E**2*r**4+(E**2*(a**2-J**2)+J**2)*r**2+2*M*(J-a)*(E**2*(J-a)-2*J)*r+4*M**2*(J-a)**2
R6=sp.expand(r*Q2*DE); Kt=r*DE*(J*(r-2*M)+2*M*a)/Delta
# dE F_t = N/R6^{3/2}
Nt=sp.expand(sp.simplify(sp.diff(Kt/sp.sqrt(R6),E)*R6**sp.Rational(3,2)))
# reduction (numeric params): 2 A' R6 - A R6' + 2 R6 sum c_k r^k = 2 Nt
R6s=sp.expand(R6.subs(sub)); Nts=sp.expand(Nt.subs(sub))
ai=[sp.Symbol(f'a{i}') for i in range(6)]; ck=[sp.Symbol(f'c{i}') for i in range(5)]
A=sum(ai[i]*r**i for i in range(6)); Mp=sum(ck[i]*r**i for i in range(5))
ident=sp.expand(2*sp.diff(A,r)*R6s-A*sp.diff(R6s,r)+2*R6s*Mp-2*Nts)
sol=sp.solve(sp.Poly(ident,r).all_coeffs(),ai+ck,dict=True)[0]
cN=[float(sol[ck[i]]) for i in range(5)]
# clock second-kind vector b (P3 coeffs) and horizon R_Delta  -- SYMBOLIC then sub
rho_t=sp.cancel(sp.together((E**2*r**3-2*M*a*Kt/r)/((r-2*M)/r)))
P3s,Rrems=sp.div(sp.Poly(sp.numer(sp.cancel(rho_t)),r),sp.Poly(sp.expand(Delta),r))
P3=sp.Poly(P3s.as_expr().subs(sub),r); Rrem=Rrems.as_expr().subs(sub)
p3=list(reversed(P3.all_coeffs()))            # low->high  b_0..b_3
b=[float(p3[i]) if i<len(p3) else 0.0 for i in range(5)]
RD=sp.lambdify(r,Rrem,'numpy')
print('c_k^t =',[round(x,5) for x in cN]); print('b =',[round(x,5) for x in b])

Acal=sum(float(sol.get(ai[i],0))*r**i for i in range(6))
R6n=sp.lambdify(r,R6s,'numpy'); Nn=sp.lambdify(r,Nts,'numpy'); Dn=sp.lambdify(r,Delta.subs(sub),'numpy')
Acaln=sp.lambdify(r,Acal,'numpy')
sq=lambda x:np.sqrt(R6n(x))
def dEF(x): return Nn(x)/R6n(x)**1.5
def U(x,k): return quad(lambda t:t**k/sq(t),r0,x,limit=200)[0]
def Uk(x): return [U(x,k) for k in range(5)]
def eta2(x): return sum(b[j]*U(x,j) for j in range(5))
def eta3(x): return quad(lambda t:RD(t)/(Dn(t)*sq(t)),r0,x,limit=200)[0]
def W(x,k,j): return quad(lambda t:(U(t,k)*t**j-U(t,j)*t**k)/sq(t),r0,x,limit=150)[0]
def Ical(x,j): return quad(lambda t:Acaln(t)*t**j/R6n(t),r0,x,limit=200)[0]

# turning point (largest root of R6 in (rp,r0)); else integrate near horizon rp
rp=1+np.sqrt(1-0.9**2)
rr=np.linspace(rp+0.05,r0,40000); vv=R6n(rr); idx=np.where(np.diff(np.sign(vv)))[0]
roots=[brentq(R6n,rr[i],rr[i+1]) for i in idx if rr[i]>rp]
if roots: rt=max(roots); xf=rt+0.35; print(f'turning R6=0 at {rt:.4f}, xf={xf:.4f}')
else:     xf=rp+0.6; print(f'no turning in range (scattering); integrate to xf={xf:.4f}')

x=xf
# --- I_direct = int dE F * eta2 ---
I_dir=quad(lambda t:dEF(t)*eta2(t),r0,xf,limit=200)[0]
# --- I_assembled ---
Uv=Uk(x); Qsum=0.0
for k in range(5):
    for j in range(k+1,5):
        Qkj=cN[k]*b[j]-cN[j]*b[k]
        if Qkj!=0.0: Qsum+=Qkj*W(x,k,j)
I_asm=-0.5*Qsum + eta2(x)*(Acaln(x)/sq(x)+0.5*sum(cN[k]*Uv[k] for k in range(5))) - sum(b[j]*Ical(x,j) for j in range(5))
print(f'\nSECOND-KIND assembly (eq:psi-split):')
print(f'  I_direct    = {I_dir:.10f}')
print(f'  I_assembled = {I_asm:.10f}')
_results.add('residual', 'second_kind_assembly', abs(I_dir-I_asm))
print(f'  diff = {abs(I_dir-I_asm):.2e}   <-- W_kj reorganization for t-branch')

# --- full clock closure: int dE F * t = I + int dE F * eta3 ---
full_dir=quad(lambda t:dEF(t)*(eta2(t)+eta3(t)),r0,xf,limit=200)[0]
dphi3=quad(lambda t:dEF(t)*eta3(t),r0,xf,limit=200)[0]
print(f'\nFULL clock t = eta2 + eta3 (horizon 3rd-kind):')
print(f'  int dEF*t direct         = {full_dir:.10f}')
print(f'  I_assembled + int dEF*eta3 = {I_asm+dphi3:.10f}')
_results.add('residual', 'full_clock_assembly', abs(full_dir-(I_asm+dphi3)))
print(f'  diff = {abs(full_dir-(I_asm+dphi3)):.2e}')
print(f'  horizon piece int dEF*eta3 = {dphi3:.6f}  (= log-sigma at r_pm, weight-2)')

_results.finish()
