# ACCRETION/EVAPORATION ASYMMETRY of the adiabatic correction delta phi_V.
# Accretion (mdot>0): INGOING Vaidya  ds^2=-f dv^2+2 dv dr, ADVANCED clock  v = t + r_*.
# Evaporation (mdot<0): OUTGOING Vaidya ds^2=-f du^2-2 du dr, RETARDED clock u = t - r_*.
# DERIVATION of the clocks on the INGOING arc (launched at r0 towards the turning point):
#   u, v are defined by the SAME Killing time t through u=t-r_*, v=t+r_*.  On the ingoing arc
#   the elapsed Schwarzschild time is T=-E U_3>0 (dt/dr=-E r^3/sqrt(S)), hence
#     v = T + r_*  (dv/dr = -E r^3/sqrt(S) + r/(r-2m)),   u = T - r_*.
#   (r_* = (r-r0) + 2m ln((r-2m)/(r0-2m)),  dr_*/dr=r/(r-2m).)  The densities E r^3/sqrt(S) +/- r/(r-2m)
#   hold on the OUTGOING arc; using them here swapped advanced and retarded (corrected on 20/10/2026).
# The same reduction dM F (it depends only on the frozen S: ingoing/outgoing Schwarzschild = the same geometry);
# ONLY the clock (v vs u) and the sign of mdot change (accretion: m up / evaporation: m down).
# => delta phi/mdot|_accr = A + B ,  delta phi/mdot|_evap = A - B  (COEFFICIENTS, /mdot with sign),
#    A = int d_m(dphi/dr) * T dr   (part at INFINITY, polylogarithm), along the motion
#    B = int d_m(dphi/dr) * r_* dr  (HORIZON part, tortoise dilogarithm)
# PHYSICS (delta phi = mdot*[.]; mdot>0 accr, mdot<0 evap): delta phi_accr=+|mdot|(A+B),
#    delta phi_evap=-|mdot|(A-B). Hence (at equal |mdot|):
#    physical DIFFERENCE = accr - evap = |mdot|*2A   (carried by A -- NOT by B!)
#    physical SUM        = accr + evap = |mdot|*2B   = FAILURE OF ANTISYMMETRY (carried by B, horizon)
# (Correction: B controls the failure of antisymmetry, not the physical difference.)
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
K=Js*Emu
dMF=(sp.diff(K,Ms)/sp.sqrt(Ssym)-K*sp.diff(Ssym,Ms)/(2*Ssym**sp.Rational(3,2)))
dMFn=sp.lambdify(r,dMF.subs(sub),'numpy')

rs=np.linspace(2.01,r0,20000); vv=Sn(rs); idx=np.where(np.diff(np.sign(vv)))[0]
rmin=max(brentq(Sn,rs[i],rs[i+1]) for i in idx if rs[i]>2.0); xf=rmin+0.4
def U3(x): return quad(lambda t:t**3/sq(t),r0,x,limit=200)[0]
def rstar(x): return (x-r0)+2*M*(np.log(x-2*M)-np.log(r0-2*M))
# The arc is INGOING: launched at r0 and run inward to xf.  On it the Schwarzschild
# time elapsed is T = -E U_3 > 0 (U_3 anchored at r0 is negative for x < r0), so
#   advanced  v = T + r_*,   retarded  u = T - r_*,
# and the shape integrand along the motion is dphi/dr = -J D_E/sqrt(S), hence -dMF.
# An earlier version used the OUTGOING densities E U_3 +/- r_* on this ingoing arc:
# that swaps the advanced and retarded coefficients and flips the sign of B.
def T_el(x): return -E*U3(x)
def v_adv(x): return T_el(x)+rstar(x)     # accretion (ingoing, elapsed advanced time)
def u_ret(x): return T_el(x)-rstar(x)     # evaporation (outgoing, u = v - 2 r_*)
dMF_mot=lambda x: -dMFn(x)                 # d_m (dphi/dr) along the ingoing motion

# --- independent check of the clock and of the orientation against the TRUE frozen
#     Hamiltonian flow of H_tau (eq:Htau-vaidya), whose flow parameter is v ---
from scipy.integrate import solve_ivp
rr_,pp_=sp.symbols('rr pp')
fv=1-2*M/rr_; wv=E**2-fv
JJ_=sp.Symbol('JJ')
HsJ=pp_*(fv-E**2)-E+sp.sqrt(wv)*sp.sqrt((E*pp_+1)**2+JJ_**2/rr_**2)
Hs=HsJ.subs(JJ_,J)
Hn=sp.lambdify((rr_,pp_),Hs,'numpy'); Hp=sp.lambdify((rr_,pp_),sp.diff(Hs,pp_),'numpy')
Hr=sp.lambdify((rr_,pp_),sp.diff(Hs,rr_),'numpy')
HJ=sp.lambdify((rr_,pp_),sp.diff(HsJ,JJ_).subs(JJ_,J),'numpy')
pg=np.linspace(-30,30,60001); hv=Hn(r0,pg)
roots=[brentq(lambda q:Hn(r0,q),pg[i],pg[i+1]) for i in range(len(pg)-1) if hv[i]*hv[i+1]<0]
p0=[q for q in roots if Hp(r0,q)<0][0]                       # ingoing root
ev=lambda t,y: y[0]-xf; ev.terminal=True
fl=solve_ivp(lambda t,y:[Hp(y[0],y[1]),-Hr(y[0],y[1]),HJ(y[0],y[1])],[0,200],[r0,p0,0.0],
             rtol=1e-12,atol=1e-14,events=ev,dense_output=True)
worst_v=worst_phi=0.0
for lam in np.linspace(0.1*fl.t[-1],fl.t[-1],6):
    x,_,ph=fl.sol(lam)
    worst_v=max(worst_v,abs(lam-v_adv(x)))
    worst_phi=max(worst_phi,abs(ph-quad(lambda t:-J*((E**2-1)*t+2*M)/sq(t),r0,x,limit=200)[0]))
print(f"clock check against the H_tau flow: max|v_flow - (T + r_*)| = {worst_v:.1e}, "
      f"max|phi_flow - int(-J D_E/sqrt S)| = {worst_phi:.1e}")
assert worst_v<1e-8 and worst_phi<1e-8, "ingoing clock or orientation inconsistent with the flow"

# --- parts A (symmetric) and B (antisymmetric = horizon) ---
A=quad(lambda x:dMF_mot(x)*T_el(x), r0,xf,limit=150)[0]
B=quad(lambda x:dMF_mot(x)*rstar(x),r0,xf,limit=150)[0]
# --- direct values with the two clocks ---
accr=quad(lambda x:dMF_mot(x)*v_adv(x),r0,xf,limit=150)[0]
evap=quad(lambda x:dMF_mot(x)*u_ret(x),r0,xf,limit=150)[0]

print("=== delta phi_V : coefficients (/mdot with sign) and PHYSICAL quantities ===")
print(f"  A (INFINITY part, polylog)   = {A:+.8f}")
print(f"  B (HORIZON part, dilog)      = {B:+.8f}")
_results.add('residual', 'accretion', abs(A+B-accr))
print(f"  coeff accr = A+B (clock v)   = {A+B:+.8f}   direct={accr:+.8f}  diff={abs(A+B-accr):.1e}")
_results.add('residual', 'evaporation', abs(A-B-evap))
print(f"  coeff evap = A-B (clock u)   = {A-B:+.8f}   direct={evap:+.8f}  diff={abs(A-B-evap):.1e}")
# --- PHYSICS: delta phi = mdot*[.], mdot>0 accr, mdot<0 evap => signs ---
dphi_accr=(A+B)     # /|mdot|, mdot>0
dphi_evap=-(A-B)    # /|mdot|, mdot<0 (delta phi_evap = mdot(A-B) = -|mdot|(A-B))
print("\n=== PHYSICS (delta phi / |mdot|) ===")
print(f"  delta phi_accr/|mdot| = +(A+B) = {dphi_accr:+.6f}")
print(f"  delta phi_evap/|mdot| = -(A-B) = {dphi_evap:+.6f}")
print(f"  DIFFERENCE  accr-evap = 2A = {dphi_accr-dphi_evap:+.6f}   (carried by A, INFINITY -- NOT by B!)")
print(f"  SUM         accr+evap = 2B = {dphi_accr+dphi_evap:+.6f}   = FAILURE OF ANTISYMMETRY (B, horizon)")
print("  (Correction: the old wording 'asymmetry=2B' confused the difference of the coefficients")
print("   (A+B)-(A-B)=2B with the PHYSICAL difference accr-evap=2A; B measures the deviation from antisymmetry.)")

# --- structure of B: bulk (second kind) + HORIZON sector (dilogarithm D_k) ---
# r_* = (r-r0)  +  2m ln((r-2m)/(r0-2m))
#        \___bulk___/     \_____horizon log_____/
B_bulk=quad(lambda x:dMF_mot(x)*(x-r0), r0,xf,limit=150)[0]
B_hor =quad(lambda x:dMF_mot(x)*2*M*(np.log(x-2*M)-np.log(r0-2*M)), r0,xf,limit=150)[0]
print("\n=== structure of the asymmetry B = B_bulk + B_hor ===")
print(f"  B_bulk (second kind, r-r0)         = {B_bulk:+.8f}")
print(f"  B_hor  (horizon log -> dilog D_k)  = {B_hor:+.8f}")
_results.add('residual', 'horizon_split', abs(B_bulk+B_hor-B))
print(f"  sum = {B_bulk+B_hor:+.8f}   direct B = {B:+.8f}   diff={abs(B_bulk+B_hor-B):.1e}")
# B_hor by parts: 2m[ A_m ln(r-2m) ]_{r0}^{xf} - 2m int A_m/(r-2m) dr,  int A_m/(r-2m) = sum c_k D_k
print("  => B_hor = 2m[A_m ln(r-2m)] - 2m*int A_m/(r-2m);  int A_m/(r-2m) = HORIZON dilogarithm D_k.")
print("     The transcendental accr/evap asymmetry (horizon log) IS the dilogarithm D_k derived above.")

_results.finish()
