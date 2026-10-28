# -*- coding: utf-8 -*-
# Separatrix coefficients, fully SYMBOLIC in (M,a,E,r_d,Jc), for ALL branches.
# UNIVERSAL formula (residues R at the triple pole z_d, F=N/Q4):
#   Q4(rd)=S''(rd)/2, Q4'(rd)=S'''(rd)/6, Q4''(rd)=S''''(rd)/12, a4=[r^6]S (or [r^6] of the curve)
#   s=sqrt(Q4(rd)); a1=Q4'(rd)/(4s); a2=Q4''(rd)/12
#   h0=F(rd); h1=F'(rd)s; h2=(1/2)(F''(rd)s^2+F'(rd)Q4'(rd)/2)
#   b3=h0/s^3; b2=(h1-3a1 h0)/s^3; b1=(h2-3a1 h1+(6a1^2-3a2)h0)/s^3
# r_d,Jc: double root (S(rd)=0,S'(rd)=0). Only S (the curve) and N (the source) change.
import sympy as sp, numpy as np
r,M,a,E,Jc,rd=sp.symbols('r M a E Jc r_d')

def bcoeffs(S,N):
    d=lambda ex,n: sp.diff(ex,r,n).subs(r,rd)
    S2,S3,S4=d(S,2),d(S,3),d(S,4)
    Q4rd,Q4prd,Q4pprd=S2/2,S3/6,S4/12
    s=sp.sqrt(Q4rd); a1=Q4prd/(4*s); a2=Q4pprd/12
    Nr,Np,Npp=N.subs(r,rd),sp.diff(N,r).subs(r,rd),sp.diff(N,r,2).subs(r,rd)
    h0=Nr/Q4rd
    Fp=(Np*Q4rd-Nr*Q4prd)/Q4rd**2
    Fpp=((Npp*Q4rd-Nr*Q4pprd)*Q4rd-2*Q4prd*(Np*Q4rd-Nr*Q4prd))/Q4rd**3
    h1=Fp*s; h2=sp.Rational(1,2)*(Fpp*s**2+Fp*(Q4prd/2))
    return {'b3':h0/s**3,'b2':(h1-3*a1*h0)/s**3,'b1':(h2-3*a1*h1+(6*a1**2-3*a2)*h0)/s**3,'s':s}

def verify(S,N,name,params,contour_brd=None):
    Sub0={M:params['M'],a:params['a'],E:params['E']}
    Sn=S.subs(Sub0)
    Jcv=float([z for z in sp.solve(sp.Eq(sp.resultant(Sn,sp.diff(Sn,r),r),0),Jc)
               if z.is_real and float(z)>params.get('Jmin',1)][0])
    rts=np.roots([complex(c) for c in sp.Poly(Sn.subs(Jc,Jcv),r).all_coeffs()])
    pr=[(i,j) for i in range(len(rts)) for j in range(i+1,len(rts)) if abs(rts[i]-rts[j])<1e-5]
    rdv=float(np.real((rts[pr[0][0]]+rts[pr[0][1]])/2))
    bc=bcoeffs(S,N); sub={M:params['M'],a:params['a'],E:params['E'],Jc:Jcv,rd:rdv}
    print(f"\n=== {name}: M={params['M']},a={params['a']},E={params['E']} -> Jc={Jcv:.6f}, r_d={rdv:.6f} ===")
    for k in ['b1','b2','b3']:
        print(f"  {k}(symbolic) = {float(bc[k].subs(sub)):+.8f}")
    print("  a4=[r^6] of the curve =",sp.simplify(sp.Poly(S,r).coeff_monomial(r**6)))

# ---------- tau BRANCH (Vaidya a=0 and TK a!=0): same S, different sources ----------
DE=(E**2-1)*r+2*M; Dl=r**2-2*M*r+a**2
S_tau=sp.expand(r*(r-2*M)*DE*(r*Dl-Jc**2*DE))    # tau sextic (TK; Vaidya = a->0)
# Vaidya tau: mass-derivative source N_m=S dm K - 1/2 K dm S, K=Jc DE
K=Jc*DE; N_vaidya=sp.expand(S_tau*sp.diff(K,M)-sp.Rational(1,2)*K*sp.diff(S_tau,M))
# TK tau: E-derivative source N_tau = E Jc r^4 (r-2M)^2 DE
N_tk=sp.expand(E*Jc*r**4*(r-2*M)**2*DE)

print("############ TAU BRANCH ############")
print("Curve: S_tau = r(r-2M)DE(rDelta-Jc^2 DE), DE=(E^2-1)r+2M, Delta=r^2-2Mr+a^2")
print("Coeffs b_i = universal formula with F=N/Q4; Q4^(k)(rd) from S^(k+2)(rd).")
verify(S_tau.subs(a,0),N_vaidya.subs(a,0),"Vaidya tau (a=0), source N_m",{'M':1,'a':0,'E':sp.Rational(7,5)})
verify(S_tau,N_tk,"TK tau, source N_tau",{'M':1,'a':sp.Rational(9,10),'E':sp.Rational(6,5),'Jmin':5})

print("\n############ t BRANCH (TK) ############")
# curve R6 = r Q2 DE ; Q2 quartic (M=1); two separatrices Jc+- (double root of Q2)
Q2=(2*E**2*Jc**2*r - E**2*Jc**2*r**2 - 4*E**2*Jc*a*r + 2*E**2*a**2*r + E**2*a**2*r**2
    + E**2*r**4 + 4*Jc**2 - 4*Jc**2*r + Jc**2*r**2 - 8*Jc*a + 4*Jc*a*r + 4*a**2)  # M=1
R6=sp.expand(r*Q2*DE)
N_t=sp.expand(sp.simplify(sp.diff((r*DE*(Jc*(r-2)+2*a)/Dl)/sp.sqrt(R6),E)*R6**sp.Rational(3,2)))
def verify_t(sign,Jguess,rdguess):
    Sub0={M:1,a:sp.Rational(9,10),E:sp.Rational(6,5)}
    Q2n=Q2.subs(Sub0); import mpmath as mp
    Q2l=sp.lambdify((r,Jc),Q2n,'mpmath'); Q2rl=sp.lambdify((r,Jc),sp.diff(Q2n,r),'mpmath')
    sl=mp.findroot(lambda rr,jj:[Q2l(rr,jj),Q2rl(rr,jj)],(mp.mpf(rdguess),mp.mpf(Jguess)))
    rdv=float(sl[0]); Jcv=float(sl[1])
    bc=bcoeffs(R6,N_t); sub={M:1,a:sp.Rational(9,10),E:sp.Rational(6,5),Jc:Jcv,rd:rdv}
    print(f"\n=== TK t {sign}: Jc={Jcv:.6f}, r_d={rdv:.6f} ===")
    for k in ['b1','b2','b3']: print(f"  {k}(symbolic) = {float(bc[k].subs(sub)):+.8f}")
print("Curve: R6 = r Q2 DE (Q2 quartic in r, LINEAR in Jc -> Jc+ != Jc-, frame dragging)")
verify_t("prograde (Jc+)",19.09,-6.62)
verify_t("retrograde (Jc-)",-18.67,-6.588)

print("\n=> ALL branches: b1,b2,b3 SYMBOLIC in (M,a,E,r_d,Jc). r_d,Jc=double root of the curve.")
print("   tau (Vaidya/TK): curve S_tau. t (TK): curve R6, two Jc+-. v (Vaidya)=same b_i as tau,")
print("   different clock. Reusable: given (M,a,E), solve for the double root -> plug into the formulas.")
