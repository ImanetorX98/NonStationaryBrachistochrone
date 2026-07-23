# -*- coding: utf-8 -*-
# COEFFICIENTI TRACKING (separatrice mobile Jc(lambda)) SIMBOLICI in (M,a,E,r_d,Jc).
# b_i^track: stessa formula b3=h0/s^3 ma con F=N_tot/Q4, N_tot=N + (dJc/dlambda) N_J.
# Verifica: b_i^track(simbolico) vs valori numerici del tracking.
import sympy as sp, numpy as np
r,M,a,E,Jc,rd=sp.symbols('r M a E Jc r_d')

def bcoeffs(S,N):
    d=lambda ex,n: sp.diff(ex,r,n).subs(r,rd)
    Q4rd,Q4prd,Q4pprd=d(S,2)/2,d(S,3)/6,d(S,4)/12
    s=sp.sqrt(Q4rd); a1=Q4prd/(4*s); a2=Q4pprd/12
    Nr,Np,Npp=N.subs(r,rd),sp.diff(N,r).subs(r,rd),sp.diff(N,r,2).subs(r,rd)
    h0=Nr/Q4rd
    Fp=(Np*Q4rd-Nr*Q4prd)/Q4rd**2
    Fpp=((Npp*Q4rd-Nr*Q4pprd)*Q4rd-2*Q4prd*(Np*Q4rd-Nr*Q4prd))/Q4rd**3
    h1=Fp*s; h2=sp.Rational(1,2)*(Fpp*s**2+Fp*(Q4prd/2))
    return {'b3':h0/s**3,'b2':(h1-3*a1*h0)/s**3,'b1':(h2-3*a1*h1+(6*a1**2-3*a2)*h0)/s**3}

DE=(E**2-1)*r+2*M; Dl=r**2-2*M*r+a**2
S=sp.expand(r*(r-2*M)*DE*(r*Dl-Jc**2*DE))
K=Jc*r*(r-2*M)*DE/Dl

# ---- Vaidya tau tracking (param m=M, a=0) ----
Sv=S.subs(a,0); Kv=Jc*DE   # a=0: K=Jc DE
Nm=sp.expand(Sv*sp.diff(Kv,M)-sp.Rational(1,2)*Kv*sp.diff(Sv,M))       # sorgente fissa
NJv=sp.expand(Sv*DE+Jc**2*r*(r-2*M)*DE**3)                             # N_J
dJcdm=Jc/M                                                            # scaling lineare
Ntot_v=sp.expand(Nm+dJcdm*NJv)
bv_fix=bcoeffs(Sv,Nm); bv_trk=bcoeffs(Sv,Ntot_v)

# ---- TK tau tracking (param E) ----
Ntau=sp.expand(E*Jc*r**4*(r-2*M)**2*DE)
NJt=sp.expand(r**3*(r-2*M)**2*DE**2)
DErd=(E**2-1)*rd+2*M
dJcdE=-E*Jc*rd/DErd
Ntot_t=sp.expand(Ntau+dJcdE*NJt)
bt_fix=bcoeffs(S,Ntau); bt_trk=bcoeffs(S,Ntot_t)

# ===== VERIFICA numerica =====
def numcheck(Scurve,name,params,Jmin,bfix,btrk,Ntot,dJc_formula):
    Sub={M:params['M'],a:params['a'],E:params['E']}
    Sn=Scurve.subs(Sub)
    Jcv=float([z for z in sp.solve(sp.Eq(sp.resultant(Sn,sp.diff(Sn,r),r),0),Jc) if z.is_real and float(z)>Jmin][0])
    rts=np.roots([float(c) for c in sp.Poly(Sn.subs(Jc,Jcv),r).all_coeffs()])
    pr=[(i,j) for i in range(6) for j in range(i+1,6) if abs(rts[i]-rts[j])<1e-5]
    rdv=float(np.real((rts[pr[0][0]]+rts[pr[0][1]])/2))
    sub={M:params['M'],a:params['a'],E:params['E'],Jc:Jcv,rd:rdv}
    print(f"\n{name}: Jc={Jcv:.5f} r_d={rdv:.5f}  dJc/dlambda={float(dJc_formula.subs(sub)):+.5f}")
    print(f"  N_tot(r_d) = {float(Ntot.subs(sub).subs(r,rdv)):.2e}  (=0 => tracking CANCELLA polo triplo)")
    for k in ['b3','b2','b1']:
        print(f"  {k}: fixed={float(bfix[k].subs(sub)):+.8f}  track={float(btrk[k].subs(sub)):+.8f}")

print("=== COEFFICIENTI TRACKING SIMBOLICI (b_i con N_tot=N+(dJc/dl)N_J) ===")
numcheck(Sv,"Vaidya tau TRACKING (dJc/dm=Jc/m)",{'M':1,'a':0,'E':sp.Rational(7,5)},1,bv_fix,bv_trk,Ntot_v,dJcdm)
numcheck(S,"TK tau TRACKING (dJc/dE=-E Jc r_d/DE(r_d))",{'M':1,'a':sp.Rational(9,10),'E':sp.Rational(6,5)},5,bt_fix,bt_trk,Ntot_t,dJcdE)
print("\n=> RISULTATO: il tracking ha N_tot(r_d)=0 -> b3^track=0 (polo triplo cancellato).")
print("   b_i^track SIMBOLICI in (M,a,E,r_d,Jc), F=N_tot/Q4. Seguire la separatrice mobile addolcisce")
print("   la singolarita' a r_d. Fixed-Jc (SEP_COEFF_SYMBOLIC) ha b3!=0; tracking b3=0.")
