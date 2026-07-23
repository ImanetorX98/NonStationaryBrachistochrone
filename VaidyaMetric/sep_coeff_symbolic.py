# -*- coding: utf-8 -*-
# COEFFICIENTI SEPARATRICE pienamente SIMBOLICI in (M,E,r_d,Jc) [Vaidya tau, a=0].
# b1,b2,b3 (residui R al polo triplo z_d) e residui clock, via DERIVATE DI S a r_d:
#   Q4(rd)=S''(rd)/2, Q4'(rd)=S'''(rd)/6, Q4''(rd)=S''''(rd)/12, a4=[r^6]S.
#   s=sqrt(Q4(rd)); a1=Q4'(rd)/(4s); a2=Q4''(rd)/12; F=N/Q4; h0,h1,h2.
#   b3=h0/s^3; b2=(h1-3a1 h0)/s^3; b1=(h2-3a1 h1+(6a1^2-3a2)h0)/s^3.
# r_d,Jc: doppia radice => S(rd)=0, S'(rd)=0 (2 vincoli in M,E). Riutilizzabile.
import sympy as sp, numpy as np

r,M,E,Jc,rd=sp.symbols('r M E Jc r_d')
# Vaidya tau sestica (a=0): S=r(r-2M)DE(r^2(r-2M)-J^2 DE), DE=(E^2-1)r+2M
DE=(E**2-1)*r+2*M
S=sp.expand(r*(r-2*M)*DE*(r**2*(r-2*M)-Jc**2*DE))
# sorgente fissa N_m = S dm K - 1/2 K dm S, K=Jc DE  (m=M parametro)
K=Jc*DE
Nm=sp.expand(S*sp.diff(K,M)-sp.Rational(1,2)*K*sp.diff(S,M))
# derivate di S e N a r_d (simboliche)
def datrd(expr): return expr.subs(r,rd)
S2=datrd(sp.diff(S,r,2)); S3=datrd(sp.diff(S,r,3)); S4=datrd(sp.diff(S,r,4))
a4=sp.Poly(S,r).coeff_monomial(r**6)
Q4rd=S2/2; Q4prd=S3/6; Q4pprd=S4/12
s=sp.sqrt(Q4rd)
a1=Q4prd/(4*s); a2=Q4pprd/12
Nrd=datrd(Nm); Nprd=datrd(sp.diff(Nm,r)); Npprd=datrd(sp.diff(Nm,r,2))
# F=N/Q4; F',F'' a rd via quoziente (Q4(rd)=Q4rd, Q4'(rd)=Q4prd, Q4''(rd)=Q4pprd)
h0=Nrd/Q4rd
Fp=(Nprd*Q4rd-Nrd*Q4prd)/Q4rd**2
Fpp=((Npprd*Q4rd-Nrd*Q4pprd)*Q4rd-2*Q4prd*(Nprd*Q4rd-Nrd*Q4prd))/Q4rd**3
h1=Fp*s; h2=sp.Rational(1,2)*(Fpp*s**2+Fp*(Q4prd/2))
b3=h0/s**3; b2=(h1-3*a1*h0)/s**3; b1=(h2-3*a1*h1+(6*a1**2-3*a2)*h0)/s**3
# clock residui (tau): eta'=(r^3-2M r^2)/(r-r_d); e1_zd=(rd^3-2M rd^2)/s
e1_zd=(rd**3-2*M*rd**2)/s
e2_zi=1/a4
print("=== COEFF SEPARATRICE Vaidya tau, SIMBOLICI in (M,E,r_d,Jc) ===")
print("s^2=Q4(rd)=S''(rd)/2 =",sp.simplify(Q4rd))
print("a4=[r^6]S =",sp.simplify(a4))
print("b3 = h0/s^3, h0 = N(rd)/Q4(rd) =",sp.simplify(h0))
print("e1_zd = (rd^3-2M rd^2)/s ; e2_zi=1/a4 =",sp.simplify(e2_zi))
print("(b1,b2 espliciti ma lunghi; struttura simbolica completa in (M,E,rd,Jc))")

# ===== VERIFICA numerica vs valori noti (M=1,E=7/5) =====
Ev=sp.Rational(7,5); Mv=1
Sn=S.subs({M:Mv,E:Ev})
Jcval=float([z for z in sp.solve(sp.Eq(sp.resultant(Sn,sp.diff(Sn,r),r),0),Jc) if z.is_real and float(z)>1][0])
Sc=sp.lambdify(r,Sn.subs(Jc,Jcval),'numpy'); rts=np.roots([float(c) for c in sp.Poly(Sn.subs(Jc,Jcval),r).all_coeffs()])
pr=[(i,j) for i in range(6) for j in range(i+1,6) if abs(rts[i]-rts[j])<1e-6]
rdval=float(np.real((rts[pr[0][0]]+rts[pr[0][1]])/2))
sub={M:Mv,E:Ev,Jc:Jcval,rd:rdval}
print(f"\n  a M=1,E=7/5: Jc={Jcval:.6f}, r_d={rdval:.6f}")
print(f"  b3(sym)={float(b3.subs(sub)):+.8f}  (atteso contorno ~ valore residui script)")
print(f"  b2(sym)={float(b2.subs(sub)):+.8f}")
print(f"  b1(sym)={float(b1.subs(sub)):+.8f}")
print(f"  e1_zd(sym)={float(e1_zd.subs(sub)):+.8f}   e2_zi=1/a4(sym)={float(e2_zi.subs(sub)):+.8f}")
print("\n=> b1,b2,b3,e_i SIMBOLICI in (M,E,r_d,Jc); r_d,Jc = doppia radice (S(rd)=S'(rd)=0). Riutilizzabile.")
