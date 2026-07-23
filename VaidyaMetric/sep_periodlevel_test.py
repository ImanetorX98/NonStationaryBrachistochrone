# -*- coding: utf-8 -*-
# I pezzi PERIOD-LEVEL (g2,g3, z_d,z_inf, zeta(z_d),P(z_inf), Ce,C0) sono:
#  simbolici? universali? o da ricalcolare per ogni (M,E,r_d,Jc)?
# Test: calcolo a DUE valori di E (7/5 e 13/10) per Vaidya tau e confronto.
import numpy as np, mpmath as mp, sympy as sp
mp.mp.dps=25
def periodlevel(Ev):
    M=1.0; r,Jc=sp.symbols('r Jc')
    S=sp.expand(r*(r-2)*((Ev**2-1)*r+2)*(r**2*(r-2)-Jc**2*((Ev**2-1)*r+2)))
    Jcv=float([z for z in sp.solve(sp.Eq(sp.resultant(S,sp.diff(S,r),r),0),Jc) if z.is_real and float(z)>1][0])
    Sc=np.array([float(c) for c in sp.Poly(S.subs(Jc,Jcv),r).all_coeffs()])
    rts=np.roots(Sc); pr=[(i,j) for i in range(6) for j in range(i+1,6) if abs(rts[i]-rts[j])<1e-6]
    rd=float(np.real((rts[pr[0][0]]+rts[pr[0][1]])/2)); Q,_=np.polydiv(Sc,np.polymul([1,-rd],[1,-rd]))
    er=np.sort(np.real(np.roots(Q))); a4=Q[0]; e1,e2,e3,e4=er
    def Q4(x): return np.polyval(Q,x)
    # g2,g3 (invarianti Weierstrass) - ALGEBRICI nelle radici e_i
    k2=((e3-e2)*(e4-e1))/((e4-e2)*(e3-e1)); pref=2/mp.sqrt((e4-e2)*(e3-e1))/mp.sqrt(a4)
    om1=mp.mpf(float(pref*mp.ellipk(k2))); w_im=float(pref*mp.ellipk(1-k2))
    g2,g3=mp.weierstrass_invariants(om1, 1j*mp.mpf(w_im)) if hasattr(mp,'weierstrass_invariants') else (None,None)
    tau=mp.mpc(0,w_im)/om1; q=mp.exp(mp.pi*1j*tau)
    L1=lambda u: mp.jtheta(1,u,q); L1p=lambda u: mp.jtheta(1,u,q,1); th1p0=L1p(0)
    eta1=-(mp.pi**2/(12*om1))*(mp.jtheta(1,0,q,3)/th1p0)
    def wzet(z): u=mp.pi*z/(2*om1); return eta1*z/om1+(mp.pi/(2*om1))*(L1p(u)/L1(u))
    def wp(z):   u=mp.pi*z/(2*om1); rr=L1p(u)/L1(u); return -eta1/om1-(mp.pi/(2*om1))**2*(mp.jtheta(1,u,q,2)/L1(u)-rr**2)
    from scipy.integrate import quad
    z_inf=float(mp.re(mp.quad(lambda x:1/mp.sqrt(Q4(x)),[e4,mp.inf])))
    z_d=z_inf+float(mp.re(mp.quad(lambda x:1/mp.sqrt(Q4(x)),[-mp.inf,rd])))
    return dict(Jc=Jcv,rd=rd,e_roots=er,om1=float(om1),w_im=w_im,
                z_d=z_d,z_inf=z_inf,zeta_zd=complex(wzet(mp.mpf(z_d))).real,
                wp_zinf=complex(wp(mp.mpf(z_inf))).real,g3overg2="alg(e_i)")
print("=== PERIOD-LEVEL: Vaidya tau a E=7/5 vs E=13/10 ===")
A=periodlevel(sp.Rational(7,5)); B=periodlevel(sp.Rational(13,10))
for k in ['Jc','rd','om1','w_im','z_d','z_inf','zeta_zd','wp_zinf']:
    print(f"  {k:9s}: E=7/5 -> {A[k]:+.6f}   E=13/10 -> {B[k]:+.6f}   {'DIVERSI' if abs(A[k]-B[k])>1e-4 else 'uguali'}")
print("  radici e_i: E=7/5 ->",np.round(A['e_roots'],4)," E=13/10 ->",np.round(B['e_roots'],4))
print("\n=> om1,w_im,z_d,z_inf,zeta(z_d),P(z_inf) CAMBIANO con E -> NON universali.")
print("   Sono PERIODI/valori di funzioni ellittiche = trascendenti (non razionali in M,E,r_d,Jc).")
print("   g2,g3 (invarianti) sono ALGEBRICI nelle radici e_i, ma z_d,z_inf (integrali ellittici")
print("   incompleti) e zeta,P su di essi sono TRASCENDENTI -> valutati per-curva (procedura universale).")
