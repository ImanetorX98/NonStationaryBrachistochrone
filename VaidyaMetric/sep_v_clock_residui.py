# -*- coding: utf-8 -*-
# Residui del CLOCK v (Vaidya v) in forma SIMBOLICA. Sorgente b_i = Vaidya tau (identica).
# v_z = E r^3/(r-r_d) + r sqrt(Q4)/(r-2M).
#   - termine 1 (E r^3/(r-r_d)): polo semplice a z_d, residuo = E r_d^3/s  (s=sqrt(Q4(r_d)))
#   - termine 2 (r sqrt(Q4)/(r-2M)): r=2M=e3 radice di Q4 -> sqrt(Q4)~kappa sqrt(r-2M),
#     z-z_h = 2 sqrt(r-2M)/kappa -> termine ~ 2r/(z-z_h) -> residuo 3a specie = 2*(2M)=4M a z_h=i w_im.
# Verifica numerica vs estrazione a contorno.
import numpy as np, mpmath as mp, sympy as sp
from scipy.integrate import quad
mp.mp.dps=30
E=1.4; M=1.0
r,Jc,rd=sp.symbols('r Jc r_d')
Es=sp.Rational(7,5)
S=sp.expand(r*(r-2)*((Es**2-1)*r+2)*(r**2*(r-2)-Jc**2*((Es**2-1)*r+2)))
Jcv=float([z for z in sp.solve(sp.Eq(sp.resultant(S,sp.diff(S,r),r),0),Jc) if z.is_real and float(z)>1][0])
Sc=np.array([float(c) for c in sp.Poly(S.subs(Jc,Jcv),r).all_coeffs()])
rts=np.roots(Sc); pr=[(i,j) for i in range(6) for j in range(i+1,6) if abs(rts[i]-rts[j])<1e-6]
rdv=float(np.real((rts[pr[0][0]]+rts[pr[0][1]])/2)); Q,_=np.polydiv(Sc,np.polymul([1,-rdv],[1,-rdv]))
er=np.sort(np.real(np.roots(Q))); a4=Q[0]; e1,e2,e3,e4=er
def Q4(x): return np.polyval(Q,x)
s=np.sqrt(Q4(rdv))
# --- RESIDUI SIMBOLICI (formule) ---
# residuo z_d (termine 1): E r_d^3 / s
res_zd_formula = E*rdv**3/s
# residuo orizzonte z_h (termine 2): 2*(2M) = 4M  (indip. dai e_i!)
res_zh_formula = 2*(2*M)
print("=== RESIDUI CLOCK v (formule simboliche) ===")
print("  residuo a z_d (da E r^3/(r-r_d)):  E r_d^3/s =",res_zd_formula)
print("  residuo a z_h=i w_im (orizzonte):  2*(2M)=4M =",res_zh_formula)
print("  (e3=2M=%.4f e' radice di Q4: Q4(2M)=%.2e ~0)"%(e3,Q4(2.0)))

# --- VERIFICA a contorno in z (Weierstrass) ---
k2=((e3-e2)*(e4-e1))/((e4-e2)*(e3-e1)); pref=2/mp.sqrt((e4-e2)*(e3-e1))/mp.sqrt(a4)
om1=mp.mpf(float(pref*mp.ellipk(k2))); w_im=float(pref*mp.ellipk(1-k2)); tau=mp.mpc(0,w_im)/om1; q=mp.exp(mp.pi*1j*tau)
L1=lambda u: mp.jtheta(1,u,q); L1p=lambda u: mp.jtheta(1,u,q,1); th1p0=L1p(0)
eta1=-(mp.pi**2/(12*om1))*(mp.jtheta(1,0,q,3)/th1p0)
def wzet(z): u=mp.pi*z/(2*om1); return eta1*z/om1+(mp.pi/(2*om1))*(L1p(u)/L1(u))
z_inf=float(mp.re(mp.quad(lambda x:1/mp.sqrt(Q4(x)),[e4,mp.inf]))); sa=float(mp.sqrt(a4))
c_r=float(mp.re(e4-(2/sa)*wzet(z_inf))); z_d=z_inf+float(mp.re(mp.quad(lambda x:1/mp.sqrt(Q4(x)),[-mp.inf,rdv])))
iw=1j*mp.mpf(w_im)
def r_of_z(z): return c_r-(1/sa)*(wzet(z-z_inf)-wzet(z+z_inf))
def sqrtQ4_z(z): return (1/sa)*(wp_(z-z_inf)-wp_(z+z_inf))
def L1pp(u): return mp.jtheta(1,u,q,2)
def wp_(z): u=mp.pi*z/(2*om1); rr=L1p(u)/L1(u); return -eta1/om1-(mp.pi/(2*om1))**2*(L1pp(u)/L1(u)-rr**2)
def vz(z):
    rr=r_of_z(z); return E*rr**3/(rr-rdv)+rr*sqrtQ4_z(z)/(rr-2*M)
def lau(fun,aa,order,eps=1e-3):
    return {n:complex(mp.quad(lambda th: fun(aa+eps*mp.exp(1j*th))*(eps*mp.exp(1j*th))**(n-1)*(1j*eps*mp.exp(1j*th)),[0,2*mp.pi])/(2j*mp.pi)) for n in range(1,order+1)}
c_zd=lau(vz,mp.mpf(z_d),1); c_zh=lau(vz,iw,1)
print("\n=== VERIFICA a contorno ===")
print(f"  res z_d:  formula={res_zd_formula:+.6f}  contorno={c_zd[1].real:+.6f}  diff={abs(res_zd_formula-c_zd[1].real):.1e}")
print(f"  res z_h:  formula={res_zh_formula:+.6f}  contorno={c_zh[1].real:+.6f}  diff={abs(res_zh_formula-c_zh[1].real):.1e}")
print("\n=> se diff~0: residui clock v SIMBOLICI: E r_d^3/s (z_d), 4M (orizzonte). Sorgente b_i = Vaidya tau.")
