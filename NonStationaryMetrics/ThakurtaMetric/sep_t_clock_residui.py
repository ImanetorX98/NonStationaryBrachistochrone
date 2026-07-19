# -*- coding: utf-8 -*-
# Residui del CLOCK t (TK t) in forma SIMBOLICA. Clock dt/dr=rho_t/sqrt(R6),
# rho_t=P3+R_Delta/Delta ; in z: etpz=rho_t/(r-r_d). Poli: z_d, orizzonti z(r_pm), z_inf.
#   - z_d:      residuo = rho_t(r_d)/s          (s=sqrt(Q4(r_d)))
#   - z(r_pm):  residuo 3a specie = R_Delta(r_pm)/((r_pm - r_mp)(r_pm - r_d) sqrt(Q4(r_pm)))
#   dove Q4(r_pm)=R6(r_pm)/(r_pm - r_d)^2, R6=r Q2 DE.
# Verifica vs estrazione a contorno (prograda Jc+).
import numpy as np, mpmath as mp, sympy as sp
from scipy.integrate import quad
mp.mp.dps=30
M=1.0; a=0.9; E=1.2
r,J=sp.symbols('r J'); Es=sp.Rational(6,5); asub={}
DE=(Es**2-1)*r+2; Delta=r**2-2*r+sp.Rational(81,100)
Q2=(2*Es**2*J**2*r-Es**2*J**2*r**2-4*Es**2*J*sp.Rational(9,10)*r+2*Es**2*sp.Rational(81,100)*r
    +Es**2*sp.Rational(81,100)*r**2+Es**2*r**4+4*J**2-4*J**2*r+J**2*r**2-8*J*sp.Rational(9,10)
    +4*J*sp.Rational(9,10)*r+4*sp.Rational(81,100))
R6=sp.expand(r*Q2*DE)
Kt=r*DE*(J*(r-2)+2*sp.Rational(9,10))/Delta
rho_t=sp.cancel(sp.together((Es**2*r**3-2*sp.Rational(9,10)*Kt/r)/((r-2)/r)))
P3poly,_=sp.div(sp.Poly(sp.numer(rho_t),r),sp.Poly(sp.denom(rho_t),r))
P3=sp.expand(P3poly.as_expr()); RD=sp.expand(sp.simplify((rho_t-P3)*Delta))
# raffina Jc+,r_d (doppia radice Q2)
Q2n=sp.lambdify((r,J),Q2,'mpmath'); Q2rn=sp.lambdify((r,J),sp.diff(Q2,r),'mpmath')
sol=mp.findroot(lambda rr,jj:[Q2n(rr,jj),Q2rn(rr,jj)],(mp.mpf('-6.62'),mp.mpf('19.089')))
rdv=float(sol[0]); Jc=float(sol[1])
R6c=np.array([float(c) for c in sp.Poly(R6.subs(J,Jc),r).all_coeffs()])
Q,_=np.polydiv(R6c,np.polymul([1,-rdv],[1,-rdv])); a4=Q[0]
def Q4(x): return np.polyval(Q,x)
rp=1+np.sqrt(1-a*a); rm=1-np.sqrt(1-a*a)
RDn=sp.lambdify(r,RD.subs(J,Jc),'numpy'); P3n=sp.lambdify(r,P3.subs(J,Jc),'numpy'); Dn=sp.lambdify(r,Delta,'numpy')
# --- RESIDUI SIMBOLICI (formule) ---
s=np.sqrt(complex(Q4(rdv)))
rho_rd = P3n(rdv)+RDn(rdv)/Dn(rdv)
res_zd = rho_rd/s
# sigma = segno del foglio di sqrt(Q4) agli orizzonti (fissato dalla parametrizzazione orbita)
sig=-1.0
res_rp = sig*RDn(rp)/((rp-rm)*(rp-rdv)*np.sqrt(complex(Q4(rp))))
res_rm = sig*RDn(rm)/((rm-rp)*(rm-rdv)*np.sqrt(complex(Q4(rm))))
print("  [invariante] res(r+)+res(r-) =",np.round((res_rp+res_rm).real,6)," (atteso 2M=%.1f)"%(2*M))
print("=== RESIDUI CLOCK t (TK, prograda Jc+) formule simboliche ===")
print(f"  Jc={Jc:.6f}  r_d={rdv:.6f}  r+={rp:.6f}  r-={rm:.6f}")
print(f"  res z_d = rho_t(r_d)/s = {res_zd:.6f}")
print(f"  res z(r+) = R_D(r+)/((r+-r-)(r+-r_d)sqrtQ4(r+)) = {res_rp:.6f}")
print(f"  res z(r-) = R_D(r-)/((r--r+)(r--r_d)sqrtQ4(r-)) = {res_rm:.6f}")
# --- VERIFICA a contorno in z ---
er=np.sort(np.real(np.roots(Q))); e1,e2,e3,e4=er
k2=((e3-e2)*(e4-e1))/((e4-e2)*(e3-e1)); pref=2/mp.sqrt((e4-e2)*(e3-e1))/mp.sqrt(a4)
om1=mp.mpf(float(pref*mp.ellipk(k2))); w_im=float(pref*mp.ellipk(1-k2)); tau=mp.mpc(0,w_im)/om1; q=mp.exp(mp.pi*1j*tau)
L1=lambda u: mp.jtheta(1,u,q); L1p=lambda u: mp.jtheta(1,u,q,1); th1p0=L1p(0)
eta1=-(mp.pi**2/(12*om1))*(mp.jtheta(1,0,q,3)/th1p0)
def wzet(z): u=mp.pi*z/(2*om1); return eta1*z/om1+(mp.pi/(2*om1))*(L1p(u)/L1(u))
z_inf=float(mp.re(mp.quad(lambda x:1/mp.sqrt(Q4(x)),[e4,mp.inf]))); sa=float(mp.sqrt(a4))
c_r=float(mp.re(e4-(2/sa)*wzet(z_inf))); z_d=z_inf+float(mp.re(mp.quad(lambda x:1/mp.sqrt(Q4(x)),[-mp.inf,rdv])))
iw=1j*mp.mpf(w_im)
def r_of_z(z): return c_r-(1/sa)*(wzet(z-z_inf)-wzet(z+z_inf))
def etpz(z): rr=r_of_z(z); return (P3n(rr)+RDn(rr)/Dn(rr))/(rr-rdv)
def find_zimg(target):
    for zg in [iw/2, iw/2+0.3, iw*0.7, om1+iw/2, 0.5+iw/2, iw/3]:
        try:
            z=mp.findroot(lambda z: r_of_z(z)-target, zg)
            if abs(complex(r_of_z(z))-target)<1e-8: return z
        except Exception: pass
    return None
zrp=find_zimg(rp); zrm=find_zimg(rm)
def lau(fun,aa,eps=5e-4):
    return complex(mp.quad(lambda th: fun(aa+eps*mp.exp(1j*th))*(1j*eps*mp.exp(1j*th)),[0,2*mp.pi])/(2j*mp.pi))
print("\n=== VERIFICA a contorno ===")
print(f"  res z_d:   formula={res_zd:+.6f}  contorno={lau(etpz,mp.mpf(z_d)):+.6f}")
print(f"  res z(r+): formula={res_rp:+.6f}  contorno={lau(etpz,zrp):+.6f}")
print(f"  res z(r-): formula={res_rm:+.6f}  contorno={lau(etpz,zrm):+.6f}")
print("\n=> residui clock t SIMBOLICI: rho_t(r_d)/s (z_d), R_D(r_pm)/((r_pm-r_mp)(r_pm-r_d)sqrtQ4(r_pm)) (orizzonti).")
