# COEFFICIENTI ANALITICI dei dilog ellittici: c_ab = Res_a(R) * Res_b(eta').
# Residui in FORMA CHIUSA dai dati della curva (r_d, Q4^{(k)}(r_d), N_m(r_d), a4, zeta/P a 2z_inf).
# Verifica vs estrazione a contorno (che erano i numeri ~1.00232, 0.28752...).
import numpy as np, mpmath as mp, sympy as sp
from scipy.integrate import quad
mp.mp.dps=25
E=1.4; m=1.0; r0=12.0
r,J=sp.symbols('r J'); Es=sp.Rational(7,5)
Ssym=sp.expand(r*(r-2)*((Es**2-1)*r+2)*(r**2*(r-2)-J**2*((Es**2-1)*r+2)))
Jc=float([s for s in sp.solve(sp.Eq(sp.resultant(Ssym,sp.diff(Ssym,r),r),0),J) if s.is_real and float(s)>1][0])
def poly_S_m(mv,Jv):
    DE=np.array([E**2-1,2*mv]); p=np.polymul(np.polymul([1,0],[1,-2*mv]),DE)
    br=np.polysub(np.polymul([1,0,0],[1,-2*mv]),Jv**2*DE); return np.polymul(p,br)
Sc=poly_S_m(m,Jc); hh=1e-6; dmS=(poly_S_m(m+hh,Jc)-poly_S_m(m-hh,Jc))/(2*hh)
K=np.polymul([Jc],[E**2-1,2.0]); Nm=np.polysub(np.polymul(Sc,np.array([2*Jc])),0.5*np.polymul(K,dmS))
rts=np.roots(Sc); pr=[(i,j) for i in range(6) for j in range(i+1,6) if abs(rts[i]-rts[j])<1e-6]
rd=float(np.real((rts[pr[0][0]]+rts[pr[0][1]])/2)); Q,_=np.polydiv(Sc,np.polymul([1,-rd],[1,-rd]))
er=np.sort(np.real(np.roots(Q))); a4=Q[0]; e1,e2,e3,e4=er
Qp=np.polyder(Q); Qpp=np.polyder(Qp); Nmp=np.polyder(Nm)
def Q4(x): return np.polyval(Q,x)
# ---- Weierstrass (per zeta(2 z_inf)) ----
k2=((e3-e2)*(e4-e1))/((e4-e2)*(e3-e1)); pref=2/mp.sqrt((e4-e2)*(e3-e1))/mp.sqrt(a4)
om1=mp.mpf(float(pref*mp.ellipk(k2))); w_im=float(pref*mp.ellipk(1-k2)); tau=mp.mpc(0,w_im)/om1; q=mp.exp(mp.pi*1j*tau)
L1=lambda u: mp.jtheta(1,u,q); L1p=lambda u: mp.jtheta(1,u,q,1); th1p0=L1p(0)
eta1=-(mp.pi**2/(12*om1))*(mp.jtheta(1,0,q,3)/th1p0)
def wzet(z): u=mp.pi*z/(2*om1); return eta1*z/om1+(mp.pi/(2*om1))*(L1p(u)/L1(u))
z_inf=float(mp.quad(lambda x:1/mp.sqrt(Q4(x)),[e4,mp.inf])); sa=float(mp.sqrt(a4))
c_r=float(mp.re(e4-(2/sa)*wzet(z_inf))); z_d=z_inf+float(mp.quad(lambda x:1/mp.sqrt(Q4(x)),[-mp.inf,rd]))

# =========== RESIDUI ANALITICI ===========
s=np.sqrt(Q4(rd))                       # = dr/dz a z_d
Q4rd=Q4(rd); Q4prd=np.polyval(Qp,rd); Q4pprd=np.polyval(Qpp,rd)
# eta' residui
P3=lambda x: x**3-2*x**2
e1_zd = P3(rd)/s
A=-1/sa; B=c_r+(1/sa)*float(mp.re(wzet(2*z_inf)))
e1_zi = A*(2*B+rd-2); e2_zi = A**2
# R residui al polo triplo z_d:  F=N_m/Q4;  a1=Q4'/(4s), a2=Q4''/12; h0,h1,h2
F  = lambda x: np.polyval(Nm,x)/np.polyval(Q,x)
Fp = lambda x: (np.polyval(Nmp,x)*np.polyval(Q,x)-np.polyval(Nm,x)*np.polyval(Qp,x))/np.polyval(Q,x)**2
Nmpp=np.polyder(Nmp)
Fpp= lambda x: ( (np.polyval(Nmpp,x)*np.polyval(Q,x)-np.polyval(Nm,x)*np.polyval(Qpp,x))*np.polyval(Q,x)
                 -2*np.polyval(Qp,x)*(np.polyval(Nmp,x)*np.polyval(Q,x)-np.polyval(Nm,x)*np.polyval(Qp,x)) )/np.polyval(Q,x)**3
a1=Q4prd/(4*s); a2=Q4pprd/12
h0=F(rd); h1=Fp(rd)*s; h2=0.5*(Fpp(rd)*s**2+Fp(rd)*(Q4prd/2))
b3_zd=h0/s**3
b2_zd=(h1-3*a1*h0)/s**3
b1_zd=(h2-3*a1*h1+(6*a1**2-3*a2)*h0)/s**3

print("=== RESIDUI ANALITICI vs CONTORNO ===")
# contorno per confronto
def r_of_z(z): return c_r-(1/sa)*(wzet(z-z_inf)-wzet(z+z_inf))
def Rf(z): rr=r_of_z(z); return np.polyval(Nm,rr)/((rr-rd)**3*np.polyval(Q,complex(rr)))
def etap(z): rr=r_of_z(z); return (rr**3-2*rr**2)/(rr-rd)
def lau(fun,a,order,eps=1e-3):
    return {n:complex(mp.quad(lambda th: fun(a+eps*mp.exp(1j*th))*(eps*mp.exp(1j*th))**(n-1)*(1j*eps*mp.exp(1j*th)),[0,2*mp.pi])/(2j*mp.pi)) for n in range(1,order+1)}
cR=lau(Rf,z_d,3); cE_zd=lau(etap,z_d,1); cE_zi=lau(etap,z_inf,2)
print(f"eta' Res z_d:  analitico={e1_zd:+.8f}  contorno={cE_zd[1].real:+.8f}  diff={abs(e1_zd-cE_zd[1].real):.1e}")
print(f"eta' Res z_inf(b1): analitico={e1_zi:+.8f}  contorno={cE_zi[1].real:+.8f}  diff={abs(e1_zi-cE_zi[1].real):.1e}")
print(f"eta' Res z_inf(b2=1/a4): analitico={e2_zi:+.8f}  contorno={cE_zi[2].real:+.8f}  diff={abs(e2_zi-cE_zi[2].real):.1e}")
print(f"R Res z_d b1: analitico={b1_zd:+.8f}  contorno={cR[1].real:+.8f}  diff={abs(b1_zd-cR[1].real):.1e}")
print(f"R Res z_d b2: analitico={b2_zd:+.8f}  contorno={cR[2].real:+.8f}  diff={abs(b2_zd-cR[2].real):.1e}")
print(f"R Res z_d b3: analitico={b3_zd:+.8f}  contorno={cR[3].real:+.8f}  diff={abs(b3_zd-cR[3].real):.1e}")
print("\n=== COEFF DILOG c_ab = b1_a(R) * e1_b(eta') (analitici) ===")
for sa_,bR in [('+z_d',b1_zd),('-z_d',-b1_zd)]:
    for sb_,eE in [('+z_d',e1_zd),('-z_d',-e1_zd),('+z_inf',e1_zi),('-z_inf',-e1_zi)]:
        print(f"  D(a={sa_}, b={sb_})  coeff = {bR*eE:+.6f}")

# =========== VALORI DEI DILOG ELLITTICI D(a,b) (le funzioni speciali stesse) ===========
import mpmath as mp
def wsig(z): u=mp.pi*z/(2*om1); return (2*om1/mp.pi)*mp.exp(eta1*z**2/(2*om1))*L1(u)/th1p0
def zr(rv): return float(quad(lambda x:1/np.sqrt(Q4(x)),e4,rv,limit=400)[0])
z0f=mp.mpf(zr(r0)); ztest=mp.mpf(zr(10.0))
def D(a,b,z): return complex(mp.quad(lambda t: mp.log(wsig(t-a))*wzet(t-b),[z0f,z]))
print("\n=== DILOG ELLITTICI D(a,b)=int lnσ(z-a)ζ(z-b)dz  (a r=10) ===")
print("  (= Brown-Levin Gamma-tilde(1,1) / length-2 Kronecker, = iterated_integral GiNaC 1e-18)")
for (na,a) in [('z_d',z_d),('-z_d',-z_d)]:
    for (nb,b) in [('z_d',z_d),('-z_d',-z_d),('z_inf',z_inf),('-z_inf',-z_inf)]:
        print(f"  D({na:5s},{nb:6s}) = {D(a,mp.mpf(b) if not isinstance(b,mp.mpf) else b,ztest).real:+.8f}")
