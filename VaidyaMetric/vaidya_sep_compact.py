# COMPATTIFICAZIONE di delta phi|_sep via simmetria +/- (involuzione z->-z, r(z) pari).
# Blocchi PARI: Z_a=ζ(z-a)-ζ(z+a), P_a=℘(z-a)+℘(z+a), Pp_a=℘'(z-a)-℘'(z+a).
# R (pari) = C0 + b1 Z_zd + b2 P_zd - (b3/2) Pp_zd + b2h0 ℘(z) + b2h2 ℘(z-iw)   [6 blocchi]
# eta'(pari)= Ce + e1d Z_zd + e2i P_zinf + e1i Z_zinf                            [4 blocchi]
# delta phi = sum_{blkR,blk eta'} coeff * J[blkR,blk eta'],  J = iterato length-2.
import numpy as np, mpmath as mp, sympy as sp
from scipy.integrate import quad
mp.mp.dps=40
E=1.4; m=1.0; r0=12.0
r,J=sp.symbols('r J'); Es=sp.Rational(7,5)
Ssym=sp.expand(r*(r-2)*((Es**2-1)*r+2)*(r**2*(r-2)-J**2*((Es**2-1)*r+2)))
Jc=float([s for s in sp.solve(sp.Eq(sp.resultant(Ssym,sp.diff(Ssym,r),r),0),J) if s.is_real and float(s)>1][0])
def pS(mv,Jv):
    DE=np.array([E**2-1,2*mv]); p=np.polymul(np.polymul([1,0],[1,-2*mv]),DE)
    return np.polymul(p,np.polysub(np.polymul([1,0,0],[1,-2*mv]),Jv**2*DE))
Sc=pS(m,Jc); dmS=(pS(m+1e-6,Jc)-pS(m-1e-6,Jc))/2e-6
K=np.polymul([Jc],[E**2-1,2.0]); Nm=np.polysub(np.polymul(Sc,np.array([2*Jc])),0.5*np.polymul(K,dmS))
rts=np.roots(Sc); pr=[(i,j) for i in range(6) for j in range(i+1,6) if abs(rts[i]-rts[j])<1e-6]
rd=float(np.real((rts[pr[0][0]]+rts[pr[0][1]])/2)); Q,_=np.polydiv(Sc,np.polymul([1,-rd],[1,-rd]))
er=np.sort(np.real(np.roots(Q))); a4=Q[0]; e1r,e2r,e3r,e4r=er
Qp=np.polyder(Q); Qpp=np.polyder(Qp); Nmp=np.polyder(Nm); Nmpp=np.polyder(Nmp)
def Q4(x): return np.polyval(Q,x)
k2=((e3r-e2r)*(e4r-e1r))/((e4r-e2r)*(e3r-e1r)); pref=2/mp.sqrt((e4r-e2r)*(e3r-e1r))/mp.sqrt(a4)
om1=mp.mpf(float(pref*mp.ellipk(k2))); w_im=float(pref*mp.ellipk(1-k2)); tau=mp.mpc(0,w_im)/om1; qn=mp.exp(mp.pi*1j*tau)
L1=lambda u: mp.jtheta(1,u,qn); L1p=lambda u: mp.jtheta(1,u,qn,1); L1pp=lambda u: mp.jtheta(1,u,qn,2); th1p0=L1p(0)
eta1=-(mp.pi**2/(12*om1))*(mp.jtheta(1,0,qn,3)/th1p0)
def wzet(z): u=mp.pi*z/(2*om1); return eta1*z/om1+(mp.pi/(2*om1))*(L1p(u)/L1(u))
def wp(z):   u=mp.pi*z/(2*om1); rr=L1p(u)/L1(u); return -eta1/om1-(mp.pi/(2*om1))**2*(L1pp(u)/L1(u)-rr**2)
def wpp(z):  return (wp(z+mp.mpf('1e-12'))-wp(z-mp.mpf('1e-12')))/mp.mpf('2e-12')
z_inf=float(mp.re(mp.quad(lambda x:1/mp.sqrt(Q4(x)),[e4r,mp.inf]))); sa=float(mp.sqrt(a4))
c_r=float(mp.re(e4r-(2/sa)*wzet(z_inf))); z_d=z_inf+float(mp.re(mp.quad(lambda x:1/mp.sqrt(Q4(x)),[-mp.inf,rd])))
def zr(rv): return float(quad(lambda x:1/np.sqrt(Q4(x)),e4r,rv,limit=400)[0]); 
z0f=mp.mpf(zr(r0))
# residui analitici
s=np.sqrt(Q4(rd)); a1=np.polyval(Qp,rd)/(4*s); a2=np.polyval(Qpp,rd)/12
F=lambda x: np.polyval(Nm,x)/np.polyval(Q,x)
Fp=lambda x:(np.polyval(Nmp,x)*np.polyval(Q,x)-np.polyval(Nm,x)*np.polyval(Qp,x))/np.polyval(Q,x)**2
Fpp=lambda x:((np.polyval(Nmpp,x)*np.polyval(Q,x)-np.polyval(Nm,x)*np.polyval(Qpp,x))*np.polyval(Q,x)
   -2*np.polyval(Qp,x)*(np.polyval(Nmp,x)*np.polyval(Q,x)-np.polyval(Nm,x)*np.polyval(Qp,x)))/np.polyval(Q,x)**3
h0=F(rd); h1=Fp(rd)*s; h2=0.5*(Fpp(rd)*s**2+Fp(rd)*(np.polyval(Qp,rd)/2))
b1zd=(h2-3*a1*h1+(6*a1**2-3*a2)*h0)/s**3; b2zd=(h1-3*a1*h0)/s**3; b3zd=h0/s**3
def b2h(ei): return np.polyval(Nm,ei)/((ei-rd)**3)*(4/np.polyval(Qp,ei)**2)
e1_zd=(rd**3-2*rd**2)/s; B=c_r+(1/sa)*float(mp.re(wzet(2*z_inf))); Aa=-1/sa
e1_zi=Aa*(2*B+rd-2); e2_zi=Aa**2
# blocchi PARI
Zf =lambda a: (lambda z: wzet(z-a)-wzet(z+a))
Pf =lambda a: (lambda z: wp(z-a)+wp(z+a))
Ppf=lambda a: (lambda z: wpp(z-a)-wpp(z+a))
one=lambda z: mp.mpf(1); wp0=lambda z: wp(z); wpiw=lambda z: wp(z-1j*mp.mpf(w_im))
# liste (coeff, funzione, nome)
zdm=mp.mpf(z_d); zim=mp.mpf(z_inf)
Rblk=[('C0',one),(b1zd,Zf(zdm)),(b2zd,Pf(zdm)),(-b3zd/2,Ppf(zdm)),(b2h(e4r),wp0),(b2h(e3r),wpiw)]
Eblk=[('Ce',one),(e1_zd,Zf(zdm)),(e2_zi,Pf(zim)),(e1_zi,Zf(zim))]
# fissa C0,Ce
def r_of_z(z): return c_r-(1/sa)*(wzet(z-z_inf)-wzet(z+z_inf))
def Rexact(z): rr=r_of_z(z); return np.polyval(Nm,rr)/((rr-rd)**3*np.polyval(Q,complex(rr)))
def Epexact(z): rr=r_of_z(z); return (rr**3-2*rr**2)/(rr-rd)
zt=mp.mpf('0.19')
C0=complex(Rexact(zt)-sum(c*f(zt) for c,f in Rblk if not isinstance(c,str)))
Ce=complex(Epexact(zt)-sum(c*f(zt) for c,f in Eblk if not isinstance(c,str)))
Rblk[0]=(C0,one); Eblk[0]=(Ce,one)
print("blocchi R:",len(Rblk)," blocchi eta':",len(Eblk)," prodotti:",len(Rblk)*len(Eblk))
print("check R (blocchi pari):",float(max(abs(complex(Rexact(mp.mpf(t))-sum(c*f(mp.mpf(t)) for c,f in Rblk))) for t in [0.15,0.22])))
print("check eta' (blocchi pari):",float(max(abs(complex(Epexact(mp.mpf(t))-sum(c*f(mp.mpf(t)) for c,f in Eblk))) for t in [0.15,0.22])))
# J[blkR,blk eta'] iterato + assemblaggio
_ic={}
def inner(f2,t):
    key=(id(f2),float(t))
    if key in _ic: return _ic[key]
    v=mp.quad(lambda ss: f2(ss),[z0f,t]); _ic[key]=v; return v
def Jblk(f1,f2,z): return mp.quad(lambda t: f1(t)*inner(f2,t),[z0f,z],maxdegree=8)
def dphi_compact(z): return sum(c1*c2*Jblk(f1,f2,z) for c1,f1 in Rblk for c2,f2 in Eblk)
def sQ(x): return np.sqrt(Q4(x))
def U(x,kk): return quad(lambda t:t**kk/((t-rd)*sQ(t)),r0,x,limit=200)[0]
def eta(x): return U(x,3)-2*U(x,2)
def dmF(x):
    def Fm(mv): return np.polyval(np.polymul([Jc],[E**2-1,2*mv]),x)/np.sqrt(np.polyval(pS(mv,Jc),x))
    return (Fm(m+1e-6)-Fm(m-1e-6))/2e-6
def dphi_direct(x): return quad(lambda t: dmF(t)*eta(t), r0, x, limit=200)[0]
print("\n=== delta phi COMPATTO (24 prodotti di blocchi pari) vs DIRETTO ===")
for rr in [11.0,10.0,9.2]:
    z=mp.mpf(zr(rr)); dc=complex(dphi_compact(z)).real; dd=dphi_direct(rr)
    print(f"  r={rr:4.1f}  compatto={dc:+.8f}  diretto={dd:+.8f}  diff={abs(dc-dd):.1e}")
nz=sum(1 for c1,f1 in Rblk for c2,f2 in Eblk if abs(c1*c2)>1e-9)
print(f"\ntermini non nulli: {nz} (era 63).  Blocchi PARI: R={len(Rblk)}, eta'={len(Eblk)}.")
