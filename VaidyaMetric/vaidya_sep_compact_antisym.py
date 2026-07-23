# COMPATTIFICAZIONE (a): forma ANTISIMMETRICA (shuffle).
# delta phi = 1/2 G~ eta + 1/2 sum_{i<j} (c_i d_j - c_j d_i) A[e_i,e_j],
#   A[g,h]=int(g*Ph - h*Pg)dz  (antisimmetrico=dilog ellittico genuino);  Pg=int_{z0}^z g.
# e = 8 funzioni PARI comuni; c=coeff in R, d=coeff in eta'. R×R e eta'×eta' -> 0 automatico.
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
def zr(rv): return float(quad(lambda x:1/np.sqrt(Q4(x)),e4r,rv,limit=400)[0])
z0f=mp.mpf(zr(r0))
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
zdm=mp.mpf(z_d); zim=mp.mpf(z_inf)
Zf =lambda a: (lambda z: wzet(z-a)-wzet(z+a))
Pf =lambda a: (lambda z: wp(z-a)+wp(z+a))
Ppf=lambda a: (lambda z: wpp(z-a)-wpp(z+a))
one=lambda z: mp.mpf(1)
# lista unificata di 8 funzioni; c=coeff R, d=coeff eta'
efun=[one, Zf(zdm), Pf(zdm), Ppf(zdm), (lambda z: wp(z)), (lambda z: wp(z-1j*mp.mpf(w_im))), Pf(zim), Zf(zim)]
enam=['1','Z_zd','P_zd','Pp_zd','wp_0','wp_iw','P_zi','Z_zi']
def r_of_z(z): return c_r-(1/sa)*(wzet(z-z_inf)-wzet(z+z_inf))
def Rexact(z): rr=r_of_z(z); return np.polyval(Nm,rr)/((rr-rd)**3*np.polyval(Q,complex(rr)))
def Epexact(z): rr=r_of_z(z); return (rr**3-2*rr**2)/(rr-rd)
zt=mp.mpf('0.19')
cvec=[None,b1zd,b2zd,-b3zd/2,b2h(e4r),b2h(e3r),0.0,0.0]        # coeff R (indice0=C0)
dvec=[None,e1_zd,0.0,0.0,0.0,0.0,e2_zi,e1_zi]                  # coeff eta' (indice0=Ce)
C0=complex(Rexact(zt)-sum(cvec[i]*efun[i](zt) for i in range(1,8))); cvec[0]=C0
Ce=complex(Epexact(zt)-sum(dvec[i]*efun[i](zt) for i in range(1,8))); dvec[0]=Ce
# primitive P_e = int_{z0}^z e
_pc={}
def Pe(i,z):
    key=(i,float(z))
    if key in _pc: return _pc[key]
    v=mp.quad(lambda t: efun[i](t),[z0f,z]); _pc[key]=v; return v
def A(i,j,z):  # antisimmetrico
    return mp.quad(lambda t: efun[i](t)*Pe(j,t)-efun[j](t)*Pe(i,t),[z0f,z],maxdegree=8)
def Gt(z): return sum(cvec[i]*Pe(i,z) for i in range(8))     # int R
def etaf(z): return sum(dvec[i]*Pe(i,z) for i in range(8))   # int eta'
# coppie sopravvissute
pairs=[(i,j) for i in range(8) for j in range(i+1,8) if abs(cvec[i]*dvec[j]-cvec[j]*dvec[i])>1e-9]
def dphi_anti(z):
    tot=0.5*Gt(z)*etaf(z)
    for i,j in pairs: tot+=0.5*(cvec[i]*dvec[j]-cvec[j]*dvec[i])*A(i,j,z)
    return tot
def sQ(x): return np.sqrt(Q4(x))
def U(x,kk): return quad(lambda t:t**kk/((t-rd)*sQ(t)),r0,x,limit=200)[0]
def eta_dir(x): return U(x,3)-2*U(x,2)
def dmF(x):
    def Fm(mv): return np.polyval(np.polymul([Jc],[E**2-1,2*mv]),x)/np.sqrt(np.polyval(pS(mv,Jc),x))
    return (Fm(m+1e-6)-Fm(m-1e-6))/2e-6
def dphi_direct(x): return quad(lambda t: dmF(t)*eta_dir(t), r0, x, limit=200)[0]
print("coppie A sopravvissute:",len(pairs))
for i,j in pairs: print("   A[%s,%s] coeff=%+.5f"%(enam[i],enam[j],(cvec[i]*dvec[j]-cvec[j]*dvec[i]).real))
print("\n=== delta phi ANTISIMMETRICO (1/2 G~eta + 1/2 sum A) vs DIRETTO ===")
for rr in [11.0,10.0,9.2]:
    z=mp.mpf(zr(rr)); dc=complex(dphi_anti(z)).real; dd=dphi_direct(rr)
    print(f"  r={rr:4.1f}  antisym={dc:+.8f}  diretto={dd:+.8f}  diff={abs(dc-dd):.1e}")
