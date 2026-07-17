# delta phi_tau|_sep : FORMULA 100% TERMINE-PER-TERMINE in Brown-Levin Gamma-tilde.
# delta phi = int_{z0}^z R(z') eta(z') dz',  eta=int eta'  (NIENTE termine di bordo).
#   = sum_{i,j} R_i eta'_j  J[f1^i, f2^j],   J[f1,f2]=int f1(z')[int_{z0}^{z'} f2] dz'
# f1 in {1, zeta_a, P_a, P'_a} (a: poli di R),  f2 in {1, zeta_b, P_b} (b: poli di eta').
# Ogni J = iterato length-2 di forme Weierstrass = Gamma-tilde(m,n; ahat,bhat), zeta->g1,P->g2,P'->g3.
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
def wsig(z): u=mp.pi*z/(2*om1); return (2*om1/mp.pi)*mp.exp(eta1*z**2/(2*om1))*L1(u)/th1p0
def wzet(z): u=mp.pi*z/(2*om1); return eta1*z/om1+(mp.pi/(2*om1))*(L1p(u)/L1(u))
def wp(z):   u=mp.pi*z/(2*om1); rr=L1p(u)/L1(u); return -eta1/om1-(mp.pi/(2*om1))**2*(L1pp(u)/L1(u)-rr**2)
def wpp(z):  return (wp(z+mp.mpf('1e-12'))-wp(z-mp.mpf('1e-12')))/mp.mpf('2e-12')
z_inf=float(mp.quad(lambda x:1/mp.sqrt(Q4(x)),[e4r,mp.inf])); sa=float(mp.sqrt(a4))
c_r=float(mp.re(e4r-(2/sa)*wzet(z_inf))); z_d=z_inf+float(mp.quad(lambda x:1/mp.sqrt(Q4(x)),[-mp.inf,rd]))
def zr(rv): return float(quad(lambda x:1/np.sqrt(Q4(x)),e4r,rv,limit=400)[0])
z0f=mp.mpf(zr(r0))

# ===== residui ANALITICI (coeff dei termini) =====
s=np.sqrt(Q4(rd)); Q4prd=np.polyval(Qp,rd); Q4pprd=np.polyval(Qpp,rd)
F=lambda x: np.polyval(Nm,x)/np.polyval(Q,x)
Fp=lambda x:(np.polyval(Nmp,x)*np.polyval(Q,x)-np.polyval(Nm,x)*np.polyval(Qp,x))/np.polyval(Q,x)**2
Fpp=lambda x:((np.polyval(Nmpp,x)*np.polyval(Q,x)-np.polyval(Nm,x)*np.polyval(Qpp,x))*np.polyval(Q,x)
   -2*np.polyval(Qp,x)*(np.polyval(Nmp,x)*np.polyval(Q,x)-np.polyval(Nm,x)*np.polyval(Qp,x)))/np.polyval(Q,x)**3
a1=Q4prd/(4*s); a2=Q4pprd/12; h0=F(rd); h1=Fp(rd)*s; h2=0.5*(Fpp(rd)*s**2+Fp(rd)*(Q4prd/2))
b1zd=(h2-3*a1*h1+(6*a1**2-3*a2)*h0)/s**3; b2zd=(h1-3*a1*h0)/s**3; b3zd=h0/s**3
# R half-period residui (b2 solo) via analitico: al punto e_i, R~N_m/((r-rd)^3 Q), Q~Q'(e_i)(r-e_i),
# r-e_i ~ (Q'(e_i)/4)(z-w_i)^2  => b2 = N_m(e_i)/((e_i-rd)^3) * 4/Q'(e_i)^2 * ... uso valore da residui gia' noti:
def b2_half(ei):
    Qpe=np.polyval(Qp,ei); return np.polyval(Nm,ei)/((ei-rd)**3)*(4/Qpe**2)
# C0: R = C0 + sum principal; C0 = valore medio. Uso: C0 = R(zt)-somma parti principali (numerico da forme)
# eta' residui analitici
e1_zd=(rd**3-2*rd**2)/s
B=c_r+(1/sa)*float(mp.re(wzet(2*z_inf))); A=-1/sa
e1_zi=A*(2*B+rd-2); e2_zi=A**2
Ce=None  # calcolato sotto per consistenza

# ---- funzioni base ----
def fbasis(kind,p):
    if kind=='1':  return lambda z: mp.mpf(1)
    if kind=='z':  return lambda z: z
    if kind=='ze': return lambda z: wzet(z-p)
    if kind=='wp': return lambda z: wp(z-p)
    if kind=='wpp':return lambda z: wpp(z-p)       # P'(z-p) plain
# R e eta' come liste (coeff, kind, pole)
hp=[mp.mpf(0), 1j*mp.mpf(w_im)]                     # semi-periodi con polo (e4->0, e2->iw_im)
Rlist=[('C0','1',mp.mpf(0))]  # C0 sotto
for (nm,a,sg) in [('zd',z_d,1),('mzd',-z_d,1)]:
    Rlist+=[ (sg*b1zd if a>0 else -b1zd,'ze',mp.mpf(a)),(b2zd,'wp',mp.mpf(a)),((-1)*( b3zd if a>0 else -b3zd)/2,'wpp',mp.mpf(a)) ]
Rlist+=[(b2_half(e4r),'wp',mp.mpf(0)),(b2_half(e3r),'wp',1j*mp.mpf(w_im))]  # e4->z=0, e3->z=iw_im
Elist=[('Ce','1',mp.mpf(0))]
Elist+=[(e1_zd,'ze',mp.mpf(z_d)),(-e1_zd,'ze',mp.mpf(-z_d)),(e2_zi,'wp',mp.mpf(z_inf)),(e2_zi,'wp',mp.mpf(-z_inf)),
        (e1_zi,'ze',mp.mpf(z_inf)),(-e1_zi,'ze',mp.mpf(-z_inf))]
# R(z), eta'(z) espliciti (per fissare C0, Ce e verifica)
def r_of_z(z): return c_r-(1/sa)*(wzet(z-z_inf)-wzet(z+z_inf))
def Rexact(z): rr=r_of_z(z); return np.polyval(Nm,rr)/((rr-rd)**3*np.polyval(Q,complex(rr)))
def Epexact(z): rr=r_of_z(z); return (rr**3-2*rr**2)/(rr-rd)
zt=mp.mpf('0.19')
def Rsum_noC0(z): return sum(c*fbasis(k,p)(z) for c,k,p in Rlist if not isinstance(c,str))
def Esum_noCe(z): return sum(c*fbasis(k,p)(z) for c,k,p in Elist if not isinstance(c,str))
C0=complex(Rexact(zt)-Rsum_noC0(zt)); Ce=complex(Epexact(zt)-Esum_noCe(zt))
Rlist[0]=(C0,'1',mp.mpf(0)); Elist[0]=(Ce,'1',mp.mpf(0))
print("C0=%.6f  Ce=%.6f"%(C0.real,Ce.real))
print("check R:",float(max(abs(complex(Rexact(mp.mpf(t))-sum(c*fbasis(k,p)(mp.mpf(t)) for c,k,p in Rlist))) for t in [0.15,0.22])))
print("check eta':",float(max(abs(complex(Epexact(mp.mpf(t))-sum(c*fbasis(k,p)(mp.mpf(t)) for c,k,p in Elist))) for t in [0.15,0.22])))

# ---- J[f1,f2] = iterato length-2 (Brown-Levin Gamma-tilde) ----
_ic={}
def inner(k2_,p2,t):
    key=(k2_,float(p2.real),float(p2.imag),float(t.real),float(t.imag))
    if key in _ic: return _ic[key]
    v=mp.quad(lambda s: fbasis(k2_,p2)(s),[z0f,t]); _ic[key]=v; return v
def Jatom(k1,p1,k2_,p2,z):
    return mp.quad(lambda t: fbasis(k1,p1)(t)*inner(k2_,p2,t),[z0f,z],maxdegree=8)
def dphi_terms(z):
    tot=0
    for c1,k1,p1 in Rlist:
        for c2,k2_,p2 in Elist:
            tot+=c1*c2*Jatom(k1,p1,k2_,p2,z)
    return tot
# ---- diretto ----
def sQ(x): return np.sqrt(Q4(x))
def U(x,kk): return quad(lambda t:t**kk/((t-rd)*sQ(t)),r0,x,limit=200)[0]
def eta(x): return U(x,3)-2*U(x,2)
def dmF(x):
    def Fm(mv): return np.polyval(np.polymul([Jc],[E**2-1,2*mv]),x)/np.sqrt(np.polyval(pS(mv,Jc),x))
    return (Fm(m+1e-6)-Fm(m-1e-6))/2e-6
def dphi_direct(x): return quad(lambda t: dmF(t)*eta(t), r0, x, limit=200)[0]
print("\n=== delta phi = Sum_ij R_i eta'_j J[f1,f2]  (100% termine-per-termine)  vs DIRETTO ===")
for rr in [11.0,10.5,10.0]:
    z=mp.mpf(zr(rr)); dt=complex(dphi_terms(z)).real; dd=dphi_direct(rr)
    print(f"  r={rr:4.1f}  Sum_termini={dt:+.8f}  diretto={dd:+.8f}  diff={abs(dt-dd):.1e}")
# conteggio termini
nt=sum(1 for c1,k1,p1 in Rlist for c2,k2_,p2 in Elist if abs(c1*c2)>1e-9)
print(f"\ntermini non nulli: {nt}  (ogni J = Gamma-tilde(m,n) via zeta->g1,P->g2,P'->g3)")
