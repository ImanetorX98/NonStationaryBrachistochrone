# VALIDAZIONE MULTI-PARAMETRO del pipeline separatrice (ramo tau).
# Gira l'INTERO metodo (coeff ANALITICI) per vari E; ogni E -> Jc, r_d, curva diversi.
# Verifica: Sum_{ij} R_i eta'_j J[f1,f2]  ==  int diretto ∂_mF·eta,  per ogni E.
import numpy as np, mpmath as mp, sympy as sp
from scipy.integrate import quad
mp.mp.dps=40
def run(Efloat):
    E=Efloat; m=1.0; r0=12.0
    r,J=sp.symbols('r J'); Es=sp.nsimplify(Efloat)
    Ssym=sp.expand(r*(r-2)*((Es**2-1)*r+2)*(r**2*(r-2)-J**2*((Es**2-1)*r+2)))
    sols=[s for s in sp.solve(sp.Eq(sp.resultant(Ssym,sp.diff(Ssym,r),r),0),J) if s.is_real and float(s)>0.5]
    # scegli Jc con separatrice fisica: turning e4 in (2,r0), r_d reale
    Jc=None
    for sc in sorted(set(float(x) for x in sols)):
        Sc0=np.polymul(np.polymul(np.polymul([1,0],[1,-2]),[E**2-1,2.0]),
                       np.polysub(np.polymul([1,0,0],[1,-2]),sc**2*np.array([E**2-1,2.0])))
        rr=np.sort(np.real(np.roots(Sc0))); 
        pr=[(i,j) for i in range(6) for j in range(i+1,6) if abs(np.roots(Sc0)[i]-np.roots(Sc0)[j])<1e-6]
        if pr and any(2.0<x<r0 for x in rr): Jc=sc
    if Jc is None: return None
    def pS(mv,Jv):
        DE=np.array([E**2-1,2*mv]); p=np.polymul(np.polymul([1,0],[1,-2*mv]),DE)
        return np.polymul(p,np.polysub(np.polymul([1,0,0],[1,-2*mv]),Jv**2*DE))
    Sc=pS(m,Jc); dmS=(pS(m+1e-6,Jc)-pS(m-1e-6,Jc))/2e-6
    K=np.polymul([Jc],[E**2-1,2.0]); Nm=np.polysub(np.polymul(Sc,np.array([2*Jc])),0.5*np.polymul(K,dmS))
    rts=np.roots(Sc); pr=[(i,j) for i in range(6) for j in range(i+1,6) if abs(rts[i]-rts[j])<1e-6]
    rd=float(np.real((rts[pr[0][0]]+rts[pr[0][1]])/2)); Q,rem=np.polydiv(Sc,np.polymul([1,-rd],[1,-rd]))
    if np.max(np.abs(rem))>1e-6: return ('badfactor',Jc,rd)
    er=np.sort(np.real(np.roots(Q))); a4=Q[0]; e1r,e2r,e3r,e4r=er
    if not (2.0<e4r<r0): return ('noorbit',Jc,e4r)
    Qp=np.polyder(Q); Qpp=np.polyder(Qp); Nmp=np.polyder(Nm); Nmpp=np.polyder(Nmp)
    def Q4(x): return np.polyval(Q,x)
    k2=((e3r-e2r)*(e4r-e1r))/((e4r-e2r)*(e3r-e1r)); pref=2/mp.sqrt((e4r-e2r)*(e3r-e1r))/mp.sqrt(a4)
    om1=mp.mpf(float(pref*mp.ellipk(k2))); w_im=float(pref*mp.ellipk(1-k2)); tau=mp.mpc(0,w_im)/om1; qn=mp.exp(mp.pi*1j*tau)
    L1=lambda u: mp.jtheta(1,u,qn); L1p=lambda u: mp.jtheta(1,u,qn,1); L1pp=lambda u: mp.jtheta(1,u,qn,2); th1p0=L1p(0)
    eta1=-(mp.pi**2/(12*om1))*(mp.jtheta(1,0,qn,3)/th1p0)
    def wzet(z): u=mp.pi*z/(2*om1); return eta1*z/om1+(mp.pi/(2*om1))*(L1p(u)/L1(u))
    def wp(z):   u=mp.pi*z/(2*om1); rr=L1p(u)/L1(u); return -eta1/om1-(mp.pi/(2*om1))**2*(L1pp(u)/L1(u)-rr**2)
    def wpp(z):  return (wp(z+mp.mpf('1e-12'))-wp(z-mp.mpf('1e-12')))/mp.mpf('2e-12')
    z_inf=float(mp.re(mp.quad(lambda x:1/mp.sqrt(Q4(x)),[e4r,mp.inf]))); sa=float(mp.re(mp.sqrt(a4)))
    c_r=float(mp.re(e4r-(2/sa)*wzet(z_inf))); z_d=z_inf+float(mp.re(mp.quad(lambda x:1/mp.sqrt(Q4(x)),[-mp.inf,rd])))
    def zr(rv): return float(quad(lambda x:1/np.sqrt(Q4(x)),e4r,rv,limit=400)[0])
    z0f=mp.mpf(zr(r0))
    # residui ANALITICI
    s=np.sqrt(Q4(rd)); Q4prd=np.polyval(Qp,rd); Q4pprd=np.polyval(Qpp,rd)
    F=lambda x: np.polyval(Nm,x)/np.polyval(Q,x)
    Fp=lambda x:(np.polyval(Nmp,x)*np.polyval(Q,x)-np.polyval(Nm,x)*np.polyval(Qp,x))/np.polyval(Q,x)**2
    Fpp=lambda x:((np.polyval(Nmpp,x)*np.polyval(Q,x)-np.polyval(Nm,x)*np.polyval(Qpp,x))*np.polyval(Q,x)
       -2*np.polyval(Qp,x)*(np.polyval(Nmp,x)*np.polyval(Q,x)-np.polyval(Nm,x)*np.polyval(Qp,x)))/np.polyval(Q,x)**3
    a1=Q4prd/(4*s); a2=Q4pprd/12; h0=F(rd); h1=Fp(rd)*s; h2=0.5*(Fpp(rd)*s**2+Fp(rd)*(Q4prd/2))
    b1zd=(h2-3*a1*h1+(6*a1**2-3*a2)*h0)/s**3; b2zd=(h1-3*a1*h0)/s**3; b3zd=h0/s**3
    def b2_half(ei): return np.polyval(Nm,ei)/((ei-rd)**3)*(4/np.polyval(Qp,ei)**2)
    e1_zd=(rd**3-2*rd**2)/s; B=c_r+(1/sa)*float(mp.re(wzet(2*z_inf))); Aa=-1/sa
    e1_zi=Aa*(2*B+rd-2); e2_zi=Aa**2
    # trova quale e_i (tra e2,e3) mappa a z=iw_im: quello con N_m(e_i)!=0
    ei_iw=None
    for ei in [e1r,e2r,e3r]:
        if abs(np.polyval(Nm,ei))>1e-3*abs(np.polyval(Nm,e4r)): ei_iw=ei
    def fbasis(kind,p):
        if kind=='1': return lambda z: mp.mpf(1)
        if kind=='ze': return lambda z: wzet(z-p)
        if kind=='wp': return lambda z: wp(z-p)
        if kind=='wpp':return lambda z: wpp(z-p)
    Rlist=[(0,'1',mp.mpf(0)),
           (b1zd,'ze',mp.mpf(z_d)),(b2zd,'wp',mp.mpf(z_d)),(-b3zd/2,'wpp',mp.mpf(z_d)),
           (-b1zd,'ze',mp.mpf(-z_d)),(b2zd,'wp',mp.mpf(-z_d)),(b3zd/2,'wpp',mp.mpf(-z_d)),
           (b2_half(e4r),'wp',mp.mpf(0)),(b2_half(ei_iw),'wp',1j*mp.mpf(w_im))]
    Elist=[(0,'1',mp.mpf(0)),(e1_zd,'ze',mp.mpf(z_d)),(-e1_zd,'ze',mp.mpf(-z_d)),
           (e2_zi,'wp',mp.mpf(z_inf)),(e2_zi,'wp',mp.mpf(-z_inf)),(e1_zi,'ze',mp.mpf(z_inf)),(-e1_zi,'ze',mp.mpf(-z_inf))]
    def r_of_z(z): return c_r-(1/sa)*(wzet(z-z_inf)-wzet(z+z_inf))
    def Rexact(z): rr=r_of_z(z); return np.polyval(Nm,rr)/((rr-rd)**3*np.polyval(Q,complex(rr)))
    def Epexact(z): rr=r_of_z(z); return (rr**3-2*rr**2)/(rr-rd)
    zt=mp.mpf('0.19')
    C0=complex(Rexact(zt)-sum(c*fbasis(k,p)(zt) for c,k,p in Rlist)); Rlist[0]=(C0,'1',mp.mpf(0))
    Ce=complex(Epexact(zt)-sum(c*fbasis(k,p)(zt) for c,k,p in Elist)); Elist[0]=(Ce,'1',mp.mpf(0))
    chkR=float(max(abs(complex(Rexact(mp.mpf(t))-sum(c*fbasis(k,p)(mp.mpf(t)) for c,k,p in Rlist))) for t in [0.15,0.22]))
    _ic={}
    def inner(k2_,p2,t):
        key=(k2_,float(p2.real),float(p2.imag),float(t)); 
        if key in _ic: return _ic[key]
        v=mp.quad(lambda ss: fbasis(k2_,p2)(ss),[z0f,t]); _ic[key]=v; return v
    def Jatom(k1,p1,k2_,p2,z): return mp.quad(lambda t: fbasis(k1,p1)(t)*inner(k2_,p2,t),[z0f,z],maxdegree=8)
    def dphi_terms(z): return sum(c1*c2*Jatom(k1,p1,k2_,p2,z) for c1,k1,p1 in Rlist for c2,k2_,p2 in Elist)
    def sQ(x): return np.sqrt(Q4(x))
    def U(x,kk): return quad(lambda t:t**kk/((t-rd)*sQ(t)),r0,x,limit=200)[0]
    def eta(x): return U(x,3)-2*U(x,2)
    def dmF(x):
        def Fm(mv): return np.polyval(np.polymul([Jc],[E**2-1,2*mv]),x)/np.sqrt(np.polyval(pS(mv,Jc),x))
        return (Fm(m+1e-6)-Fm(m-1e-6))/2e-6
    def dphi_direct(x): return quad(lambda t: dmF(t)*eta(t), r0, x, limit=200)[0]
    xtest=0.5*(e4r+r0)  # meta' orbita
    z=mp.mpf(zr(xtest)); dt=complex(dphi_terms(z)).real; dd=dphi_direct(xtest)
    return dict(E=E,Jc=Jc,rd=rd,e4=e4r,chkR=chkR,xtest=xtest,dterms=dt,ddir=dd,diff=abs(dt-dd))

print("=== VALIDAZIONE MULTI-PARAMETRO (M=1, vari E) ===")
print(" E      Jc        r_d       e4(turn)  checkR    delta_terms   delta_direct   diff")
for Ev in [1.30,1.40,1.50,1.60,1.25]:
    res=run(Ev)
    if res is None or 'dterms' not in res: print(f" {Ev}: {res}"); continue
    print(f" {res['E']:.2f}  {res['Jc']:.5f}  {res['rd']:+.5f}  {res['e4']:.4f}  {res['chkR']:.1e}  {res['dterms']:+.8f}  {res['ddir']:+.8f}  {res['diff']:.1e}")
