# SEPARATRICE Thakurta-Kerr ramo t (tempo COORDINATA) con TRACKING di Jc(E).
# Segue la separatrice mobile: J=Jc(E_eff). Sorgente totale:
#   N_tot = N_t + (dJc/dE) N_J^t,  N_t=dE(K/sqrtR6)R6^{3/2}, N_J^t=dJ(K/sqrtR6)R6^{3/2}.
#   dJc/dE = -Q2_E/Q2_J  a (r_d,Jc)  (doppia radice di Q2: Q2_r=0 -> Q2_J dJc+Q2_E dE=0).
# Clock t = eta_t (2a+3a specie) INVARIATO. Due separatrici: Jc+ (prograda), Jc- (retrograda).
import numpy as np, mpmath as mp, sympy as sp, sys, logging
from scipy.integrate import quad
try: from tqdm import tqdm
except ImportError:
    def tqdm(it,**k): return it
logging.basicConfig(level=logging.INFO,format='[%(asctime)s] %(message)s',datefmt='%H:%M:%S',stream=sys.stdout)
log=logging.getLogger(); mp.mp.dps=30
M=1.0; a=0.9; E=1.2; r0=20.0
r,J,Esy,Msy,asy=sp.symbols('r J E M a'); asub={Msy:1,asy:sp.Rational(9,10)}
DE=(Esy**2-1)*r+2*Msy; Delta=r**2-2*Msy*r+asy**2
Q2=(2*Esy**2*J**2*Msy*r-Esy**2*J**2*r**2-4*Esy**2*J*Msy*asy*r+2*Esy**2*Msy*asy**2*r+Esy**2*asy**2*r**2+Esy**2*r**4
    +4*J**2*Msy**2-4*J**2*Msy*r+J**2*r**2-8*J*Msy**2*asy+4*J*Msy*asy*r+4*Msy**2*asy**2)
R6f=sp.expand(r*Q2*DE); Ktf=r*DE*(J*(r-2*Msy)+2*Msy*asy)/Delta
log.info("derivo sorgenti N_t e N_J^t (simbolico, puo' richiedere ~1min)")
Ntf=sp.expand(sp.simplify(sp.diff(Ktf/sp.sqrt(R6f),Esy)*R6f**sp.Rational(3,2)))
NJtf=sp.expand(sp.simplify(sp.diff(Ktf/sp.sqrt(R6f),J)*R6f**sp.Rational(3,2)))
rho_tf=sp.cancel(sp.together((Esy**2*r**3-2*Msy*asy*Ktf/r)/((r-2*Msy)/r)))
# sostituisco M,a numerici (E,J restano simbolici)
R6=sp.expand(R6f.subs(asub)); Nt=sp.expand(Ntf.subs(asub)); NJt=sp.expand(NJtf.subs(asub))
rho_t=sp.cancel(rho_tf.subs(asub)); DeltaN=Delta.subs(asub)
# clock: rho_t = P3 + RD/Delta (denom vero = 2500*Delta)
P3poly,_=sp.div(sp.Poly(sp.numer(rho_t),r),sp.Poly(sp.denom(rho_t),r))
P3=sp.expand(P3poly.as_expr()); RD=sp.expand(sp.simplify((rho_t-P3)*DeltaN)); Delta=DeltaN
# derivate di Q2 per dJc/dE (E,J simbolici, poi valuto a r_d,Jc,E=6/5)
Q2n_a=Q2.subs(asub); Q2E=sp.diff(Q2n_a,Esy); Q2J=sp.diff(Q2n_a,J)
Es56=sp.Rational(6,5)
Q2E_l=sp.lambdify((r,J),Q2E.subs(Esy,Es56),'mpmath'); Q2J_l=sp.lambdify((r,J),Q2J.subs(Esy,Es56),'mpmath')

def run(Jc0, rd_guess, SIGN, pyanchor_fixed):
    log.info("=== %s: Jc0=%.4f, raffino (Jc,r_d) ==="%(SIGN,Jc0))
    Q2E5=Q2.subs(asub).subs(Esy,Es56)
    Q2l=sp.lambdify((r,J),Q2E5,'mpmath'); Q2rl=sp.lambdify((r,J),sp.diff(Q2E5,r),'mpmath')
    sol=mp.findroot(lambda rd_,Jc_:[Q2l(rd_,Jc_),Q2rl(rd_,Jc_)],(mp.mpf(str(rd_guess)),mp.mpf(str(Jc0))))
    rd_hp=float(sol[0]); Jc=float(sol[1]); rd=rd_hp
    dJcdE=float(-Q2E_l(rd,Jc)/Q2J_l(rd,Jc))
    log.info("%s raffinato: Jc=%.10f r_d=%.10f  dJc/dE=%.6f"%(SIGN,Jc,rd,dJcdE))
    # coeff polinomi a Jc
    R6c=np.array([float(c) for c in sp.Poly(R6.subs(Esy,Es56).subs(J,Jc),r).all_coeffs()])
    Ntc=np.array([float(c) for c in sp.Poly(sp.expand(Nt.subs(Esy,Es56)).subs(J,Jc),r).all_coeffs()])
    NJc=np.array([float(c) for c in sp.Poly(sp.expand(NJt.subs(Esy,Es56)).subs(J,Jc),r).all_coeffs()])
    # N_tot = N_t + dJc/dE * N_J^t (allineo gradi)
    L=max(len(Ntc),len(NJc)); Ntot=np.r_[np.zeros(L-len(Ntc)),Ntc]+dJcdE*np.r_[np.zeros(L-len(NJc)),NJc]
    P3n=sp.lambdify(r,P3.subs(Esy,Es56).subs(J,Jc),'numpy'); RDn=sp.lambdify(r,RD.subs(Esy,Es56).subs(J,Jc),'numpy'); Dn=sp.lambdify(r,Delta,'numpy')
    Q,rem=np.polydiv(R6c,np.polymul([1,-rd],[1,-rd]))
    log.info("resto polydiv R6/(r-r_d)^2 = %.1e"%np.max(np.abs(rem)))
    er=np.sort(np.real(np.roots(Q))); a4=Q[0]; e1r,e2r,e3r,e4r=er
    def Q4(x): return np.polyval(Q,x)
    k2=((e3r-e2r)*(e4r-e1r))/((e4r-e2r)*(e3r-e1r)); pref=2/mp.sqrt((e4r-e2r)*(e3r-e1r))/mp.sqrt(a4)
    om1=mp.mpf(float(pref*mp.ellipk(k2))); w_im=float(pref*mp.ellipk(1-k2)); tau=mp.mpc(0,w_im)/om1; qn=mp.exp(mp.pi*1j*tau)
    L1=lambda u: mp.jtheta(1,u,qn); L1p=lambda u: mp.jtheta(1,u,qn,1); L1pp=lambda u: mp.jtheta(1,u,qn,2); th1p0=L1p(0)
    eta1=-(mp.pi**2/(12*om1))*(mp.jtheta(1,0,qn,3)/th1p0)
    def wsig(z): u=mp.pi*z/(2*om1); return (2*om1/mp.pi)*mp.exp(eta1*z**2/(2*om1))*L1(u)/th1p0
    def wzet(z): u=mp.pi*z/(2*om1); return eta1*z/om1+(mp.pi/(2*om1))*(L1p(u)/L1(u))
    def wp(z):   u=mp.pi*z/(2*om1); rr=L1p(u)/L1(u); return -eta1/om1-(mp.pi/(2*om1))**2*(L1pp(u)/L1(u)-rr**2)
    def wpp(z):
        u=mp.pi*z/(2*om1); T1=L1p(u)/L1(u); T2=L1pp(u)/L1(u); T3=mp.jtheta(1,u,qn,3)/L1(u)
        return -(mp.pi/(2*om1))**3*(T3-3*T1*T2+2*T1**3)
    z_inf=float(mp.re(mp.quad(lambda x:1/mp.sqrt(Q4(x)),[e4r,mp.inf]))); sa=float(mp.sqrt(a4))
    c_r=float(mp.re(e4r-(2/sa)*wzet(z_inf))); z_d=z_inf+float(mp.re(mp.quad(lambda x:1/mp.sqrt(Q4(x)),[-mp.inf,rd])))
    def zr(rv): return float(quad(lambda x:1/np.sqrt(Q4(x)),e4r,rv,limit=400)[0])
    z0=mp.mpf(zr(r0)); iw=1j*mp.mpf(w_im)
    def r_of_z(z): return c_r-(1/sa)*(wzet(z-z_inf)-wzet(z+z_inf))
    def Rf(z): rr=r_of_z(z); return np.polyval(Ntot,rr)/((rr-rd)**3*np.polyval(Q,complex(rr)))  # sorgente N_tot
    log.info("check r(z): %.1e"%max(abs(float(np.real(r_of_z(mp.mpf(zr(rv)))))-rv) for rv in [19.,18.]))
    def lau(fun,aa,order,eps=1e-3):
        return {n:complex(mp.quad(lambda th: fun(aa+eps*mp.exp(1j*th))*(eps*mp.exp(1j*th))**(n-1)*(1j*eps*mp.exp(1j*th)),[0,2*mp.pi])/(2j*mp.pi)) for n in range(1,order+1)}
    PPr={'zd':(mp.mpf(z_d),lau(Rf,mp.mpf(z_d),3)),'mzd':(mp.mpf(-z_d),lau(Rf,mp.mpf(-z_d),3)),'h0':(mp.mpf(0),lau(Rf,mp.mpf(0),2)),'h2':(iw,lau(Rf,iw,2))}
    C0=complex(Rf(mp.mpf('0.05'))-sum(b.get(1,0)*wzet(mp.mpf('0.05')-aa)+b.get(2,0)*wp(mp.mpf('0.05')-aa)+b.get(3,0)*(-wpp(mp.mpf('0.05')-aa)/2) for aa,b in PPr.values()))
    def Gt(z):
        s=C0*z
        for nm,(aa,b) in PPr.items(): s+= b.get(1,0)*mp.log(wsig(z-aa))-b.get(2,0)*wzet(z-aa)-(b.get(3,0)/2)*wp(z-aa)
        return s
    log.info("check dG/dz=R_tot: %.1e"%max(abs(complex((Gt(mp.mpf(t)+mp.mpf('1e-7'))-Gt(mp.mpf(t)-mp.mpf('1e-7')))/mp.mpf('2e-7')-Rf(mp.mpf(t)))) for t in [0.05,0.09]))
    def sQ(x): return np.sqrt(Q4(x))
    # orizzonti z(r_pm) per la decomposizione della clock
    rp=1+np.sqrt(1-a*a); rm=1-np.sqrt(1-a*a)
    def r_solve(target,zg): return mp.findroot(lambda z: r_of_z(z)-target, zg)
    def find_zimg(target):
        for zg in [iw/2, iw/2+0.3, iw*0.7, om1+iw/2, 0.5+iw/2, iw/3]:
            try:
                z=r_solve(target,zg)
                if abs(complex(r_of_z(z))-target)<1e-8: return z
            except Exception: pass
        return None
    zrp=find_zimg(rp); zrm=find_zimg(rm)
    log.info("z(r+)=%s z(r-)=%s"%(mp.nstr(zrp,6) if zrp else None,mp.nstr(zrm,6) if zrm else None))
    def etpz(z): rr=r_of_z(z); return (P3n(rr)+RDn(rr)/Dn(rr))/(rr-rd)
    def lauv(fun,aa,order,eps=5e-4):
        return {n:complex(mp.quad(lambda th: fun(aa+eps*mp.exp(1j*th))*(eps*mp.exp(1j*th))**(n-1)*(1j*eps*mp.exp(1j*th)),[0,2*mp.pi])/(2j*mp.pi)) for n in range(1,order+1)}
    poles=[('zd',mp.mpf(z_d)),('mzd',mp.mpf(-z_d)),('zi',mp.mpf(z_inf),2),('mzi',mp.mpf(-z_inf),2),
           ('rp',zrp),('mrp',-zrp),('rm',zrm),('mrm',-zrm)]
    PPe={}
    for p in poles:
        nm=p[0]; aa=p[1]; order=p[2] if len(p)>2 else 1
        PPe[nm]=(aa,lauv(etpz,aa,order))
    ztv=mp.mpf('0.03')
    def prcv(z): return sum(b.get(1,0)*wzet(z-aa)+b.get(2,0)*wp(z-aa) for aa,b in PPe.values())
    Cv=complex(etpz(ztv)-prcv(ztv))
    recerr=max(abs(complex(etpz(mp.mpf(t))-(Cv+prcv(mp.mpf(t))))) for t in [0.02,0.04])
    log.info("check eta_t' ricostruito: %.1e"%recerr)
    def etat(x): return quad(lambda t:(P3n(t)+RDn(t)/Dn(t))/((t-rd)*sQ(t)), r0, x, limit=200)[0]
    def dEF(x): return np.polyval(Ntot,x)/np.polyval(R6c,x)**1.5  # sorgente TOTALE
    def dphi_direct(x): return quad(lambda t: dEF(t)*etat(t), r0, x, limit=200)[0]
    # BLOCK: delta phi_t = G~(z) eta_t - sum c_i c_j int gf_i gf_j
    Gterms=[(C0,'z',mp.mpf(0))]
    for nm,(aa,b) in PPr.items(): Gterms+=[(b.get(1,0),'lns',aa),(-b.get(2,0),'zet',aa),(-(b.get(3,0))/2,'wp',aa)]
    Vterms=[(Cv,'one',mp.mpf(0))]
    for nm,(aa,b) in PPe.items():
        if abs(b.get(1,0))>1e-11: Vterms.append((b.get(1,0),'zet',aa))
        if abs(b.get(2,0))>1e-11: Vterms.append((b.get(2,0),'wp',aa))
    def gf(k,p,z):
        if k=='z': return z
        if k=='one': return mp.mpf(1)
        if k=='lns': return mp.log(wsig(z-p))
        if k=='zet': return wzet(z-p)
        if k=='wp': return wp(z-p)
    def atomP(ki,pi,kj,pj,z): return mp.quad(lambda t: gf(ki,pi,t)*gf(kj,pj,t),[z0,z],maxdegree=8)
    def intGeta(z): return sum(ci*cj*atomP(ki,pi,kj,pj,z) for ci,ki,pi in Gterms for cj,kj,pj in Vterms)
    def dphi_block(x):
        z=mp.mpf(zr(x)); return complex(Gt(z)).real*etat(x)-complex(intGeta(z)).real
    ndl=sum(1 for ci,ki,pi in Gterms for cj,kj,pj in Vterms if ki=='lns' and kj=='zet' and abs(ci*cj)>1e-9)
    print("=== TK t %s TRACKING BLOCK: G~ eta_t - sum atomi vs DIRETTO (N_tot) ==="%SIGN)
    for x in tqdm([19.0,18.0,17.0],desc="TK-t-track-%s"%SIGN):
        db=dphi_block(x); dd=dphi_direct(x)
        print(f"  r={x:4.1f}  block={db:+.8e}  dir={dd:+.8e}  diff={abs(db-dd):.1e}")
    print("%s: dJc/dE=%.6f ; dilog (lnsigma x zeta)=%d ; prodotti=%d ; N_tot=N_t+(dJc/dE)N_J^t."%(
        SIGN,dJcdE,ndl,len(Gterms)*len(Vterms)))
    return dJcdE

run(19.088436,-6.62,'prograda',None)
run(-18.671,-6.588,'retrograda',None)
log.info("FINE tracking ramo t (entrambe le separatrici)")
