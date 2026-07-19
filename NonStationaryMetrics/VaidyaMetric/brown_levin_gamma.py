# BROWN-LEVIN Gamma-tilde(1,1) SCRITTA ESPLICITAMENTE + serie-q TABULABILE.
# g^(1)(x,tau) = d/dx log theta1(pi x, q),  q=e^{i pi tau}, reticolo (1,tau).
#   serie-q (Fourier):  g^(1)(x) = pi cot(pi x) + 4 pi sum_{n>=1} q^{2n}/(1-q^{2n}) sin(2 pi n x)
# Gamma-tilde(1,1; x1,x2; y) = int_0^y g1(t-x1) [ int_0^t g1(s-x2) ds ] dt   (length-2 Kronecker)
# elliptic dilog E2(x;y) = int_0^y Ltheta(t-x) dt,  Ltheta=log theta1(pi.,q)  (peso-2, tabulabile)
import mpmath as mp
mp.mp.dps=30
def make(tau):
    q=mp.exp(mp.pi*1j*tau)
    Lth =lambda x: mp.log(mp.jtheta(1,mp.pi*x,q))
    g1_th=lambda x: mp.pi*mp.jtheta(1,mp.pi*x,q,1)/mp.jtheta(1,mp.pi*x,q)   # via theta1
    def g1_q(x,N=80):                                                       # via serie-q (TABULABILE)
        s=mp.pi/mp.tan(mp.pi*x)
        for n in range(1,N+1): s+=4*mp.pi*(q**(2*n)/(1-q**(2*n)))*mp.sin(2*mp.pi*n*x)
        return s
    return q,Lth,g1_th,g1_q
tau=1j*mp.mpf('0.9059733802550'); q,Lth,g1_th,g1_q=make(tau)
print("=== g^(1): serie-q (tabulabile) vs theta1 ===")
for x in [0.13,0.31,0.47]:
    x=mp.mpf(x); print(f"  x={float(x):.2f}  g1_q={complex(g1_q(x)).real:+.10f}  g1_th={complex(g1_th(x)).real:+.10f}  diff={float(abs(g1_q(x)-g1_th(x))):.1e}")
# --- Gamma-tilde(1,1) definizione esplicita (iterato) ---
def Gamma11(x1,x2,y,y0=mp.mpf(0)):
    inner=lambda t: mp.quad(lambda s: g1_th(s-x2),[y0,t])
    return mp.quad(lambda t: g1_th(t-x1)*inner(t),[y0,y])
# verifica vs GiNaC-style nested (gia' 1e-18) e antisimmetria/shuffle: Gamma11(x1,x2)+Gamma11(x2,x1)=Lth-prod
x1,x2,y=mp.mpf('0.30')*1j,mp.mpf('0.42')*1j,mp.mpf('0.15')
G12=Gamma11(x1,x2,y); G21=Gamma11(x2,x1,y)
shuffle=(Lth(y-x1)-Lth(-x1))*(Lth(y-x2)-Lth(-x2))
print("\n=== Gamma-tilde(1,1) shuffle:  G(x1,x2)+G(x2,x1) = [int g1(-x1)][int g1(-x2)] ===")
print(f"  G12+G21={complex(G12+G21).real:+.10f}  shuffle={complex(shuffle).real:+.10f}  diff={float(abs(G12+G21-shuffle)):.1e}")

# ============ MAPPA D(a,b) -> Gamma-tilde  (dati veri della separatrice Vaidya) ============
import numpy as np, sympy as sp
from scipy.integrate import quad
E=1.4; m=1.0; r0=12.0
rS,JS=sp.symbols('r J'); Es=sp.Rational(7,5)
Ssym=sp.expand(rS*(rS-2)*((Es**2-1)*rS+2)*(rS**2*(rS-2)-JS**2*((Es**2-1)*rS+2)))
Jc=float([s for s in sp.solve(sp.Eq(sp.resultant(Ssym,sp.diff(Ssym,rS),rS),0),JS) if s.is_real and float(s)>1][0])
def pS(mv,Jv):
    DE=np.array([E**2-1,2*mv]); p=np.polymul(np.polymul([1,0],[1,-2*mv]),DE)
    return np.polymul(p,np.polysub(np.polymul([1,0,0],[1,-2*mv]),Jv**2*DE))
Sc=pS(m,Jc); rts=np.roots(Sc); pr=[(i,j) for i in range(6) for j in range(i+1,6) if abs(rts[i]-rts[j])<1e-6]
rd=float(np.real((rts[pr[0][0]]+rts[pr[0][1]])/2)); Q,_=np.polydiv(Sc,np.polymul([1,-rd],[1,-rd]))
er=np.sort(np.real(np.roots(Q))); a4=Q[0]; e1r,e2r,e3r,e4r=er
def Q4(x): return np.polyval(Q,x)
k2=((e3r-e2r)*(e4r-e1r))/((e4r-e2r)*(e3r-e1r)); pref=2/mp.sqrt((e4r-e2r)*(e3r-e1r))/mp.sqrt(a4)
om1=mp.mpf(float(pref*mp.ellipk(k2))); w_im=float(pref*mp.ellipk(1-k2)); tau2=mp.mpc(0,w_im)/om1; q2=mp.exp(mp.pi*1j*tau2)
L1=lambda u: mp.jtheta(1,u,q2); L1p=lambda u: mp.jtheta(1,u,q2,1); th1p0=L1p(0)
eta1=-(mp.pi**2/(12*om1))*(mp.jtheta(1,0,q2,3)/th1p0)
def wsig(z): u=mp.pi*z/(2*om1); return (2*om1/mp.pi)*mp.exp(eta1*z**2/(2*om1))*L1(u)/th1p0
def wzet(z): u=mp.pi*z/(2*om1); return eta1*z/om1+(mp.pi/(2*om1))*(L1p(u)/L1(u))
z_inf=float(mp.quad(lambda x:1/mp.sqrt(Q4(x)),[e4r,mp.inf]))
z_d=z_inf+float(mp.quad(lambda x:1/mp.sqrt(Q4(x)),[-mp.inf,rd]))
def zr(rv): return float(quad(lambda x:1/np.sqrt(Q4(x)),e4r,rv,limit=400)[0])
# theta/g1 in QUESTO reticolo:
q,Lth,g1_th,g1_q=make(tau2)
# (1) verifica identita' Weierstrass <-> theta1/g1:
zt=mp.mpf('0.19'); bb=z_d
print("\n=== identita' Weierstrass <-> theta1/g1 (reticolo Vaidya) ===")
lhsz=wzet(zt-bb); rhsz=(eta1/om1)*(zt-bb)+(1/(2*om1))*g1_th((zt-bb)/(2*om1))
print(f"  zeta_W(z-b): W={complex(lhsz).real:+.8f}  theta={complex(rhsz).real:+.8f}  diff={float(abs(lhsz-rhsz)):.1e}")
aa=z_d
lhss=mp.log(wsig(zt-aa)); rhss=eta1*(zt-aa)**2/(2*om1)+Lth((zt-aa)/(2*om1))+mp.log(2*om1/mp.pi)-mp.log(th1p0)
print(f"  lnσ_W(z-a): W={complex(lhss).real:+.8f}  theta={complex(rhss).real:+.8f}  diff={float(abs(lhss-rhss)):.1e}")
# (2) NUCLEO TRASCENDENTE di D:  int Lth(u'-ahat) g1(u'-bhat) du' = Gamma11 + peso-1
z0f=mp.mpf(zr(r0)); zx=mp.mpf(zr(10.0))
for (na,a),(nb,b) in [(('z_d',z_d),('z_inf',z_inf)),(('z_d',z_d),('z_d',z_d))]:
    ah=a/(2*om1); bh=b/(2*om1); u0=z0f/(2*om1); uz=zx/(2*om1)
    core=mp.quad(lambda t: Lth(t-ah)*g1_th(t-bh),[u0,uz])
    G11=Gamma11(ah,bh,uz,u0)
    w1=Lth(u0-ah)*(Lth(uz-bh)-Lth(u0-bh))
    print(f"  D-core({na},{nb}): ∫Lθ·g1={complex(core).real:+.8f}  =Γ̃11+w1={complex(G11+w1).real:+.8f}  diff={float(abs(core-(G11+w1))):.1e}")
