# CHIUSURA weight-2 di delta phi|_sep: int G~ eta' dz = somma di ATOMI dilog ellittico.
# G~ e eta' scomposti in parti principali (coeff ANALITICI da contorno). Prodotto costruito
# PROGRAMMATICAMENTE (nessuna riduzione a mano) -> ogni coppia (kind_i x kind_j) = 1 atomo.
# atomi lnσ×ζ = dilog ellittico D; ζ×ζ = C; ecc. Verifica: somma atomi = int diretto.
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
def Q4(x): return np.polyval(Q,x)
k2=((e3-e2)*(e4-e1))/((e4-e2)*(e3-e1)); pref=2/mp.sqrt((e4-e2)*(e3-e1))/mp.sqrt(a4)
om1=mp.mpf(float(pref*mp.ellipk(k2))); w_im=float(pref*mp.ellipk(1-k2)); tau=mp.mpc(0,w_im)/om1; q=mp.exp(mp.pi*1j*tau)
L1=lambda u: mp.jtheta(1,u,q); L1p=lambda u: mp.jtheta(1,u,q,1); L1pp=lambda u: mp.jtheta(1,u,q,2); th1p0=L1p(0)
eta1=-(mp.pi**2/(12*om1))*(mp.jtheta(1,0,q,3)/th1p0)
def wsig(z): u=mp.pi*z/(2*om1); return (2*om1/mp.pi)*mp.exp(eta1*z**2/(2*om1))*L1(u)/th1p0
def wzet(z): u=mp.pi*z/(2*om1); return eta1*z/om1+(mp.pi/(2*om1))*(L1p(u)/L1(u))
def wp(z):   u=mp.pi*z/(2*om1); rr=L1p(u)/L1(u); return -eta1/om1-(mp.pi/(2*om1))**2*(L1pp(u)/L1(u)-rr**2)
def wpp(z):  return (wp(z+1e-8)-wp(z-1e-8))/2e-8
z_inf=float(mp.quad(lambda x:1/mp.sqrt(Q4(x)),[e4,mp.inf])); sa=float(mp.sqrt(a4))
c_r=float(mp.re(e4-(2/sa)*wzet(z_inf))); z_d=z_inf+float(mp.quad(lambda x:1/mp.sqrt(Q4(x)),[-mp.inf,rd]))
def r_of_z(z): return c_r-(1/sa)*(wzet(z-z_inf)-wzet(z+z_inf))
def zr(rv): return float(quad(lambda x:1/np.sqrt(Q4(x)),e4,rv,limit=400)[0])
z0f=mp.mpf(zr(r0))
# --- R=G~' e eta' ---
def R(z): rr=r_of_z(z); return np.polyval(Nm,rr)/((rr-rd)**3*np.polyval(Q,complex(rr)))
def etap(z): rr=r_of_z(z); return (rr**3-2*rr**2)/(rr-rd)
def laurent(fun,a,order,eps=1e-3):
    return {n:complex(mp.quad(lambda th: fun(a+eps*mp.exp(1j*th))*(eps*mp.exp(1j*th))**(n-1)*(1j*eps*mp.exp(1j*th)),[0,2*mp.pi])/(2j*mp.pi)) for n in range(1,order+1)}
PPr={'zd':(z_d,laurent(R,z_d,3)),'mzd':(-z_d,laurent(R,-z_d,3)),'h0':(mp.mpf(0),laurent(R,mp.mpf(0),2)),'h2':(1j*mp.mpf(w_im),laurent(R,1j*mp.mpf(w_im),2))}
PPe={'zd':(z_d,laurent(etap,z_d,1)),'mzd':(-z_d,laurent(etap,-z_d,1)),'zi':(z_inf,laurent(etap,z_inf,2)),'mzi':(-z_inf,laurent(etap,-z_inf,2))}
def const_of(fun,PP,zt):
    return complex(fun(zt)-sum(b.get(1,0)*wzet(zt-a)+b.get(2,0)*wp(zt-a)+b.get(3,0)*(-wpp(zt-a)/2) for a,b in PP.values()))
C0=const_of(R,PPr,mp.mpf('0.19')); Ce=const_of(etap,PPe,mp.mpf('0.19'))
# --- G~ e eta' come LISTE di termini (coeff, kind, pole);  G~ = int R ---
# int di R:  ζ->lnσ(coeff b1), ℘->-ζ(coeff b2), -℘'/2 -> +℘/2? uso G~ principal: b1 lnσ - b2 ζ - (b3/2)℘, +C0 z
G_terms=[(C0,'z',mp.mpf(0))]
for a,b in PPr.values():
    G_terms+=[(b.get(1,0),'lns',a),(-b.get(2,0),'zet',a),(-b.get(3,0)/2,'wp',a)]
Ep_terms=[(Ce,'one',mp.mpf(0))]
for a,b in PPe.items():
    aa=b[0]; bb=b[1]; Ep_terms+=[(bb.get(1,0),'zet',aa),(bb.get(2,0),'wp',aa)]
def gfun(kind,p,z):
    if kind=='z': return z
    if kind=='one': return mp.mpf(1)
    if kind=='lns': return mp.log(wsig(z-p))
    if kind=='zet': return wzet(z-p)
    if kind=='wp':  return wp(z-p)
# verifica ricostruzioni
def Gt_rec(z): return sum(c*gfun(k,p,z) for c,k,p in G_terms)-sum(c*gfun(k,p,z0f) for c,k,p in G_terms if k=='z')*0
def Ep_rec(z): return sum(c*gfun(k,p,z) for c,k,p in Ep_terms)
print("check G~'(=R): ", max(abs(complex((Gt_rec(mp.mpf(t)+mp.mpf('1e-7'))-Gt_rec(mp.mpf(t)-mp.mpf('1e-7')))/mp.mpf('2e-7')-R(mp.mpf(t)))) for t in [0.15,0.22]))
print("check eta'_rec:", max(abs(complex(Ep_rec(mp.mpf(t))-etap(mp.mpf(t)))) for t in [0.15,0.22,0.30]))
# --- ATOMI: coppia (kind_i x kind_j), integrale ∫_{z0}^z prod dz ---
_cache={}
def atom(ki,pi,kj,pj,z):
    key=(ki,float(pi.real),float(pi.imag),kj,float(pj.real),float(pj.imag),float(z))
    if key in _cache: return _cache[key]
    v=mp.quad(lambda t: gfun(ki,pi,t)*gfun(kj,pj,t),[z0f,z]); _cache[key]=v; return v
def E2_closed(z):
    return sum(ci*cj*atom(ki,pi,kj,pj,z) for ci,ki,pi in G_terms for cj,kj,pj in Ep_terms)
def E2_direct(z): return mp.quad(lambda t: Gt_rec(t)*etap(t),[z0f,z])
print("\n=== delta phi weight-2 = int G~ eta' dz : SOMMA ATOMI vs DIRETTO ===")
for rr in [11.0,10.5,10.0]:
    z=mp.mpf(zr(rr)); ec=complex(E2_closed(z)); ed=complex(E2_direct(z))
    print(f"  r={rr:4.1f}  atomi={ec.real:+.8f}  diretto={ed.real:+.8f}  diff={abs(ec-ed):.1e}")
# classifica atomi dilog (lns x zet)
print("\natomi irriducibili (dilog ellittico) presenti = coppie lnσ×ζ:")
seen=set()
for ci,ki,pi in G_terms:
    for cj,kj,pj in Ep_terms:
        if ki=='lns' and kj=='zet' and abs(ci*cj)>1e-6:
            print(f"  D(a={complex(pi):+.4f}, b={complex(pj):+.4f})  coeff={ (ci*cj):+.5f}")

# ============ delta phi COMPLETO = G~(z) eta(r) - E2_closed(z)  vs DIRETTO VERO ============
def sQ(x): return np.sqrt(Q4(x))
def U(x,k): return quad(lambda t:t**k/((t-rd)*sQ(t)),r0,x,limit=200)[0]
def eta(x): return U(x,3)-2*U(x,2)
def dmF(x):
    def Fm(mv): return np.polyval(np.polymul([Jc],[E**2-1,2*mv]),x)/np.sqrt(np.polyval(poly_S_m(mv,Jc),x))
    return (Fm(m+1e-6)-Fm(m-1e-6))/2e-6
def dphi_direct(x): return quad(lambda t: dmF(t)*eta(t), r0, x, limit=200)[0]
def Gt_at(x):  # G~(z(r)) GREZZO (IBP indipendente da costante additiva)
    z=mp.mpf(zr(x)); return complex(Gt_rec(z)).real
print("\n=== delta phi_tau|_sep COMPLETO: [G~ η - Σ atomi] vs diretto ===")
for rr in [11.0,10.5,10.0]:
    z=mp.mpf(zr(rr))
    dclosed=Gt_at(rr)*eta(rr)-complex(E2_closed(z)).real
    ddir=dphi_direct(rr)
    print(f"  r={rr:4.1f}  chiuso(G~η-atomi)={dclosed:+.8f}  diretto={ddir:+.8f}  diff={abs(dclosed-ddir):.1e}")
