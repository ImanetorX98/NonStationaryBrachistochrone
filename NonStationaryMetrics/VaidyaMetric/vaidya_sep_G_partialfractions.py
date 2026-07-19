# ALGORITMO (Hermite-Ostrogradsky via parti principali sul toro) - NIENTE ansatz/gradi.
# G(z)=int dm F = C0 z + sum_a [ b1_a lnσ(z-a) - b2_a ζ(z-a) - (b3_a/2) ℘(z-a) ]
# b_n^a = (1/2πi) ∮ R(z)(z-a)^{n-1} dz,  R(z)=N_m(r(z))/((r(z)-r_d)^3 Q(r(z))).
import numpy as np, mpmath as mp, sympy as sp
from scipy.integrate import quad
mp.mp.dps=30
E=1.4; m=1.0; r0=12.0
# --- Jc, r_d, S, Q, N_m (numerico robusto) ---
r,J=sp.symbols('r J'); Es=sp.Rational(7,5)
Ssym=sp.expand(r*(r-2)*((Es**2-1)*r+2)*(r**2*(r-2)-J**2*((Es**2-1)*r+2)))
_Js=sp.solve(sp.Eq(sp.resultant(Ssym,sp.diff(Ssym,r),r),0),J)
Jc=float([s for s in _Js if s.is_real and float(s)>1][0])
def poly_S_m(mv,Jv):
    DE=np.array([E**2-1,2*mv]); p=np.polymul(np.polymul([1,0],[1,-2*mv]),DE)
    br=np.polysub(np.polymul([1,0,0],[1,-2*mv]),Jv**2*DE); return np.polymul(p,br)
Sc=poly_S_m(m,Jc)
hh=1e-6; dmS=(poly_S_m(m+hh,Jc)-poly_S_m(m-hh,Jc))/(2*hh)
K=np.polymul([Jc],[E**2-1,2.0]); dmK=np.array([2*Jc])
Nm=np.polysub(np.polymul(Sc,dmK),0.5*np.polymul(K,dmS))
rts=np.roots(Sc); pr=[(i,j) for i in range(6) for j in range(i+1,6) if abs(rts[i]-rts[j])<1e-6]
rd=float(np.real((rts[pr[0][0]]+rts[pr[0][1]])/2))
Q,_=np.polydiv(Sc,np.polymul([1,-rd],[1,-rd]))          # deg4, = a4c*Q4monic
er=np.sort(np.real(np.roots(Q))); a4=Q[0]; e1,e2,e3,e4=er
def Q4(x): return np.polyval(Q,x)
print(f"Jc={Jc:.8f} r_d={rd:.8f} radici Q={er} a4={a4:.5f}")
# --- Weierstrass da theta1 (periodi 2om1,2w_im) ---
k2=((e3-e2)*(e4-e1))/((e4-e2)*(e3-e1)); pref=2/mp.sqrt((e4-e2)*(e3-e1))/mp.sqrt(a4)
om1=mp.mpf(float(pref*mp.ellipk(k2))); w_im=float(pref*mp.ellipk(1-k2)); tau=mp.mpc(0,w_im)/om1; q=mp.exp(mp.pi*1j*tau)
L1=lambda u: mp.jtheta(1,u,q); L1p=lambda u: mp.jtheta(1,u,q,1); L1pp=lambda u: mp.jtheta(1,u,q,2); th1p0=L1p(0)
eta1=-(mp.pi**2/(12*om1))*(mp.jtheta(1,0,q,3)/th1p0)
def wsig(z): u=mp.pi*z/(2*om1); return (2*om1/mp.pi)*mp.exp(eta1*z**2/(2*om1))*L1(u)/th1p0
def wzet(z): u=mp.pi*z/(2*om1); return eta1*z/om1+(mp.pi/(2*om1))*(L1p(u)/L1(u))
def wp(z):   u=mp.pi*z/(2*om1); rr=L1p(u)/L1(u); return -eta1/om1-(mp.pi/(2*om1))**2*(L1pp(u)/L1(u)-rr**2)
def wpp(z):  return (wp(z+1e-8)-wp(z-1e-8))/2e-8
z_inf=float(mp.quad(lambda x:1/mp.sqrt(Q4(x)),[e4,mp.inf]))
sa=float(mp.sqrt(a4)); c_r=float(mp.re(e4-(2/sa)*wzet(z_inf)))
z_d=z_inf+float(mp.quad(lambda x:1/mp.sqrt(Q4(x)),[-mp.inf,rd]))
def r_of_z(z): return c_r-(1/sa)*(wzet(z-z_inf)-wzet(z+z_inf))    # r(z) esplicito
# --- R(z) elliptic ---
def Nmv(x): return np.polyval(Nm,x)
def R(z):
    rr=r_of_z(z); return Nmv(rr)/((rr-rd)**3*np.polyval(Q,complex(rr)))
# --- estrazione parti principali via contorno ---
def laurent(a,order,eps=1e-3,Npt=4000):
    b={}
    for n in range(1,order+1):
        f=lambda th: R(a+eps*mp.exp(1j*th))*(eps*mp.exp(1j*th))**(n-1)*(1j*eps*mp.exp(1j*th))
        val=mp.quad(f,[0,2*mp.pi])/(2j*mp.pi)
        b[n]=complex(val)
    return b
poles_o3=[('zd',z_d),('mzd',-z_d)]
halfp=[('h0',mp.mpf(0)),('h1',om1),('h2',1j*mp.mpf(w_im)),('h3',om1+1j*mp.mpf(w_im))]
PP={}
for nm,a in poles_o3: PP[nm]=(a,laurent(a,3))
for nm,a in halfp:    PP[nm]=(a,laurent(a,2))
print("\nparti principali (b1,b2,b3):")
for nm in PP:
    a,b=PP[nm]; print(f"  {nm:4s} a={complex(a):+.4f}  b1={b.get(1,0):+.4e} b2={b.get(2,0):+.4e} b3={b.get(3,0):+.4e}")
# --- C0 da match a z regolare ---
zt=mp.mpf('0.19')
def princ(z):
    s=0
    for nm,(a,b) in PP.items():
        s+= b.get(1,0)*wzet(z-a)+b.get(2,0)*wp(z-a)-(b.get(3,0)/2)*wpp(z-a)
    return s
C0=complex(R(zt)-princ(zt))
print("\nC0=",C0)
# --- G(z) esplicito ---
def G(z):
    s=C0*z
    for nm,(a,b) in PP.items():
        s+= b.get(1,0)*mp.log(wsig(z-a)) - b.get(2,0)*wzet(z-a) - (b.get(3,0)/2)*wp(z-a)
    return s
# --- verifica dG/dz = R ---
print("\nverifica dG/dz = R(z):")
for zt in [0.19,0.30,0.12]:
    zt=mp.mpf(zt); dG=(G(zt+mp.mpf('1e-7'))-G(zt-mp.mpf('1e-7')))/mp.mpf('2e-7')
    print(f"  z={float(zt):.2f}  dG={complex(dG):+.6f}  R={complex(R(zt)):+.6f}  diff={abs(complex(dG-R(zt))):.1e}")

# ================= ASSEMBLAGGIO FINALE delta phi|_sep =================
# delta phi = G~(x) eta(x) - int_{r0}^x G~(r) eta'(r) dr,  G~=G-G(z0),  eta=U3-2U2,
#   eta'(r)=(r^3-2r^2)/sqrtS = (r^3-2r^2)/((r-r_d) sqrtQ);  in z: eta'dz-form=(r^3-2r^2)/(r-r_d).
def zr(rv): return float(quad(lambda x:1/np.sqrt(Q4(x)), e4, rv, limit=400)[0])
z0=mp.mpf(zr(r0))
sQ=lambda x: np.sqrt(Q4(x))
def U(x,k): return quad(lambda t:t**k/((t-rd)*sQ(t)),r0,x,limit=200)[0]
def eta(x): return U(x,3)-2*U(x,2)
Gt=lambda z: G(z)-G(z0)                                  # G~ con base r0
# dm F diretto
def dmF(x):
    def Fm(mv): return np.polyval(np.polymul([Jc],[E**2-1,2*mv]),x)/np.sqrt(np.polyval(poly_S_m(mv,Jc),x))
    return (Fm(m+1e-6)-Fm(m-1e-6))/2e-6
def dphi_direct(x): return quad(lambda t: dmF(t)*eta(t), r0, x, limit=200)[0]
def dphi_asm(x):
    zx=mp.mpf(zr(x))
    bnd=float(mp.re(Gt(zx)))*eta(x)
    intg=quad(lambda t: float(mp.re(Gt(mp.mpf(zr(t)))))*(t**3-2*t**2)/((t-rd)*sQ(t)), r0, x, limit=120)[0]
    return bnd-intg
print("\n=== delta phi_tau|_sep : ASSEMBLATO (G esplicito) vs DIRETTO ===")
for x in [11.0,10.0,9.2,8.9,8.8]:
    da=dphi_asm(x); dd=dphi_direct(x)
    print(f"  r={x:4.2f}  asm={da:+.8f}  dir={dd:+.8f}  diff={abs(da-dd):.1e}")
# valore separatrice al turning (limite)
import numpy as _np
xt=e4-1e-4
print(f"\n  delta phi_tau|_sep (r->turning e4={e4:.4f}) ~ dir={dphi_direct(e4-1e-3):.8f}")

# ============ CHIUSURA WEIGHT-2: int G~ eta' dz  ->  atomi dilog ellittico ============
# eta'-forma = (r^3-2r^2)/(r-r_d) = funzione ellittica meromorfa. Parti principali (contorno):
def etap(z):
    rr=r_of_z(z); return (rr**3-2*rr**2)/(rr-rd)
def laurent_f(fun,a,order,eps=1e-3):
    b={}
    for n in range(1,order+1):
        f=lambda th: fun(a+eps*mp.exp(1j*th))*(eps*mp.exp(1j*th))**(n-1)*(1j*eps*mp.exp(1j*th))
        b[n]=complex(mp.quad(f,[0,2*mp.pi])/(2j*mp.pi))
    return b
# poli eta': z_d,-z_d (ord1 da 1/(r-rd)); z_inf,-z_inf (ord2 da r^3)
PPe={'zd':(z_d,laurent_f(etap,z_d,1)), 'mzd':(-z_d,laurent_f(etap,-z_d,1)),
     'zi':(z_inf,laurent_f(etap,z_inf,2)), 'mzi':(-z_inf,laurent_f(etap,-z_inf,2))}
# costante media di eta'
zt=mp.mpf('0.19')
def prcE(z):
    s=0
    for nm,(a,b) in PPe.items(): s+= b.get(1,0)*wzet(z-a)+b.get(2,0)*wp(z-a)
    return s
Ce=complex(etap(zt)-prcE(zt))
print("\n=== eta' parti principali ===")
for nm,(a,b) in PPe.items(): print(f"  {nm}: a={complex(a):+.4f} b1={b.get(1,0):+.4e} b2={b.get(2,0):+.4e}")
print("  Ce=",Ce)
# verifica ricostruzione eta' = Ce + sum[b1 zeta + b2 P]
def etap_rec(z):
    s=Ce
    for nm,(a,b) in PPe.items(): s+= b.get(1,0)*wzet(z-a)+b.get(2,0)*wp(z-a)
    return s
print("  check eta'_rec vs eta':", max(abs(complex(etap_rec(mp.mpf(t))-etap(mp.mpf(t)))) for t in [0.12,0.19,0.30]))

# --- atomi weight-2 (definizione = valore funzione speciale, via quadratura) ---
z0f=mp.mpf(zr(r0))
def D_at(a,b):   # int_{z0}^{z} lnσ(z'-a) ζ(z'-b) dz'   (dilog ellittico)
    return lambda z: mp.quad(lambda t: mp.log(wsig(t-a))*wzet(t-b),[z0f,z])
def C_at(a,b):   # int ζ(z'-a) ζ(z'-b) dz'
    return lambda z: mp.quad(lambda t: wzet(t-a)*wzet(t-b),[z0f,z])
def Pz_at(b):    # int z' ζ(z'-b) dz'
    return lambda z: mp.quad(lambda t: t*wzet(t-b),[z0f,z])
def Lns_at(a):   # int lnσ(z'-a) dz'
    return lambda z: mp.quad(lambda t: mp.log(wsig(t-a)),[z0f,z])

# G~ = C0 z + sum_a [ B1 lnσ(z-a) - B2 ζ(z-a) - (B3/2) P(z-a) ]   (dai coeff PP di R)
# Assemblo  int G~ eta' dz  usando eta' = Ce + sum_p[ e1_p ζ(z-zp) + e2_p P(z-zp) ].
# Prodotti -> atomi:
#   [B1 lnσ(z-a)]*[e1_p ζ(z-zp)]     -> B1 e1_p  D(a,zp)
#   [B1 lnσ(z-a)]*[Ce]               -> B1 Ce    Lns(a)
#   [B1 lnσ(z-a)]*[e2_p P(z-zp)]     -> IBP: B1 e2_p [ -D'?]  (P=-ζ') -> ridotto sotto
#   [-B2 ζ(z-a)]*[e1_p ζ(z-zp)]      -> -B2 e1_p  C(a,zp)
#   [-B2 ζ(z-a)]*[Ce]                -> -B2 Ce lnσ(z-a)  (weight-1, = ∫ζ)
#   [C0 z]*[e1_p ζ(z-zp)]            -> C0 e1_p  Pz(zp)
#   ... (termini con P e con -B3/2 P: ridotti via IBP a D,C + weight-1; qui li includo numericamente come atomi P)
def PP_at(a,b):  # int P(z-a) ζ(z-b) dz   (atomo con Weierstrass P)
    return lambda z: mp.quad(lambda t: wp(t-a)*wzet(t-b),[z0f,z])
def zP_at(b):    # int z P(z-b) dz
    return lambda z: mp.quad(lambda t: t*wp(t-b),[z0f,z])
def LnP_at(a,b): # int lnσ(z-a) P(z-b) dz
    return lambda z: mp.quad(lambda t: mp.log(wsig(t-a))*wp(t-b),[z0f,z])
def PzP_at(a,b): return lambda z: mp.quad(lambda t: wp(t-a)*wp(t-b),[z0f,z])

# coeff di R (=G~'): B1,B2,B3 per polo (da PP di R)
GR=[(PP[nm][0], PP[nm][1]) for nm in PP]   # (a, {1,2,3})
def E2_closed(z):
    tot=0
    for (a,B) in GR:
        B1=B.get(1,0); B2=B.get(2,0); B3=B.get(3,0)
        for nm,(zp,e) in PPe.items():
            e1=e.get(1,0); e2=e.get(2,0)
            tot+= B1*e1*D_at(a,zp)(z)   - B2*e1*C_at(a,zp)(z)
            tot+= B1*e2*LnP_at(a,zp)(z) - B2*e2*PP_at(a,zp)(z)
            tot+= -(B3/2)*e1*PP_at(a,zp)(z) - (B3/2)*e2*PzP_at(a,zp)(z)
        # termini * Ce
        tot+= B1*Ce*Lns_at(a)(z) - B2*Ce*mp.log(wsig(z-a)) - (B3/2)*Ce*wzet(z-a)*0 - (B3/2)*Ce*(-wzet(z-a))
    # C0 z * eta'
    for nm,(zp,e) in PPe.items():
        e1=e.get(1,0); e2=e.get(2,0)
        tot+= C0*e1*Pz_at(zp)(z) + C0*e2*zP_at(zp)(z)
    tot+= C0*Ce*(z**2/2 - z0f**2/2)
    return tot

# --- VERIFICA CHIUSURA: E2_closed (somma atomi nominati) vs E2_direct (quadratura) ---
print("\n=== CHIUSURA weight-2: int G~ eta' dz  = somma atomi dilog ellittico ? ===")
def E2_direct(z):
    return mp.quad(lambda t: mp.re(Gt(t))*etap(t),[z0f,z])
for rr in [11.0,10.5,10.0]:
    z=mp.mpf(zr(rr))
    ec=complex(E2_closed(z)); ed=complex(E2_direct(z))
    print(f"  r={rr:4.1f}  E2_closed={ec.real:+.8f}  E2_direct={ed.real:+.8f}  diff={abs(ec-ed):.1e}")
