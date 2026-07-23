# -*- coding: utf-8 -*-
# NAMING peso-1 a livello INTEGRALE (TK tau, genus-2): U_k = funzioni tabulate theta[delta]
# agli e_pm. Le primitive sono valutate DIRETTAMENTE (log theta-ratio, zeta_delta), non
# integrando le derivate -> verifica vera. Poi assemblo il clock eta=U3-2M U2 nominato.
#   U_k(r) = a_k u1 + b_k u2 + c_k L_rat + d_k Zp1 + e_k Zp2 + f_k Zm1 + g_k Zm2 + const
#   L_rat = log[theta[d](w-e+)/theta[d](w-e-)] (3a specie), Zpm = zeta_delta(w-e_pm) (2a).
from sage.all import QQ
from sage.schemes.riemann_surfaces.riemann_surface import RiemannSurface as SageRS
from abelfunctions import RiemannTheta
import numpy as np
from scipy.integrate import quad
import itertools

M,a,J,Ehat = 1.0,0.9,2.5,1.4
Rs = PolynomialRing(QQ, ['s','Y']); s,Y = Rs.gens()
lam=[QQ(1200),QQ(-2300),QQ(-11428),QQ(-5519),QQ(24700),QQ(62500)]
qpoly=sum(lam[i]*s**i for i in range(6))
X=SageRS(Y**2-qpoly, prec=80)
omega=np.array(X.matrix_of_integral_values([Rs(1),s]),dtype=complex)[:,:2]
tau=np.array(X.riemann_matrix(),dtype=complex); ominv=np.linalg.inv(omega)
def qn(sv): return sum(float(lam[i])*sv**i for i in range(6))
rmin=4.046197656444178; s_b=1.0/rmin
def Iu(s_to):
    U=np.sqrt(s_b-s_to)
    g0=lambda u:2*u*1.0/np.sqrt(qn(s_b-u**2)); g1=lambda u:2*u*(s_b-u**2)/np.sqrt(qn(s_b-u**2))
    return -np.array([quad(g0,0,U,limit=200)[0],quad(g1,0,U,limit=200)[0]])
def w_of(s_to): return ominv@Iu(s_to)
e_plus=w_of(0.0); e_minus=-e_plus
def theta_d(z,av,bv,order=0):
    av=np.array(av); bv=np.array(bv); zz=z+tau@av+bv
    pref=np.exp(1j*np.pi*(av@tau@av)+2j*np.pi*(av@(z+bv))); th=complex(RiemannTheta(zz,tau))
    if order==0: return pref*th
    g=np.array([complex(RiemannTheta(zz,tau,derivs=[e])) for e in ([1,0],[0,1])])
    if order==1: return 2j*np.pi*av+g/th
    H=np.array([[complex(RiemannTheta(zz,tau,derivs=[p,q])) for q in ([1,0],[0,1])] for p in ([1,0],[0,1])])
    return H/th-np.outer(g/th,g/th)
halfs=[np.array(v)/2 for v in itertools.product([0,1],repeat=2)]
odd=[(av,bv) for av in halfs for bv in halfs if abs(theta_d(np.zeros(2),av,bv))<1e-6]
av,bv=odd[1]   # delta #1 (migliore)
print("uso delta #1 a=",av," b=",bv)

def Sn(rv):
    Dl=rv**2-2*M*rv+a**2; Emu=(Ehat**2-1)*rv+2*M
    return rv*(rv-2*M)*Emu*(rv*Dl-J**2*Emu)
# --- (1) coefficienti dal DIFFERENZIALE: r^k/sqrtS = sum coeff * basis_diff ---
rg=np.linspace(11.5,rmin+0.4,22); sr=1.0/rg; dsdr=-1.0/rg**2
du1dr=(1/np.sqrt(qn(sr)))*dsdr; du2dr=(sr/np.sqrt(qn(sr)))*dsdr
W=np.array([w_of(sv) for sv in sr]); dwdr=(ominv@np.vstack([du1dr,du2dr])).T; sqrtS=np.sqrt(Sn(rg))
D3=np.array([ (theta_d(W[i]-e_plus,av,bv,1)-theta_d(W[i]-e_minus,av,bv,1))@dwdr[i] for i in range(len(rg))])
Gp=np.array([ theta_d(W[i]-e_plus,av,bv,2)@dwdr[i] for i in range(len(rg))])
Gm=np.array([ theta_d(W[i]-e_minus,av,bv,2)@dwdr[i] for i in range(len(rg))])
Bdiff=np.column_stack([du1dr,du2dr,D3,Gp[:,0],Gp[:,1],Gm[:,0],Gm[:,1]])
coeffs={}
for k in [2,3,4]:
    coef,_,_,_=np.linalg.lstsq(Bdiff,rg**k/sqrtS,rcond=None)
    coeffs[k]=coef
    print(f"  k={k}: coeff(u1,u2,Lrat,Zp1,Zp2,Zm1,Zm2)=",np.round(coef,4))

# --- (2) PRIMITIVE nominate valutate DIRETTAMENTE (non integrando derivate) ---
def named_prims(rv):
    sv=1.0/rv; w=w_of(sv)
    Lrat=np.log(theta_d(w-e_plus,av,bv,0)/theta_d(w-e_minus,av,bv,0))
    Zp=theta_d(w-e_plus,av,bv,1); Zm=theta_d(w-e_minus,av,bv,1)
    # u1,u2 (olomorfe) = integrale abeliano diretto (coordinate Jacobiane, "nominate" = Abel map)
    return w,Lrat,Zp,Zm
def u_hol(rv):  # (u1,u2) unnormalizzate base r0
    s_to=1.0/rv
    I0=quad(lambda t:1.0/np.sqrt(qn(t)),1.0/11.5,s_to,limit=200)[0]
    I1=quad(lambda t:t/np.sqrt(qn(t)),1.0/11.5,s_to,limit=200)[0]
    return np.array([I0,I1])
def named_vec(rv):   # vettore [u1,u2,Lrat,Zp1,Zp2,Zm1,Zm2]
    _,Lrat,Zp,Zm=named_prims(rv); uh=u_hol(rv)
    return np.array([uh[0],uh[1],Lrat,Zp[0],Zp[1],Zm[0],Zm[1]])
def Uk_direct(k,rv): return quad(lambda x:x**k/np.sqrt(Sn(x)),11.5,rv,limit=200)[0]

r0=11.5; base=named_vec(r0)
print("\n=== VERIFICA INTEGRALE: U_k(r) = coeff . [named(r)-named(r0)]  (floor theta ~1e-4) ===")
rtest=[10.0,8.0,6.0,4.6]
for k in [2,3,4]:
    print(f"  k={k}:")
    for rv in rtest:
        rec=complex(coeffs[k]@(named_vec(rv)-base)).real
        dr=Uk_direct(k,rv); print(f"    r={rv}: U_k(named)={rec:+.6f}  U_k(dir)={dr:+.6f}  diff={abs(rec-dr):.1e}")
# --- (3) clock eta = U3 - 2M U2 nominato ---
print("\n=== clock eta = U3 - 2M U2 in forma nominata theta[delta] ===")
ceta=coeffs[3]-2*M*coeffs[2]
print("  coeff eta (u1,u2,Lrat,Zp1,Zp2,Zm1,Zm2)=",np.round(ceta,4))
for rv in rtest:
    rec=complex(ceta@(named_vec(rv)-base)).real
    dr=Uk_direct(3,rv)-2*M*Uk_direct(2,rv); print(f"    r={rv}: eta(named)={rec:+.6f} eta(dir)={dr:+.6f} diff={abs(rec-dr):.1e}")
print("\n=> se diff~1e-4: clock e U_k NOMINATI in theta[delta] agli e_pm (funzioni tabulate).")
print("   3a specie=log theta-ratio; 2a specie=zeta_delta. phi0 e clock CHIUSI. Coeff dal differenziale.")
