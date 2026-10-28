# -*- coding: utf-8 -*-
# Genus-2 dilog Q-SERIES -- STEP 3: Omega = log theta-ratio as a convergent NOME SERIES.
# tau has q12=6.7 (>1). Use S=Fricke: tau'=-tau^{-1} (nomes <1). Modular transformation of theta:
#   theta[delta](z;tau) = pref(z) * theta[S delta](tau^{-1} z; tau') * exp(-i pi z^T tau^{-1} z)
# In the RATIO Omega the z-independent prefactors cancel:
#   Omega = [log theta'(t^{-1}(w-e+)) - log theta'(t^{-1}(w-e-))] - i pi[(w-e+)^T t^{-1}(w-e+)
#            - (w-e-)^T t^{-1}(w-e-)]      (theta' = nome series in tau', converges)
# Check: Omega(nome series) vs Omega(direct, RiemannTheta at tau).
from sage.all import QQ
from sage.schemes.riemann_surfaces.riemann_surface import RiemannSurface as SageRS
from abelfunctions import RiemannTheta
import numpy as np
from scipy.integrate import quad
import itertools

M,a,J,Ehat = 1.0,0.9,2.5,1.4
Rs = PolynomialRing(QQ,['s','Y']); s,Y=Rs.gens()
lam=[QQ(1200),QQ(-2300),QQ(-11428),QQ(-5519),QQ(24700),QQ(62500)]
qpoly=sum(lam[i]*s**i for i in range(6)); X=SageRS(Y**2-qpoly,prec=80)
omega=np.array(X.matrix_of_integral_values([Rs(1),s]),dtype=complex)[:,:2]
tau=np.array(X.riemann_matrix(),dtype=complex); ominv=np.linalg.inv(omega)
taup=-np.linalg.inv(tau); tinv=np.linalg.inv(tau)
print("nomes of tau':",[f"{abs(np.exp(1j*np.pi*taup[i,j])):.4f}" for i,j in [(0,0),(1,1),(0,1)]])

halfs=[np.array(v)/2 for v in itertools.product([0,1],repeat=2)]
def theta_af(z,av,bv,t):
    av=np.array(av); bv=np.array(bv); zz=z+t@av+bv
    pref=np.exp(1j*np.pi*(av@t@av)+2j*np.pi*(av@(z+bv)))
    return pref*complex(RiemannTheta(zz,t))
def theta_latt(z,av,bv,t,N=7):  # NOME series (lattice sum, in tau')
    tot=0j
    for n in itertools.product(range(-N,N+1),repeat=2):
        nn=np.array(n)+np.array(av)
        tot+=np.exp(1j*np.pi*(nn@t@nn)+2j*np.pi*(nn@(np.array(z)+np.array(bv))))
    return tot
odd=[(av,bv) for av in halfs for bv in halfs if abs(theta_af(np.zeros(2),av,bv,tau))<1e-6]
av,bv=odd[1]
# S-transformed characteristic: try them all, pick the one that matches the ratio
def qn(sv): return sum(float(lam[i])*sv**i for i in range(6))
rmin=4.046197656444178; s_b=1.0/rmin
def Iu(s_to):
    U=np.sqrt(s_b-s_to)
    g0=lambda u:2*u/np.sqrt(qn(s_b-u**2)); g1=lambda u:2*u*(s_b-u**2)/np.sqrt(qn(s_b-u**2))
    return -np.array([quad(g0,0,U,limit=200)[0],quad(g1,0,U,limit=200)[0]])
def w_of(sv): return ominv@Iu(sv)
e_plus=w_of(0.0); e_minus=-e_plus
def Omega_direct(rv):
    w=w_of(1.0/rv); return np.log(theta_af(w-e_plus,av,bv,tau)/theta_af(w-e_minus,av,bv,tau))
# theta with ZERO characteristic (lattice sum in tau'), absorb delta into the argument:
#  theta[delta](z;tau) = theta(z+tau*a+b; tau) * exp(i pi a^T tau a + 2 pi i a^T (z+b))
#  -> Omega = [log theta0(zeta1;tau)-log theta0(zeta2;tau)] + 2 pi i a^T (e_minus - e_plus)
#  zeta_i = (w-e_pm) + tau*a + b ; then S on theta0 (zero char, clean):
#  theta0(zeta;tau) = C exp(-i pi zeta^T tinv zeta) theta0(tinv*zeta; tau')
def theta0_latt(z,t,N=8):
    tot=0j
    for n in itertools.product(range(-N,N+1),repeat=2):
        nn=np.array(n); tot+=np.exp(1j*np.pi*(nn@t@nn)+2j*np.pi*(nn@np.array(z)))
    return tot
def Omega_nome(rv):
    w=w_of(1.0/rv); z1=w-e_plus; z2=w-e_minus
    zeta1=z1+tau@av+bv; zeta2=z2+tau@av+bv
    # S-transform theta0 (char zero): log ratio = -i pi(zeta1 tinv zeta1 - zeta2 tinv zeta2)
    #   + log[theta0(tinv zeta1;tau')/theta0(tinv zeta2;tau')]
    Lth=np.log(theta0_latt(tinv@zeta1,taup)/theta0_latt(tinv@zeta2,taup))
    quad_=-1j*np.pi*(zeta1@tinv@zeta1 - zeta2@tinv@zeta2)
    phase=2j*np.pi*(av@(e_minus-e_plus))    # characteristic phase (elementary)
    return Lth+quad_+phase
print("\n=== STEP 3: Omega nome series (tau', zero char + S) vs direct (mod 2pi i, log branch) ===")
def modbranch(d):  # reduce the imaginary part mod 2pi (branch of the complex log)
    return d.real + 1j*((d.imag+np.pi)%(2*np.pi)-np.pi)
errs=[]
for rv in [10.0,8.0,6.5,5.0]:
    on=Omega_nome(rv); od=Omega_direct(rv); d=modbranch(on-od); errs.append(abs(d))
    print(f"  r={rv}: Omega_nome={on:+.5f}  Omega_dir={od:+.5f}  diff(mod 2pi i)={abs(d):.2e}")
if max(errs)<1e-5:
    print("\n=> Omega = log theta-ratio HAS a convergent NOME series (tau'=-tau^{-1}, zero char+S).")
    print("   third kind of the dilog in nomes. STEP 4: Kronecker-Eisenstein kernel for Lambda.")
else:
    print(f"\n=> residual diff {max(errs):.1e}: check the quadratic term/log branch/det.")
