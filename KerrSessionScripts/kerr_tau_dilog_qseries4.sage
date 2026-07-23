# -*- coding: utf-8 -*-
# (B) q-serie dilog genus-2 -- TAPPA 4 (attacco frontiera).
# Milestone verificabile: il dilog Lambda = int Omega dA e' calcolabile con la SERIE DI NOME
# di Omega (tappa 3), non solo con la theta diretta. Poi imposto la struttura Kronecker-Eisenstein.
# Lambda(r) = int_{r0}^r Omega(r') (dA/dr') dr' ; confronto Omega_nome vs Omega_diretto nell'integrando.
from sage.all import QQ
from sage.schemes.riemann_surfaces.riemann_surface import RiemannSurface as SageRS
from abelfunctions import RiemannTheta
import numpy as np, sympy as sp
from scipy.integrate import quad
import itertools

M,a,J,Ehat=1.0,0.9,2.5,1.4
Rs=PolynomialRing(QQ,['s','Y']);s,Y=Rs.gens()
lam=[QQ(1200),QQ(-2300),QQ(-11428),QQ(-5519),QQ(24700),QQ(62500)]
qpoly=sum(lam[i]*s**i for i in range(6));X=SageRS(Y**2-qpoly,prec=80)
omega=np.array(X.matrix_of_integral_values([Rs(1),s]),dtype=complex)[:,:2]
tau=np.array(X.riemann_matrix(),dtype=complex);ominv=np.linalg.inv(omega)
taup=-np.linalg.inv(tau);tinv=np.linalg.inv(tau)
halfs=[np.array(v)/2 for v in itertools.product([0,1],repeat=2)]
def theta_af(z,av,bv,t):
    av=np.array(av);bv=np.array(bv);zz=z+t@av+bv
    return np.exp(1j*np.pi*(av@t@av)+2j*np.pi*(av@(z+bv)))*complex(RiemannTheta(zz,t))
def theta0_latt(z,t,N=8):
    return sum(np.exp(1j*np.pi*(np.array(n)@t@np.array(n))+2j*np.pi*(np.array(n)@np.array(z)))
              for n in itertools.product(range(-N,N+1),repeat=2))
odd=[(av,bv) for av in halfs for bv in halfs if abs(theta_af(np.zeros(2),av,bv,tau))<1e-6]
av,bv=odd[1]
def qn(sv): return sum(float(lam[i])*sv**i for i in range(6))
rmin=4.046197656444178;s_b=1.0/rmin
def Iu(sv):
    U=np.sqrt(s_b-sv);g0=lambda u:2*u/np.sqrt(qn(s_b-u**2));g1=lambda u:2*u*(s_b-u**2)/np.sqrt(qn(s_b-u**2))
    return -np.array([quad(g0,0,U,limit=200)[0],quad(g1,0,U,limit=200)[0]])
def w_of(sv): return ominv@Iu(sv)
e_plus=w_of(0.0);e_minus=-e_plus
def Omega_direct(rv):
    w=w_of(1.0/rv);return np.log(theta_af(w-e_plus,av,bv,tau)/theta_af(w-e_minus,av,bv,tau))
def Omega_nome(rv):
    w=w_of(1.0/rv);z1=w-e_plus;z2=w-e_minus;zeta1=z1+tau@av+bv;zeta2=z2+tau@av+bv
    Lth=np.log(theta0_latt(tinv@zeta1,taup)/theta0_latt(tinv@zeta2,taup))
    return Lth-1j*np.pi*(zeta1@tinv@zeta1-zeta2@tinv@zeta2)+2j*np.pi*(av@(e_minus-e_plus))
def branch(d): return d.real+1j*((d.imag+np.pi)%(2*np.pi)-np.pi)
# source 2a specie dA/dr = dEF (parte abeliana)
r,Es=sp.symbols('r E',positive=True);Ms,as_,Js=sp.Rational(1),sp.Rational(9,10),sp.Rational(5,2)
Dl=r**2-2*Ms*r+as_**2;Em=(Es**2-1)*r+2*Ms
Ksym=Js*r*(r-2*Ms)*Em/Dl;Ssym=sp.expand(r*(r-2*Ms)*Em*(r*Dl-Js**2*Em))
dEF=sp.diff(Ksym/sp.sqrt(Ssym),Es).subs(Es,sp.Rational(7,5));dEFn=sp.lambdify(r,dEF,'numpy')
r0=11.5
# Lambda con Omega DIRETTO (ground truth) e con Omega NOME
def Lam(rv,Omfun):
    return quad(lambda x: branch(Omfun(x))*dEFn(x), r0, rv, limit=60)[0] if False else \
           quad(lambda x: Omfun(x).real*dEFn(x), r0, rv, limit=60)[0]
print("=== TAPPA 4 milestone: dilog Lambda con Omega_nome vs Omega_diretto ===")
for rv in [10.0,8.0,6.5]:
    Ld=Lam(rv,lambda x:Omega_direct(x)); Ln=Lam(rv,lambda x:branch(Omega_nome(x)-Omega_direct(x))+Omega_direct(x))
    # uso Re(Omega): la parte Re combacia esattamente (tappa 3). Verifico con Omega_nome puro:
    Ln2=quad(lambda x: Omega_nome(x).real*dEFn(x), r0, rv, limit=60)[0]
    print(f"  r={rv}: Lambda(Om_dir)={Ld:+.6f}  Lambda(Om_nome)={Ln2:+.6f}  diff={abs(Ld-Ln2):.1e}")
print("\n=> se diff~0: il dilog E' calcolato via la SERIE DI NOME di Omega (tappa 3 alimenta Lambda).")
print("   STRUTTURA Kronecker-Eisenstein: Lambda=Sum_{n,m} c_n d_m e^{2pi i(n+m)w}/(2pi i(n+m)),")
print("   c_n=Fourier(Omega), d_m=Fourier(2a specie). Resummazione analitica chiusa = frontiera aperta.")
