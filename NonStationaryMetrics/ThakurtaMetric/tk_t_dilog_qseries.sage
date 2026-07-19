# -*- coding: utf-8 -*-
# q-serie dilog genus-2 TK-t generico (M=1,a=9/10,E=6/5,J=5 scattering). Template TK-tau:
# naming theta[delta] + Omega=log theta-ratio nome-serie + dilog split Q(elem)+L(dilog).
# NB: il ramo t ha ANCHE dilog agli orizzonti (da rho_t); qui il pezzo Omega (infinito).
from sage.all import QQ
from sage.schemes.riemann_surfaces.riemann_surface import RiemannSurface as SageRS
from abelfunctions import RiemannTheta
import numpy as np, sympy as sp
from scipy.integrate import quad
import itertools
Rs=PolynomialRing(QQ,['s','Y']);s,Y=Rs.gens()
lam=[QQ(9900),QQ(45000),QQ(-67606),QQ(-538212),QQ(-587325),QQ(2101250)]  # TK-t J=5 odd model
qpoly=sum(lam[i]*s**i for i in range(6)); X=SageRS(Y**2-qpoly,prec=80)
omega=np.array(X.matrix_of_integral_values([Rs(1),s]),dtype=complex)[:,:2]
tau=np.array(X.riemann_matrix(),dtype=complex);ominv=np.linalg.inv(omega)
taup=-np.linalg.inv(tau);tinv=np.linalg.inv(tau)
print("nomi tau':",[f"{abs(np.exp(1j*np.pi*taup[i,j])):.4f}" for i,j in [(0,0),(1,1),(0,1)]])
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
# curva even + sorgente dE F_t
r,a,E,J=sp.symbols('r a E J'); M=1
DE=(E**2-1)*r+2*M; Delta=r**2-2*M*r+a**2
Q2=(2*E**2*J**2*r-E**2*J**2*r**2-4*E**2*J*a*r+2*E**2*a**2*r+E**2*a**2*r**2+E**2*r**4+4*J**2-4*J**2*r+J**2*r**2-8*J*a+4*J*a*r+4*a**2)
R6=sp.expand(r*Q2*DE); Kt=r*DE*(J*(r-2*M)+2*M*a)/Delta
sub={a:sp.Rational(9,10),E:sp.Rational(6,5),J:sp.Rational(5)}
dEF=sp.diff(Kt/sp.sqrt(R6),E).subs(sub); dEFn=sp.lambdify(r,dEF,'numpy')
R6c=[float(c) for c in sp.Poly(R6.subs(sub),r).all_coeffs()]
def R6n(x): return np.polyval(R6c,x)
rmin=3.079; s_b=1.0/rmin
def Iu(sv):
    U=np.sqrt(abs(s_b-sv));g0=lambda u:2*u/np.sqrt(abs(qn(s_b-u**2)));g1=lambda u:2*u*(s_b-u**2)/np.sqrt(abs(qn(s_b-u**2)))
    sg=-1.0 if sv<s_b else 1.0
    return sg*np.array([quad(g0,0,U,limit=200)[0],quad(g1,0,U,limit=200)[0]])
def w_of(sv): return ominv@Iu(sv)
e_plus=w_of(0.0);e_minus=-e_plus
def Q_part(rv):
    w=w_of(1.0/rv);z1=w-e_plus;z2=w-e_minus;zeta1=z1+tau@av+bv;zeta2=z2+tau@av+bv
    return (-1j*np.pi*(zeta1@tinv@zeta1-zeta2@tinv@zeta2)+2j*np.pi*(av@(e_minus-e_plus))).real
def L_partN(rv,N=8):
    w=w_of(1.0/rv);zeta1=(w-e_plus)+tau@av+bv;zeta2=(w-e_minus)+tau@av+bv
    return np.log(theta0_latt(tinv@zeta1,taup,N)/theta0_latt(tinv@zeta2,taup,N)).real
r0=20.0
print("\n=== TK-t dilog (pezzo Omega, infinito): split Q(elem)+L(dilog) ===")
for rv in [15.,10.,5.]:
    LQ=quad(lambda x:Q_part(x)*dEFn(x),r0,rv,limit=60)[0]
    LL=quad(lambda x:L_partN(x)*dEFn(x),r0,rv,limit=60)[0]
    print(f"  r={rv}: int Q dA={LQ:+.6f}  int L dA(dilog)={LL:+.6f}")
print("\n=== convergenza q-serie del dilog in ordine di nome N (r=10) ===")
rv=10.0; ref=quad(lambda x:L_partN(x,8)*dEFn(x),r0,rv,limit=60)[0]
for N in [1,2,3,4,6]:
    val=quad(lambda x:L_partN(x,N)*dEFn(x),r0,rv,limit=60)[0]
    print(f"  N={N}: int L dA = {val:+.7f}  |diff da N=8| = {abs(val-ref):.1e}")
print("\n=> TK-t: dilog(infinito) = serie di nome KE convergente. + dilog ORIZZONTE (da rho_t, a z(r_pm)).")
