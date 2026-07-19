# -*- coding: utf-8 -*-
# q-serie dilog genus-2 VAIDYA generico (a=0, param m). Template TK-tau: naming theta[delta]
# + Omega=log theta-ratio nome-serie + dilog Lambda split Q(elem)+L(dilog), convergenza nome.
from sage.all import QQ
from sage.schemes.riemann_surfaces.riemann_surface import RiemannSurface as SageRS
from abelfunctions import RiemannTheta
import numpy as np, sympy as sp
from scipy.integrate import quad
import itertools
mfloat,Efloat,Jfloat=1.0,1.4,2.5
Rs=PolynomialRing(QQ,['s','Y']);s,Y=Rs.gens()
lam=[QQ(24),QQ(-46),QQ(-248),QQ(-112),QQ(575),QQ(1250)]   # Vaidya odd model
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
# curva even + sorgente dm F (Vaidya)
r,m,E,J=sp.symbols('r m E J')
DE=(E**2-1)*r+2*m; S=sp.expand(r*(r-2*m)*DE*(r**2*(r-2*m)-J**2*DE)); K=J*DE
dmF=sp.diff(K/sp.sqrt(S),m).subs({m:1,E:sp.Rational(7,5),J:sp.Rational(5,2)})
dmFn=sp.lambdify(r,dmF,'numpy')
def Sn(x):
    E7=1.4;DEn=(E7**2-1)*x+2; return float(x*(x-2)*DEn*(x**2*(x-2)-6.25*DEn))
# turning e branch points (per w_of base)
Scoef=[float(c) for c in sp.Poly(S.subs({m:1,E:sp.Rational(7,5),J:sp.Rational(5,2)}),r).all_coeffs()]
rts=np.sort(np.real(np.roots(Scoef))); rmin=rts[rts>3][0] if any(rts>3) else rts[-1]
s_b=1.0/rmin
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
r0=12.0
print("\n=== VAIDYA dilog: split Lambda = int Q dA (elem) + int L dA (dilog) ===")
for rv in [10.,8.,6.5]:
    LQ=quad(lambda x:Q_part(x)*dmFn(x),r0,rv,limit=60)[0]
    LL=quad(lambda x:L_partN(x)*dmFn(x),r0,rv,limit=60)[0]
    print(f"  r={rv}: int Q dA(elem)={LQ:+.6f}  int L dA(dilog)={LL:+.6f}")
print("\n=== convergenza q-serie del dilog in ordine di nome N (r=8) ===")
rv=8.0; ref=quad(lambda x:L_partN(x,8)*dmFn(x),r0,rv,limit=60)[0]
for N in [1,2,3,4,6]:
    val=quad(lambda x:L_partN(x,N)*dmFn(x),r0,rv,limit=60)[0]
    print(f"  N={N}: int L dA = {val:+.7f}  |diff da N=8| = {abs(val-ref):.1e}")
print("\n=> VAIDYA: dilog = elem + serie di nome KE genus-2 convergente. Naming theta[delta] agli e_pm.")
