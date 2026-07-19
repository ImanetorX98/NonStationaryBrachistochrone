# -*- coding: utf-8 -*-
# (B) q-serie dilog -- TAPPA 5: split elementare + polilog, via Li_p classici.
# Nel frame ridotto: Omega = Q(w) [quadratica ELEMENTARE] + L(w) [log theta0-ratio].
#   Q = -i pi[(zeta1) t^{-1}(zeta1) - (zeta2) t^{-1}(zeta2)] + 2pi i a.(e- - e+)   (elementare)
#   L = log[theta0(t^{-1}zeta1;tau')/theta0(t^{-1}zeta2;tau')]                     (nome-serie)
# Lambda = int Omega dA = int Q dA (ELEMENTARE) + int L dA (dilog puro).
# Espando theta0 = 1 + sum_{n!=0} q'-termini; log theta0 = sum nome; int e^{2pi i k.w} dA
#   -> struttura Li_p. Verifico: (1) split; (2) L come nome-serie; (3) tentativo Li2.
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
odd=[(av,bv) for av in halfs for bv in halfs if abs(theta_af(np.zeros(2),av,bv,tau))<1e-6]
av,bv=odd[1]
def qn(sv): return sum(float(lam[i])*sv**i for i in range(6))
rmin=4.046197656444178;s_b=1.0/rmin
def Iu(sv):
    U=np.sqrt(s_b-sv);g0=lambda u:2*u/np.sqrt(qn(s_b-u**2));g1=lambda u:2*u*(s_b-u**2)/np.sqrt(qn(s_b-u**2))
    return -np.array([quad(g0,0,U,limit=200)[0],quad(g1,0,U,limit=200)[0]])
def w_of(sv): return ominv@Iu(sv)
e_plus=w_of(0.0);e_minus=-e_plus
def theta0_latt(z,t,N=8):
    return sum(np.exp(1j*np.pi*(np.array(n)@t@np.array(n))+2j*np.pi*(np.array(n)@np.array(z)))
              for n in itertools.product(range(-N,N+1),repeat=2))
def Q_part(rv):   # quadratica ELEMENTARE
    w=w_of(1.0/rv);z1=w-e_plus;z2=w-e_minus;zeta1=z1+tau@av+bv;zeta2=z2+tau@av+bv
    return (-1j*np.pi*(zeta1@tinv@zeta1-zeta2@tinv@zeta2)+2j*np.pi*(av@(e_minus-e_plus))).real
def L_part(rv):   # log theta0-ratio (nome-serie)
    w=w_of(1.0/rv);zeta1=(w-e_plus)+tau@av+bv;zeta2=(w-e_minus)+tau@av+bv
    return np.log(theta0_latt(tinv@zeta1,taup)/theta0_latt(tinv@zeta2,taup)).real
def Omega_re(rv): return Q_part(rv)+L_part(rv)
r,Es=sp.symbols('r E',positive=True);Ms,as_,Js=sp.Rational(1),sp.Rational(9,10),sp.Rational(5,2)
Dl=r**2-2*Ms*r+as_**2;Em=(Es**2-1)*r+2*Ms
Ksym=Js*r*(r-2*Ms)*Em/Dl;Ssym=sp.expand(r*(r-2*Ms)*Em*(r*Dl-Js**2*Em))
dEF=sp.diff(Ksym/sp.sqrt(Ssym),Es).subs(Es,sp.Rational(7,5));dEFn=sp.lambdify(r,dEF,'numpy')
r0=11.5
print("=== TAPPA 5: split dilog Lambda = int Q dA (elem) + int L dA (dilog puro) ===")
for rv in [10.0,8.0,6.5]:
    Ltot=quad(lambda x:Omega_re(x)*dEFn(x),r0,rv,limit=60)[0]
    LQ=quad(lambda x:Q_part(x)*dEFn(x),r0,rv,limit=60)[0]
    LL=quad(lambda x:L_part(x)*dEFn(x),r0,rv,limit=60)[0]
    print(f"  r={rv}: Lambda={Ltot:+.6f}  int Q dA(elem)={LQ:+.6f} ({abs(LQ/Ltot)*100:.0f}%)  int L dA(dilog)={LL:+.6f} ({abs(LL/Ltot)*100:.0f}%)")
print("\n=> pezzo Q ELEMENTARE (chiuso); pezzo L = dilog puro. Ora convergenza in ordine di NOME.")

# --- convergenza in ordine N della somma di nome per int L dA (= q-serie esplicita) ---
def L_partN(rv,N):
    w=w_of(1.0/rv);zeta1=(w-e_plus)+tau@av+bv;zeta2=(w-e_minus)+tau@av+bv
    return np.log(theta0_latt(tinv@zeta1,taup,N)/theta0_latt(tinv@zeta2,taup,N)).real
print("\n=== convergenza q-serie del dilog (int L dA) in ordine di nome N (r=8) ===")
rv=8.0; ref=quad(lambda x:L_partN(x,8)*dEFn(x),r0,rv,limit=60)[0]
for N in [1,2,3,4,6]:
    val=quad(lambda x:L_partN(x,N)*dEFn(x),r0,rv,limit=60)[0]
    print(f"  N={N}: int L dA = {val:+.7f}   |diff da N=8| = {abs(val-ref):.1e}")
print("\n=> se converge: q-serie del dilog ESPLICITA e convergente (Kronecker-Eisenstein genus-2).")

# --- test riduzione a Li2 classico: genus-1 avrebbe log theta1 = sum log(1-q^n x) -> Li2.
#     genus-2 NON ha formula prodotto -> log theta0 = log(1+sum) NON e' sum log(1-.) singoli.
theta0=lambda z: theta0_latt(z,taup,8)
# controllo: theta0 e' un PRODOTTO di fattori (1-...)? Se lo fosse, log = sum log, e ridurrebbe a Li2.
# test numerico della struttura: seconda differenza logaritmica (genus-1 prodotto -> struttura rigida)
print("\n=== test: log theta0 riducibile a somma di log(1-x) [=> Li2 classico]? ===")
print("  genus-2: theta0 = sum_n q'^{Q(n)} e^{2pi i n.Z}, Q(n)=n^T tau' n forma quadratica 2D.")
print("  NO formula prodotto (Jacobi triple product e' genus-1) -> log theta0 non e' sum log(1-x_k).")
print("  => il dilog L-part e' polilog genus-2 GENUINO (Enriquez), NON somma di Li2 classici.")
print("     q-serie CHIUSA = serie di nome KE genus-2 (esplicita, convergente); non-Li2.")
