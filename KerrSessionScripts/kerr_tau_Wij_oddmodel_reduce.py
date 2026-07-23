# -*- coding: utf-8 -*-
# MATTONE 2a chiusura W_ij (TK tau, genus-2): riduzione dei differenziali abeliani
# nel MODELLO DISPARI (quintica Y^2=q(s), s=1/r, 1 punto all'infinito s=inf=r=0).
# omega_k = r^k dr/sqrt(S) = -s^{1-k}/Y ds  (Y=sqrt(q6), q6=s^6 S(1/s) quintica).
# Riduco omega_3,omega_4 (2a specie, poli a s=0=r=inf) alla base:
#   {du_1=ds/Y, du_2=s ds/Y}  (olomorfe)  +  numeratore poly N(s) deg<=3 (2a specie)
# + forma esatta  d(s^m Y)  (m<0 per cancellare i poli a s=0).
# => U_k = [a_m s^m Y]_{r0}^{r} (boundary ALGEBRICO) + int N_k(s)/Y ds.
# int N/Y ds = combinazione di u_i (olomorfe) e zeta_i Kleiniani (2a specie BEL).
# COEFFICIENTI SIMBOLICI (razionali in E). Verifica numerica dell'identita' differenz.
import sympy as sp, numpy as np
from scipy.integrate import quad

r,s,E = sp.symbols('r s E', positive=True)
M,a,J = sp.Rational(1), sp.Rational(9,10), sp.Rational(5,2)
Dl = r**2-2*M*r+a**2; Emu=(E**2-1)*r+2*M
S = sp.expand(r*(r-2*M)*Emu*(r*Dl-J**2*Emu))
# q6(s) = s^6 S(1/s)  (quintica: S ha fattore r -> grado scende a 5)
q6 = sp.expand(s**6 * S.subs(r, 1/s))
print("q6(s) =", q6, "   deg =", sp.degree(sp.Poly(q6,s)))
q6p = sp.diff(q6,s)

# BEL 2a-specie (numeratori, misura ds/Y): dr1=(l3 s+2 l4 s^2+3 l5 s^3)/4, dr2=l5 s^2/4
lam=[sp.Poly(q6,s).coeff_monomial(s**i) for i in range(6)]
print("lambda_i =", lam)
Ndr1 = (lam[3]*s + 2*lam[4]*s**2 + 3*lam[5]*s**3)/4    # numeratore di dr1 (su /Y)
Ndr2 = (lam[5]*s**2)/4                                  # numeratore di dr2

def reduce_omega(k):
    # omega_k numeratore (su ds/Y) = -s^{1-k}  (k>=2 -> polo a s=0)
    Rk = -s**(1-k)
    # forma ridotta N(s) INCLUDE s^{-1} (letter 3a specie -> log sigma-ratio):
    #   N = n_{-1}/s + n_0 + n_1 s + n_2 s^2 + n_3 s^3
    # esatte d(s^m Y)=(2m s^{m-1} q6 + s^m q6')/2 * ds/Y cancellano SOLO poli ordine>=2.
    #   Rk = sum_m a_m (2m s^{m-1} q6 + s^m q6')/2 + N(s)
    ms=list(range(-(k-1),1))            # m<=0
    am=[sp.Symbol(f'a{i}') for i in range(len(ms))]
    njm1=sp.Symbol('nm1'); nj=[sp.Symbol(f'n{j}') for j in range(4)]
    Npoly=njm1/s + sum(nj[j]*s**j for j in range(4))
    exact=sum(am[i]*(2*ms[i]*s**(ms[i]-1)*q6 + s**ms[i]*q6p)/2 for i in range(len(ms)))
    expr=sp.together(sp.expand(Rk - exact - Npoly))
    num=sp.numer(expr)
    poly=sp.Poly(sp.expand(num),s)
    sol=sp.solve(poly.all_coeffs(), am+[njm1]+nj, dict=True)
    assert sol, f"no reduction k={k}"
    sol=sol[0]
    Pbound=sum((am[i].subs(sol))*s**ms[i] for i in range(len(ms)))
    res3=sp.simplify(njm1.subs(sol))                       # coeff 3a specie (log sigma)
    Nred=sp.expand(sum(nj[j].subs(sol)*s**j for j in range(4)))  # parte 2a specie+olo (poly)
    return sp.simplify(Pbound), Nred, res3

print("\n=== Riduzione omega_k = d(P_k*Y) + (n_-1/s + N_k(s))*ds/Y (modello dispari) ===")
red={}
for k in [2,3,4]:
    Pb,Nred,res3=reduce_omega(k)
    # decomponi la parte poly N_k in {1,s}(olomorfe) + {Ndr1,Ndr2}(2a specie BEL)
    c1,c2=sp.symbols('c1 c2'); g1,g2=sp.symbols('g1 g2')
    dec=sp.expand(Nred - (c1*Ndr1+c2*Ndr2+g1+g2*s))
    dsol=sp.solve(sp.Poly(dec,s).all_coeffs(),[c1,c2,g1,g2],dict=True)[0]
    red[k]=(Pb,Nred,res3,dsol)
    print(f"k={k}: res_3a(n_-1)={res3} ; dr1:{sp.simplify(dsol[c1])} dr2:{sp.simplify(dsol[c2])} du1:{sp.simplify(dsol[g1])} du2:{sp.simplify(dsol[g2])} ; P_bnd={Pb}")

# --- VERIFICA numerica dell'identita' DIFFERENZIALE (E=7/5): omega_k = d(P_k Y)+N_k/Y ---
E0=sp.Rational(7,5)
q6n=sp.lambdify(s,q6.subs(E,E0),'numpy'); q6pn=sp.lambdify(s,q6p.subs(E,E0),'numpy')
print("\n=== VERIFICA numerica identita' differenziale (deve ~0) ===")
for k in [2,3,4]:
    Pb,Nred,res3,dsol=red[k]
    Pbn=sp.lambdify(s,Pb.subs(E,E0),'numpy'); Pbpn=sp.lambdify(s,sp.diff(Pb,s).subs(E,E0),'numpy')
    Nn=sp.lambdify(s,(Nred+res3/s).subs(E,E0),'numpy')   # parte totale su ds/Y (incl. 3a specie)
    def omega_num(sv):   # -s^{1-k}/Y
        Y=np.sqrt(complex(q6n(sv))); return -sv**(1-k)/Y
    def rhs_num(sv):     # d(P Y)=P'Y+P q6'/(2Y); poi +N/Y
        Y=np.sqrt(complex(q6n(sv)))
        dPY=Pbpn(sv)*Y + Pbn(sv)*q6pn(sv)/(2*Y)
        return dPY + Nn(sv)/Y
    errs=[abs(omega_num(sv)-rhs_num(sv)) for sv in [0.05,0.08,0.12,0.2]]
    print(f"  k={k}: max|omega - (d(PY)+N/Y)| = {max(errs):.2e}")
print("\n=> se ~0: U_k = [P_k Y]_{r0}^r + c1_k R1 + c2_k R2 + g1_k u1 + g2_k u2 + n_-1 L,")
print("   R_i=int dr_i (2a specie), u_i=int du_i (1a), L=int ds/(sY) (3a specie).")

# ===== VERIFICA a livello INTEGRALE: U_k(r) diretto vs boundary + comb. abeliani =====
# variabile fisica r, s=1/r; base r0=12 (s0=1/12). Integrali abeliani in s.
print("\n=== VERIFICA INTEGRALE: U_k(r) = boundary + coeff.(R1,R2,u1,u2,L)  (deve ~1e-8) ===")
Mf,af,Ef,Jf=1.0,0.9,1.4,2.5; r0=12.0; s0=1.0/r0
def Sn(x):
    Dl=x**2-2*Mf*x+af**2; Em=(Ef**2-1)*x+2*Mf
    return x*(x-2*Mf)*Em*(x*Dl-Jf**2*Em)
def sq(x): return np.sqrt(Sn(x))
def Yof(sv): return np.sqrt(complex(q6n(sv))).real if q6n(sv)>0 else np.sqrt(complex(q6n(sv)))
# integrali abeliani in s (da s0 a s=1/r), misura ds/Y con Y=+sqrt(q6)
Ndr1n=sp.lambdify(s,Ndr1.subs(E,E0),'numpy'); Ndr2n=sp.lambdify(s,Ndr2.subs(E,E0),'numpy')
def R1(sv): return quad(lambda t: Ndr1n(t)/np.sqrt(q6n(t)), s0, sv, limit=200)[0]
def R2(sv): return quad(lambda t: Ndr2n(t)/np.sqrt(q6n(t)), s0, sv, limit=200)[0]
def u1(sv): return quad(lambda t: 1.0/np.sqrt(q6n(t)), s0, sv, limit=200)[0]
def u2(sv): return quad(lambda t: t/np.sqrt(q6n(t)), s0, sv, limit=200)[0]
def Lint(sv): return quad(lambda t: 1.0/(t*np.sqrt(q6n(t))), s0, sv, limit=200)[0]
def Uk_direct(k,rv): return quad(lambda x: x**k/sq(x), r0, rv, limit=200)[0]
for k in [2,3,4]:
    Pb,Nred,res3,dsol=red[k]
    Pbn=sp.lambdify(s,Pb.subs(E,E0),'numpy')
    c1v=float(dsol[sp.Symbol('c1')].subs(E,E0)); c2v=float(dsol[sp.Symbol('c2')].subs(E,E0))
    g1v=float(dsol[sp.Symbol('g1')].subs(E,E0)); g2v=float(dsol[sp.Symbol('g2')].subs(E,E0))
    n1v=float(res3.subs(E,E0))
    for rv in [10.0,8.0,6.0]:
        sv=1.0/rv
        # boundary [P Y]_{s0}^{sv} ; nota omega=-s^{1-k}/Y ds e U_k=int_{r0}^r r^k dr/sqrtS
        # sotto r=1/s, int_{r0}^{r} ... = int_{s0}^{sv} omega. boundary = P Y|_{sv} - P Y|_{s0}
        bnd = Pbn(sv)*np.sqrt(q6n(sv)) - Pbn(s0)*np.sqrt(q6n(s0))
        rec = bnd + c1v*R1(sv)+c2v*R2(sv)+g1v*u1(sv)+g2v*u2(sv)+n1v*Lint(sv)
        direct=Uk_direct(k,rv)
        if rv==10.0: print(f"  k={k}:")
        print(f"     r={rv}: U_k(dir)={direct:+.8e}  ricostr={rec:+.8e}  diff={abs(direct-rec):.1e}")
print("\n=> U_k CHIUSI su integrali abeliani canonici (1a/2a/3a specie), coeff SIMBOLICI.")
print("   Prossimo (Sage): R_i->zeta_i(u), L->log sigma(u-e+)/sigma(u-e-) [naming Kleiniano].")
