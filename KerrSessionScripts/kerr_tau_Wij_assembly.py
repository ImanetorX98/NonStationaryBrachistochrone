# -*- coding: utf-8 -*-
# MATTONE 2b chiusura W_ij (TK tau, genus-2): MONTAGGIO nella base canonica.
# psi = 1/2 Ehat (rho-rho~) = 1/2 Ehat int(A dη - η dA),
#   A = source primitive = sum_k c_k U_k ,  η = clock = sum_k b_k U_k,  b=(0,0,-2M,1,0).
# Ogni U_k = sum_alpha M_{k,alpha} V_alpha (+ boundary alg), V in {u1,u2,R1,R2,L}
#   (mattone 2a). => A=sum_a a_a V_a, η=sum_a h_a V_a, a_a=sum_k c_k M_ka, h_a=sum_k b_k M_ka.
# psi = 1/2 Ehat sum_{a<b} P_ab w_ab,  P_ab=a_a h_b - a_b h_a (SIMBOLICO),
#   w_ab = int(V_a dV_b - V_b dV_a) integrali iterati canonici.
# Coppie con L = dilog genus-2 (peso 2 irriducibile); resto = Kleinian peso-1.
import sympy as sp, numpy as np
from scipy.integrate import quad

r,s,E = sp.symbols('r s E', positive=True)
M,a,J = sp.Rational(1), sp.Rational(9,10), sp.Rational(5,2)
Dl = r**2-2*M*r+a**2; Emu=(E**2-1)*r+2*M
S = sp.expand(r*(r-2*M)*Emu*(r*Dl-J**2*Emu)); Sp=sp.diff(S,r)
K = J*r*(r-2*M)*Emu/Dl; F=K/sp.sqrt(S)

# --- c_k: riduzione dE F = d(Acal/sqrtS) + sum_k c_k r^k/sqrtS (2N=2S A' - A S' + 2 S Mpoly) ---
dEF=sp.diff(F,E); N=E*J*r**4*(r-2*M)**2*Emu
assert sp.simplify(dEF-N/S**sp.Rational(3,2))==0
Ac=[sp.Symbol(f'A{i}') for i in range(6)]; ck=[sp.Symbol(f'c{i}') for i in range(5)]
Acal=sum(Ac[i]*r**i for i in range(6)); Mp=sum(ck[i]*r**i for i in range(5))
eq=sp.expand(2*N-(2*S*sp.diff(Acal,r)-Acal*Sp+2*S*Mp))
sol=sp.solve(sp.Poly(eq,r).all_coeffs(),Ac+ck,dict=True)[0]
cvec=[sp.simplify(sol[ck[i]]) for i in range(5)]
print("c_k (source, razionali in E):"); [print(f"  c{i}=",cvec[i]) for i in range(5)]
bvec=[0,0,-2*M,1,0]   # clock

# --- M_{k,alpha}: U_k in base {u1,u2,R1,R2,L} (da mattone 2a) ---
# ricostruisco la riduzione odd-model qui (stesse formule)
q6=sp.expand(s**6*S.subs(r,1/s)); q6p=sp.diff(q6,s)
lam=[sp.Poly(q6,s).coeff_monomial(s**i) for i in range(6)]
Ndr1=(lam[3]*s+2*lam[4]*s**2+3*lam[5]*s**3)/4; Ndr2=lam[5]*s**2/4
Pbnd={0:sp.Integer(0),1:sp.Integer(0)}   # boundary di U_k: U_k=[P_k Y]+sum M V
def reduce_omega(k):
    Rk=-s**(1-k); ms=list(range(-(k-1),1))
    am=[sp.Symbol(f'a{i}') for i in range(len(ms))]
    njm1=sp.Symbol('nm1'); nj=[sp.Symbol(f'n{j}') for j in range(4)]
    Npoly=njm1/s+sum(nj[j]*s**j for j in range(4))
    exact=sum(am[i]*(2*ms[i]*s**(ms[i]-1)*q6+s**ms[i]*q6p)/2 for i in range(len(ms)))
    expr=sp.together(sp.expand(Rk-exact-Npoly)); poly=sp.Poly(sp.expand(sp.numer(expr)),s)
    ss=sp.solve(poly.all_coeffs(),am+[njm1]+nj,dict=True)[0]
    res3=sp.simplify(njm1.subs(ss)); Nred=sp.expand(sum(nj[j].subs(ss)*s**j for j in range(4)))
    Pbnd[k]=sp.simplify(sum(am[i].subs(ss)*s**ms[i] for i in range(len(ms))))
    c1,c2,g1,g2=sp.symbols('c1 c2 g1 g2')
    dec=sp.solve(sp.Poly(sp.expand(Nred-(c1*Ndr1+c2*Ndr2+g1+g2*s)),s).all_coeffs(),[c1,c2,g1,g2],dict=True)[0]
    return {'R1':sp.simplify(dec[c1]),'R2':sp.simplify(dec[c2]),'u1':sp.simplify(dec[g1]),
            'u2':sp.simplify(dec[g2]),'L':res3}
alphas=['u1','u2','R1','R2','L']
Mmat={0:{'u2':sp.Integer(-1)},1:{'u1':sp.Integer(-1)}}  # U0=-u2, U1=-u1
for k in [2,3,4]: Mmat[k]=reduce_omega(k)
Acal=Acal.subs(sol)   # 𝒜(r): parte algebrica del source (A=[Acal/sqrtS]+sum c_k U_k)
def Mka(k,al): return Mmat[k].get(al,sp.Integer(0))
# a_alpha, h_alpha
a_al={al:sp.simplify(sum(cvec[k]*Mka(k,al) for k in range(5))) for al in alphas}
h_al={al:sp.simplify(sum(bvec[k]*Mka(k,al) for k in range(5))) for al in alphas}
print("\nsource in base canonica a_alpha:"); [print(f"  {al}:",a_al[al]) for al in alphas]
print("clock in base canonica h_alpha:");   [print(f"  {al}:",h_al[al]) for al in alphas]
# P_ab
P={}
for i in range(5):
    for j in range(i+1,5):
        ai,aj=alphas[i],alphas[j]; Pab=sp.simplify(a_al[ai]*h_al[aj]-a_al[aj]*h_al[ai])
        if Pab!=0: P[(ai,aj)]=Pab
print("\nP_ab NONZERO (coeff simbolici degli iterati canonici w_ab):")
for key,val in P.items():
    tag="DILOG (peso2, ha L)" if 'L' in key else "peso1 Kleinian"
    print(f"  {key}: {val}     [{tag}]")
print("\nCoppia olomorfa (u1,u2):", P.get(('u1','u2'),0), " (deve 0: niente olo x olo, cf Q_01=0)")

# ===== VERIFICA end-to-end: psi_montaggio vs psi_diretto (E=7/5) =====
E0=sp.Rational(7,5); Ehat=1.4; Mf,af,Ef,Jf=1.0,0.9,1.4,2.5; r0=12.0; s0=1.0/r0
q6n=sp.lambdify(s,q6.subs(E,E0),'numpy')
Nnum={'u1':lambda t:1.0,'u2':lambda t:t,
      'R1':sp.lambdify(s,Ndr1.subs(E,E0),'numpy'),'R2':sp.lambdify(s,Ndr2.subs(E,E0),'numpy'),
      'L':lambda t:1.0/t}
def Vint(al,rv):  # V_alpha(r)=int_{s0}^{s} N_al/Y ds,  s=1/r
    sv=1.0/rv; return quad(lambda t: Nnum[al](t)/np.sqrt(q6n(t)), s0, sv, limit=200)[0]
def dVdr(al,rv):  # dV/dr = -s^2 N_al(s)/Y,  s=1/r
    sv=1.0/rv; return -sv**2*Nnum[al](sv)/np.sqrt(q6n(sv))
def w_ab(a_,b_,rv):  # int_{r0}^{r}(V_a dV_b - V_b dV_a)
    return quad(lambda x: Vint(a_,x)*dVdr(b_,x)-Vint(b_,x)*dVdr(a_,x), r0, rv, limit=120)[0]
# psi diretto (ground truth) = 1/2 Ehat (rho-rho~)
def Sn(x):
    Dl=x**2-2*Mf*x+af**2; Em=(Ef**2-1)*x+2*Mf; return x*(x-2*Mf)*Em*(x*Dl-Jf**2*Em)
def sqn(x): return np.sqrt(Sn(x))
dEFn=sp.lambdify(r,dEF.subs(E,E0),'numpy')
def A_of(x): return quad(dEFn,r0,x,limit=200)[0]
def Uk(x,k): return quad(lambda t:t**k/sqn(t),r0,x,limit=200)[0]
def B_of(x): return Uk(x,3)-2*Mf*Uk(x,2)
def psi_direct(rv):
    rho=quad(lambda x:A_of(x)*(x**3-2*Mf*x**2)/sqn(x),r0,rv,limit=100)[0]
    rhot=quad(lambda x:B_of(x)*dEFn(x),r0,rv,limit=100)[0]
    return 0.5*Ehat*(rho-rhot)
Pn={key:float(val.subs(E,E0)) for key,val in P.items()}
# --- parti ALGEBRICHE (boundary di U_k, e Acal/sqrtS del source) ---
Yof=lambda sv: np.sqrt(q6n(sv))
Pbn={k:sp.lambdify(s,Pbnd[k].subs(E,E0),'numpy') for k in range(5)}
def bd(k,rv):   # boundary di U_k = [P_k Y]_{s0}^{s}
    sv=1.0/rv; return Pbn[k](sv)*Yof(sv)-Pbn[k](s0)*Yof(s0)
Acaln=sp.lambdify(r,Acal.subs(E,E0),'numpy')
cval=[float(cvec[k].subs(E,E0)) for k in range(5)]; bval=[0,0,-2*Mf,1,0]
a_aln={al:float(a_al[al].subs(E,E0)) for al in alphas}; h_aln={al:float(h_al[al].subs(E,E0)) for al in alphas}
def A_alg(rv): return (Acaln(rv)/sqn(rv)-Acaln(r0)/sqn(r0))+sum(cval[k]*bd(k,rv) for k in range(5))
def eta_alg(rv): return sum(bval[k]*bd(k,rv) for k in range(5))
def A_ab(rv): return sum(a_aln[al]*Vint(al,rv) for al in alphas)
def eta_ab(rv): return sum(h_aln[al]*Vint(al,rv) for al in alphas)
print("\n=== VERIFICA decomposizione A,eta = algebrico + abeliano (coeff simbolici) ===")
for rv in [10.0,8.0,6.5]:
    Ad=A_alg(rv)+A_ab(rv); Ao=A_of(rv); ed=eta_alg(rv)+eta_ab(rv); eo=B_of(rv)
    print(f"  r={rv}: A_dec-A_dir={abs(Ad-Ao):.1e}  eta_dec-eta_dir={abs(ed-eo):.1e}")
# T_alg INDIPENDENTE: integrali con almeno una parte algebrica
etatot_p=lambda x:(x**3-2*Mf*x**2)/sqn(x)   # dη/dr totale
def dA_ab(rv): return sum(a_aln[al]*dVdr(al,rv) for al in alphas)
def deta_ab(rv): return sum(h_aln[al]*dVdr(al,rv) for al in alphas)
def dA_alg(rv): return dEFn(rv)-dA_ab(rv)
def deta_alg(rv): return etatot_p(rv)-deta_ab(rv)
def T_alg_formula(rv):
    integrand=lambda x:(A_alg(x)*deta_alg(x)-eta_alg(x)*dA_alg(x)
                        +A_alg(x)*deta_ab(x)-eta_alg(x)*dA_ab(x)
                        +A_ab(x)*deta_alg(x)-eta_ab(x)*dA_alg(x))
    return 0.5*Ehat*quad(integrand,r0,rv,limit=120)[0]
print("\n=== psi = 1/2 Ehat[ Sum P_ab w_ab (trascend.) + T_alg (algebrico) ] end-to-end ===")
for rv in [10.0,8.0,6.5]:
    trans=0.5*Ehat*sum(Pn[key]*w_ab(key[0],key[1],rv) for key in Pn)
    dil=0.5*Ehat*sum(Pn[key]*w_ab(key[0],key[1],rv) for key in Pn if 'L' in key)
    pd=psi_direct(rv); Talg_f=T_alg_formula(rv); tot=trans+Talg_f
    print(f"  r={rv}: psi_dir={pd:+.7f}  trascend+T_alg={tot:+.7f}  diff={abs(pd-tot):.1e}  (dilog={dil:+.6f})")
print("\n=> psi = 1/2 Ehat [ P(u1,R1)w + P(u2,R1)w + P(R1,R2)w  (KLEINIAN peso1) ]")
print("        + 1/2 Ehat [ P(.,L) w(.,L)  (UN dilog genus-2, sorgente = 3a specie clock) ]")
print("        + termini ALGEBRICI (boundary elementari, gia' chiusi).")
print("   A,eta decomposti in alg+abeliano con coeff SIMBOLICI (razionali E), verificato.")
