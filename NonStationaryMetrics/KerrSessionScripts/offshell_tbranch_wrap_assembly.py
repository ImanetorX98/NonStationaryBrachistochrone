# -*- coding: utf-8 -*-
# GENERIC off-shell wrap (t-branch, main11 "open" piece) -- SECOND-KIND block closed in the
# canonical genus-2 basis {u1,u2,R1,R2,L}, REUSING the verified tau-branch machinery
# (kerr_tau_Wij_assembly.py). Only the curve and the two input letters change:
#   curve  S = r*Emu*Q2  (t-branch genuine quartic Q2),  Emu=(E^2-1)r+2M
#   letter Phi_p = sum_j p_j U_j  (kernel A_poly=p1 r+p0),  U_j=int r^j/sqrt(S)
#   letter Phi_a = sum_k a_k U_k  (dilation radial-action polynomial part)
#   d phi_extra^(2nd) = -int (A_poly/sqrt S) Phi_a dr
#                     = 1/2 sum_{al<be} P_ab w_ab  + 1/2 T_alg  - 1/2 [Phi_p Phi_a]  (products)
#   P_ab = a_al h_be - a_be h_al ,  a_al=sum_j p_j M_ja , h_al=sum_k a_k M_ka  (SYMBOLIC in E).
# Ground truth: direct2 (single nested integral) from adiabatic_offshell_coefficients._verify.
import sympy as sp, numpy as np
from scipy.integrate import quad

r,s,E = sp.symbols('r s E', positive=True)
Mf,af,Jf = 1.0, 0.9, 6.0
M,a,J = sp.Rational(1), sp.Rational(9,10), sp.Integer(6)
E0=sp.Rational(7,5)                                   # E=1.4
Emu=(E**2-1)*r+2*M; Dl=r**2-2*M*r+a**2
Q2=(2*E**2*J**2*M*r-E**2*J**2*r**2-4*E**2*J*M*a*r+2*E**2*M*a**2*r+E**2*a**2*r**2
    +E**2*r**4+4*J**2*M**2-4*J**2*M*r+J**2*r**2-8*J*M**2*a+4*J*M*a*r+4*M**2*a**2)
S=sp.expand(r*Emu*Q2)

# ---- off-shell coefficient vectors (symbolic in E) ----
Qp,N3=sp.div(sp.Poly(sp.expand(E**2*J*r**4*Emu),r), sp.Poly(sp.expand(Q2),r))
p_j=[sp.simplify(c) for c in Qp.all_coeffs()[::-1]]            # p_0,p_1
deng=sp.expand(Dl*Emu); Qd,R=sp.div(sp.Poly(sp.expand(S),r),sp.Poly(deng,r))
a_k=[sp.simplify(c) for c in Qd.all_coeffs()[::-1]]           # a_0..a_3
print("kernel p_j =",p_j)
print("dilation a_k =",a_k)

# ---- odd model + U_k reduction (identical machinery) ----
q6=sp.expand(s**6*S.subs(r,1/s)); q6p=sp.diff(q6,s)
lam=[sp.Poly(q6,s).coeff_monomial(s**i) for i in range(7)]
Ndr1=(lam[3]*s+2*lam[4]*s**2+3*lam[5]*s**3)/4; Ndr2=lam[5]*s**2/4
Pbnd={0:sp.Integer(0),1:sp.Integer(0)}
def reduce_omega(k):
    Rk=-s**(1-k); ms=list(range(-(k-1),1))
    am=[sp.Symbol(f'am{i}') for i in range(len(ms))]
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
Mmat={0:{'u2':sp.Integer(-1)},1:{'u1':sp.Integer(-1)}}
for k in [2,3]: Mmat[k]=reduce_omega(k)
def Mka(k,al): return Mmat[k].get(al,sp.Integer(0))
# source Phi_p (p_j, j=0,1) and clock Phi_a (a_k, k=0..3) in canonical basis
a_al={al:sp.simplify(sum(p_j[j]*Mka(j,al) for j in range(2))) for al in alphas}   # a_alpha (from kernel)
h_al={al:sp.simplify(sum(a_k[k]*Mka(k,al) for k in range(4))) for al in alphas}   # h_alpha (from dilation)
P={}
for i in range(5):
    for j in range(i+1,5):
        ai,aj=alphas[i],alphas[j]; Pab=sp.simplify(a_al[ai]*h_al[aj]-a_al[aj]*h_al[ai])
        if Pab!=0: P[(ai,aj)]=Pab
print("\nP_ab NONZERO (symbolic coeff of canonical iterated integrals w_ab):")
for key,val in P.items():
    tag="DILOG (weight-2, has L)" if 'L' in key else "weight-1 Kleinian"
    print(f"  {key}: [{tag}]")

# ================= numeric end-to-end verification =================
q6n=sp.lambdify(s,q6.subs(E,E0),'numpy')
Nnum={'u1':lambda t:1.0,'u2':lambda t:t,
      'R1':sp.lambdify(s,Ndr1.subs(E,E0),'numpy'),'R2':sp.lambdify(s,Ndr2.subs(E,E0),'numpy'),
      'L':lambda t:1.0/t}
r0,rf=11.0,6.0; s0=1.0/r0
def Sn(x):
    Em=(1.96-1)*x+2*Mf; Dln=x**2-2*Mf*x+af**2
    Q2n=(2*1.96*Jf**2*Mf*x-1.96*Jf**2*x**2-4*1.96*Jf*Mf*af*x+2*1.96*Mf*af**2*x+1.96*af**2*x**2
         +1.96*x**4+4*Jf**2*Mf**2-4*Jf**2*Mf*x+Jf**2*x**2-8*Jf*Mf**2*af+4*Jf*Mf*af*x+4*Mf**2*af**2)
    return x*Em*Q2n
def sqn(x): return np.sqrt(max(Sn(x),0.0))
def Vint(al,rv):
    sv=1.0/rv; return quad(lambda t: Nnum[al](t)/np.sqrt(q6n(t)), s0, sv, limit=200)[0]
def dVdr(al,rv):
    sv=1.0/rv; return -sv**2*Nnum[al](sv)/np.sqrt(q6n(sv))
def w_ab(a_,b_,rv):
    return quad(lambda x: Vint(a_,x)*dVdr(b_,x)-Vint(b_,x)*dVdr(a_,x), r0, rv, limit=120)[0]
def Uk(x,k): return quad(lambda t:t**k/sqn(t),r0,x,limit=200)[0]
pjn=[float(x.subs(E,E0)) for x in p_j]; akn=[float(x.subs(E,E0)) for x in a_k]
# ground truth: direct2 = -int (A_poly/sqrtS) * (sum a_k U_k) dr
def Phi_a_of(x): return sum(akn[k]*Uk(x,k) for k in range(4))
def direct2(rv): return -quad(lambda x:((pjn[0]+pjn[1]*x)/sqn(x))*Phi_a_of(x), r0, rv, limit=150)[0]
# products term: -1/2 Phi_p(rf) Phi_a(rf), Phi_p=sum p_j U_j
def Phi_p_of(x): return sum(pjn[j]*Uk(x,j) for j in range(2))
# boundary of U_k = [P_k Y]
Yof=lambda sv: np.sqrt(q6n(sv))
Pbn={k:sp.lambdify(s,Pbnd[k].subs(E,E0),'numpy') for k in range(4)}
def bd(k,rv): sv=1.0/rv; return Pbn[k](sv)*Yof(sv)-Pbn[k](s0)*Yof(s0)
a_aln={al:float(a_al[al].subs(E,E0)) for al in alphas}; h_aln={al:float(h_al[al].subs(E,E0)) for al in alphas}
Pn={key:float(val.subs(E,E0)) for key,val in P.items()}
# alg/ab split of each letter (a_al->Phi_p, h_al->Phi_a)
def Phip_ab(rv): return sum(a_aln[al]*Vint(al,rv) for al in alphas)
def Phia_ab(rv): return sum(h_aln[al]*Vint(al,rv) for al in alphas)
def Phip_alg(rv): return sum(pjn[j]*bd(j,rv) for j in range(2))
def Phia_alg(rv): return sum(akn[k]*bd(k,rv) for k in range(4))
def dPhip_ab(rv): return sum(a_aln[al]*dVdr(al,rv) for al in alphas)
def dPhia_ab(rv): return sum(h_aln[al]*dVdr(al,rv) for al in alphas)
dPhip_tot=lambda x:(pjn[0]+pjn[1]*x)/sqn(x); dPhia_tot=lambda x:sum(akn[k]*x**k for k in range(4))/sqn(x)
def dPhip_alg(rv): return dPhip_tot(rv)-dPhip_ab(rv)
def dPhia_alg(rv): return dPhia_tot(rv)-dPhia_ab(rv)
def T_alg(rv):  # boundary-involving part of 1/2 int(Phi_p dPhi_a - Phi_a dPhi_p)
    integrand=lambda x:(Phip_alg(x)*dPhia_alg(x)-Phia_alg(x)*dPhip_alg(x)
                        +Phip_alg(x)*dPhia_ab(x)-Phia_alg(x)*dPhip_ab(x)
                        +Phip_ab(x)*dPhia_alg(x)-Phia_ab(x)*dPhip_alg(x))
    return 0.5*quad(integrand,r0,rv,limit=120)[0]
print("\n=== decomposition Phi_p,Phi_a = alg + abelian (symbolic coeff) ===")
for rv in [10.0,8.0,6.5]:
    print(f"  r={rv}: Phi_p diff={abs(Phip_alg(rv)+Phip_ab(rv)-Phi_p_of(rv)):.1e}  "
          f"Phi_a diff={abs(Phia_alg(rv)+Phia_ab(rv)-Phi_a_of(rv)):.1e}")
print("\n=== d phi_extra^(2nd): assembly (Kleinian w_ab + dilog + T_alg + products) vs direct2 ===")
for rv in [10.0,8.0,6.5]:
    trans=0.5*sum(Pn[key]*w_ab(key[0],key[1],rv) for key in Pn)
    dil=0.5*sum(Pn[key]*w_ab(key[0],key[1],rv) for key in Pn if 'L' in key)
    prod=-0.5*Phi_p_of(rv)*Phi_a_of(rv)
    asm=trans+T_alg(rv)+prod; dd=direct2(rv)
    print(f"  r={rv}: direct2={dd:+.8f}  assembly={asm:+.8f}  diff={abs(dd-asm):.1e}  (dilog L-part={dil:+.6f})")
print("\n=> t-branch generic off-shell SECOND-KIND block = 1/2 sum P_ab w_ab (Kleinian weight-1)")
print("   + 1/2 P(.,L) w(.,L) (ONE genus-2 dilog) + T_alg + products, ALL symbolic coeff in E.")
