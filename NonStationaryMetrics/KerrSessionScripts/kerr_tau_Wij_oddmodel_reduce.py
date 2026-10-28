# -*- coding: utf-8 -*-
# BUILDING BLOCK 2a of the W_ij closure (TK tau, genus-2): reduction of the abelian differentials
# in the ODD MODEL (quintic Y^2=q(s), s=1/r, one point at infinity s=inf=r=0).
# omega_k = r^k dr/sqrt(S) = -s^{1-k}/Y ds  (Y=sqrt(q6), q6=s^6 S(1/s) quintic).
# Reduce omega_3,omega_4 (second kind, poles at s=0=r=inf) to the basis:
#   {du_1=ds/Y, du_2=s ds/Y}  (holomorphic)  +  polynomial numerator N(s) deg<=3 (second kind)
# + an exact form  d(s^m Y)  (m<0 to cancel the poles at s=0).
# => U_k = [a_m s^m Y]_{r0}^{r} (ALGEBRAIC boundary term) + int N_k(s)/Y ds.
# int N/Y ds = a combination of u_i (holomorphic) and Kleinian zeta_i (BEL second kind).
# SYMBOLIC COEFFICIENTS (rational in E). Numerical check of the differential identity.
import sympy as sp, numpy as np
from scipy.integrate import quad

r,s,E = sp.symbols('r s E', positive=True)
M,a,J = sp.Rational(1), sp.Rational(9,10), sp.Rational(5,2)
Dl = r**2-2*M*r+a**2; Emu=(E**2-1)*r+2*M
S = sp.expand(r*(r-2*M)*Emu*(r*Dl-J**2*Emu))
# q6(s) = s^6 S(1/s)  (quintic: S has a factor r -> the degree drops to 5)
q6 = sp.expand(s**6 * S.subs(r, 1/s))
print("q6(s) =", q6, "   deg =", sp.degree(sp.Poly(q6,s)))
q6p = sp.diff(q6,s)

# BEL second kind (numerators, measure ds/Y): dr1=(l3 s+2 l4 s^2+3 l5 s^3)/4, dr2=l5 s^2/4
lam=[sp.Poly(q6,s).coeff_monomial(s**i) for i in range(6)]
print("lambda_i =", lam)
Ndr1 = (lam[3]*s + 2*lam[4]*s**2 + 3*lam[5]*s**3)/4    # numerator of dr1 (over /Y)
Ndr2 = (lam[5]*s**2)/4                                  # numerator of dr2

def reduce_omega(k):
    # numerator of omega_k (over ds/Y) = -s^{1-k}  (k>=2 -> pole at s=0)
    Rk = -s**(1-k)
    # the reduced form N(s) INCLUDES s^{-1} (third-kind letter -> log sigma-ratio):
    #   N = n_{-1}/s + n_0 + n_1 s + n_2 s^2 + n_3 s^3
    # the exact forms d(s^m Y)=(2m s^{m-1} q6 + s^m q6')/2 * ds/Y cancel ONLY poles of order>=2.
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
    res3=sp.simplify(njm1.subs(sol))                       # third-kind coeff (log sigma)
    Nred=sp.expand(sum(nj[j].subs(sol)*s**j for j in range(4)))  # second-kind + holomorphic part (poly)
    return sp.simplify(Pbound), Nred, res3

print("\n=== Reduction omega_k = d(P_k*Y) + (n_-1/s + N_k(s))*ds/Y (odd model) ===")
red={}
for k in [2,3,4]:
    Pb,Nred,res3=reduce_omega(k)
    # split the polynomial part N_k into {1,s} (holomorphic) + {Ndr1,Ndr2} (BEL second kind)
    c1,c2=sp.symbols('c1 c2'); g1,g2=sp.symbols('g1 g2')
    dec=sp.expand(Nred - (c1*Ndr1+c2*Ndr2+g1+g2*s))
    dsol=sp.solve(sp.Poly(dec,s).all_coeffs(),[c1,c2,g1,g2],dict=True)[0]
    red[k]=(Pb,Nred,res3,dsol)
    print(f"k={k}: res_3rd(n_-1)={res3} ; dr1:{sp.simplify(dsol[c1])} dr2:{sp.simplify(dsol[c2])} du1:{sp.simplify(dsol[g1])} du2:{sp.simplify(dsol[g2])} ; P_bnd={Pb}")

# --- NUMERICAL CHECK of the DIFFERENTIAL identity (E=7/5): omega_k = d(P_k Y)+N_k/Y ---
E0=sp.Rational(7,5)
q6n=sp.lambdify(s,q6.subs(E,E0),'numpy'); q6pn=sp.lambdify(s,q6p.subs(E,E0),'numpy')
print("\n=== NUMERICAL CHECK of the differential identity (must be ~0) ===")
for k in [2,3,4]:
    Pb,Nred,res3,dsol=red[k]
    Pbn=sp.lambdify(s,Pb.subs(E,E0),'numpy'); Pbpn=sp.lambdify(s,sp.diff(Pb,s).subs(E,E0),'numpy')
    Nn=sp.lambdify(s,(Nred+res3/s).subs(E,E0),'numpy')   # total part over ds/Y (incl. third kind)
    def omega_num(sv):   # -s^{1-k}/Y
        Y=np.sqrt(complex(q6n(sv))); return -sv**(1-k)/Y
    def rhs_num(sv):     # d(P Y)=P'Y+P q6'/(2Y); then +N/Y
        Y=np.sqrt(complex(q6n(sv)))
        dPY=Pbpn(sv)*Y + Pbn(sv)*q6pn(sv)/(2*Y)
        return dPY + Nn(sv)/Y
    errs=[abs(omega_num(sv)-rhs_num(sv)) for sv in [0.05,0.08,0.12,0.2]]
    print(f"  k={k}: max|omega - (d(PY)+N/Y)| = {max(errs):.2e}")
print("\n=> if ~0: U_k = [P_k Y]_{r0}^r + c1_k R1 + c2_k R2 + g1_k u1 + g2_k u2 + n_-1 L,")
print("   R_i=int dr_i (second kind), u_i=int du_i (first), L=int ds/(sY) (third kind).")

# ===== CHECK at the INTEGRAL level: U_k(r) direct vs boundary + abelian combination =====
# physical variable r, s=1/r; base point r0=12 (s0=1/12). Abelian integrals in s.
print("\n=== INTEGRAL CHECK: U_k(r) = boundary + coeff.(R1,R2,u1,u2,L)  (must be ~1e-8) ===")
Mf,af,Ef,Jf=1.0,0.9,1.4,2.5; r0=12.0; s0=1.0/r0
def Sn(x):
    Dl=x**2-2*Mf*x+af**2; Em=(Ef**2-1)*x+2*Mf
    return x*(x-2*Mf)*Em*(x*Dl-Jf**2*Em)
def sq(x): return np.sqrt(Sn(x))
def Yof(sv): return np.sqrt(complex(q6n(sv))).real if q6n(sv)>0 else np.sqrt(complex(q6n(sv)))
# abelian integrals in s (from s0 to s=1/r), measure ds/Y with Y=+sqrt(q6)
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
        # boundary [P Y]_{s0}^{sv} ; note omega=-s^{1-k}/Y ds and U_k=int_{r0}^r r^k dr/sqrtS
        # under r=1/s, int_{r0}^{r} ... = int_{s0}^{sv} omega. boundary = P Y|_{sv} - P Y|_{s0}
        bnd = Pbn(sv)*np.sqrt(q6n(sv)) - Pbn(s0)*np.sqrt(q6n(s0))
        rec = bnd + c1v*R1(sv)+c2v*R2(sv)+g1v*u1(sv)+g2v*u2(sv)+n1v*Lint(sv)
        direct=Uk_direct(k,rv)
        if rv==10.0: print(f"  k={k}:")
        print(f"     r={rv}: U_k(dir)={direct:+.8e}  reconstr={rec:+.8e}  diff={abs(direct-rec):.1e}")
print("\n=> U_k CLOSED on canonical abelian integrals (first/second/third kind), SYMBOLIC coeffs.")
print("   Next (Sage): R_i->zeta_i(u), L->log sigma(u-e+)/sigma(u-e-) [Kleinian naming].")
