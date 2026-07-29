# -*- coding: utf-8 -*-
# Vaidya v-branch off-shell: EXPLICIT closed-form assembly (block1 A+B+C + block2 dilog+M2m), verified.
# kernel=A_kernel_v/sqrt(S_v), A_kernel_v=E^2 J DE r^4/Q2v = A_poly+N3/Q2v (Hermite).
# inner 2nd-kind: Sigma_2nd=int P_inner_v/sqrt(S_v);  P_inner_v = poly + a2/(r-2m)^2 + a1/(r-2m) + b/DE.
#   double pole 1/(r-2m)^2 reduced: int dr/((r-2m)^2 sqrtS) via d[sqrtS/(r-2m)].
# Sigma_elem = 4m^2/(r-2m) - 2m log(r-2m)  (closed).
# wrap = -int(A_kernel_v/sqrtS)(Sigma_2nd+Sigma_elem) = block1 + block2. Verify == direct.
import sympy as sp, numpy as np
from scipy.integrate import quad
from scipy.optimize import brentq
r,pr,m,E,J=sp.symbols('r p_r m E J',positive=True)
f=1-2*m/r; w=E**2-f; DE=(E**2-1)*r+2*m; K1=sp.simplify(f-E**2); K2=sp.simplify(w*E**2-K1**2)
Q2v=(2*E**2*J**2*m*r-E**2*J**2*r**2+E**2*r**4+4*J**2*m**2-4*J**2*m*r+J**2*r**2)
S_v=sp.expand(r*DE*Q2v); D_v=sp.simplify(K1**2-K2*(w*J**2/r**2-1)); P=sp.simplify(K2+K1**2)
SECp=sp.simplify(2*w*K1*D_v+2*w*K1*P+P**2+K1**2*D_v)
A_kernel_v=sp.simplify(E**2*J*DE*r**4/Q2v)
P_inner_v=sp.cancel(-m*SECp*r**2/(K2**2*w))
# --- partial fractions of P_inner_v over (r-2m)^2 * DE ---
den=sp.expand((r-2*m)**2*DE)
polyq,rem=sp.div(sp.Poly(sp.expand(P_inner_v*den),r),sp.Poly(den,r))   # P_inner_v = polyq + rem/den
pf=sp.apart(rem.as_expr()/den, r, full=False)
print("P_inner_v poly part coeffs:",[sp.simplify(c) for c in polyq.all_coeffs()[::-1]])
print("P_inner_v pole part (apart):",pf)
a2=sp.simplify((rem.as_expr()/DE).subs(r,2*m))                     # coeff of 1/(r-2m)^2
a1=sp.simplify(sp.diff(rem.as_expr()/DE,r).subs(r,2*m))           # coeff of 1/(r-2m)
bDE=sp.simplify((rem.as_expr()/(r-2*m)**2).subs(r,sp.solve(DE,r)[0]))  # coeff of 1/DE
print(f"a2(1/(r-2m)^2)={a2}  a1(1/(r-2m))={a1}  b(1/DE)={bDE}")

# ===== numeric verification of the full closed form =====
sub={m:1.0,E:sp.Rational(7,5),J:sp.Rational(5,2)}; r0=12.0
Svn=sp.lambdify(r,S_v.subs(sub),'numpy'); sq=lambda x:np.sqrt(max(Svn(x),0.0))
Akv=sp.lambdify(r,A_kernel_v.subs(sub),'numpy'); Pinv=sp.lambdify(r,P_inner_v.subs(sub),'numpy')
mm=1.0
# direct ground truth
Hv=pr*(f-E**2)-1+sp.sqrt(w)*sp.sqrt(E**2*pr**2+J**2/r**2); la=lambda ex:sp.lambdify((r,pr,m,E,J),ex,'numpy')
Hn=la(Hv); Hp=la(sp.diff(Hv,pr)); Hm=la(sp.diff(Hv,m)); Gpr=la(sp.diff(sp.diff(Hv,J)/sp.diff(Hv,pr),pr))
mv,Ev,Jv=1.0,1.4,2.5
def pr0(R):
    g=lambda p:Hn(R,p,mv,Ev,Jv); ps=np.linspace(-40,-1e-4,4000); vs=[g(p) for p in ps]
    s=np.where(np.diff(np.sign(vs)))[0]; return brentq(g,ps[s[0]],ps[s[0]+1])
def inner(x): p=pr0(x); return mv*Hm(x,p,mv,Ev,Jv)/Hp(x,p,mv,Ev,Jv)
def Sig(x): return quad(inner,r0,x,limit=200)[0]
def wrap_direct(x): return -quad(lambda t:(Gpr(t,pr0(t),mv,Ev,Jv)/Hp(t,pr0(t),mv,Ev,Jv))*Sig(t),r0,x,limit=120)[0]
# explicit letters (numeric nested integrals)
def U(x,k): return quad(lambda t:t**k/sq(t),r0,x,limit=200)[0]
def Pi_pole(x,rp): return quad(lambda t:1.0/((t-rp)*sq(t)),r0,x,limit=200)[0]      # third-kind
def Pi2_pole(x,rp): return quad(lambda t:1.0/((t-rp)**2*sq(t)),r0,x,limit=200)[0]  # double pole (2nd-kind)
def PiDE(x): return quad(lambda t:1.0/((Ev**2-1)*t+2*mm)/sq(t),r0,x,limit=200)[0]
# Sigma_2nd via letters: poly part + a2 Pi2 + a1 Pi + b PiDE
polyc=[float(c.subs(sub)) for c in polyq.all_coeffs()[::-1]]; a2n=float(a2.subs(sub)); a1n=float(a1.subs(sub)); bn=float(bDE.subs(sub))
def Sig2nd(x): return sum(polyc[k]*U(x,k) for k in range(len(polyc)))+a2n*Pi2_pole(x,2*mm)+a1n*Pi_pole(x,2*mm)+bn*PiDE(x)
def Sigelem(x): return 4*mm**2/(x-2*mm)-2*mm*np.log(x-2*mm) - (4*mm**2/(r0-2*mm)-2*mm*np.log(r0-2*mm))
# wrap via explicit Sigma = Sig2nd + Sigelem
def wrap_explicit(x): return -quad(lambda t:(Akv(t)/sq(t))*(Sig2nd(t)+Sigelem(t)),r0,x,limit=120)[0]
print("\nverify Sigma_2nd (letters) == numeric int P_inner_v/sqrtS, and full wrap:")
for x in [11.0,10.0,9.0]:
    s2_let=Sig2nd(x); s2_num=quad(lambda t:Pinv(t)/sq(t),r0,x,limit=200)[0]
    d=wrap_direct(x); we=wrap_explicit(x)
    print(f"  r={x}: Sig2nd d={abs(s2_let-s2_num):.1e} | wrap direct={d:+.7f} explicit={we:+.7f} diff={abs(d-we):.1e}")
print("\n=> if diffs ~0: P_inner_v reduced to U_k+double/simple third-kind; Sigma_elem closed; wrap explicit.")
