# (b) fast FAY: relations among the A[e_i,e_j] modulo weight 1, through DERIVATIVES (closed form).
# A'_k = e_i Pe_j - e_j Pe_i  (weight 1, closed: Pe = zeta/lnsigma/z).  W_basis' = derivatives of the products.
# rank of [A' | W'] -> null space with a non-zero A part = Fay relations.
import numpy as np, mpmath as mp, sympy as sp
from scipy.integrate import quad
mp.mp.dps=30
E=1.4; m=1.0; r0=12.0
r,J=sp.symbols('r J'); Es=sp.Rational(7,5)
Ssym=sp.expand(r*(r-2)*((Es**2-1)*r+2)*(r**2*(r-2)-J**2*((Es**2-1)*r+2)))
Jc=float([s for s in sp.solve(sp.Eq(sp.resultant(Ssym,sp.diff(Ssym,r),r),0),J) if s.is_real and float(s)>1][0])
def pS(mv,Jv):
    DE=np.array([E**2-1,2*mv]); p=np.polymul(np.polymul([1,0],[1,-2*mv]),DE)
    return np.polymul(p,np.polysub(np.polymul([1,0,0],[1,-2*mv]),Jv**2*DE))
Sc=pS(m,Jc); dmS=(pS(m+1e-6,Jc)-pS(m-1e-6,Jc))/2e-6
K=np.polymul([Jc],[E**2-1,2.0]); Nm=np.polysub(np.polymul(Sc,np.array([2*Jc])),0.5*np.polymul(K,dmS))
rts=np.roots(Sc); pr=[(i,j) for i in range(6) for j in range(i+1,6) if abs(rts[i]-rts[j])<1e-6]
rd=float(np.real((rts[pr[0][0]]+rts[pr[0][1]])/2)); Q,_=np.polydiv(Sc,np.polymul([1,-rd],[1,-rd]))
er=np.sort(np.real(np.roots(Q))); a4=Q[0]; e1r,e2r,e3r,e4r=er
def Q4(x): return np.polyval(Q,x)
k2=((e3r-e2r)*(e4r-e1r))/((e4r-e2r)*(e3r-e1r)); pref=2/mp.sqrt((e4r-e2r)*(e3r-e1r))/mp.sqrt(a4)
om1=mp.mpf(float(pref*mp.ellipk(k2))); w_im=float(pref*mp.ellipk(1-k2)); tau=mp.mpc(0,w_im)/om1; qn=mp.exp(mp.pi*1j*tau)
L1=lambda u: mp.jtheta(1,u,qn); L1p=lambda u: mp.jtheta(1,u,qn,1); L1pp=lambda u: mp.jtheta(1,u,qn,2); th1p0=L1p(0)
eta1=-(mp.pi**2/(12*om1))*(mp.jtheta(1,0,qn,3)/th1p0)
def wsig(z): u=mp.pi*z/(2*om1); return (2*om1/mp.pi)*mp.exp(eta1*z**2/(2*om1))*L1(u)/th1p0
def wzet(z): u=mp.pi*z/(2*om1); return eta1*z/om1+(mp.pi/(2*om1))*(L1p(u)/L1(u))
def wp(z):   u=mp.pi*z/(2*om1); rr=L1p(u)/L1(u); return -eta1/om1-(mp.pi/(2*om1))**2*(L1pp(u)/L1(u)-rr**2)
def wpp(z):
    u=mp.pi*z/(2*om1); T1=L1p(u)/L1(u); T2=L1pp(u)/L1(u); T3=mp.jtheta(1,u,qn,3)/L1(u)
    return -(mp.pi/(2*om1))**3*(T3-3*T1*T2+2*T1**3)
z_inf=float(mp.re(mp.quad(lambda x:1/mp.sqrt(Q4(x)),[e4r,mp.inf]))); sa=float(mp.sqrt(a4))
z_d=z_inf+float(mp.re(mp.quad(lambda x:1/mp.sqrt(Q4(x)),[-mp.inf,rd])))
def zr(rv): return float(quad(lambda x:1/np.sqrt(Q4(x)),e4r,rv,limit=400)[0])
z0=mp.mpf(zr(r0)); zdm=mp.mpf(z_d); zim=mp.mpf(z_inf); iw=1j*mp.mpf(w_im)
# functions e_i and their antiderivatives Pe_i (closed), based at z0
ef=[lambda z: mp.mpf(1),
    lambda z: wzet(z-zdm)-wzet(z+zdm),
    lambda z: wp(z-zdm)+wp(z+zdm),
    lambda z: wpp(z-zdm)-wpp(z+zdm),
    lambda z: wp(z),
    lambda z: wp(z-iw),
    lambda z: wp(z-zim)+wp(z+zim),
    lambda z: wzet(z-zim)-wzet(z+zim)]
def Pe(i,z):
    if i==0: return z-z0
    if i==1: return (mp.log(wsig(z-zdm))-mp.log(wsig(z+zdm)))-(mp.log(wsig(z0-zdm))-mp.log(wsig(z0+zdm)))
    if i==2: return (-(wzet(z-zdm)+wzet(z+zdm)))-(-(wzet(z0-zdm)+wzet(z0+zdm)))
    if i==3: return ((wp(z-zdm)-wp(z+zdm)))-((wp(z0-zdm)-wp(z0+zdm)))
    if i==4: return (-wzet(z))-(-wzet(z0))
    if i==5: return (-wzet(z-iw))-(-wzet(z0-iw))
    if i==6: return (-(wzet(z-zim)+wzet(z+zim)))-(-(wzet(z0-zim)+wzet(z0+zim)))
    if i==7: return (mp.log(wsig(z-zim))-mp.log(wsig(z+zim)))-(mp.log(wsig(z0-zim))-mp.log(wsig(z0+zim)))
Rall=[0,1,2,3,4,5]; Eall=[0,1,6,7]
pairs=[(i,j) for i in range(8) for j in range(i+1,8) if (i in Rall and j in Eall) or (i in Eall and j in Rall)]
enam=['1','Z_zd','P_zd','Pp_zd','wp0','wpiw','P_zi','Z_zi']
# A'_k(z) = e_i Pe_j - e_j Pe_i (closed, weight 1)
def dA(i,j,z): return ef[i](z)*Pe(j,z)-ef[j](z)*Pe(i,z)
# weight-1 basis W and its derivatives W' = d/dz(Pe_i Pe_j) = e_i Pe_j + e_j Pe_i ; d/dz(z Pe_i)=Pe_i+z e_i; d/dz Pe_i=e_i
wbasis=[]  # list of functions z -> derivative
for i in range(8):
    for j in range(i,8): wbasis.append((lambda i,j: (lambda z: ef[i](z)*Pe(j,z)+ef[j](z)*Pe(i,z)))(i,j))
for i in range(8): wbasis.append((lambda i:(lambda z: Pe(i,z)+z*ef[i](z)))(i))
for i in range(8): wbasis.append((lambda i:(lambda z: ef[i](z)))(i))
wbasis.append(lambda z: 2*z); wbasis.append(lambda z: mp.mpf(1))
zs=[mp.mpf(zr(rr)) for rr in np.linspace(11.9,8.9,120)]
M=np.array([[complex(dA(i,j,z)) for (i,j) in pairs]+[complex(w(z)) for w in wbasis] for z in zs])
nA=len(pairs)
# SVD; count the null space and how many null vectors have a non-trivial A part
U,S,Vt=np.linalg.svd(M,full_matrices=True)
S=S/S[0]; nz=int(np.sum(S>1e-9)); nnull=M.shape[1]-nz
# project: rank of the A part alone modulo the W part
MA=M[:,:nA]; MW=M[:,nA:]
MAperp=MA-MW@np.linalg.lstsq(MW,MA,rcond=None)[0]
sA=np.linalg.svd(MAperp,compute_uv=False); sA=sA/sA[0]
rankA=int(np.sum(sA>1e-7))
print(f"pairs A = {nA}")
print(f"rank of A modulo weight 1 = {rankA}")
print(f"=> independent Fay relations = {nA-rankA}")
print("singular values of A_perp (normalised):",np.array2string(sA,precision=2,max_line_width=220))

# ===== CHECK that the relations are REAL (on the nested A, not only on derivatives) =====
# null space of A_perp -> combinations lambda with sum lambda_k A_k = weight 1.
Uu,Ss,Vt2=np.linalg.svd(MAperp,full_matrices=True)
null=Vt2[rankA:]   # rows = null vectors (relations), 16 of them
print("\ncheck of %d relations on true NESTED integrals:"%null.shape[0])
_pc={}
def Penum(i,z):
    k=(i,float(z.real),float(z.imag))
    if k in _pc: return _pc[k]
    v=mp.quad(lambda t: ef[i](t),[z0,z]); _pc[k]=v; return v
def Atrue(i,j,z): return mp.quad(lambda t: ef[i](t)*Penum(j,t)-ef[j](t)*Penum(i,t),[z0,z],maxdegree=7)
ztest=[mp.mpf(zr(11.0)),mp.mpf(zr(10.0)),mp.mpf(zr(9.5))]
Avals=np.array([[complex(Atrue(i,j,z)) for (i,j) in pairs] for z in ztest])  # 3 x 21
Wvals=np.array([[complex(w(z)) for w in wbasis]+[1.0] for z in ztest])       # weight 1 (adding a constant)
# for each relation: R(z)=sum lambda A_k(z); it must lie in span(W) (weight 1)
for ri in range(min(4,null.shape[0])):
    lam=null[ri]
    Rz=Avals@lam                       # 3 values
    # residual of Rz with respect to span(Wvals)
    mu,res,rk,sv=np.linalg.lstsq(Wvals,Rz,rcond=None)
    resid=np.max(np.abs(Wvals@mu-Rz))
    print(f"  relation {ri}: |sum lambda*A_true - (weight-1 fit)| = {resid:.1e}   (|lambda|max={np.max(np.abs(lam)):.2f})")
