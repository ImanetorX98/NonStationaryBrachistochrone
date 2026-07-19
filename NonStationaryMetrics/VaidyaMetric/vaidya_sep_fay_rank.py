# (b) Riduzione di FAY: rango dei 21 A[e_i,e_j] MODULO peso-1, via SVD su molti z.
# rango < 21 => esistono relazioni (riducibile). rango = 21 => shuffle e' gia' minimale.
import numpy as np, mpmath as mp, sympy as sp
from scipy.integrate import quad
mp.mp.dps=25
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
def wzet(z): u=mp.pi*z/(2*om1); return eta1*z/om1+(mp.pi/(2*om1))*(L1p(u)/L1(u))
def wp(z):   u=mp.pi*z/(2*om1); rr=L1p(u)/L1(u); return -eta1/om1-(mp.pi/(2*om1))**2*(L1pp(u)/L1(u)-rr**2)
def wpp(z):  return (wp(z+mp.mpf('1e-10'))-wp(z-mp.mpf('1e-10')))/mp.mpf('2e-10')
def wsig(z): u=mp.pi*z/(2*om1); return (2*om1/mp.pi)*mp.exp(eta1*z**2/(2*om1))*L1(u)/th1p0
z_inf=float(mp.re(mp.quad(lambda x:1/mp.sqrt(Q4(x)),[e4r,mp.inf]))); sa=float(mp.sqrt(a4))
z_d=z_inf+float(mp.re(mp.quad(lambda x:1/mp.sqrt(Q4(x)),[-mp.inf,rd])))
def zr(rv): return float(quad(lambda x:1/np.sqrt(Q4(x)),e4r,rv,limit=400)[0])
z0f=mp.mpf(zr(r0))
zdm=mp.mpf(z_d); zim=mp.mpf(z_inf); iw=1j*mp.mpf(w_im)
efun=[lambda z: mp.mpf(1), lambda z: wzet(z-zdm)-wzet(z+zdm), lambda z: wp(z-zdm)+wp(z+zdm),
      lambda z: wpp(z-zdm)-wpp(z+zdm), lambda z: wp(z), lambda z: wp(z-iw),
      lambda z: wp(z-zim)+wp(z+zim), lambda z: wzet(z-zim)-wzet(z+zim)]
enam=['1','Z_zd','P_zd','Pp_zd','wp0','wpiw','P_zi','Z_zi']
Rset=[1,2,3,4,5]; Eset=[1,6,7]  # indici non-costanti in R e eta' (piu' l'1 condiviso)
Rall=[0,1,2,3,4,5]; Eall=[0,1,6,7]
pairs=[(i,j) for i in range(8) for j in range(i+1,8)
       if (i in Rall and j in Eall) or (j in Rall and i in Eall)]
# togli coppie con coeff sempre nullo (R×R o E×E): tienile se i in R&j in E o viceversa
pairs=[(i,j) for (i,j) in pairs if (i in Rall and j in Eall) or (i in Eall and j in Rall)]
_pc={}
def Pe(i,z):
    key=(i,float(z.real),float(z.imag))
    if key in _pc: return _pc[key]
    v=mp.quad(lambda t: efun[i](t),[z0f,z]); _pc[key]=v; return v
def Aij(i,j,z): return mp.quad(lambda t: efun[i](t)*Pe(j,t)-efun[j](t)*Pe(i,t),[z0f,z],maxdegree=6)
# campiona z su molti punti fisici
zs=[mp.mpf(zr(rr)) for rr in np.linspace(11.8,9.0,26)]
print("valuto 21 A e base peso-1 su",len(zs),"punti...")
# matrice A: righe=z, colonne=coppie
Amat=np.array([[complex(Aij(i,j,z)) for (i,j) in pairs] for z in zs])
# base peso-1: prodotti Pe_i*Pe_j (i<=j), z*Pe_i, Pe_i, z^2,z,1  (span del peso-1)
Wcols=[]
for z in zs:
    row=[]
    P=[complex(Pe(i,z)) for i in range(8)]; zz=complex(z)
    for i in range(8):
        for j in range(i,8): row.append(P[i]*P[j])
    for i in range(8): row.append(zz*P[i])
    for i in range(8): row.append(P[i])
    row+=[zz*zz,zz,1.0]
    Wcols.append(row)
Wmat=np.array(Wcols)
# proietta A sull'ortocomplemento di span(W): A_perp = A - W (W^+ A)
Wp=np.linalg.pinv(Wmat)
Aperp=Amat-Wmat@(Wp@Amat)
# rango di Aperp
sv=np.linalg.svd(Aperp,compute_uv=False)
sv=sv/sv[0]
rank=int(np.sum(sv>1e-8))
print("valori singolari (norm.):",np.array2string(sv,precision=2,max_line_width=200))
print(f"\nRANGO dei {len(pairs)} A modulo peso-1 = {rank}")
print(f"=> relazioni di Fay indipendenti: {len(pairs)-rank}")
