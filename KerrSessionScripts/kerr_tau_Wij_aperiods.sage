# -*- coding: utf-8 -*-
# (A) CHIUSURA parte OLOMORFA via a-PERIODI (TK tau, genus-2).
# I coeff olomorfi alpha_k,beta_k di U_k = (poli, coeff simbolici g_i) + alpha_k u1 + beta_k u2
# sono determinati dagli a-PERIODI (principio, non fit; base pole normalizzata a a-periodo zero):
#   (alpha_k, beta_k) = omega^{-1} . [oint_{a_j} omega_k],   omega_k = r^k dr/y (modello PARI).
# Sage matrix_of_integral_values integra i differenziali polinomiali r^k dr/y sui cicli.
# Sanity: k=0,1 (olomorfi) -> vettori unita' (U_0=u2,U_1=u1). Period-level come C0/Ce separatrice.
from sage.all import QQ
from sage.schemes.riemann_surfaces.riemann_surface import RiemannSurface as SageRS
import numpy as np

M,a,E,J = QQ(1),QQ(9)/10,QQ(7)/5,QQ(5)/2
Rr=PolynomialRing(QQ,['x','y']); x,y=Rr.gens()
Dl=x**2-2*M*x+a**2; Em=(E**2-1)*x+2*M
Sint=(1250*x*(x-2*M)*Em*(x*Dl-J**2*Em)).change_ring(QQ)   # modello intero (come psi_forward_abel)
Xrs=SageRS(y**2-Sint, prec=80)
# a-periodi (2 cicli) dei differenziali x^k dx/y, k=0..4
diffs=[Rr(x**k) for k in range(5)]
MIV=np.array(Xrs.matrix_of_integral_values(diffs),dtype=complex)   # 5 x 4 (diff x cicli)
aper=MIV[:,:2]   # a-periodi (prime 2 colonne)
print("a-periodi oint_{a_j} x^k dx/y  (righe k=0..4, colonne a1,a2):")
for k in range(5): print(f"  k={k}: {np.round(aper[k],5)}")
omega=aper[:2,:]      # 2x2: a-periodi olomorfi (k=0,1)
ominv=np.linalg.inv(omega)
print("\nomega (a-per olomorfi k=0,1) =\n",np.round(omega,5))
print("\n=== coeff olomorfi (alpha_k,beta_k) = omega^{-1} . a-per(omega_k) ===")
for k in range(5):
    ab=aper[k]@ominv   # aper[k] = [alpha,beta] . omega  =>  [alpha,beta] = aper[k] . omega^{-1}
    tag=" (olomorfo: deve ~ e_k)" if k<2 else ""
    print(f"  k={k}: (alpha,beta)=({ab[0]:+.5f},{ab[1]:+.5f}){tag}")
print("\n=> k=0,1 danno vettori unita' (sanity). k=2,3,4: coeff olomorfi period-level di U_k")
print("   (determinati dagli a-periodi = principio, NON fit; dipendono dai periodi come C0/Ce).")
print("   Combinati coi coeff di polo simbolici g_i -> decomposizione canonica COMPLETA.")
