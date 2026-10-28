# -*- coding: utf-8 -*-
# q-SERIES STEP 2: Siegel reduction of the Riemann matrix tau (genus-2).
# Minkowski-reduce Im(tau) through U in GL(2,Z) (tau -> U tau U^T) + shift Re(tau) mod 1
# (tau -> tau - B, B integer symmetric) -> all nomes |q_ij|<1 -> the q-series converges.
import numpy as np

tau=np.array([[0.13979+0.99018j, 0.30247-0.60597j],
              [0.30247-0.60597j,-0.57313+1.49301j]])
print("initial tau=\n",np.round(tau,4))
def nomes(t): return [abs(np.exp(1j*np.pi*t[i,j])) for i,j in [(0,0),(1,1),(0,1)]]
print("initial nomes |q11|,|q22|,|q12| =",np.round(nomes(tau),4))

# --- Minkowski reduce Im(tau) through U in GL(2,Z), apply tau->U tau U^T ---
U=np.eye(2,dtype=int)
for _ in range(50):
    Y=tau.imag
    if Y[0,0]>Y[1,1]:
        Sw=np.array([[0,1],[1,0]]); U=Sw@U; tau=Sw@tau@Sw.T; continue
    n=int(round(Y[0,1]/Y[0,0]))
    if n==0: break
    T=np.array([[1,-n],[0,1]]); U=T@U; tau=T@tau@T.T
print("\nU (GL(2,Z)) =",U.tolist()," det=",int(round(np.linalg.det(U))))
print("reduced Im(tau)=\n",np.round(tau.imag,4)," (Minkowski: 2|Y12|<=Y11<=Y22)")
# --- shift Re(tau) mod 1 (B integer symmetric) ---
B=np.round(tau.real).astype(int); B=np.round((B+B.T)/2).astype(int)  # symmetric
tau=tau-B
print("B (shift Re) =",B.tolist())
print("reduced tau=\n",np.round(tau,4))
print("reduced nomes |q11|,|q22|,|q12| =",np.round(nomes(tau),4))
red=nomes(tau)
if max(red)<1: print("\n=> ALL nomes <1: the q-series converges. Siegel reduction OK.")
else:
    print("\n=> q12 still >=1: the S transformation (tau->-tau^{-1}) or a full modular one is needed.")
    # try S on the block: tau -> -tau^{-1} (Fricke), then reduce again
    tS=-np.linalg.inv(tau); print("   tau after S=-tau^{-1}:\n",np.round(tS,4)," nomes:",np.round(nomes(tS),4))
print("\nNote: a positive definite Im(tau) always keeps the lattice sum convergent;")
print("the reduction is needed only for the POWER SERIES in the nomes (Kronecker-Eisenstein).")
