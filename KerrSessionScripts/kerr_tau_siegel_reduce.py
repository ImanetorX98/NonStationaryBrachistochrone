# -*- coding: utf-8 -*-
# TAPPA 2 q-serie: riduzione di Siegel della matrice di Riemann tau (genus-2).
# Minkowski-riduce Im(tau) via U in GL(2,Z) (tau -> U tau U^T) + shift Re(tau) mod 1
# (tau -> tau - B, B intero simmetrico) -> tutti i nomi |q_ij|<1 -> q-serie converge.
import numpy as np

tau=np.array([[0.13979+0.99018j, 0.30247-0.60597j],
              [0.30247-0.60597j,-0.57313+1.49301j]])
print("tau iniziale=\n",np.round(tau,4))
def nomes(t): return [abs(np.exp(1j*np.pi*t[i,j])) for i,j in [(0,0),(1,1),(0,1)]]
print("nomi iniziali |q11|,|q22|,|q12| =",np.round(nomes(tau),4))

# --- Minkowski reduce Im(tau) tramite U in GL(2,Z), applica tau->U tau U^T ---
U=np.eye(2,dtype=int)
for _ in range(50):
    Y=tau.imag
    if Y[0,0]>Y[1,1]:
        Sw=np.array([[0,1],[1,0]]); U=Sw@U; tau=Sw@tau@Sw.T; continue
    n=int(round(Y[0,1]/Y[0,0]))
    if n==0: break
    T=np.array([[1,-n],[0,1]]); U=T@U; tau=T@tau@T.T
print("\nU (GL(2,Z)) =",U.tolist()," det=",int(round(np.linalg.det(U))))
print("Im(tau) ridotta=\n",np.round(tau.imag,4)," (Minkowski: 2|Y12|<=Y11<=Y22)")
# --- shift Re(tau) mod 1 (B intero simmetrico) ---
B=np.round(tau.real).astype(int); B=np.round((B+B.T)/2).astype(int)  # simmetrico
tau=tau-B
print("B (shift Re) =",B.tolist())
print("tau ridotta=\n",np.round(tau,4))
print("nomi ridotti |q11|,|q22|,|q12| =",np.round(nomes(tau),4))
red=nomes(tau)
if max(red)<1: print("\n=> TUTTI i nomi <1: q-serie converge. Riduzione Siegel OK.")
else:
    print("\n=> q12 ancora >=1: serve trasformazione S (tau->-tau^{-1}) o modulare completa.")
    # prova S sul blocco: tau -> -tau^{-1} (Fricke), poi rireduci
    tS=-np.linalg.inv(tau); print("   tau dopo S=-tau^{-1}:\n",np.round(tS,4)," nomi:",np.round(nomes(tS),4))
print("\nInv: Im(tau) def positiva mantiene la somma reticolare convergente sempre;")
print("la riduzione serve solo per la SERIE DI POTENZE nei nomi (Kronecker-Eisenstein).")
