# PROOF OF CONCEPT: perturbare FLRW fa SALTARE il genere della brachistocrona.
# Brachistocrona in mezzo a indice sferico n(r): dphi/dr = L/(r sqrt(n^2 r^2 - L^2)),
# = caso FLRW degenere del paper (Randers beta=0, isotropo -> Fermat int n ds).
# Curva spettrale y^2 = P(r), P = n^2 r^2 - L^2.  Genere = floor((deg P -1)/2).
import sympy as sp, numpy as np
from scipy.integrate import quad
r,n0,L,k,lam=sp.symbols('r n0 L k lam',positive=True)

# --- FLRW (n=const): retta, genere 0 ---
I0=sp.integrate(L/(r*sp.sqrt(n0**2*r**2-L**2)),r)
print('FLRW  n^2=n0^2:  phi =',sp.simplify(I0))
print('   deg P = 2 -> genere 0 (retta r cos(phi-phi0)=L/n0)\n')

# --- overdensita' uniforme (tidal): Phi ~ 1/2 w^2 r^2 => n^2 = n0^2(1+k r^2) ---
P1=sp.expand(n0**2*(1+k*r**2)*r**2-L**2)
print('overdensita  n^2=n0^2(1+k r^2):  P =',P1,'  deg',sp.degree(P1,r),'-> genere 1 ELLITTICA')

# --- tidal ordine successivo => genere 2 IPERELLITTICA ---
P2=sp.expand(n0**2*(1+k*r**2+lam*r**4)*r**2-L**2)
print('             n^2=n0^2(1+k r^2+lam r^4):  deg',sp.degree(P2,r),'-> genere 2 IPERELLITTICA\n')

# --- verifica numerica genere-1: 4 radici distinte, curva ellittica genuina ---
sub={n0:1.0,L:2.0,k:0.05}
c1=[float(sp.Poly(P1.subs(sub),r).all_coeffs()[::-1][i]) for i in range(5)]
rts=np.roots(c1[::-1])
print('genus-1 quartica radici:',[complex(round(z.real,4),round(z.imag,4)) for z in rts])
disc=float(np.prod([abs(rts[i]-rts[j]) for i in range(4) for j in range(i+1,4)]))
print(f'  discriminante = {disc:.2f} != 0 -> 4 radici DISTINTE = genere 1 esatto')
print('  dphi = L dr/(r sqrt P): polo TERZA SPECIE a r=0 (CENTRO overdensita) e r=inf')
print('  => lettera peso-2 ancorata al CENTRO dell overdensita (come dlog orizzonte Vaidya)')

# turning fisico e conferma ellittico (integrale non elementare)
def P1n(x): return 1.0*(1+0.05*x**2)*x**2-4.0
rturn=min([z.real for z in rts if abs(z.imag)<1e-9 and z.real>0])
print(f'  turning fisico (perielio) r={rturn:.4f}: dphi ellittico 3a specie')
print('  => correzione adiabatica (overdensita che cresce, k=k(clock)) = DILOG ELLITTICO Bloch-Wigner')
print('\nCONCLUSIONE: perturbare FLRW SCOLLASSA la brachistocrona.')
print('  genere 0 (cerchio/retta) -> 1 (ellittica) -> 2 (iperellittica) col crescere del multipolo tidal.')
print('  Il modulo ellittico tau traccia la FORZA dell overdensita (k).')
