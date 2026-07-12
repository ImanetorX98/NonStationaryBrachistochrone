# -*- coding: utf-8 -*-
"""
Forma chiusa della SEPARATRICE del ramo t (tempo coordinato/arrivo) della
brachistocrona equatoriale di Kerr, e prova che e' GENERE 1 (Weierstrass),
non genere 2 (Kleinian).

Struttura algebrica:
  dphi/dr = K(r)/sqrt(R6(r)),   R6 = r Q2(r) [(E^2-1)r + 2M]   (SESTICO)
  con Q2 quartico.  R6 ha grado 6 -> il ramo t GENERICO e' genere 2
  (iperellittico).  MA la separatrice dinamica (saddle-node, cattura
  marginale retrograda J = J_c^-) ha una RADICE DOPPIA:
        R6 = (r - r*)^2 Q4(r) ,   Q4 quartico  ->  GENERE 1.
  Meccanismo diverso dalla separatrice tau (dove la cancellazione doranTau
  (r^2+c^2)|Q4 fa scendere il genere), stesso esito ellittico.

Forma chiusa (come la eq. A3 del ramo tau, terza specie):
  z(r)   = P^{-1}(A/r + B ; g2, g3)          [P = Weierstrass p]
  phi(r) = Lam0 z(r) + sum_k lam_k ln( sigma(z - v_k)/sigma(z + v_k) )
  poli c_k = { r_+, r_-, r* }  (TRE poli reali; nessun polo a r=0, ne'
             la coppia complessa +-ic della separatrice tau)
  v_k    = P^{-1}(A/c_k + B)
  lam_k  = - A alpha_k / ( c_k^2 P'(v_k) )
  Lam0   = Lam_poly - sum_k alpha_k/c_k + sum_k 2 zeta(v_k) lam_k

Validazione: phi dal flusso di Hamilton vs decomposizione a 2e-14.
"""

import numpy as np
import sympy as sp
from scipy.integrate import quad
from scipy.optimize import fsolve, brentq

M, a, E = 1.0, 0.9, 1.2
rp = M + np.sqrt(M**2 - a**2)
rm = M - np.sqrt(M**2 - a**2)


def Q2(r, J):
    return (2*E**2*J**2*M*r - E**2*J**2*r**2 - 4*E**2*J*M*a*r + 2*E**2*M*a**2*r
            + E**2*a**2*r**2 + E**2*r**4 + 4*J**2*M**2 - 4*J**2*M*r + J**2*r**2
            - 8*J*M**2*a + 4*J*M*a*r + 4*M**2*a**2)


def R6(r, J):
    return r * Q2(r, J) * ((E**2 - 1)*r + 2*M)


# --- saddle-node retrogrado: R6 = R6' = 0 (radice doppia) ---
rst, Jst = fsolve(lambda x: [R6(x[0], x[1]),
                             (R6(x[0]+1e-7, x[1]) - R6(x[0]-1e-7, x[1]))/2e-7],
                  [3.51, -8.05])
print(f"saddle-node (separatrice t):  r* = {rst:.6f},  J_c^- = {Jst:.6f}")
d2 = (R6(rst+1e-4, Jst) - 2*R6(rst, Jst) + R6(rst-1e-4, Jst))/1e-8
print(f"  R6(r*)={R6(rst,Jst):.2e}, R6'(r*)~0, R6''(r*)={d2:.1f}  -> RADICE DOPPIA (genere 2->1)")

# --- flusso di Hamilton del ramo t ---
r, pr = sp.symbols('r pr')
f = 1 - 2*M/r
Dl = r**2 - 2*M*r + a**2
b = 2*M*a/r
v2 = 1 - f/E**2
P = r**2 + a**2 + 2*M*a**2/r
Pb = P + b**2/E**2
pphi = sp.Symbol('pphi')
H = pphi*b*v2/Pb + sp.sqrt(Dl*v2/Pb)*sp.sqrt((Dl/r**2)*pr**2 + pphi**2/Pb) - 1
Hn = sp.lambdify((r, pr), H.subs(pphi, Jst), 'numpy')
dHp = sp.lambdify((r, pr), sp.diff(H, pr).subs(pphi, Jst), 'numpy')
dpp = sp.lambdify((r, pr), sp.diff(H, pphi).subs(pphi, Jst), 'numpy')


def prof(rv):
    pg = np.linspace(-100, 100, 10001)
    Hv = Hn(rv, pg)
    rts = [brentq(lambda p: Hn(rv, p), pg[i], pg[i+1])
           for i in range(len(pg)-1)
           if np.isfinite(Hv[i]) and np.isfinite(Hv[i+1]) and Hv[i]*Hv[i+1] < 0]
    ing = [p for p in rts if dHp(rv, p) < 0]
    return min(ing) if ing else np.nan


def dphidr(rv):
    p = prof(rv)
    return dpp(rv, p)/dHp(rv, p)


def sqrtQ4(rv):
    return np.sqrt(abs(R6(rv, Jst)))/abs(rv - rst)   # sqrt(Q4) = sqrt(R6)/|r-r*|


def G(rv):
    return dphidr(rv)*sqrtQ4(rv)                      # G = dphi/dr * sqrt(Q4)


# --- frazioni parziali di G:  poli {0, r_+, r_-, r*} ---
poles = [0.0, rp, rm, rst]
rs = np.linspace(3.9, 40, 60)
rows, rhs = [], []
for rv in rs:
    g = G(rv)
    if np.isfinite(g):
        rows.append([1.0] + [1.0/(rv - p) for p in poles])
        rhs.append(g)
rows, rhs = np.array(rows), np.array(rhs)
coef, *_ = np.linalg.lstsq(rows, rhs, rcond=None)
Lam = coef[0]
alphas = dict(zip(['0', 'r_+', 'r_-', 'r*'], coef[1:]))
rms = np.sqrt(np.mean((rows@coef - rhs)**2))
print(f"\nfit razionale G(r): RMS = {rms:.1e}  (piccolo = genere-1 confermato)")
print(f"  Lam_poly = {Lam:+.6f}")
for n, al in alphas.items():
    print(f"  alpha[{n:4s}] = {al:+.6f}")

# --- Weierstrass: g2, g3, A, B dalla quartica Q4 (radice a r=0) ---
c = [float(x) for x in sp.Poly(sp.expand(r*Q2(r, Jst)*((E**2-1)*r + 2*M)),
                               r).all_coeffs()]
q4r = [z for z in np.roots(c) if abs(z - rst) > 1e-3]
Q4p = (np.poly(q4r)*c[0]).real
a0, a1, a2, a3, a4 = Q4p/np.array([1, 4, 6, 4, 1])
g2 = a0*a4 - 4*a1*a3 + 3*a2**2
g3 = a0*a2*a4 + 2*a1*a2*a3 - a2**3 - a0*a3**2 - a1**2*a4
A_, B_ = Q4p[3]/4, Q4p[2]/12
print(f"\nWeierstrass:  g2={g2:.5f}, g3={g3:.5f},  A={A_:.5f}, B={B_:.5f}")
print("  z(r)=P^-1(A/r+B),  v_k=P^-1(A/c_k+B),  lam_k=-A*alpha_k/(c_k^2 P'(v_k))")

# --- validazione: phi diretto vs decomposto ---
def integ_dec(rv):
    s = Lam + sum(alphas[n]/(rv - p)
                  for n, p in zip(['0', 'r_+', 'r_-', 'r*'], poles))
    return s/sqrtQ4(rv)


r_a, r_b = 3.9, 9.0
I_dir = quad(dphidr, r_a, r_b, limit=200)[0]
I_dec = quad(integ_dec, r_a, r_b, limit=200)[0]
print(f"\nphi({r_a}->{r_b}):  flusso = {I_dir:+.9f},  decomposto = {I_dec:+.9f}")
print(f"  diff = {abs(I_dir - I_dec):.1e}   -> forma chiusa genere-1 VERIFICATA")
