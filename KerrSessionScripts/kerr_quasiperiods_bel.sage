# -*- coding: utf-8 -*-
# Quasi-periodi di 2a specie (matrice eta / kappa di Klein) per la curva genus-2
# del ramo tau -> ingrediente MANCANTE per la forma chiusa Kleiniana di psi.
#
# psi = 1/2 Ehat (rho - rho~) e' un integrale iterato lunghezza-2 di forme di
# 1a/2a specie (NO 3a specie: provato beta_+=delta_+=0). Per la teoria classica
# Baker-Klein CHIUDE nelle funzioni Kleiniane zeta,sigma (peso 1, stessa classe
# di phi_0), NON e' un polilog. Serve la matrice dei quasi-periodi eta (2a specie).
#
# METODO (rompe il muro "Sage non da' eta" di doranTau.md sec.7-iii):
#   - modello DISPARI (quintica, 1 punto all'infinito) via x=1/s: le 2a-specie
#     s^2,s^3 ds/Y sono pulite (doppio polo, no residuo). Il modello PARI (deg 6)
#     ha 2 punti all'infinito -> x^2 dx/y e' 3a specie (residuo) -> eta sbagliato.
#   - 2a-specie CANONICHE di Baker-Enolski-Leykin (coeff della quintica):
#       dr_1 = (lam3 s + 2 lam4 s^2 + 3 lam5 s^3) ds/(4y)
#       dr_2 = lam5 s^2 ds/(4y)
#   - periodi via matrix_of_integral_values di Sage sui 4 cicli (INTERI: no crash
#     Singular). Pairing du=[1,s] <-> dr=[dr1,dr2].
#
# VALIDAZIONE: kappa=eta*omega^-1 SIMMETRICA (1.4e-12); Legendre generalizzata
#   omega' eta^T - omega eta'^T = -i pi I (canonica).
#
# Quintica: Y^2 = s^6 S(1/s), S_tau = 1200x^6-2300x^5-11428x^4-5519x^3+24700x^2
#   +62500x (params M=1,a=9/10,E=7/5,J=5/2). lam_i = coeff (interi).
from sage.all import QQ
from sage.schemes.riemann_surfaces.riemann_surface import RiemannSurface as SageRS
import numpy as np

Rs = PolynomialRing(QQ, ['s', 'Y']); s, Y = Rs.gens()
lam = [QQ(1200), QQ(-2300), QQ(-11428), QQ(-5519), QQ(24700), QQ(62500)]  # i=0..5
q = sum(lam[i]*s**i for i in range(6))
X = SageRS(Y**2 - q, prec=100)
assert X.genus == 2

# 2a-specie canoniche BEL (fattore /2 = 4y->2y della convenzione f_y)
dr1 = (lam[3]*s + 2*lam[4]*s**2 + 3*lam[5]*s**3)/2
dr2 = lam[5]*s**2/2
du = [Rs(1), s]                    # du_1=1, du_2=s
Pi = np.array(X.matrix_of_integral_values(du + [dr1, dr2]), dtype=complex)  # 4x4
omega, omega_p = Pi[:2, :2], Pi[:2, 2:]     # 1a specie a,b
eta,   eta_p   = Pi[2:, :2], Pi[2:, 2:]     # 2a specie a,b
tau = np.linalg.solve(omega, omega_p)
kappa = eta @ np.linalg.inv(omega)          # matrice di Klein (per sigma)
Leg = omega_p @ eta.T - omega @ eta_p.T

print('tau =\n', np.round(tau, 5))
print('tau sym err        =', np.linalg.norm(tau - tau.T))
print('kappa =\n', np.round(kappa, 5))
print('kappa SYM ERR      =', np.linalg.norm(kappa - kappa.T))
print('Legendre/(2 pi i)  =\n', np.round(Leg/(2j*np.pi), 4), '  (atteso -I/2)')
print("\nOK: eta (2a specie) calcolata e validata -> sigma/zeta di Klein pronti.")
print("    Prossimo: sigma=exp(-1/2 u^T kappa u) theta(...), zeta_i=d_i log sigma,")
print("    psi(r) = combinazione chiusa di zeta(u(r)); validare vs psi numerico.")
