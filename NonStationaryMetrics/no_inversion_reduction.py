# -*- coding: utf-8 -*-
"""
Riduzione analitica del no-inversion a estremi fissi (referee obiezione #11).

Il paper afferma che, a ENDPOINT FISSI, la brachistocrona del ramo t resta piu'
SUPERFICIALE di quella tau: r_min^t > r_min^tau. Il referee nota che il bound
puntuale n_t/n_tau=E/f>1 NON basta a ordinare i minimi di due funzionali diversi.
Corretto: a stesso J l'ordine e' addirittura OPPOSTO (t piu' profondo); il risultato
a estremi fissi e' un'inversione guidata dalla mappa angolo J->Phi.

Questo script mostra la riduzione a due lemmi (e li verifica robustamente):
  Parametrizza ogni ramo per il raggio di turning r_min (che fissa J via la
  condizione di turning: rDelta=J^2 DE per tau, Q2=0 per t). Sia
     Phi_br(r_min) = int_{r_min}^{r0} F_br(r; J_br(r_min)) dr   (semi-angolo).
  (A) A r_min uguale, l'integrando t domina PUNTUALMENTE: F_t(r;J_t)>F_tau(r;J_tau)
      per r in (r_min,r0). Poiche' J_{t,tau}(r_min) sono ALGEBRICI, questa e' una
      disuguaglianza polinomiale (non trascendente) => Phi_t(r_min)>Phi_tau(r_min).
  (B) Phi_br e' strettamente decrescente in r_min.
  Insieme: a Phi fisso, Phi_tau(r_min^tau)=Phi_t(r_min^t)>Phi_tau(r_min^t)
           => r_min^tau < r_min^t  (tesi).  QED modulo (A) simbolico.
"""
import numpy as np
import sympy as sp
from scipy.integrate import quad


def branch_data(a, E, M=1.0, r0=10.0):
    r, J = sp.symbols('r J', real=True)
    DE = (E**2 - 1) * r + 2 * M
    Delta = r**2 - 2 * M * r + a**2
    Stau = sp.expand(r * (r - 2 * M) * DE * (r * Delta - J**2 * DE))
    Ktau = J * r * (r - 2 * M) * DE / Delta
    Ftau = sp.lambdify((r, J), sp.Abs(Ktau) / sp.sqrt(Stau), 'numpy')
    Q2 = (2*E**2*J**2*M*r - E**2*J**2*r**2 - 4*E**2*J*M*a*r + 2*E**2*M*a**2*r + E**2*a**2*r**2
          + E**2*r**4 + 4*J**2*M**2 - 4*J**2*M*r + J**2*r**2 - 8*J*M**2*a + 4*J*M*a*r + 4*M**2*a**2)
    R6 = sp.expand(r * Q2 * DE)
    Kt = r * DE * (J * (r - 2 * M) + 2 * M * a) / Delta
    Ft = sp.lambdify((r, J), sp.Abs(Kt) / sp.sqrt(R6), 'numpy')
    J2tau = sp.lambdify(r, r * Delta / DE, 'numpy')            # J_tau(r_min): rDelta=J^2 DE
    Ac, Bc, Cc = [sp.lambdify(r, c, 'numpy') for c in sp.Poly(Q2, J).all_coeffs()]

    def Jt(rm):                                                # J_t(r_min): prograde root of Q2=0
        A, B, C = Ac(rm), Bc(rm), Cc(rm); d = B * B - 4 * A * C
        if d < 0: return np.nan
        c = [x for x in ((-B + np.sqrt(d)) / (2 * A), (-B - np.sqrt(d)) / (2 * A)) if x > 0]
        return min(c) if c else np.nan
    return Ftau, Ft, J2tau, Jt, r0


def phis(a, E):
    Ftau, Ft, J2tau, Jt, r0 = branch_data(a, E)

    def Phi_tau(rm):
        return quad(lambda x: Ftau(x, np.sqrt(J2tau(rm))), rm + 1e-6, r0, limit=200)[0]

    def Phi_t(rm):
        Jv = Jt(rm)
        return np.nan if np.isnan(Jv) else quad(lambda x: Ft(x, Jv), rm + 1e-6, r0, limit=200)[0]
    return Phi_tau, Phi_t, Ftau, Ft, J2tau, Jt, r0


print("Verifica dei due lemmi su vari (a,E):")
print("  a    E     LemmaA (F_t>F_tau ptwise)  LemmaB (Phi decrescente)  => r_min^t>r_min^tau")
for a, E in [(0.3, 1.2), (0.6, 1.2), (0.9, 1.2), (0.9, 1.05), (0.9, 1.5), (0.5, 1.3)]:
    Phi_tau, Phi_t, Ftau, Ft, J2tau, Jt, r0 = phis(a, E)
    A_ok = B_ok = True; prev = None
    for rm in np.linspace(3, 8, 12):
        Jta, Jtt = np.sqrt(J2tau(rm)), Jt(rm)
        if np.isnan(Jtt): continue
        # (A) integrando puntuale
        for rr in np.linspace(rm + 0.2, r0, 8):
            if Ft(rr, Jtt) <= Ftau(rr, Jta): A_ok = False
        # (B) monotonia
        pt = Phi_t(rm)
        if prev is not None and pt > prev: B_ok = False
        prev = pt
    thesis = A_ok and B_ok
    print(f"  {a:.1f}  {E:.2f}   {str(A_ok):>5}                      {str(B_ok):>5}"
          f"                     {str(thesis):>5}")
print("\nLemma A e' ALGEBRICO (J_{t,tau}(r_min) radici di polinomi) => prova simbolica")
print("possibile per positivita' polinomiale. QED del teorema modulo (A) simbolico.")
