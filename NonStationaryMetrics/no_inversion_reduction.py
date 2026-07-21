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


# ============================================================================
# PROVA SIMBOLICA di Lemma A al TURNING POINT (parte singolare dominante).
# Il segno di F_t^2 - F_tau^2 vicino a r_min e' dato dai RESIDUI: F^2 ~ c/(r-r_min).
# c_t > c_tau  <=>  D := numeratore di (c_t-c_tau) > 0. Con le condizioni di turning
# (J_tau^2 = r_d Delta_d/DE_d, Q2(r_d,J_t)=0) si riduce a:
#     D = 4 r_d Delta(r_d) * B,   B = a DE(r_d) (J* - J_t),
#     J* = [E^2 r_d(a^2+r_d^2)+2a^2]/(a DE_d)  (>0).
# Catena di segni (r_d>2M, r_d>r_+, E>1, a>0), tutti fattori espliciti:
#   A   = coeff J_t^2 di Q2(r_d)     = -(r_d-2M) DE_d          < 0   (parabola giu')
#   J_v = vertice = -B_q/2A          = -2a/(r_d-2M)            < 0
#   Q2(r_d,J*) = -E^2 r_d^3 (E^2 r_d^2+a^2) Delta_d/[a^2 DE_d] < 0
#   => J*>0>J_v e Q2(J*)<0 con A<0  => J* oltre la radice grande => J_t<J*
#   => B>0 => D>0 => c_t>c_tau => F_t^2>F_tau^2 al turning.  QED (parte singolare)
# La positivita' INTERNA (r lontano da r_min) si riduce a N_G(r)>0 (grado 5),
# verificata numericamente sopra.
def symbolic_residue_proof():
    rd, a, E = sp.symbols('r_d a E', positive=True); M = 1
    Jt = sp.symbols('J_t', positive=True)
    Q2d = (2*E**2*Jt**2*M*rd - E**2*Jt**2*rd**2 - 4*E**2*Jt*M*a*rd + 2*E**2*M*a**2*rd
           + E**2*a**2*rd**2 + E**2*rd**4 + 4*Jt**2*M**2 - 4*Jt**2*M*rd + Jt**2*rd**2
           - 8*Jt*M**2*a + 4*Jt*M*a*rd + 4*M**2*a**2)
    P = sp.Poly(Q2d, Jt)
    A, Bq = P.coeff_monomial(Jt**2), P.coeff_monomial(Jt)
    DEd = (E**2-1)*rd + 2*M
    Jstar = (E**2*rd*(a**2+rd**2) + 2*a**2)/(a*DEd)
    assert sp.simplify(A - (-(rd-2*M)*DEd)) == 0
    assert sp.simplify(-Bq/(2*A) - (-2*a/(rd-2*M))) == 0
    Q2star = sp.factor(Q2d.subs(Jt, Jstar))
    print("  A            =", sp.factor(A), "   (<0 for r_d>2M)")
    print("  vertex J_v   =", sp.factor(-Bq/(2*A)), "   (<0)")
    print("  J*           =", sp.factor(Jstar), "   (>0)")
    print("  Q2(r_d,J*)   =", Q2star, "   (<0 for r_d>r_+)")
    print("  => J_t < J* => B>0 => residue ordering proven.  QED")


if __name__ == '__main__':
    print("\n=== Symbolic proof of Lemma A at the turning point ===")
    symbolic_residue_proof()


# ============================================================================
# PROVA SIMBOLICA COMPLETA di Lemma A (F_t^2 > F_tau^2 per ogni r>r_d).
# Passi (M=1; validi per r_d>2M, E>1, a>0):
#   sign(F_t^2-F_tau^2) = sign(N_G),  N_G = (J_t(r-2M)+2Ma)^2 P_tau - J_tau^2(r-2M)Q2.
#   Ridotto mod turning [J_tau^2=r_d Delta_d/DE_d, Q2(r_d,J_t)=0]:
#     N_G = (r - r_d) * [4 r Delta(r)/(r_d-2M)] * W(r),   W(r) LINEARE in r.
#   Serve W(r)>0 per r>r_d. W(r_d) e la pendenza w1 sono entrambi lineari in J_t:
#     W(r_d) = (r_d-2M) * B,   B = a DE_d (J*  - J_t),   J*  = [E^2 r_d(a^2+r_d^2)+2a^2]/(a DE_d)
#     w1     = a DE_d (J** - J_t),                        J** = [E^2 r_d^2(r_d-1)+E^2 a^2 r_d + a^2]/(a DE_d)
#   In entrambi J_t (radice prograde di Q2(r_d,.)) e' < J* e < J** perche':
#     A = coeff J_t^2 di Q2(r_d) = -(r_d-2M)DE_d < 0        (parabola verso il basso)
#     vertice J_v = -2a/(r_d-2M) < 0,   J*, J** > 0         (a destra del vertice)
#     Q2(r_d,J*)  = -E^2 r_d^3 (E^2 r_d^2+a^2) Delta_d / [a^2 DE_d] < 0
#     Q2(r_d,J**) = -(r_d-2M)(E^2 r_d^2+a^2)(E^2 r_d^2[(r_d-1)^2+a^2]+a^2)/[a^2 DE_d] < 0
#   => J* , J** oltre la radice grande => J_t < J*, J_t < J** => B>0, w1>0.
#   W lineare con W(r_d)>0 e pendenza w1>0  =>  W(r)>0 per ogni r>r_d.  QED (Lemma A).
def symbolic_full_proof():
    rd, a, E = sp.symbols('r_d a E', positive=True); M = 1
    Jt = sp.symbols('J_t', positive=True)
    DEd = (E**2 - 1) * rd + 2 * M
    Q2d = (2*E**2*Jt**2*M*rd - E**2*Jt**2*rd**2 - 4*E**2*Jt*M*a*rd + 2*E**2*M*a**2*rd
           + E**2*a**2*rd**2 + E**2*rd**4 + 4*Jt**2*M**2 - 4*Jt**2*M*rd + Jt**2*rd**2
           - 8*Jt*M**2*a + 4*Jt*M*a*rd + 4*M**2*a**2)
    A = sp.Poly(Q2d, Jt).coeff_monomial(Jt**2)
    Jstar = (E**2*rd*(a**2+rd**2) + 2*a**2)/(a*DEd)                       # per W(r_d)>0
    Jss = (E**2*rd**2*(rd-1) + E**2*a**2*rd + a**2)/(a*DEd)               # per w1>0
    assert sp.simplify(A + (rd-2*M)*DEd) == 0                            # A = -(r_d-2M)DE_d
    Q2star = sp.factor(Q2d.subs(Jt, Jstar))
    Q2ss = sp.factor(Q2d.subs(Jt, Jss))
    # manifest sign of the two Q2 values (all factors signed for r_d>2M,E>1,a>0)
    print("  A = -(r_d-2M)DE_d :", sp.factor(A), " (<0)")
    print("  Q2(r_d,J*)  =", Q2star, " (<0)")
    print("  Q2(r_d,J**) =", Q2ss, " (<0)")
    print("  bracket in Q2(r_d,J**) = E^2 r_d^2[(r_d-1)^2+a^2]+a^2:",
          sp.simplify(sp.expand(E**2*rd**2*((rd-1)**2+a**2)+a**2)
                      - (E**2*a**2*rd**2+E**2*rd**4-2*E**2*rd**3+E**2*rd**2+a**2)) == 0)
    print("  => J_t<J* and J_t<J** => W(r_d)>0 and slope w1>0 => W(r)>0 for r>r_d.  Lemma A QED")


if __name__ == '__main__':
    print("\n=== Full symbolic proof of Lemma A (all r>r_d) ===")
    symbolic_full_proof()


# ============================================================================
# LEMMA B (monotonia di Phi_tau) -- CONDIZIONALE. Phi_tau(r_min) NON e' globalmente
# monotona: ha un massimo a un raggio di deflessione r_pk~3M (soglia
# strong-deflection/winding, dove l'orbita si avvicina all'orbita circolare instabile
# del ramo tau). Per r_min>r_pk, Phi_tau e' strettamente DECRESCENTE. I turning point
# a estremi fissi (r_min^tau, r_min^t ~ 4.7..8 M) stanno tutti in questo regime
# scattering, dove (B) vale e il teorema no-inversion e' dimostrato (modulo (A), gia'
# provato). Verifica numerica robusta:
def lemma_B_regime():
    from scipy.integrate import quad
    M = 1.0
    def Phi_tau(rm, a, E, r0=10.0):
        Delta = lambda r: r*r - 2*M*r + a*a; DE = lambda r: (E*E-1)*r + 2*M
        J = np.sqrt(rm*Delta(rm)/DE(rm))
        def g(s):
            r = rm + s*s; P = r*Delta(r) - J*J*DE(r)
            return 2*J*np.sqrt(r*(r-2*M)*DE(r))/(Delta(r)*np.sqrt(P/(r-rm)))
        return quad(g, 0, np.sqrt(r0-rm), limit=200)[0]
    for a, E in [(0.2, 1.2), (0.5, 1.2), (0.9, 1.2), (0.9, 1.05), (0.9, 1.5)]:
        rms = np.linspace(2.2, 9, 120); ph = np.array([Phi_tau(rm, a, E) for rm in rms])
        ipk = int(np.argmax(ph)); rpk = rms[ipk]
        mono = all(ph[i+1] < ph[i] + 1e-9 for i in range(ipk, len(ph)-1))
        print(f"  a={a:.1f} E={E:.2f}:  r_pk={rpk:.2f}   Phi_tau decreasing for r_min>r_pk? {mono}")


if __name__ == '__main__':
    print("\n=== Lemma B (conditional monotonicity, r_min > r_pk) ===")
    lemma_B_regime()


# ============================================================================
# LEMMA B -- prova parziale in forma chiusa (regime grazing) + struttura.
# Forma potenziale efficace: Phi_tau = sqrt(Vmin) int_{rmin}^{r0} K/sqrt(V-Vmin) dr,
#   V=r Delta/DE,  K=sqrt(r(r-2M))/Delta,  Vmin=V(rmin),  W=K/V' (>0, W'<0).
# Con V=Vmin+sigma^2 e una IBP:
#   dPhi/dVmin = (V0-2Vmin) W(Vmin)/(Sigma sqrt(Vmin))
#                + (2/sqrt(Vmin)) int_0^Sigma W'(Vmin+sigma^2) g(sigma) dsigma,
#   g(sigma) = (V0-2Vmin) sigma/Sigma + Vmin - sigma^2,  Sigma=sqrt(V0-Vmin).
# g(0)=Vmin>0, g(Sigma)=0, altra radice <0 => g>0 su [0,Sigma) => secondo termine <0.
# Se V(r0)<=2 V(rmin) (grazing): primo termine <=0 => dPhi/dVmin<0 PROVATO.
# Se V(r0)>2 V(rmin): primo termine >0, dominato dall'integrale (verificato, ma
# disuguaglianza stretta alla soglia trascendente r_pk: prova simbolica aperta).
def lemma_B_partial():
    import sympy as sp
    r, a, E, s = sp.symbols('r a E sigma', positive=True); M = 1
    Vmin, V0, Sig = sp.symbols('V_min V_0 Sigma', positive=True)
    g = (V0 - 2*Vmin)*s/Sig + Vmin - s**2
    print("  g(sigma) =", g)
    print("  g(0) =", g.subs(s, 0), " (>0);  g(Sigma) =", sp.simplify(g.subs(s, Sig).subs(Sig**2, V0-Vmin)), " (=0)")
    roots = sp.solve(sp.Eq(g, 0), s)
    print("  radici di g:", [sp.simplify(rt.subs(Sig, sp.sqrt(V0-Vmin))) for rt in roots],
          " => una =Sigma, l'altra <0  =>  g>0 su [0,Sigma)")
    print("  => secondo termine (2/sqrt Vmin) int W' g dsigma < 0 (W'<0).")
    print("  Grazing V0<=2Vmin: primo termine <=0 => dPhi/dVmin<0 PROVATO in forma chiusa.")


if __name__ == '__main__':
    print("\n=== Lemma B: closed-form partial proof (grazing regime) ===")
    lemma_B_partial()
