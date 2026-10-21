# -*- coding: utf-8 -*-
"""
Termine mancante d_J F nella sorgente adiabatica conforme di Thakurta-Kerr.

I pesi conformi E_eff = Ehat/A e J_eff = J/A sono quelli della derivazione canonica
del manoscritto (sorgente di Eulero (E d_E + J d_J)F); qui non vengono ridedotti.
Poiche' la shape frozen e' quella di Kerr, F = F(r; E_eff, J_eff), la variazione al
prim'ordine (A = A_0(1+eps eta), eps = A'/A) contiene DUE termini:

    delta phi = -int eta [ Ehat d_E F + J d_J F ] dr        (E_0=Ehat, J_0=J a A_0=1)

Questo script verifica che il termine J d_J F si riduce nella stessa base di d_E F e
tabula i coefficienti c_k e b_i con e senza di esso.
Una versione precedente usava solo -Ehat int eta d_E F dr, senza il termine J d_J F.
Struttura: (E d_E + J d_J) F e' l'operatore di scala; il suo numeratore si riduce
nella STESSA base di building block di d_E F (nessun mattone nuovo), quindi la
correzione entra come coefficienti c_k^J (generico) o b_i^J (separatrice), additivi
a quelli in E.

NB: Vaidya (parametro m, non conforme) NON e' toccato: E,J entrambi conservati,
sorgente d_m F, nessun termine d_J F.

Verifica a livello di ASSEMBLY completo (non solo coefficienti): sostituendo la
sorgente di Eulero (E d_E + J d_J)F in tk_t_sep_blockassembly.py, il block assembly
con TUTTE le funzioni speciali (Weierstrass zeta,wp via theta + dilog genus-2,
143 prodotti) riproduce l'integrale diretto a ~1e-8 sulla separatrice TK-t --
stessa precisione del solo termine E. Quindi la forma chiusa e' completa anche col
termine J: nessuna funzione speciale nuova, solo coefficienti c_k^tot=c_k^E+c_k^J.

Uscita: c_k^J e b_i corretti per TK ramo t e tau, con verifica di chiusura.
Numeratori di ramo: t  K_t = r D_E (J(r-2M)+2Ma)/Delta;  tau  K_tau = Jc r (r-2M) D_E/Delta
(il pattern sferico J D_E vale solo per a = 0).
I parametri numerici sono sostituiti subito dopo le derivate in E e J, e il J critico e'
ottenuto dalle radici reali dei fattori del risultante: stessa matematica, pochi secondi
invece di ore di semplificazione simbolica.
"""
import sympy as sp
import numpy as np

r, M, a, E, J, Jc, rd = sp.symbols('r M a E J Jc r_d', real=True)


# ---------------------------------------------------------------- generico
def reduce_2nd(F, R6, sub, label):
    """Riduce d_E F, d_J F e (E d_E + J d_J)F nella base di 2a specie r^k/sqrt(R6)."""
    out = {}
    for tag, src in (('E', E * sp.diff(F, E)), ('J', J * sp.diff(F, J)),
                     ('tot', E * sp.diff(F, E) + J * sp.diff(F, J))):
        N = sp.expand(sp.cancel((src * R6**sp.Rational(3, 2)).subs(sub)))
        R6s = sp.expand(R6.subs(sub))
        Ac = [sp.Symbol(f'A{i}') for i in range(6)]
        ck = [sp.Symbol(f'c{i}') for i in range(5)]
        A5 = sum(Ac[i] * r**i for i in range(6))
        Mp = sum(ck[i] * r**i for i in range(5))
        eq = sp.expand(2 * N - (2 * R6s * sp.diff(A5, r) - A5 * sp.diff(R6s, r) + 2 * R6s * Mp))
        sol = sp.solve(sp.Poly(eq, r).all_coeffs(), Ac + ck, dict=True)
        assert sol, f"{label}/{tag}: riduzione 2a specie NON chiude"
        out[tag] = [float(sol[0][ck[i]]) for i in range(5)]
    return out


print("=" * 70)
print("GENERICO ramo t (J qualsiasi): c_k^E, c_k^J, c_k^tot  [chiudono nella stessa base]")
print("=" * 70)
DE = (E**2 - 1) * r + 2 * M
Dl = r**2 - 2 * M * r + a**2
Q2 = (2*E**2*J**2*M*r - E**2*J**2*r**2 - 4*E**2*J*M*a*r + 2*E**2*M*a**2*r + E**2*a**2*r**2
      + E**2*r**4 + 4*J**2*M**2 - 4*J**2*M*r + J**2*r**2 - 8*J*M**2*a + 4*J*M*a*r + 4*M**2*a**2)
R6g = sp.expand(r * Q2 * DE)
Ktg = r * DE * (J * (r - 2 * M) + 2 * M * a) / Dl
subg = {M: 1, a: sp.Rational(9, 10), E: sp.Rational(6, 5), J: sp.Rational(5, 2)}
red = reduce_2nd(Ktg / sp.sqrt(R6g), R6g, subg, "t-gen")
print("  c_k^E   =", [round(x, 5) for x in red['E']], "  (riduzione di E d_E F)")
print("  c_k^J   =", [round(x, 5) for x in red['J']], "  <-- NUOVO (riduzione di J d_J F)")
print("  c_k^tot =", [round(x, 5) for x in red['tot']], "  (riduzione di (E d_E + J d_J)F)")
chk = [red['E'][i] + red['J'][i] for i in range(5)]
print("  c^E + c^J =", [round(float(x), 5) for x in chk], "  (deve = c_k^tot, per linearita')")


# ---------------------------------------------------------------- separatrice
def bcoeffs(S, N):
    d = lambda ex, n: sp.diff(ex, r, n).subs(r, rd)
    S2, S3, S4 = d(S, 2), d(S, 3), d(S, 4)
    Q4rd, Q4prd, Q4pprd = S2 / 2, S3 / 6, S4 / 12
    s = sp.sqrt(Q4rd); a1 = Q4prd / (4 * s); a2 = Q4pprd / 12
    Nr, Np, Npp = N.subs(r, rd), sp.diff(N, r).subs(r, rd), sp.diff(N, r, 2).subs(r, rd)
    h0 = Nr / Q4rd
    Fp = (Np * Q4rd - Nr * Q4prd) / Q4rd**2
    Fpp = ((Npp * Q4rd - Nr * Q4pprd) * Q4rd - 2 * Q4prd * (Np * Q4rd - Nr * Q4prd)) / Q4rd**3
    h1 = Fp * s; h2 = sp.Rational(1, 2) * (Fpp * s**2 + Fp * (Q4prd / 2))
    return {'b3': h0 / s**3, 'b2': (h1 - 3 * a1 * h0) / s**3,
            'b1': (h2 - 3 * a1 * h1 + (6 * a1**2 - 3 * a2) * h0) / s**3}


def sep_bi(S, K, params, name, Jmin=5):
    Jv = sp.Symbol('Jvar')
    Sj, Kj = S.subs(Jc, Jv), K.subs(Jc, Jv)
    F = Kj / sp.sqrt(Sj)
    Sub0 = {M: params['M'], a: params['a'], E: params['E']}
    N_E = sp.expand(sp.cancel((sp.diff(F, E) * Sj**sp.Rational(3, 2)).subs(Sub0))).subs(Jv, Jc)
    N_tot = sp.expand(sp.cancel(((E * sp.diff(F, E) + Jv * sp.diff(F, Jv)) * Sj**sp.Rational(3, 2)).subs(Sub0))).subs(Jv, Jc)
    S = sp.expand(S.subs(Sub0))
    Sn = S.subs(Sub0)
    res = sp.Poly(sp.resultant(Sn, sp.diff(Sn, r), r), Jc)
    reals = []
    for fac, _ in res.factor_list()[1]:
        reals += [float(z) for z in fac.real_roots()]
    Jcv = min(z for z in reals if z > Jmin)
    rts = np.roots([complex(c) for c in sp.Poly(Sn.subs(Jc, Jcv), r).all_coeffs()])
    prs = [(i, j) for i in range(len(rts)) for j in range(i + 1, len(rts)) if abs(rts[i] - rts[j]) < 1e-5]
    rdv = float(np.real((rts[prs[0][0]] + rts[prs[0][1]]) / 2))
    sub = {M: params['M'], a: params['a'], E: params['E'], Jc: Jcv, rd: rdv}
    bE, bT = bcoeffs(S, N_E), bcoeffs(S, N_tot)
    print(f"\n  {name}:  Jc={Jcv:.5f}  r_d={rdv:.5f}   (NB r_d<0: cfr. separatrix_classification.py)")
    E0 = float(params['E'])
    for k in ['b1', 'b2', 'b3']:
        v0, v1 = float(bE[k].subs(sub)), float(bT[k].subs(sub))
        en = E0 * v0                    # energy-only source Ehat d_E F
        print(f"    {k}:  d_E F={v0:+.6f}   Ehat d_E F={en:+.6f}   totale={v1:+.6f}   "
              f"assiale J d_J F={v1 - en:+.6f}  ({100 * (v1 - en) / en:+.2f}% di Ehat d_E F)")


print()
print("=" * 70)
print("SEPARATRICE: b_i paper (solo d_E) vs corretto (E d_E + J d_J)")
print("=" * 70)
S_tau = sp.expand(r * (r - 2 * M) * DE * (r * Dl - Jc**2 * DE))
sep_bi(S_tau, Jc * r * (r - 2 * M) * DE / Dl, {'M': 1, 'a': sp.Rational(9, 10), 'E': sp.Rational(6, 5)}, "TK tau")
Q2c = Q2.subs(J, Jc)
sep_bi(sp.expand(r * Q2c * DE), r * DE * (Jc * (r - 2 * M) + 2 * M * a) / Dl,
       {'M': 1, 'a': sp.Rational(9, 10), 'E': sp.Rational(6, 5)}, "TK t")
print()
print("Vaidya: NON toccato (parametro m, non conforme; E,J conservati; sorgente d_m F).")
print("FATTO.")
