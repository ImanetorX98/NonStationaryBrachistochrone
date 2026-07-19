# Stato della verifica — brachistocrone adiabatiche (separatrici + J generico)

Documento riepilogativo di TUTTE le verifiche. Metodi:
- **SYM** = identità algebrica esatta (sympy/Mathematica Solve/Simplify)
- **NUM** = verifica numerica Python (riduzione, contorno, ODE)
- **SAGE** = periodi/θ genus-2 via Sage + abelfunctions (RiemannTheta)
- **WL** = cross-check INDIPENDENTE Mathematica (tool diverso)

═══════════════════════════════════════════════════════════════════
## 1. SEPARATRICI (genus-1) — coefficienti sorgente b1,b2,b3
═══════════════════════════════════════════════════════════════════
Formula universale b_i = residui R al polo triplo z_d, F=N/Q4, Q4^(k)(r_d) da S^(k+2)(r_d).
Simbolici in (M,a,E,r_d,Jc). Script: `SEP_COEFF_SYMBOLIC.py`, `sep_coeff_symbolic.py`.

| Ramo | Jc, r_d | b1,b2,b3 | SYM | NUM | WL |
|---|---|---|---|---|---|
| Vaidya τ | 7.0266, −3.3637 | 0.2704, 0.0326, 0.00987 | ✓ | contorno 1e-7 | Laurent 1e-16 |
| Vaidya v | = τ | = τ (sorgente uguale) | ✓ | ✓ | (=τ) |
| TK τ | 20.328, −7.130 | −1.836, −0.044, −0.048 | ✓ | ✓ | Laurent 1e-16 |
| TK t+ | 19.089, −6.621 | −1.617, −0.0737, −0.0353 | ✓ | block 2e-8 | Laurent 1e-38 |
| TK t− | −18.671, −6.588 | +1.617, +0.0708, +0.0343 | ✓ | ✓ | Laurent 1e-38 |
Cross-check WL: `crosscheck_sep_bi.wl` (τ), `crosscheck_tkt_bi.wl` (t±). Metodo INDIPENDENTE
(r(t) da ODE dr/dt=√Q4 + estrazione Laurent) vs formula h0/s³.

═══════════════════════════════════════════════════════════════════
## 2. SEPARATRICI — residui del CLOCK
═══════════════════════════════════════════════════════════════════
| Ramo | Residui | SYM | NUM | WL |
|---|---|---|---|---|
| τ (Vaidya/TK) | e1_zd=(r_d³−2M r_d²)/s, e2_zi=1/(E²−1) | ✓ | ✓ | (formula) |
| v (Vaidya) | z_d: E r_d³/s ; orizzonte: **4M** | ✓ | contorno 0 | Weierstrass nativo: z_d 1e-15, 4M esatto + PROVA SIMBOLICA |
| t (TK) | z_d: ρ_t(r_d)/s ; z(r±): R_Δ(r±)/[(r±−r∓)(r±−r_d)√Q4(r±)] | ✓ | contorno 1e-6 | formula + invariante 2M |
Invarianti puliti: v orizzonte=4M; t: res(r+)+res(r−)=**2M**; a4=E²−1.
Script: `sep_v_clock_residui.py`, `sep_t_clock_residui.py`, `crosscheck_clock_res.wl`.

═══════════════════════════════════════════════════════════════════
## 3. SEPARATRICI — costanti additive Ce, C0
═══════════════════════════════════════════════════════════════════
Ce = η'(0)+2e1_zd ζ(z_d)−2e2_zi ℘(z_∞)+2e1_zi ζ(z_∞) ; C0 = −Σ_a[b1 ζ+b2 ℘−b3/2 ℘'](z_∞−a).
Coeff (e_i,b_i) SIMBOLICI; valori ζ,℘ ai punti = **period-level**. Verif: `vaidya_sep_C0Ce_closed.py` 1e-8.

═══════════════════════════════════════════════════════════════════
## 4. TRACKING della separatrice (Jc mobile)
═══════════════════════════════════════════════════════════════════
N_tot = N + (dJc/dλ)N_J. TEOREMA: N_tot(r_d)=0 ⟹ b3^track=0 (polo triplo cancellato),
perché N_tot(r_d)=½K(r_d)S'(r_d)(dr_d/dλ)=0 dato S'(r_d)=0 (doppia radice).
- dJc/dλ: Vaidya dJc/dm=Jc/m ; TK dJc/dE=−E Jc r_d/DE(r_d).
- b_i^track SIMBOLICI in (M,a,E,r_d,Jc).
| | SYM | WL |
|---|---|---|
| N_tot(r_d)=0 mod {S,S'} | ✓ | **Mathematica: 0 esatto (Vaidya + TK)** |
Script: `sep_tracking_coeff.py`, `crosscheck_tracking.wl`.

═══════════════════════════════════════════════════════════════════
## 5. J GENERICO (genus-2) — coefficienti simbolici
═══════════════════════════════════════════════════════════════════
c_k (riduzione 2ª specie di ∂F), Q_kj=c_k b_j−c_j b_k, g_i (Taylor q6^{−1/2}), P(r) (T_alg).
| Ramo | Simbolici in | SYM | NUM | WL |
|---|---|---|---|---|
| TK-τ | (a,E,J) M=1 | ✓ | riduz 1e-15 | ✓ (g_i,c_k,Q_kj,T_alg) |
| Vaidya | (m,E,J) | ✓ | riduz 1e-15 | ✓ identità=0, c_k match Python |
| TK-t | (a,E,J) M=1 | ✓ | riduz 2.6e-17 | ✓ N_t poly, identità=0, g_i |
g_i notevoli: TK-τ/Vaidya g0=1/√(E²−1); TK-t g0=1/(E√(E²−1)) [a4=E²(E²−1)].
Script: `kerr_psi_explicit_verified.py`, `vaidya_generic_coeff.py`, `tk_t_generic_coeff.py`,
`crosscheck_generic.wl`, `crosscheck_fullsymbolic.wl`.

═══════════════════════════════════════════════════════════════════
## 6. J GENERICO — struttura ψ, naming θ, q-serie dilog
═══════════════════════════════════════════════════════════════════
ψ = ½Ê Σ Q_kj W_kj + ½Ê G_alg. G_alg elementare (P(r)+Σ log, verif 0). W_kj → θ[δ] agli e_±.
Dilog Λ = serie di nome Kronecker-Eisenstein genus-2 (NON Li₂; genus-2 senza formula prodotto).
| Ramo | montaggio ψ | naming θ | q-serie dilog (conv. geometrica) |
|---|---|---|---|
| TK-τ | end-to-end 1e-15 (NUM) | integrale 1e-6 (SAGE) | N=4→1e-13 (SAGE) |
| Vaidya | (stesso schema) | (SAGE) | N=4→2e-14 (SAGE) |
| TK-t | (stesso schema) | (SAGE) | N=4→4e-16 (SAGE) |
NB: θ genus-2 NON è nativa in Mathematica ⟹ questi restano verificati solo via Sage/abelfunctions.
TK-t ha ANCHE dilog agli orizzonti (da ρ_t, a z(r±)): stessa struttura, punti shiftati.
Script: `kerr_tau_Wij_*.py/.sage`, `*_dilog_qseries*.sage`.

═══════════════════════════════════════════════════════════════════
## 7. GERARCHIA period-level (per separatrici e J generico)
═══════════════════════════════════════════════════════════════════
- 🟢 RAZIONALE-simbolico (formule universali): residui b_i,e_i,g_i,c_k,Q_kj,P,res, coeff clock.
- 🟡 ALGEBRICO: radici e_i, invarianti g2,g3 (simmetrici nelle radici).
- 🔴 TRASCENDENTE (period-level): τ, punti marcati z_d,z_∞/e_±, ζ,℘,θ ai punti, Ce,C0, α,β.
  NON razionali, NON universali (cambiano coi parametri — dimostrato), valutati per-curva via
  procedura universale (radici→periodi→ζ,℘/θ). Come K(m): formula universale, valore per-modulo.

═══════════════════════════════════════════════════════════════════
## SINTESI
═══════════════════════════════════════════════════════════════════
- Ogni coefficiente ALGEBRICO è simbolico in tutte le variabili naturali (M,a,E + r_d,Jc per
  separatrice; M,a,E,J per generico) E cross-checkato con Mathematica (tool indipendente).
- I pezzi TRASCENDENTI (period-level) sono funzioni speciali valutate per-curva (via Sage per
  genus-2, Weierstrass nativo per genus-1).
- Risultati notevoli verificati: tracking cancella il polo triplo (b3=0); invarianti 4M, 2M;
  dilog genus-2 = serie di nome convergente non-Li₂.
- Documenti: `SEP_COEFF_SYMBOLIC.md`, `GENUS2_CLOSED_FORM.md`, `progress.md`, questo file.
