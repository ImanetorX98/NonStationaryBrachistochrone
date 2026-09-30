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
- dJc/dλ: Vaidya dJc/dm=Jc/m ; TK, **ramo τ**: dJc/dE=−E Jc r_d/DE(r_d). Per il ramo t vale invece la differenziazione implicita dJ/dE = −Q_{2,E}/Q_{2,J} sulla radice doppia.
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
## 8. AUDIT DELLE DIMOSTRAZIONI (24–30 settembre 2026)
═══════════════════════════════════════════════════════════════════
Conti espliciti fatti o rifatti durante gli audit (Claude, gpt-oss, qwen, GPT-6 Astra).
Metodo: SYM = identità esatta; INT = intervalli outward-rounded; HP = alta precisione
(non prova); CAS2 = secondo sistema indipendente.

### Paper I
| Risultato | Conto | Metodo | Script |
|---|---|---|---|
| Hamiltoniana di ramo H = h − ℓ (Thm I.1); H_v = h − 1 | identità | SYM | `paper1/verification/verify_paper1_core.py` |
| H_τ = h(p_r + 1/Ê, J) − f/Ê | identità | SYM | idem |
| Eulero con costo: (J∂_J + p_r∂_{p_r})H = H + ℓ_T; per τ, ℓ_τ = (f − ∂_{p_r}H_τ)/Ê | identità | SYM | idem |
| Termine di bordo S_D = [r p_r] − ∫ℓ_T dλ: −λ (ramo v), −Δτ (ramo τ) | identità d(r p_r)/dλ = H + ℓ + ΘH | SYM | idem |
| Controesempio: ℓ_τ = 7/5 + √291/15 ≠ 1 a m=1, r=6, Ê=7/5, J=0, p_r = −5√291/97 (ingoing) | valutazione esatta | SYM | idem |
| J = 0: S = r³(r−2m)²D_E, rango M = 8; disc_r S ∝ Ê¹²J¹⁰m¹⁶ q(J²) | rango, fattorizzazione | SYM | idem |
| Residui all'infinito di r^k dr/√S: 0, 0, −1/√s6, s5/(2 s6^{3/2}); dV₂: −1/√a4 | serie in t = 1/r | SYM | idem |
| Rotaia circolare Schwarzschild M=1, r=6, Ê=7/5: a^r = −241/1800, g(a,∂_t) = 0 | Christoffel | SYM | idem |
| Lemma I.A: velocità coordinate u^i/u^χ, mappa proiettiva; controesempio χ = t + x/2 | a mano | — | testo |
| rem:maxwell: 𝒥′ > 0 ⇒ periastro unico, niente multi-escursioni; sweep < π per raggi diversi | a mano | — | testo |
| a = e^{Ht} su sezioni chiuse: R = 12H² + 6e^{−2Ht}/R_c² (non de Sitter); de Sitter chiuso: ∫H dt/cosh Ht = π/2 − arctan sinh Ht₀ < π | sympy | SYM | scratch, formula nel testo |
| Vaidya ottica: ∂_v a_rr ∝ m′(2Ê² − 3f), ∂_v a_φφ ∝ m′(Ê² − 2f), mai nulli insieme su f > 0 | sympy | SYM | testo |
| Controesempio discesa: g = −f(t)dt² + [f(E₀²−f)/E₀²]δ ⇒ a_ij = δ_ij, W non CKV | sympy | SYM | testo |

### Paper II
| Risultato | Conto | Metodo | Script |
|---|---|---|---|
| Prova A in M simbolico: J* = [E²r(a²+r²)+2Ma²]/(aD_E), J** = […+Ma²]/(aD_E), vertice −2Ma/(r−2M), Q₂(J**) con M²a² | identità | SYM | `no_inversion_reduction.py` |
| Fattorizzazione N_G = (r−r_min)·4MrΔ/[(r_min−2M)DE_min]·W(r) (manca 1/DE_min nella versione vecchia, anche a M=1) | resto mod Q₂ | SYM + numerico M=1,2 | idem |
| (A) vale anche a a=0: W(r_min) = E²r_min³(r_min−2M), pendenza E²r_min²(r_min−M) | sostituzione | SYM | testo |
| Numeratori per ramo: τ: N_E = EJ r⁴(r−2M)²D_E, N_J = r³(r−2M)²D_E²; t: N_E = E r⁵ D_E[J(r−2M)+2Ma], N_J = E²r⁵D_E²; grado 7 | derivazione | SYM | testo |
| |B′|/√P = 2Ma[3(E²−1)r² − (5E²−13)Mr − 14M²]√(ΔD_E)/(E³r^{13/2}(r−2M)): W_sup algebrica | confronto | numerico | testo |
| V = rΔ/D_E crescente: numeratore di V′ in r = 2M+s a coefficienti positivi | polinomio | SYM | testo |
| Quarter regime: controesempio astratto W = 10 − V + V²/100, Φ′(1) = 1528√3/375 > 0 | integrale | SYM | testo |
| Lemma B: g = VW crescente e concava per b ≥ 1/2 (certificati A_c, P_c); α > 1/2 (r ≥ 5/2), α > 2/3 (r ≥ 3), α > 2/3 statico (r ≥ 9/4); c* = sech²U*, U* tanh U* = 1; C_c² = 16(b+1)⁴/[4(b+1)+c]³ | certificati polinomiali | SYM | `paper2/verification/verify_lemma_B.py` |
| Controesempi: quarto a b=1, r₀=2.01 (xΦ′ ∈ [0.005438, 0.005460]); unimodalità a b=1/100, r₀=10⁴ (segni +,−,+,− a r_min = 3, 6, 30, 9000) | quadratura | INT | idem |
| CAP RUN 6: 7 configurazioni, finestre fino a r_{1/2}, raccordo V(r₀)/2 (minimo 1.85e−10) | certificato | INT | `no_inversion_schwarzschild_CAP_grid.py`, `CAP_grid_certificates.log` |
| Figura BVP r₀=6, E=1.2: r_min^t = 4.735 vs r_{1/2} = 4.784; V/V₀ = 0.484 (istanza di (i),(ii)) | BVP | numerico | `KerrScripts/bvp_estremi_fissi.py` |
| Nessuna orbita circolare τ esterna: ∂_r[rΔ − J²D_E] = D_E V′ > 0 sulla shell | identità | SYM | testo |
| Costati: (J_τ/A)/J_η = f/E a a=0 (10/21 a f=2/3, E=7/5); J_η = J_t/A | derivata del funzionale | SYM | testo |
| TK: ∂_η a = 2Ê²AA′h₀/(Ê² − A²f₀)² ≠ 0 se A′ ≠ 0 | derivata | SYM | testo (`eq:aT-running`) |
| TK non Einstein: Ric_{ηr} = 2(A′/A)Γ̄^η_{ηr}, Γ̄^η_{ηr} = M(r²+a²)/(r²Δ) all'equatore | Ricci completo (a=0) | SYM | testo |
| Spettro dell'operatore di curvatura ottico (Schwarzschild congelato): K_t, K_t, K_tan; R_opt = 2(2K_t + K_tan); fuori piano cos²χ K_t + sin²χ K_tan | identità | SYM + CAS2 | `verify_myers_optical.wls` |
| Metrica nulla: K = −M(2r−3M)/r⁴ < 0 | identità | SYM | testo |
| Taş Prop. 4 con L = F_T²/2 (Legendre forte); lancio comune = W(0)=0; controesempio alla planarità affine (cono) e B_L − B_Lᵀ ≠ db | confronto | SYM | verifica indipendente nell'audit di GPT-6 Astra; testo |
| Bound W ≤ W_sup(r) < 0 sulle finestre r/M ∈ [6, 9.275] a a=0.9, E = 1.3 e 1.4 | certificato | INT | `paper2/verification/verify_wsup_window.py` |

Stato dei verificatori al 30/9: Paper I 2 Python + 9 WLS, Paper II 11 Python + 23 WLS, tutti
eseguiti senza FAIL sulla macchina dell'autore (alcuni WLS di Paper I non legano l'exit code ai
FAIL: fa fede l'output stampato); `make_provenance.py --check` passa (32 voci, digest 42ed5ab1779c5813).
Tabelle affermazione → script: `paper1/verification/README.md`, `paper2/verification/README.md`.

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
