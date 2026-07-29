# TODO MASTER — stato consolidato (agg. 28 lug 2026)

Fonte unica dei punti aperti. Consolida: referee report 2026-07-22 (13 issue), propagazione testo
main11, chiusura off-shell, tracking separatrice TK, e la nota derivazioni_post_referaggio.
Legenda: [x] fatto · [~] parziale · [ ] aperto · **BLOCKER/MAJOR/MINOR** severità referee.

---
## A. Fondamenti control-theory (referee Issues 1,2,6) — coperti dalla nota derivazioni
- [x] **Esistenza** minimo assoluto (dominio compatto regolare, indicatrice unif. convessa) — Teorema 5.1, VALIDO
- [x] **Normalità** PMP (ogni minimo regolare = estremale normale) — Prop 3.1, VALIDO
- [x] **Limite stazionario** `H≡0`, recupero energia fissata — Prop 3.2, VALIDO
- [x] **Certificato globale HJB** (sub-soluzione Lipschitz ⟹ ottimalità) — Teorema 7.1, VALIDO
- [x] **Framework coniugati/Maxwell/cut locus** — §8, VALIDO
- [ ] **MAJOR — Equivalenza esplicita con Perlick** (`F_rail=F_Perlick` o `X_{H_rail}∥X_{H_Perlick}` su `H=0`):
      §3.3 la lascia da mostrare. Serve come proposizione autonoma nel paper.
- [~] **Issue 2 — `p_φ` costato vs momento meccanico `L_mech`** — DERIVAZIONE FATTA (commit ad3a372,
      `costate_finsler_derivation.py`); resta solo la riscrittura del testo §4.1. NON invalida i risultati.
      Derivazione rigorosa di Finsler (Euler-verificata):
      * il `(p_r0,J)` del paper **È** il costato di Pontryagin, NON il meccanico geodetico (p_r0=1.29475
        combacia col costato Finsler 1.29487; geodetico u_r=1.049 sbagliato del 20%);
      * costato `p_φ` **≠** meccanico `u_φ=g(u,∂_φ)` (rapporto 1.125; brachistocrona non-geodetica) →
        il paper fa BENE a usare il costato;
      * il costato grezzo `J` è **conservato (peso 0)**; `J_eff=J/A` è la NORMALIZZAZIONE ottica `P=p/A`
        (confermato dall'oracolo `test_adiabatic_noreg.py:98`: `E0,J0` FISSI, `J_eff=J0/A` passata a H̄).
      Nessuna tensione: `J` (costato, peso 0) e `J_eff` (normalizzato, peso −1) legati da `/A`.
      NEXT (solo TESTO §4.1): sostituire la giustificazione "peso conforme delle cariche di Kerr"
      (argomento meccanico) con: "`J` = costato di Pontryagin conservato, distinto da `g(u,∂_φ)` perché
      la brachistocrona è non-geodetica; `J_eff=J/A` è la normalizzazione ottica". Unificare
      `rail_conservation.py` (nomenclatura costato vs meccanico).
- [ ] Integrare il "Testo proposto" (§10 nota) come blocco control-theory nel paper (main.tex + PRD).
- [ ] **MAJOR(Issue 6) — esclusioni**: minimalità dentro l'ergosfera (selettore spacelike, continuazione
      analitica ≠ ottimo fisico); PMP necessario ≠ esistenza/minimalità globale (già chiarito dalla nota,
      da riflettere nel testo); freezing≠turning≠separatrice (notazione).

## B. Correzione adiabatica off-shell (referee Issue 1, main11)
- [x] Termine di dilatazione `−(Ȧ/A)P_r`, Euler completo `D=Θ+P_r∂_Pr`, sorgente `S_D` — verificato vs vero flusso (slope~2), oracolo `test_adiabatic_noreg.py`
- [x] Riduzione livello-B `∫p_r dr = Σa_k U_k + terza specie`, coeff simbolici (M,a,E,J)
- [x] **Wrap off-shell GENERICO TK-t CHIUSO** in funzioni speciali: A(2ª specie→Kleiniane ζ,σ + dilog s=0)
      + B(Δ=0 dilog genus-2) + C(elementare via Hermite). Verificato 1e-14 (2 config) + ancora fisica (slope~2).
      Scripts: `KerrSessionScripts/offshell_tbranch_FULL_assembly.py`, `physics_anchor_offshell_closed.py`
- [x] **Wrap off-shell GENERICO TK-τ CHIUSO** — stessa struttura A+B+C, kernel `A=r²D(EJr−2Ma)/Q3`, azione
      `(−2J²M,−J²(E²−1),−2M,1)`, terza specie solo a Δ=0. Verificato 1e-14 (2 config) + ancora fisica
      (link 7.8e-9 vs sotto-pezzo del flusso). Scripts: `KerrSessionScripts/offshell_taubranch_closed_form_codex.py`,
      `offshell_taubranch_physics_anchor_codex.py` (verificati girare/passare in questa sessione)
- [ ] **FRONTIERA — theta-nome naming** dei dilog genus-2 (`D_{j,root}` a Δ=0, dilog s=0) nella classe
      tabulata Baune/D'Hoker (referee Issue 7). Serve rappresentazione theta + teorema di identificazione.
      Vale per TK-t e TK-τ (entrambi lasciano i letters genus-2 come trascendenti non nominati).
- [x] **Vaidya-τ generico off-shell CHIUSO** (`VaidyaMetric/vaidya_offshell_FULL_assembly.py`,
      `vaidya_offshell_shellpoly.py`). Macchina a polinomio-di-shell (no sqrt hell). Semplificazioni
      analitiche sulla shell pulita: `G=J/(Δp_r)`, `H_pr=(r−2m)DE p_r/(Er²)`, `H_m=C0+C2 p_r²`,
      kernel `A_V=JEr³DE/Q3`, inner `P_inner=−mN4/((r−2m)DE)` — verificati vs numerico 1e-16.
      FULL A+B+C: A(2ª specie W_jk) + B(terza specie a **r=2m orizzonte** e DE=0, dilogs genus-2) +
      C(Hermite elementare, remainder=0 SIMBOLICO). A+B+C == diretto a 1e-16 (config 1) e 1e-15 (config 2).
      **Coeff simbolici in TUTTI (m,E,J)** (a=0 toglie un parametro → all-param passa, a differenza di TK
      E-simbolico). Sorgente `Θ=m∂_m` (no dilatazione, corretto). ρ(r=2m)=−2E²J²m², ρ(DE)=8E²m⁴/…
      **ANCORATO ALLA FISICA**: closed==sub-pezzo off-shell del vero flusso τ a 6.9e-10; totale slope=2.03
      (~2). Il physics anchor ha beccato un bug di segno (Vaidya `m` cresce → `δp_r=(S_D−λΘH)/H_pr`),
      corretto → slope 2.03. Doc: `vaidya_offshell_closure.md`. Anchor: `vaidya_tau_physics_anchor.py`.
- [x] **Vaidya-v generico off-shell CHIUSO** (esplicito, verificato, fisica-ancorato).
      BUG FIXATO (sessione dedicata): avevo derivato `d_pr G` dalla `G` GIÀ on-shell; serve la derivata
      OFF-shell valutata on-shell. Estratto numericamente e verificato:
      * curva `S_v=r·DE·Q2v` (genus-2, `Q2v`=quartica a=0), diversa dal τ;
      * **kernel PURO 2ª specie** `= A_kernel_v/√S_v`, `A_kernel_v=+E²J·DE·r⁴/Q2v` (verificato 1e-16;
        niente parte elementare/D_v^{3/2} — erano artefatti del bug);
      * inner = `elem_inn`(razionale) + 2ª specie (`A_inn·sign(y)/√D_v`), verificato 1e-17;
      * **wrap = block1 + block2**, verificato a 3.5e-18:
        - block1 = kernel×Σ_2nd = **genus-2 A+B+C** (stessa macchina del τ, dilog a r=2m e DE=0);
        - block2 = kernel×Σ_elem = **DOMINANTE** (585%, grande cancellazione).
      block2 via IBP: `∫(rᵏ/√S_v)log(r−2m)dr = Uₖ log(r−2m) − ∫Uₖ/(r−2m)dr` → `∫Uₖ/(r−2m)dr` è una
      **classe weight-2 DISTINTA** dai dilog puri (Abeliano/fattore-lineare, senza √S; tabulata ma diverso
      sotto-tipo). Il v-branch è più ricco del τ (costo "−1" → parte additiva → block2). RIMANE: Hermite
      kernel (Q2v) → g_k; riduzione Σ_2nd e Σ_elem; assemblaggio block1(A+B+C)+block2(IBP → classe Uₖ/(r-2m));
      verifica + physics-anchor. **PHYSICS-ANCHORED**: closed wrap == sub-pezzo off-shell del vero flusso v a
      7.6e-11, totale slope=1.96 (`vaidya_v_physics_anchor.py`). Σ_elem CHIUSO=4m²/(r−2m)−2m·log(r−2m),
      una sola lettera log-Abeliana M^{2m} (coeff 2m). RIMANE (meccanico): block1 A+B+C esplicito
      (P_inner_v ha DOPPIO polo a r=2m) + naming M^{2m}. Scripts: `vaidya_v_FULL_closure.py`,
      `vaidya_v_physics_anchor.py`, `vaidya_v_offshell_structure.py`.
- [~] **Coeff simbolici all-(M,a,J)**: tabelle già simboliche; Hermite `rem_k`/`rho` mostrati E-simbolici;
      inverso modulare all-param = muro perf SymPy → usare Singular / tower QQ(a,E,J)[M]
- [ ] Cosmetico: cancellazione grande A≈−75 vs C≈+73 nella decomposizione di Hermite (decomp più naturale)
- [ ] Paper Tab.1: upgrade "off-shell closed form: open" → "assemblato+ancorato su TUTTE le 4 branch (TK t/τ, Vaidya τ/v); solo theta-naming genus-2 deferred"

## C. Separatrice TK t-branch — tracking counterterm
- [x] Diagnosi: polo di POTENZA `1/(r−r_d)²` on-path (Jc=2.9364, r_d=1.5123 in ergosfera), residuo ∝ ΔS(r_d)=34.75
- [x] Charge tracking `dJc/dE=−0.115` muove il root (`dr_d/dE=+0.0515`) ma NON cancella il polo del kernel
- [x] Non-uniformità LOCALIZZATA a r_d (arco pulito matcha flusso a 99.4% Jc; overshoot satura vicino r_d)
- [ ] **Costruire il controtermine `dr_d/dλ`**: derivare il moto del doppio root lungo il raggio di
      dilatazione (E,J)→(E/A,J/A), verificare se sottrae il polo 34.75; OPPURE analisi di strato limite
      (matched asymptotics inner/outer alla separatrice). Vaidya immune (r_d off-path).

## D. Fisica/classificazione (referee Issues 3,4,5)
- [ ] **MAJOR(3)** — notazione soglie unica: `J_ergo/J_deg/J_sep^phys/J_pen`; nessun root a raggio negativo
      chiamato separatrice fisica; `7.0266` etichettato `J_deg/m` (righe 737, 1829 ancora `\Jc`)
- [ ] **MAJOR(4)** — Vaidya USCENTE: derivare l'Hamiltoniana `u`-branch da `ds²=−f du²−2 du dr`; correggere
      il segno `δ_accr−δ_evap=2qA` (non `2qB`); no "evaporazione" per `m'<0` in chart ingoing
- [ ] **MAJOR(5)** — Lemma B: rendere condizionale fuori dal sotto-regime provato; togliere "never"/"theorem
      in scattering regime"; togliere il ratio puntuale `n_t/n_τ=E/f` dalla tabella come "prova"
- [ ] **MAJOR(7)** — genus-2 Kronecker–Eisenstein: definire la classe o declassare a "integrali abeliani
      iterati depth-2"; rank-5 come evidenza numerica; allineare Mathematica/Sage con VERIFICATION_STATUS.md
- [ ] Nomenclatura orizzonte/ergosfera: TK = oggetto compatto non BH; `Δ=0` = superficie null seme Kerr
      (vedi memoria thakurta-kerr-not-black-hole)

## E. Editoriale/bibliografia/riproducibilità (main11 sec.8-9) — MINOR
- [ ] Rinominare Eq (32)/(38)/(51), titolo App. B (on-shell), caption Fig. 10/14
- [ ] Abstract/Tab.1: split "on-shell closed" / "complete first-order (S_D)" / "off-shell assembly"
- [ ] Bib: [15][31][43][44][45][46][49][54 Nario→Natario][58][59][64 Zenodo DOI] — titoli/DOI/date
- [ ] Overfull boxes (166/106/68 pt); sync PRD; regression runner; environment lockfile
- [ ] Rilascio taggato/Zenodo DOI legato al PDF; coefficienti simbolici machine-readable

---
### Riferimenti
- Referee: `REFEREE_REPORT_CURRENT_2026-07-22.md` (13 issue)
- Nota derivazioni: `DERIVAZIONI_POST_REFERAGGIO_review.md`
- Off-shell/tracking dettaglio: `TODO_offshell_closure_and_tracking.md`
- Propagazione testo: `TODO_main11_text_propagation.md`
- Scripts chiusura: `KerrSessionScripts/offshell_tbranch_*.py`, `physics_anchor_offshell_closed.py`
- Tracking: `ThakurtaMetric/tk_sep_offshell_divergence.py`, `tk_sep_tracking_vs_trueflow.py`
