# TODO MASTER — stato consolidato (agg. 30 lug 2026)

Fonte unica dei punti aperti. Consolida: referee report 2026-07-22 (13 issue), propagazione testo
main11, chiusura off-shell, tracking separatrice TK, e la nota derivazioni_post_referaggio.
Legenda: [x] fatto · [~] parziale · [ ] aperto · **BLOCKER/MAJOR/MINOR** severità referee.

---
## A. Fondamenti control-theory (referee Issues 1,2,6) — ✅ CHIUSO (nota derivazioni + propagato ai paper)
- [x] **Esistenza** minimo assoluto (dominio compatto regolare, indicatrice unif. convessa) — Teorema 5.1, VALIDO
- [x] **Normalità** PMP (ogni minimo regolare = estremale normale) — Prop 3.1, VALIDO
- [x] **Limite stazionario** `H≡0`, recupero energia fissata — Prop 3.2, VALIDO
- [x] **Certificato globale HJB** (sub-soluzione Lipschitz ⟹ ottimalità) — Teorema 7.1, VALIDO
- [x] **Framework coniugati/Maxwell/cut locus** — §8, VALIDO
- [x] **MAJOR — Equivalenza esplicita con Perlick DIMOSTRATA** (`F_rail=F_Perlick`): la Finsler tempo-
      d'arrivo È la metrica di Randers di Perlick a energia fissata, `β_a=ω_a`, `a_ab=Ê²h_ab/(f(Ê²−f))`;
      limite null Ê→∞ = Fermat. Verificato 0–4.4e-16. **PROPAGATO al paper** (main.tex §W-hierarchy(a) + PRD, eq:perlick-randers, compilano).
      Scripts: `perlick_equivalence.py/.md`.
- [x] **Issue 2 — `p_φ` costato vs momento meccanico `L_mech`** — CHIUSO. Derivazione (commit ad3a372,
      `costate_finsler_derivation.py`) + testo §4.1 COMPLETO in entrambi i paper: paragrafo "Branch costate
      vs mechanical angular momentum" (main line 586-607, PRD 549-555): `J=p_φ` = costato di Noether NON a
      priori meccanico `g(u,∂_φ)`; coincidono qui perché `∂_φ` Killing di metrica+problema di controllo;
      Legendre branch-specific; `J_eff=J/A` = rescaling ottico `P_φ=p_φ/A`, non riscalatura del costato.
      Ancorato alla Finsler stazionaria (eq:perlick-randers). NON invalida i risultati.
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
- [x] Integrare il "Testo proposto" (§10 nota) come blocco control-theory nel paper (main.tex + PRD).
      FATTO — nuova subsection "Existence, normality, and global optimality" (sec:control-foundations)
      dopo indicatrice/Pontryagin in ENTRAMBI: Esistenza (Tonelli, bound ρ|v|≤F≤R|v|, coercività+
      Arzelà-Ascoli+semicont. convessa, s-coupling Grönwall) + Normalità (p₀=0 ⟹ contraddizione
      h(p)≥ρ|p|>0, costato mai nullo) + HJB (subsoluzione Lipschitz saturata ⟹ ottimo globale) +
      Conjugate/Maxwell/cut-locus. main 71→72pp, PRD 33pp, compila pulito. Commit 6813fbe.
- [x] **MAJOR(Issue 6) — esclusioni**: minimalità dentro l'ergosfera (selettore spacelike, continuazione
      analitica ≠ ottimo fisico); PMP necessario ≠ esistenza/minimalità globale; freezing≠turning≠
      separatrice. FATTO — paragrafo "What is not claimed inside the ergosphere" nel blocco control-theory
      (sec:control-foundations, main+PRD): ξ spacelike ⟹ φ(r) continuato = continuazione analitica NON
      minimo fisico certificato; esistenza/HJB asserite solo fuori ergosfera; le tre degenerazioni distinte.
      Commit 6813fbe.

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
      STATO (sessione): abelfunctions OK (`~/.sage/local`, usare sage di /Applications); convergenza
      theta-nome VERIFICATA (N=4→1e-13, `kerr_tau_dilog_qseries5.sage`); MA naming Kleiniano NON chiude
      (residui ~1e-3, ostruzione theta-divisor Θ=W_{g-1}, `kerr_tau_Wij_diffform.sage`). Teorema completo
      = ricerca dedicata / collaborazione con esperti (framework higher-genus esso stesso frontiera aperta).
- [ ] **theta-naming: MATCHING STRUTTURALE NUMERICO** (via di mezzo FATTIBILE): esprimere le lettere
      W_jk/dilog nella forma D'Hoker–Schlotterer per la nostra curva, verificare numericamente il matching
      struttura-per-struttura (abelfunctions), dichiarare "structurally matched, single-valued completion
      deferred" — alza la claim da "conjectural" a "structurally matched". NON è il teorema completo.
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
- [x] **Controtermine `dr_d/dλ` — RISOLTO (natura identificata)**: NON esiste. Sotto l'evoluzione il raggio
      dilatazione lascia la separatrice → il doppio root si DISSOLVE, l'orbita ATTRAVERSA la separatrice.
      È SEPARATRIX CROSSING (Neishtadt), non un controtermine algebrico. `tk_separatrix_crossing.md`,
      `tk_sep_tracking_counterterm.py`, propagato ad App C (main.tex+PRD, +ref Neishtadt). Commit c9bd978.
      [ex] derivare il moto del doppio root lungo il raggio di
      dilatazione (E,J)→(E/A,J/A), verificare se sottrae il polo 34.75; OPPURE analisi di strato limite
      (matched asymptotics inner/outer alla separatrice). Vaidya immune (r_d off-path).

## D. Fisica/classificazione (referee Issues 3,4,5,7) — ✅ CHIUSO (propagato ai paper; frontiere deferite esplicite)
- [x] **MAJOR(3)** — notazione soglie: tabella autoritativa (`separatrix_notation_table.md`,
      `separatrix_classification.py`) + PROPAGATO. `J_deg=7.0266, r_d=−3.3637` (radice NEGATIVA =
      degenerazione algebrica) etichettato J_deg ovunque in main (850,1948) + PRD (763), disambiguato dalla
      soglia dinamica `J_c(v₀)` ("must not be conflated"), interpretazione capture/escape rimossa. Nessun
      `\Jc=7`/"separatrice 7.0266" residuo. `\Jc` restante = solo la marginale τ `a/E` (significato unico).
- [x] **MAJOR(4)** — Vaidya USCENTE: segno FIXATO + PROPAGATO (main+PRD): `2A_∞=17.455` (differenza fisica),
      `2B_hor=12.225` (deviazione antisimmetria); clock uscente DERIVATO (`du/dr=dv/dr−2dr_*/dr`, 2.2e-16);
      nomenclatura accr=ingoing/evap=outgoing. `vaidya_asymmetry_fix.md`, commit aa56992. Inversione
      evaporativa esplicitamente SCOPED al livello frozen-clock; BVP outgoing completo + scan fisico
      no-inversion outgoing DEFERITI in-testo (§685-712 main / 643-665 PRD), non rivendicati (onesto).
- [x] **MAJOR(5)** — Lemma B: PROPAGATO (`issue5_lemmaB_status.md`). Lemma A provato in chiuso;
      Lemma B **condizionale/tiered** (chiuso per grazing/quarter `r₀≤R_*(E)` + asintotico statico +
      interval-arithmetic certificato a r₀=10M + congettura in generale, aperta a `r_pk` trascendente);
      rapporto `n_t/n_τ=E/f` "necessary but not sufficient" (no più "prova"); no "never" non qualificato.
      main 1792-1837 / PRD 1538-1590; Table protocollo con footnote necessary-not-sufficient.
- [x] **MAJOR(7)** — genus-2: DECLASSATO a "depth-two iterated Abelian integrals" + "hyperelliptic-dilog
      class conjecturally irreducible"; rank-5 = evidenza numerica; Table 1 tiered proved/machine-certified/
      conjectural; higher-genus polylog framework "deferred". main (68,75,841,159) + PRD (82,1107,1206).
      VERIFICATION_STATUS.md allineato (SYM/NUM/SAGE + WL Mathematica cross-check indipendente).
- [x] Nomenclatura orizzonte/ergosfera: TK = "conformal rotating compact object" (§977 main / 819 PRD);
      `Δ=0` = "seed Kerr null surface" (non orizzonte); linguaggio BH/orizzonte solo per Vaidya.
      Coerente con memoria thakurta-kerr-not-black-hole.

## E. Editoriale/bibliografia/riproducibilità (main11 sec.8-9) — MINOR — ✅ CHIUSO (salvo Zenodo DOI, blocco esterno)
- [x] Rinominare Eq (32)/(38)/(51), titolo App. B (on-shell), caption Fig. 10/14 — AUDIT COMPLETO.
      Linguaggio on-shell/off-shell pervasivo e CORRETTO in eq/caption/titoli/teorema in ENTRAMBI i paper
      (main.tex + PRD, 38 occorrenze PRD). Ogni claim di forma chiusa qualificata come "on-shell component",
      con "complete first-order term = extended-Hamiltonian (off-shell)" esplicito: teorema unificato
      (main 1820 / PRD 1571), Vaidya subsec (main 655 / PRD 627), Table 1 (line 150), abstract (55-60),
      App validation (PRD 2194 "complete first-order (on-shell+off-shell)"). Nessun overclaim residuo.
- [x] Abstract/Tab.1: split "on-shell closed" / "complete first-order (S_D)" / "off-shell assembly" —
      FATTO. Abstract main.tex line 55-60 ("on-shell component reduces in closed form ... off-shell costate
      piece, of the same class"); Table 1 line 150 "On-shell first-order correction"; propagato a PRD.
- [x] Bib: [54 Nario→Natario]✓ + DLMF 1.2.7✓ + SageMath 2026-05-04✓ + Neishtadt1987/2013✓ + repo SHA✓
      (commit bbdc385). AUDIT: ZERO citazioni/ref indefinite in entrambi i compile; igiene genus-1/genus-2
      CORRETTA (separatrice=ellittica usa Brown-Levin/Beilinson-Levin propriamente; disclaimer esplicito
      main line 1260 "that reference does not by itself establish a higher-genus construction"); entries
      complete (autori/pagine/DOI dove disponibili). UNICO aperto → Zenodo DOI (bloccato esterno, sotto).
- [x] Overfull boxes: main.tex script-map table 165pt→6.7pt (p{} cols); PRD 92.7pt→8.1pt (split eqs
      t-Q2/unified/radial-action/vaidya-deltaphi), 38/30/17pt eliminati. Commit 5c0120e.
- [x] Regression runner + environment lockfile: `run_regression.py` (12/12 pass, residui 0–1.2e-12,
      slope 2.00–2.21; --quick, -k) + `requirements-lock.txt` (stack pinnato + nota Sage/abelfunctions).
      Commit fb45761.
- [x] Rilascio taggato/Zenodo DOI — FATTO. `10.5281/zenodo.21707378` (release v1.0) nel bib
      (`Rosignoli2026Code`) + Data availability di entrambi i paper. Commit 42aebb3. **Ultimo blocco chiuso.**

## F. Follow-up referee main12 — ✅ CHIUSO (8 azioni)
Dettaglio in `REFEREE_FOLLOWUP_main12.md`. Verdetto: "responded successfully to the principal
technical referee criticism". Tutte 8 le azioni indirizzate (commit bc2b85c, 7c57bc8, b225b9e, 7daa912):
- [x] (1) Contraddizione abstract/Table 1 su assembly off-shell → opzione 1 (assemblato 4 branch, solo naming deferito)
- [x] (2) Validità solo su sub-archi `H_Pr≠0` (già presente) + (6) istantanea vs Neishtadt crossing (già presente)
- [x] (3) Convenzioni di curva (sheet/cuts/base point/orientamento/costanti) — paragrafo App. main+PRD
- [x] (4) PDF-solo vs repository (reconstructibility) — paragrafo main+PRD
- [x] (5) Fixed-charge vs fixed-endpoint + `J_1=−δφ_dyn/∂_J φ_0` (eq:J1-fixedendpoint)
- [x] (7) η = rappresentazione ramo di arrivo (non 3ª famiglia) + (8) Vaidya clock v,τ canonici; t frozen; u separato

## G. Follow-up referee main13 — obiezioni ✅ già coperte; estensioni = futuro
Dettaglio in `REFEREE_FOLLOWUP_main12.md` (stesso lavoro) + piano in `PLAN_split_and_extensions.md`.
Verdetto: "major conceptual objections addressed; remaining = rigor/consistency/extension, not fundamental flaws".
- [x] 4 "remaining objections" (regular sub-arcs/turning points; fixed-launch vs fixed-endpoint; off-shell
      consistency; explicit quadratures vs pipeline) — TUTTE già chiuse dal follow-up main12 (report scritto
      contro versione precedente ai fix).
- [ ] **ESTENSIONI (lavoro futuro, NON richiesto per submission corrente)**:
  - [x] Forma **generale/universale** della sorgente off-shell `S_D` — **TROVATA+VERIFICATA+PROPAGATA AL PAPER**
        (chiude estensione 3 + critica GPT). Identità universale di Finsler `(J∂_J+P_r∂_Pr)H=H+1` (provata per
        β,A_rr,A_φφ ARBITRARIE). DUE forme universali per tipo di deformazione (non metrica-specifica):
        * **Conforme (TK)**: `S_D = λ + ∫E_eff H_{E_eff} dλ` (λ universale, E H_E unico pezzo metrica-specifico).
        * **Mass-function (Vaidya)**: `S_D = [r p_r] − λ` (termine di bordo, via auto-similarità `m H_m=−(rH_r+JH_J)`).
        * Scheletro comune ±λ. Simbolico esatto (anche astratto) + numerico macchina su TK e Schwarzschild.
        Doc `UNIVERSAL_SD_source.md`; script `ThakurtaMetric/universal_SD_source_check.py`,
        `VaidyaMetric/universal_SD_source_vaidya.py`; commit 39afd74/7618675. Propagato: paragrafo
        "Universal form of the off-shell source" (main eq:finsler-euler/SD-conformal/SD-vaidya + PRD), commit 06bc3ca.
  - [x] **η-brachistocrona in TK** — FATTA+PROPAGATA (chiude estensioni 4 E 5). `F_t=A F_η` (frozen) +
        monotonia (`t_f=∫A dη` crescente in η_f) ⟹ frozen/adiabatica/separatrice η ≡ t (stessa curva,
        carica rietichettata `J_t=A J_η`, verificato 7e-9); DUE famiglie ottimali (arrivo t≡η, proprio τ),
        non tre; η = gauge conforme naturale (H_η costo unitario, `S_D=η+∫E H_E dη`). Doc
        `eta_brachistochrone.md`, script `ThakurtaMetric/eta_brachistochrone_check.py`. Paragrafo "The
        conformal-time (η) brachistochrone" main+PRD, commit 4aa16ea.
  - [ ] **Vaidya terzo clock** SOLO se legato a osservatori privilegiati (Kodama), non per analogia:
        esiste un tempo proprio di Kodama-osservatore con brachistocrona distinta da v e τ? [estensione 6, condizionale]
  - [x] Framework unificato t/τ/η-brachistocrone [estensione 5] — chiuso insieme all'η-brachistocrona:
        il framework unificato È "due famiglie (arrivo t≡η, proprio τ)". Vedi sopra, commit 4aa16ea.
  - [ ] Teoria generale dei **turning point in movimento** [estensione 1] + delle **separatrici adiabatiche**
        e correzioni first-order [estensione 2] (research-heavy)
- [ ] **DECISIONE STRATEGICA — split in 2 paper** (FLRW-Vaidya sferico | TK assisimmetrico): confine fisico
      pulito (frame dragging). Raccomandazione: prima pubblicare il combinato (verdetto positivo), POI split
      sulle estensioni. Alternativa: split subito se il combinato (74pp) è troppo lungo per il journal.
      Piano dettagliato: `PLAN_split_and_extensions.md`.

## H. Technical Revision Memorandum main14 (GPT) — change-list completa in `REFEREE_MEMO_main14.md`
Verdetto: core pubblicabile, major revision, NON splittare finché i fondamentali non sono corretti nel
master. Ordine: P0 fondamentali → P1 armonizzazione → P2 editoriale → split in 2 paper.
- [x] **P0 (4) — TUTTI FATTI**: 4.1 problema terminale formale (2ae215f); 4.2 HJB non-autonoma/verifica
      sufficiente (22f994a); 4.3 costato≠meccanico, `L_mech=((Ê²−f)/Ê)J` verificato (5c07823); 4.4 App C
      ri-derivata da variabili canoniche, `-αP_r` da regola catena, verificato (3102013). Dettaglio in
      `REFEREE_MEMO_main14.md`.
- [x] **P1 (5) — TUTTI FATTI** (commit 6ac905f): 4.5 R_*(E,a) esplicito; 4.6 tracking→derivata famiglia;
      4.7 J_deg vs separatrice (App retitolata); 4.8 evaporazione→negative-rate ingoing; 4.9 genus-2 peso-due.
- [~] **P2 editoriale**: TESTO FATTO (8b3e421, 15052b2) — tabella notazione, conclusione 3 blocchi,
      "App. Appendix" duplicati, global attractor, J_c^+(A) numerata, Fig 9 forward-ref. RESTA (fase split):
      rigenerare IMMAGINI figure (Fig 10 label/slope, Fig 14 due residui, assi, scala pannelli), spostare
      thrust/fuel in remark, rimuovere §5.1, Fig 20 cross-ref.
- [x] **Review main15 (GPT)** — ricalca main14 (GPT ha revisionato la versione post-main14 ri-segnalando
      molto di già fatto). AUDIT del master: 4 residui testuali reali che main14 aveva mancato, sistemati
      (commit 87de08c): (A) App C "Extended system" scriveva ancora `H=p_s+H̄(r,P_r)` come canonica →
      corretto a (r,p_r); (B) conclusione aveva Open problems + Future directions ridondanti → un solo
      blocco Outlook (Proved/Conditional/Outlook); (C) Vaidya sferico `δφ_sep`→`δφ_deg` + titolo
      "degeneration-family"; (D) thrust/fuel "photon rocket 45%" → remark compatto. Resto main15 = figure
      da rigenerare + §5.1 + distribuzione split → fase split (sotto).
## I. SPLIT IN DUE PAPER — piano d'esecuzione (da main15 §5.3)
Prerequisito: master coerente ✅ (P0/P1/P2-testo + figure P0/P1 fatti). Le figure P0/P1 sono già
rigenerate/etichettate correttamente (commit cb0bbc7/d076374/ff05ede) e valgono per entrambi i paper.

### I.0 Preparazione — ✅ FATTO
- [x] Creati `paper1/paper1.tex` (rail/Kodama/Vaidya) e `paper2/paper2.tex` (Thakurta-Kerr): scheletri con
      preambolo iopart, `\graphicspath{{../paper/Immagini/}}` (figure condivise), `\bibliography{../paper/refs}`
      (bib condivisa), classe iopart copiata localmente, placeholder di sezione con `\todo` di migrazione.
      Entrambi compilano (5pp, bibliografia completa via `\nocite{*}` finché non ci sono `\cite` reali).
- [x] Change-log/mappa sezioni: `SPLIT_MAP.md` (quale sezione del master va in quale paper).

### I.1 Paper I — "Controlled-rail brachistochrones..." (sferico) — ✅ CORPO+App C FATTI
- [x] Corpo migrato (commit bb72d00): intro sferica (rinvia TK al companion); fondamenta condivise
      (controlled-rail, esistenza/normalità/HJB non-auton., gerarchia W); FLRW; Vaidya (Kodama energy,
      fenomenologia, plunge law + assenza inversione evaporativa, adiabatica+dilog); 3 figure Vaidya;
      conclusioni 3 blocchi (Vaidya scope, BVP outgoing aperta). Cross-ref TK → "companion paper".
- [x] App C migrata (commit 6310ce6): first-order extended-Hamiltonian correction (eq:adiab-exact
      self-contained, derivazione canonica eq:normflow, sorgente universale, branch Vaidya). Ref TK ammorbiditi.
- [x] **Paper I compila a 33pp con ZERO undefined references.** `paper1/paper1.tex`.
- [x] App A ("Reproducibility and consistency") + App B ("genus-degeneration loci", Vaidya J_deg)
      migrate (commit 261709e): App A residui riscritti a scope Vaidya (Kodama energy, limite a=0,
      minimizzazione branch), colormap TK + figure Vaidya duplicate rimosse; App B Vaidya J_deg, branch TK
      condensate a rimando companion.
- [x] Abstract reale **227 parole** (≤300 CQG); Data availability (DOI Zenodo) + ack riempiti; zero \todo.
- [x] **Paper I COMPLETO: 38pp, zero undefined/multiply refs, self-contained.** `paper1/paper1.tex`.

### I.2 Paper II — "Brachistochrones in conformal Kerr spacetimes: ergosphere trichotomy, separatrices,
      and adiabatic response" (assisimmetrico, ~30-40pp)
- [ ] Contenuti: indicatrice TK + equivalenza t≡η (arrival) + Fig 10; ramo arrivo vs τ + dizionario
      costato/meccanico; ergosfera cuspide/tricotomia + diagramma penetrazione; forme Weierstrass/Kleiniane
      frozen + continuazione Doran; inversione + protocol dependence (same-launch vs fixed-endpoint tiered,
      R_*(E,a)); first-order on-shell/off-shell CANONICAMENTE ri-derivato + true-flow validation;
      rappresentazione iterati Abeliani lunghezza-2 con terminologia cauta (peso-1/peso-2).
- [ ] Abstract dedicato ≤300.

### I.3 Interventi strutturali durante lo split (main15 §5.1-5.2, tabella §4)
- [ ] **Rimuovere/isolare §5.1** (quasi-costanti off-equatoriali, theta/p_theta/Carter, sorgente O(a²)):
      NON trasferirla in Paper I; conservare gli script nel repo come materiale per articolo futuro 3D.
- [x] Ridurre remark thrust/fuel (già fatto nel master, commit 87de08c) — verificare che passi nei due paper.
- [ ] Fondamenta condivise (rail formalism, gerarchia W, esistenza/normalità/HJB, equivalenza Perlick):
      collocate in Paper I come framework generale, richiamate da Paper II.

### I.4 Polish figure/caption residuo (durante lo split)
- [ ] Caption A6-A10: "numerical evidence / robustness scan" non "proof"; A9 specificare norma residuo
      (assoluta/relativa); A10 = local variational check, non prova di minimo globale.
- [ ] Fig 20: caption "finite numerical scan" (già ok via label sec:inversion) — verificare dopo renumber.
- [ ] Scala pannelli galleria appendice troppo piccoli: split/ingrandire se finiscono nei paper.

### I.5 Finalizzazione
- [ ] Rigenerare TUTTI i cross-reference dopo la separazione (numeri eq/fig/sez cambiano).
- [ ] Verificare zero undefined ref/citation, zero overfull, in entrambi i paper separati.
- [ ] Aggiornare Data availability + DOI Zenodo in entrambi (già presente nel master).
- [ ] Recuperare eventuali risultati riusati tra i due paper con citazione incrociata (Paper I ↔ Paper II).
- Nota: 4.2/4.3/4.4 = lacune reali in pezzi aggiunti da me (rigore, non invalidano i numerici).

---
### Riferimenti
- Referee: `REFEREE_REPORT_CURRENT_2026-07-22.md` (13 issue)
- Nota derivazioni: `DERIVAZIONI_POST_REFERAGGIO_review.md`
- Off-shell/tracking dettaglio: `TODO_offshell_closure_and_tracking.md`
- Propagazione testo: `TODO_main11_text_propagation.md`
- Scripts chiusura: `KerrSessionScripts/offshell_tbranch_*.py`, `physics_anchor_offshell_closed.py`
- Tracking: `ThakurtaMetric/tk_sep_offshell_divergence.py`, `tk_sep_tracking_vs_trueflow.py`
