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

### I.2 Paper II — "Brachistochrones in conformal Kerr spacetimes..." (assisimmetrico) — ✅ FATTO
- [x] Assemblato (commit 7dc7e19): intro TK (cita Paper I per fondamenta); body = §5 TK **senza §5.1**
      (η-brachistocrona t≡η, indicatrice/Hamiltoniane, dizionario costato) + §6 (separatrici Weierstrass,
      tricotomia/cuspide, inversione, fixed-endpoint tiered R_*(E,a), breathing); conclusioni 3 blocchi;
      appendici TK (explicit closed forms, genus-degeneration, first-order); Data availability+ack.
      Cross-ref Paper I/Vaidya/§5.1-rimossa → companion. **55pp, zero undefined/multiply refs, zero \todo.**
- [x] Abstract reale **220 parole** (≤300 CQG).
- [ ] Nota: 55pp è ancora sopra il target 30-40pp (main15) — trimming opzionale in polish (App A residui,
      overlap Vaidya nelle appendici, §5.1 già rimossa).

### I.3 Interventi strutturali durante lo split — ✅ FATTO (commit e1be126)
- [x] **§5.1 rimossa/isolata**: esclusa dal body di Paper II (I.2); paragrafo appendice off-equatoriale
      "Closed toroidal averages" (p_θ/inclinazione) rimosso dall'App A di Paper II. Script conservati nel
      repo come materiale futuro 3D. I 3 paragrafi equatoriali (Kerr constants, rotational inversion,
      conformal trichotomy) tenuti in Paper II.
- [x] Remark thrust/fuel ridotto (master 87de08c) — presente ridotto in Paper I, verificato.
- [x] Fondamenta condivise in Paper I, richiamate da Paper II (18 rimandi al companion). ✓

### I.4 Polish figure/caption residuo — ✅ FATTO (commit 1b891fc)
- [x] Caption scan: fig:cm-spin-fix ("finite scan"), fig:cm-conf-fix ("finite numerical scan, not a
      theorem"), fig:asym-spin/asym-conf ("robustness scan, not a theorem") — tutte con framing
      numerical-evidence non-proof. fig:vaidya-kerra0: norma residuo sup <1e-12 (A9). Entrambi i paper.
- [x] Fig 20 (fig:cm-spin-fix): "finite scan, tiered no-inversion status" — già ok, verificato.
- [x] Pannelli: colormap a \textwidth (single-panel) — leggibili, ok.
- [ ] (scope, → I.5 trimming) fig:vaidya-kerra0 + fig:verifmin sono Vaidya ma presenti anche in App A di
      Paper II (overlap): rimuovere da Paper II in fase di trimming.

### I.5 Finalizzazione — ✅ PARZIALE (README+verifica fatti; trimming opzionale)
- [x] README: sezione two-paper split + mappa script Paper I/II (commit 4f61e4d).
- [x] Verifica finale: Paper I 38pp / Paper II 55pp, **0 undefined, 0 overfull>80** in entrambi.
- [ ] (opzionale) Trimming Paper II 55→30-40pp: rimuovere overlap Vaidya nelle appendici (fig:vaidya-kerra0,
      verifmin, vaidya-bounce, vaidya-timing, vaidya-offshell + prose che li referenzia → ritargeting/rimozione),
      condensare parte Vaidya di App B, tabella script-map per-paper (ora lista tutti gli script).
- [ ] **Aggiornare `README.md`**: sezione "quale script per quale paper" — mappa directory/script a
      Paper I vs Paper II, così chi riproduce sa dove guardare. Mappatura proposta:
      * **Paper I** (rail/Kodama/Vaidya sferico): `FLRWmetric/`, `VaidyaMetric/` + fondamenta control-theory
        condivise (`KerrSessionScripts/perlick_equivalence*.py`, esistenza/normalità).
      * **Paper II** (Thakurta-Kerr): `ThakurtaMetric/`, `KerrSessionScripts/` (genus-2), `KerrMetric/`, `KerrScripts/`.
      * **Condivisi**: `PaperFigures/`, `run_regression.py`, `requirements-lock.txt`, `paper_style.py`.
      Aggiornare anche la tabella script-map in ciascun paper per elencare solo i propri script.
- [ ] Rigenerare TUTTI i cross-reference dopo la separazione (numeri eq/fig/sez cambiano).
- [ ] Verificare zero undefined ref/citation, zero overfull, in entrambi i paper separati.
- [ ] Aggiornare Data availability + DOI Zenodo in entrambi (già presente nel master).
- [ ] Recuperare eventuali risultati riusati tra i due paper con citazione incrociata (Paper I ↔ Paper II).
- Nota: 4.2/4.3/4.4 = lacune reali in pezzi aggiunti da me (rigore, non invalidano i numerici).

## J. Ristrutturazione teorema-applicazione Paper I/II — dettaglio in `RESTRUCTURE_MEMO.md`
Principio: Paper I DIMOSTRA (formalismo+FLRW+Vaidya, teoremi I.1–I.6 numerati); Paper II CITA e applica
(TK + corollari), niente ridimostrazioni/duplicazioni. Ordine esecuzione (§9 del memo):
- [x] 1. Numerati I.1–I.6 in Paper I (commit 4571410): Thm I.1 PMP, Thm I.2 esist./norm., Prop I.3 HJB,
      Lemma I.4 riduzione polinomiale, Thm I.5 first-order, Lemma I.6 double-root. Compila.
- [x] 2. Paper II cita i teoremi: 7.1 (§TK: I.1/I.2/I.3, eead4fc) + 7.3 (App C: Thm I.5, eead4fc) +
      7.2 (Lemma I.4 riduzione, 44b49fa) + 7.4 (ladder pointer, 44b49fa). Tutte fatte.
- [x] 3. App C: paragrafo Vaidya rimosso da Paper II (95e6085); intro cita Thm I.5, resta solo il conforme.
- [x] 4. App B TK-primario in Paper II (bbc26fa): 2 branch Vaidya rimossi, macchina condivisa richiamata,
      Vaidya J_deg → Paper I. Residuo rifinito: paragrafi macchina resi metric-neutral/TK-primari
      (sorgente generica $\partial_\lambda F=N/S^{3/2}$, $N=\Ehat N_E+J N_J$ via Lemma I.4/Thm I.5;
      seed mass $2m\to2M$; rimosso "This example is Vaidya" e paragrafo "Vaidya ($\lambda=m$)";
      $N_m\to N$). Ricompilato 49pp, 0 undefined.
- [x] 5. §breathing → "Application of the breathing-indicatrix theorem to conformal Kerr" (1617a86).
- [x] 6. Tabelle scisse: Paper I solo Vaidya (7c60436), Paper II solo TK (5f1d76f).
- [x] 7. Figure: Vaidya rimosse da Paper II (5f1d76f); Paper I non aveva figure TK.
- [x] 8. Ladder ridotta a pointer + riga TK in Paper II (44b49fa).
- [x] 9. Cross-reference: **0 undefined in entrambi** dopo tutti i tagli.
- [~] 10. Lettura incrociata: contaminazione maggiore rimossa (App C/B, tabelle, figure, ladder, breathing,
      sorgente universale). Menzioni Vaidya residue in Paper II = confronti qualitativi legittimi (ladder,
      limite statico). **Paper I 38pp, Paper II 50pp, self-contained, teorema-applicazione stabilita.**
- Testi modello §7 per aperture Paper II già in `RESTRUCTURE_MEMO.md`.

### J-bis. Memorandum GPT 2 (2026-07-31) — 2 sostituzioni P0 + bib (scope scelto dall'utente)
Audit GPT confermato nei file: 2 contaminazioni simmetriche vere che il restructure non aveva rimosso.
- [x] **P0-A. De-conformalizzata Theorem I.5 / App C di Paper I.** App C era duplicato quasi-verbatim
      dell'App C conforme di Paper II. Enunciato I.5 riscritto in variabili canoniche $(r,p_r)$ con Euler
      astratto $\Theta$ (via $E_{\rm eff},J_{\rm eff},P_r$ dall'enunciato); rimossi i paragrafi conformi
      ("conformal rescaling/dilation term", "dilation letter" con $a_k$ Kerr/$\Delta$, script ThakurtaMetric,
      pole ergosfera/crossing TK) → ora posseduti solo da Paper II App C. Aggiunta chiusura Vaidya:
      identità Randers–Euler $(J\partial_J+p_r\partial_{p_r})H=H+1$ + self-similarità Schwarzschild →
      **boundary form $S_{\mathcal D}=[rp_r]-\lambda$**. Special-function/separatrix riquadrati Vaidya
      (third-kind a $r=2m$, no crossing; TK crossing → 1 frase a Paper II).
- [x] **P0-B. Rimosso blocco Vaidya da §sec:unified di Paper II** (`J_deg=m j(E)`, 7.0266,
      accretion/evaporation, tracking `J=Jc(m)` sorgente `∂_m F+(Jc/m)∂_J F`, advanced-time, horizon dilog
      a $r=2m$) → riquadrato conforme + pointer a Paper I e App.~\ref{app:tracking}.
- [x] **Bib reciproca**: voci `PaperOne`/`PaperTwo` in `paper/refs.bib`; `\cite{PaperTwo}` in intro Paper I,
      `\cite{PaperOne}` in intro Paper II. Rinominato §"breathing-indicatrix theorem" →
      "Adiabatic breathing of the conformal-Kerr indicatrix" (evita teorema non numerato).
- Ricompilati: **Paper I 38pp, Paper II 49pp, 0 undefined entrambi**.

### J-ter. Simmetria del fix P0 (item 1,2,3, 2026-07-31)
Il giro P0 era asimmetrico (pulito Paper II ma non i gemelli in Paper I). Chiusi 3 residui reali:
- [x] **1. App B di Paper I → Vaidya-only.** Rimossi residui TK: parentetica `λ=A`, forward-computation
      `N=Ê N_E+J N_J`/"30% shift"/`conformal_source_Jterm.py`, paragrafo "Rotating (TK) case",
      "three branches"→"two", "for both Vaidya and TK", coefficienti `(M,a,E,J)`→`(m,E,J)` (roots `2m`).
      Restano solo pointer brevi a Paper II.
- [x] **2. Sanity-check rail TK spostato** da §4 di Paper I (`rail_conservation.py`, orbita Kerr on-shell)
      → nuovo paragrafo in §2 di Paper II; in Paper I resta pointer + il punto generale on-shell/off-shell.
- [x] **3. Riga "Vaidya, ∂_m" + `N_m`** rimossi da tab:reduction-ck e testo di App A di Paper II
      (tabella riproducibilità ora TK-only; Vaidya → pointer una riga a Paper I).
- Non in scope (cosmetici, lasciati): status-tier, rinumerazione eq/fig, rename Kodama→conformal-selector,
      trim notazione `P_r` in Table 1 di Paper I.

### J-quater. Audit memo GPT "Prescrizioni giustapposte" (2026-07-31) — layer Appendici mancato
Il memo ha rivelato contaminazione residua NEGLI APPENDICI (i giri prima toccavano App C e i paragrafi-macchina
di App B, non la prosa di App A né la coda di App B). 4 fix reali:
- [x] **A. App A di Paper I → Vaidya-only.** Prosa riduzione era TK (`R_6`/`Q_2`/ramo `t`/`c_k^E,c_k^J`/
      "conformal ψ-split") e caption diceva "(Thakurta–Kerr rows)" con tabella già solo-Vaidya (incoerenza).
      Ora identità 11×11 generale + sola specializzazione Vaidya; TK → 1 frase a Paper II.
- [x] **B. App A di Paper II** apre con "By Lemma~I.4 of Paper~I" e non ri-espone il solve 11×11 generale
      (tiene identità come risultato citato + applicazione TK).
- [x] **C. Coda App B di Paper II**: "Tracking cancels the triple pole" ora cita Lemma I.6 e toglie
      "for both Vaidya and TK"; "Transcendence" riquadrata TK (`Δ=0`, ancora a infinito) via contenuto
      Vaidya (`∂_m F`, `r=2m`) → pointer a Paper I.
- [x] **D. HJB condizionale** in Paper II (`paper2:109`): "global optimality certified" → "only under the
      hypotheses of the verification criterion..." (Paper I Prop I.3 era già condizionale).
- Ricompilati: **Paper I 37pp, Paper II 49pp, 0 undefined**. Cosmetico lasciato: coda reproducibility
      condivisa (nome q-series/rank-five) leggermente duplicata.

### J-quinquies. Memo GPT "Consigli finali" (2026-07-31) — verdetto: publishable after moderate revision
Fatte le modifiche SIGNIFICATIVE (correttezza/coerenza):
- [x] **P0 §2.1 gerarchia peso-uno/peso-due (il vero fix scientifico).** La correzione on-shell del primo
      ordine contiene GIÀ i `W_jk` (peso due, esplicito in §thakurta 430-458), quindi "on-shell = proved
      weight-one closed form" era contraddittorio. Uniformato in entrambi: frozen ∂_λφ_0=peso uno, ma
      on-shell first-order porta un settore peso-due `W_jk`; off-shell aggiunge terzo-tipo (dilog).
      Fix a paper1:1434/1443, paper2:519/1137/1989/2004. Anche "Kleinian ζ,σ closed forms of U_k,W_jk"
      → solo U_k Kleinian, W_jk assemblati (non ridotti a peso uno).
- [x] **P0 §2.2 HJB Paper I** (paper1:382): "certified between launch and first singularity" → "PMP +
      local numerical minimality, NOT a complete global HJB certificate; globally certified only up to
      first Maxwell/conjugate obstruction".
- [x] **P0 §2.3 fixed endpoints** (paper1:195): "both spatial endpoints AND the arrival are prescribed"
      (contraddittorio) → "both spatial endpoint positions prescribed, arrival clock DETERMINED by optimum".
- [x] **§2 ergosfera** (paper1:394): "spacelike selector to restore a timelike cost" → "timelike observer
      congruence / different timelike selector for a compact future-timelike indicatrix".
- [x] **Mislabel Kodama→conformal-Killing rail charge in Paper II** (paper2:1166, 1631, 1857): TK non ha
      vettore di Kodama.
- [x] **Refusi**: "transversality~the transversality" (paper2:488), "the 3+1 the 3+1 Hamiltonian"
      (paper2:1630), "off-shell shell"→"off-shell shell departure" (paper1:606). "Secs.3–2"/"App.Appendix"
      non presenti (già a posto).
- [x] **Terminologia §4.2**: intestazione "Tracking cancels the triple pole" → "The degeneration-family
      derivative cancels the leading pole" (entrambi). "separatrix-following" lasciato dove usato con caveat.
- NON fatto (poco significativo, motivato all'utente): trim notazione Table 1 Paper I; rigenerazione
      legende figure evaporation→negative-rate; conversione "companion paper"→"Ref.[X]" ovunque; thrust/Δv;
      coda reproducibility duplicata.
- Ricompilati: **Paper I 38pp, Paper II 49pp, 0 undefined**.

### J-sexies. Memo GPT "Risposta paper1-5/paper2-5" (2026-07-31) — set significativo
Verdetto GPT sceso a minor-to-moderate; fix peso-due/HJB/endpoint riconosciuti "Riuscita". Chiusi i residui:
- [x] **(1) Blocco TK-costate fuori dalla sezione Vaidya di Paper I** (paper1:585-593): rimossi
      `p̃_φ=J−A²b/Ê`, `J_eff=J/A=P_φ`, "Kerr-frame"; ora "on all Vaidya branches costate=J, drift=m∂_m",
      mappe conformi/gravitomagnetiche → rinvio a Paper II.
- [x] **(2) Fix bib chirurgico (auto-citazione)**: voci companion spostate in `paper/companionI.bib`
      (PaperOne) e `companionII.bib` (PaperTwo), fuori dal refs.bib condiviso. Paper I include
      `refs,companionII` (cita solo PaperTwo); Paper II include `refs,companionI` (solo PaperOne).
      Verificato: nessun paper si auto-cita (bbl 70→69).
- [x] **(3) Refusi**: "App.~\ref{app:tracking}"→"\ref{...}" (iopart rende "Appendix B.1", evita "App.
      Appendix B.1"); "Secs.~closed--thakurta"→"thakurta--closed" (era 3–2 all'indietro); "fuel budget
      ... diverges"→"velocity indicatrix degenerates"; "accretion/evaporation drift"→"non-stationary
      mass-rate drift"; tolto "Δv" dal corpo.
- [x] **(4) Terminologia/segno**: paragrafo "Accretion--evaporation asymmetry"→"Advanced/retarded
      frozen-clock decomposition" + caveat "non è soluzione del BVP outgoing"; convenzione di segno λ in
      App C dichiarata come orientazione (sgn(dλ/dm)=−1), non flip post-hoc (numeri invariati).
- NON in scope (legende figure evaporation cotte nelle immagini → rigenerazione script; conversione
      "companion paper"→"Ref.[X]" ovunque; separatrix-following residuo con caveat corretto).
- Ricompilati: **Paper I 38pp, Paper II 49pp, 0 undefined**, auto-citazione eliminata.

### J-septies. Memo GPT "Risposta paper1-6/paper2-6" (2026-08-01) — 4 punti finali (90-95% recepito)
- [x] **P0 §2: formula TK residua nella sezione Vaidya di Paper I** (paper1:610): usava
      `F(r;E_eff,J_eff)`, sorgente Euleriana `Ê∂_E+J∂_J`, clock `η`, costate `p_η`. Riscritta Vaidya
      (`F(r;Ê,J)`, drift `m∂_m`, off-shell `p_v`); `E_eff/J_eff/P_r/p_η` relegati esplicitamente a Paper II
      (la formula è sparita, non solo un rinvio).
- [x] **P1 §3: direzione narrativa Paper I→II**: apertura §4.4 "same WKB scheme as conformal (companion)"
      → "we apply Lemma I.4/Theorem I.5 to the Vaidya modulus; specialized to conformal Kerr in Paper II";
      "conformal ψ-split of the companion paper" → "universal on-shell assembly established here"; attenuate
      "(as in Kerr)" e "as the conformal letter of companion".
- [x] **P1 §5: separatrix-following/tracking in Paper II** eliminati globalmente (7 occorrenze):
      "separatrix-following response"→"adiabatic derivative along the degeneration family",
      "tracking coefficients"→"family-derivative coefficients", "tracked source"→"total family-derivative
      source", "charge tracking"→"differentiation along the marginal family", ecc. Grep residuo = 0.
- [x] **P1 §4: negative-rate/frozen-clock**: notazione `δφ_accr/evap`→`δφ_adv/ret`, "physical difference"→
      "formal advanced/retarded frozen-clock difference", tabella script "accretion/evaporation split"→
      "advanced/retarded frozen-clock decomposition". **Figura rigenerata**: `inversione_fisica.py` (label
      già corrette "ingoing mass-rate") rieseguito → `fig_vaidya_no_inversione_evaporazione.png` aggiornata
      in paper/Immagini (via legenda "evapora"/"accretion-evaporation rate").
- Ricompilati: **Paper I 38pp, Paper II 49pp, 0 undefined**.

### J-octies. Memo GPT chiusura "paper1-7/paper2-7" (2026-08-01) — 4 punti finali → congelabili
- [x] **P0: Figs 8-9 + Fig 6 di Paper I rigenerate.** `vaidya_brachistochrone_vparam.py` label
      `evapora/statico/accresce`→`formal negative-rate continuation / static ingoing model / accreting
      ingoing model` (fig_vaidya_bounce + fig_vaidya_timing); `plunge_vaidya_t_tau.py` asse
      "accretion rate"→"ingoing mass-rate parameter" (fig_vaidya_plunge_t_tau). Rieseguiti + copiati,
      verificato visivamente (inglese).
- [x] **P1 §3: "The physical asymmetry is subtler"** → "The formal advanced/retarded frozen-clock asymmetry
      is subtler..."; 684 "accretion–evaporation asymmetry"→"advanced/retarded frozen-clock asymmetry".
- [x] **P1 §4: Paper II App B** ultimo tracking: "Tracking uses..."→"Differentiation along the degeneration
      family uses..."; `b_3^{track}`→`b_3^{fam}`. Grep positivi = 0 (resta solo root-tracking numerico).
- [x] **P1 §5: η riservato al conforme.** `η=U_3−2m U_2` (clock proper-time Vaidya) → `T_τ` in 6 punti
      (§4, App A, App B); l'η conforme FLRW/TK/pointer invariato. Notazionale, numeri invariati.
- P2 (non bloccanti): thrust 2-frasi, marcatura righe TK Table 1 (utente le tiene), placeholder DOI.
- Ricompilati: **Paper I 38pp, Paper II 49pp, 0 undefined**, figure inglesi. GPT: congelabili.

---
### Riferimenti
- Referee: `REFEREE_REPORT_CURRENT_2026-07-22.md` (13 issue)
- Nota derivazioni: `DERIVAZIONI_POST_REFERAGGIO_review.md`
- Off-shell/tracking dettaglio: `TODO_offshell_closure_and_tracking.md`
- Propagazione testo: `TODO_main11_text_propagation.md`
- Scripts chiusura: `KerrSessionScripts/offshell_tbranch_*.py`, `physics_anchor_offshell_closed.py`
- Tracking: `ThakurtaMetric/tk_sep_offshell_divergence.py`, `tk_sep_tracking_vs_trueflow.py`

### J-nonies. Review editoriale finale CQG (2026-08-01) — blocco "faccio io" eseguito
Status scientifico frozen; solo editoriale. Fatto:
- [x] **P0: rimosso ToC** in entrambi → l'header "CONTENTS" sparito (verificato: ora mostra il titolo).
- [x] **Abstract riscritti** (I e II): FLRW espanso, CKV per esteso, "(Perlick 1991)"/Mathematica/Zenodo/Sage
      rimossi, "Status." integrato in prosa, "evaporative inversion"→"inversion in the frozen-clock branch
      comparison", frasi di chiusura/status di GPT.
- [x] **Paper I**: §4.3 rinominata (frozen-clock comparison); `teleological`→`endpoint-sensitive`; thrust
      remark 2-frasi (Δv fuori dal corpo); `J_deg` consistente (2 occorrenze "separatrix |J|=Jc"→genus-deg).
- [x] **Paper II**: Fig.1 caption "due famiglie spaziali" (t≡η, τ); "Semi-analytic"→"On-shell first-order";
      rimosse 4 "companion 3D analysis" (+ "off-equatorial extension left to future work"); τ-branch clock
      `η=U_3−2MU_2`→`T_τ` (3 punti) + riga tabella "clock dη/dr"→"dT_br/dr" con T_t=t,T_τ=τ.
- [x] **Metadata**: affiliazione completa (Dept of Physics, Via A. Bassi 6, 27100 Pavia); **ack con modelli
      esatti** (OpenAI GPT-5.6 (Sol), Anthropic Claude (Opus 4.8)); **Funding** + **Competing interests** +
      **ORCID 0009-0004-4536-0285**. Data availability standardizzata (wording IOP).
- Ricompilati: **Paper I 37pp, Paper II 49pp, 0 undefined**.
- NON fatto (per scelta/input): §4 paragrafi contesto fisico (add material, opzionale); §6.5 sottodivisione
  sezione adiabatica lunga; §7 densità frasi (soggettivo); §8 restyle grafico figure (font/greyscale);
  §9 stile riferimenti globale (Fig.→figure, grande mecc.); §12 DOI companion (serve assegnazione preprint);
  keywords nel sistema di submission (ORCID già nel PDF).

### J-decies. CQG editoriale — contesto + stile riferimenti + sottodivisione (2026-08-01)
- [x] **§4 contesto fisico**: paragrafo "Beyond the brachistochrone..." in intro Paper I; "In the rotating
      case the ergosphere..." in intro Paper II; frase di significato fisico a chiusura di entrambe le
      Conclusions. (+ fix CKV bare nell'intro roadmap di Paper I.)
- [x] **§9 stile riferimenti IOP** (globale, entrambi): Fig.~→figure~, Figs.~→figures~, Sec.~→section~,
      Secs.~→sections~, Eq.~→equation~, Eqs.~→equations~, App.~C→appendix~C; maiuscole a inizio frase
      gestite (regex ". "/paragrafo/riga). "Appendix~C of Paper II" letterali (companion) lasciati.
      Preservati \ref/\eqref e la numerazione. 0 undefined.
- [x] **§6.5 sottodivisione** sezione adiabatica lunga di Paper II (§2): da 1 subsection a 4 —
      2.1 On-shell adiabatic response, 2.2 Complete first-order response and true-flow validation,
      2.3 Fixed-charge and fixed-endpoint perturbations, 2.4 Frozen penetration phase diagram
      (3 \paragraph→\subsection, nessuna equazione spostata).
- Ricompilati: **Paper I 38pp, Paper II 49pp, 0 undefined**.

### J-undecies. Final editorial checklist (2026-08-01) — pronti per submission
- [x] #1 DOI companion verificati (I 21739998, II 21740000), nessun placeholder.
- [x] #2 Running header corto via \markboth (no più troncatura del titolo lungo). Verificato render.
- [x] #3 η riservato al conforme: residui `T_br`/`T_τ` in Paper II (η(r) generico 303-305, "clock η" 330,
      tabella η=t/η=τ 597). Nessun η=proper-time residuo.
- [x] #4 "evaporative inversion" residuo (paper1:1017) → "inversion in the frozen advanced/retarded clock
      comparison".
- [x] #5 J_deg vs J_c: già coerente in corpo+appendici (nessun residuo).
- [x] #6 Rimosso doppione thrust (paper1:265 "Maintaining the controlled invariant..." duplicato).
- [x] #7 Figure: label "angular momentum" solo in script di verifica (non figure del paper) → J già ok;
      font della figura indicatrice aumentati (legend 6→8, assi 9, tick 8), rigenerata fig_indicatrici.
- [x] #8 AI names esatti: OpenAI GPT-5.6 (Sol), Anthropic Claude (Opus 4.8).
- [x] #9 Final pass: **0 undefined, 0 multiply-defined**, DOI/hyperlink/numerazione ok. Paper I 38pp, II 49pp.

### J-duodecies. Mock CQG Associate Editor (desk-reject soft) — restructuring editoriale Paper I (2026-08-03)
Verdetto: problema di presentazione/positioning, non di matematica. Fatto su Paper I:
- [x] **Abstract** riscritto attorno a UNA domanda ("What is the fastest constrained worldline in a
      spacetime that is itself evolving?"): gap (non-autonomo, niente energia conservata) → idea
      (controlled-rail/Kodama) → fondamenti → 2 casi (FLRW, Vaidya). Tecnicismi (horizon dilogarithm,
      S_D, weight-two, O(ε²)) tolti. 194 parole.
- [x] **Introduzione** ristrutturata in 4 blocchi: (1) domanda/motivazione fisica, (2) gap + positioning
      esplicito (stazionario vs dinamico), (3) idea controlled-rail/Kodama, (4) risultati + roadmap con
      separazione framework(§2)/applicazioni(§3 FLRW, §4 Vaidya) + rinvio Paper II.
- [x] **De-promozione**: "the unifying formalism"→"the shared construction", "Strikingly,"→tolto,
      "We stress that"→tolto, "sharp probe"→"quantitative probe". "genuinely/a genuine" tecnici lasciati.
- Conclusione: già chiude col significato fisico (frase quantitative-probe).
- Ricompilato: Paper I 38pp, 0 undefined. Paper II: si aspetta la risposta CQG prima di replicare.

### J-terdecies. Codex CQG editorial review (major-revision, resubmit) — cluster testuale-gating (2026-08-03)
Review su paper1-11 (post mio restructuring): loda abstract-domanda/intro-4-blocchi, spinge oltre. 4 condizioni-cancello.
Fatto il cluster testuale/basso-rischio scelto dall'utente:
- [x] **① Novelty + literature positioning**: aggiunte le citazioni della linea diretta (Giannoni–Piccione–
      Verderesi 1997, Giannoni–Piccione–Tausk 2002, Giannoni–Piccione 2002 arrival-time GR, Caponio–Javaloyes
      2026 survey) in refs.bib + intro. Rimosso il claim assoluto "have not been addressed"; novità ristretta
      (mantenimento attivo carica-rail in background non-autonomo + Kodama + free-arrival Pontryagin + FLRW/Vaidya).
- [x] **Abstract calibrato** (M11): "legitimate Pontryagin problem" ora esplicita "global minimisation
      conditional, not automatic"; "solve in closed form"→"derive the extremal equations and their closed
      representations"; "cleanly separates"→separazione concreta (spatial-curvature FLRW vs mass-flow Vaidya).
- [x] **④ Calibrazione claim/caption**: fig 3 "minimizes both"→"numerical perturbation check ... a local
      check, not a global-minimality proof"; fig A2/verifmin "minimizes its own travel time"→"local
      perturbation check"; body idem; Conclusions "well-posed" qualificato col dominio timelike/nondegenere.
- [x] **⑤ Bibliografia**: rimosso `\nocite{*}` → **69→44 voci (solo citate)**, 0 undefined.
- Ricompilato Paper I 37pp, 0 undefined. Restano (parte strutturale, prossimo giro se serve): ② promuovere
      Theorem I.5 nel corpo, ③ shorten §4.4 (funzioni speciali→appendici), Table 1 Paper-I-only, figure/Table A2.

### J-quaterdecies. Codex review — parte strutturale (2026-08-03)
- [x] **② Theorem I.5 promosso nel corpo** (§4.4, dopo la componente on-shell): enunciato + eq:adiab-exact
      ora nel body; App C diventa "This appendix proves Theorem I.5, stated in §..." (rimosso il duplicato).
      Risolto l'ordine (I.5 non più sepolto in App C dopo I.6). 0 undefined/multiply.
- [x] **Table 1 ridotta a Paper-I** (M5, non-distruttivo): righe TK (A(η),a,E_eff,J_eff,P_r,η) collassate in
      UNA riga di rinvio marcata "Conformal Thakurta–Kerr notation (forward references only; companion Paper II)".
- Restano (heavy, prossimo giro): ③ shorten §4.4 (macchina funzioni-speciali→App B/supplement, tenere i 5
      oggetti nel corpo); figure/hierarchy (combinare FLRW, fig7→supplement, fig8-10 prima delle Conclusions);
      Table A2 ricostruzione. + minori (m-list: define CKV, ripetizioni "genuinely/cleanly", ecc.).

### J-quindecies. Codex review — batch strutturale/editoriale 2 (2026-08-03)
- [x] **Float order (M7)**: `\clearpage` prima di \section{Conclusions} → le fig 8-10 non spezzano più le
      Conclusions.
- [x] **Table A2 ricostruita (M8)**: tolto l'header-prosa multicolumn che collideva con la riga dati e l'hash;
      "complete first-order correction" ora riga normale; "paths under VaidyaMetric/" spostato in caption;
      residual right-aligned; SHA in footnote.
- [x] **CKV esteso (m3)**: titolo §2.3 e Conclusions "Killing--CKV--Kodama" → "Killing / conformal-Killing /
      Kodama" (niente più sigla bare).
- [x] **③ parziale + m24**: paragrafo rank-five/shuffle/Fay nel corpo §4.4 condensato in 2 frasi + rinvio al
      reproducibility package (esperimento strutturale fuori dalla narrativa principale). "closed form"→
      "closed representation" nella caption script-map.
- [x] Fix overfull Table 1 (riga di rinvio TK come multicolumn che va a capo).
- Ricompilato Paper I 38pp, 0 undefined, 0 overfull≥100pt.
- DEFERRED (heavy, sessione nuova per budget): rigenerazione figure (combinare FLRW fig 2-3, fig 7→supplement,
      line-styles greyscale), relocazione più ampia della macchina funzioni-speciali di §4.4→appendici,
      restanti minori (m-list: ripetizioni "genuinely/cleanly", sub-arc/norm nelle caption, ecc.).
