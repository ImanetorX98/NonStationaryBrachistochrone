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
- [x] **Wrap off-shell GENERICO t-branch CHIUSO** in funzioni speciali: A(2ª specie→Kleiniane ζ,σ + dilog s=0)
      + B(Δ=0 dilog genus-2) + C(elementare via Hermite). Verificato 1e-14 (2 config) + ancora fisica (slope~2)
- [ ] **FRONTIERA — theta-nome naming** dei dilog genus-2 (`D_{j,root}` a Δ=0, dilog s=0) nella classe
      tabulata Baune/D'Hoker (referee Issue 7). Serve rappresentazione theta + teorema di identificazione.
- [ ] **Vaidya generico off-shell** — curva `S_V=r·Emu·Q2(a=0)` genus-2, terza specie a `r=2m`, sorgente
      `Θ=m∂_m` (no dilatazione). Port meccanico del t-branch.
- [~] **Coeff simbolici all-(M,a,J)**: tabelle già simboliche; Hermite `rem_k`/`rho` mostrati E-simbolici;
      inverso modulare all-param = muro perf SymPy → usare Singular / tower QQ(a,E,J)[M]
- [ ] Cosmetico: cancellazione grande A≈−75 vs C≈+73 nella decomposizione di Hermite (decomp più naturale)
- [ ] Paper Tab.1: upgrade "off-shell closed form: open" → "t-branch assemblato+ancorato; naming e Vaidya deferred"

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
