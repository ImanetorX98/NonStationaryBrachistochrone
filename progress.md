# Progress — brachistocrone in spaziotempi non stazionari (sessione)

Traccia completa dei progressi. Riferimenti: `sumUp.md` (indice figure),
`KerrMetric/doranTau.md` (ramo τ), `KerrMetric/doranT.md` (ramo t),
`paper/main.tex` (CQG), `paper/main_prd_revtex.tex` (backup PRD).

Parametri di riferimento: `M=1, a=0.9, E=1.2`, `r_e=2M=2`, `r_±=M±√(M²−a²)`
(`r_+=1.4359, r_-=0.5641`), `J_c=a/E=0.75`, `Ω_H=a/(2Mr_+)=0.3134`.

---

## 1. φ(r) generale ramo t (J_+^t) — genere 2 Kleinian
- `dφ/dr = K(r)/√R6`, `K = r[(E²−1)r+2M](J(r−2M)+2Ma)/Δ` (GENERALE in M,a,E,J).
- `R6 = r Q2(r)[(E²−1)r+2M]` sestico; 6 radici semplici → genere 2.
- `J_+^t = (4M²+2a²+a²/E²)/(2a) = 2M²/a+a+a/(2E²) = 3.4347` (turning su r_e).
- Decomposizione 1ª+3ª specie, matrice periodi τ (Sage, Riemann OK), verifica
  vs flusso di Hamilton 4×10⁻¹².
- Script: `KerrMetric/kerr_jpt_genus2_kleinian.py` (Sage periodi),
  `kerr_jpt_genus2_reduction_check.py`, `kerr_jpt_genus2_figure.py`,
  `kerr_jcm_capture_figure.py` (cattura 0.95 J_c^-, saddle-node).
- Paper: eq. t-K, t-genus2; figure fig:jpt-g2, fig:jcm-cap.

## 2. φ(r) generale ramo τ — genere 2 Kleinian (NUOVO, calcolato esplicito)
- `dφ/dr|_τ = J r√(wf)/(Δ√(Δ−J²w))`, `w=E²−f`, `f=1−2M/r` (forma BL, eq.56).
- `dφ/dr = K_τ/√S`, `K_τ = J r(r−2M)[(E²−1)r+2M]/Δ`.
- **Sestico** `S(r) = r(r−2M)[(E²−1)r+2M][rΔ − J²((E²−1)r+2M)]`,
  6 radici semplici (per J generico) → genere 2. Per J=−0.9J_c:
  branch points {−4.55, 0, 0.036±0.686i, 1.93, 2.00}.
- **Decomposizione esplicita** (poli 3ª specie agli ORIZZONTI r_±):
  `φ(r) = c1∫r dr/√S + c0∫dr/√S + α_+∫dr/((r−r_+)√S) + α_-∫dr/((r−r_-)√S)`.
- **Coefficienti generali (M,a,E,J)** — verificati diff=0:
  - `c1 = J(E²−1)`
  - `c0 = 2MJ`
  - `α_± = ∓ J a²[M(E²+1) ± (E²−1)√(M²−a²)] / (2√(M²−a²))`
  - (num: c1=0.44J, c0=2J, α_+=−2.4453J, α_-=+2.0889J).
- Matrice periodi τ (Sage, J=−0.9J_c): `[[0.782+0.833i, 0.068+0.771i],
  [0.068+0.771i, 0.257+1.114i]]`, Riemann OK (‖τ−τᵀ‖=7e-12, Im τ≻0).
- Verifiche: decomposizione vs diretta 10⁻¹⁵; traiettoria speciale vs diretta
  4.7×10⁻¹³ (fig_tau_phi_special_vs_direct).
- **Riduzione a |J|=J_c**: la f² rende r_e radice DOPPIA; la (r−2M) di K_τ ne
  cancella una → √Q4 quartica = genere 1 Weierstrass = separatrice (eq.59),
  verificato 1.8×10⁻¹⁴.
- Script: `KerrMetric/kerr_tau_general_genus2.py` (+ `_periods.sage.py`).
- Paper: eq. tau-sextic, tau-genus2, tau-coeffs (appendice, dopo separatrice τ).

## 3. Tricotomia τ vs dicotomia t (CORRETTA — era sbagliata nella mia tabella)
- **τ TRICOTOMIA** (muro √(wf): reale solo r≥r_e):
  - `|J|>J_c`: rimbalzo LISCIO fuori r_e (periapside, r_min>r_e).
  - `J=+J_c`: UNICO che PENETRA (attraversa, spirala su r_+); √f si cancella.
  - `|J|<J_c` e `J=−J_c`: rimbalzo con CUSPIDE/CORNER a r_e.
  - insieme penetrante τ = SINGOLO PUNTO {+J_c}.
- **t DICOTOMIA**: penetra INTERVALLO `(J_c^-, J_c^+)=(−8.05, 3.43)`
  (retrograde inclusi, frame dragging auto-sintonizza 𝒦(r_e)=a/E); mai cuspidi.
  - `J_c^+ = 2M²/a+a+a/(2E²)` (alto spin) o saddle-node (basso spin, a<a*≈0.65).
  - `J_c^- = −8.05` (saddle-node retrogrado, spirala su r_*=3.514).
- Paper: Table A1 (tab:penetration) riscritta; §VI trichotomy affinato;
  Conformal trichotomy e didascalie corrette (solo +J_c attraversa).
- Immagine `cuspide_ergosfera.py` rigenerata (etichette "+J_c crosses").
- Figure: fig_master_penetration_taut (asse J), fig_atlas_tau, fig_atlas_t
  (5 regimi, flusso Hamilton adattivo solve_ivp), fig_tau_jc_pm (±J_c).

## 4. Scan fine regimi τ retrogradi (ogni 0.05 J_c)
- `k=−0.05…−0.95`: CUSPIDE a r_e (tutto |J|<J_c, intervallo continuo).
- `k=−1.00` (=−J_c): GRAZING/CORNER (marginale, r_min=r_e, dφ/dr finito=−1.04).
- `k=−1.05…−1.50`: SCATTERING (periapside liscio fuori r_e, r_min cresce).
- boundary layer cuspide `ε* ≈ 0.947(J_c−|J|) → 0` avvicinandosi a −J_c.
- Script/figure: `kerr_tau_scan_regimes_table.py`, `fig_tau_scan_regimes.py`,
  `fig_tau_shooting_mJc.py`.

## 5. Chiarimenti fisici (verificati)
- **Muro √f**: τ minimizza tempo proprio, `√(wf)` reale solo r≥r_e (dentro r_e
  niente frame statico). Fa cuspare TUTTE le |J|<J_c (prograde e retrograde):
  la svolta centrifuga `Δ−J²w=0` è DENTRO r_e (es. r_t=1.65 per +0.5J_c) ma
  IRRAGGIUNGIBILE; l'orbita si ferma al muro r_e (dφ/dr→0, azimut congelato).
- **Co-rotazione**: dentro r_e `dφ/dλ=(w/Δ)J`, segno=segno(J). Conta SOLO al
  marginale: blocca l'attraversamento di −J_c (che +J_c ottiene) → corner.
- **La forma chiusa NON distingue ±J_c** (riduzione dipende da J²): dφ/dr(±J_c,
  r_e)=±1.04 entrambi finiti. La non-penetrazione di −J_c è FISICA (co-rot),
  non geometrica.
- **Doran vs BL**: le φ(r) sono in BL. r_e NON è singolare in BL (Δ(r_e)=a²≠0):
  BL attraversa già l'ergosfera. Doran serve per (a) frame regolare del tempo
  proprio a r_e, (b) rimuovere il log-winding all'ORIZZONTE r_+ (poli 3ª specie).
- **Particella libera retrograda**: PUÒ penetrare l'ergosfera; L<0 resta
  possibile (Penrose) ma dφ/dt è forzato >0 (co-rotazione). Il vincolo del
  binario (non la geometria) ferma la brachistocrona τ retrograda.
- **Spirale sull'orizzonte UNIVERSALE**: ogni geodetica → `dφ/dt=Ω_H=a/(2Mr_+)
  =0.3134`, qualunque L (verificato L=+4,−4,−8,0). Il segno di L cambia la
  traiettoria PRIMA di r_+ (una retrograda inverte a r_e), non la spirale finale.

## 6. Layout/paper (fatto)
- Reformat CQG (iopart.cls) + backup PRD (revtex), parità di contenuto.
- ToC: fix sovrapposizione "Appendix A" (\renewcommand\numberline naturale).
- Referenze contigue (\clearpage prima della bibliografia).
- Tabella ottica Randers (τ/t/v/η per FLRW, Vaidya, Thakurta-Kerr).
- Metriche ottiche di ramo (Table 1): α, β, n per ogni ramo/metrica.

---

## 7. (B) Inversione Kleinian r(φ) — TENTATA, PARZIALE
Genere 2: r(φ) NON è un singolo quoziente-θ. Due vie analitiche:
- **Klein 2-punti**: `r1+r2=℘22(u)`, `r1 r2=−℘12(u)`, `℘ij=−∂²log σ/∂ui∂uj`.
- **Enolski–Hackmann–Lämmerzahl (integrale singolo, divisore-θ)**:
  `dr/dλ=√S`, `r(λ)=−σ1(u)/σ2(u)|_{σ(u)=0}`, `u1=λ`; poi `φ(λ)=∫K_τ dλ`.
- È FORMULA ANALITICA (σ = somma reticolare θ, convergenza geometrica), non
  iterazione; unico root-finding = localizzare il punto sul divisore-θ.

### Esito (mattoni)
- **(i) θ genere-2** `kerr_genus2_theta.py`: somma reticolare, quasi-periodicità
  verificata (intera 1e-14, τ-periodo 1e-11), convergenza geometrica (stabile N=6).
  **FATTO ✓**
- **(ii) Mappa di Abel + periodi** `kerr_genus2_abel_klein.py`,
  `kerr_tau_periods_export.sage.py` (esporta τ, period matrix A|B a prec=100):
  `τ=A⁻¹B` a 1e-16, `u(r)` calcolata per integrazione diretta dei differenziali
  `dx/y, x dx/y`. **FATTO ✓**
  (Nota: τ a prec=100 differisce da prec=40 per una trasformazione Sp(4,Z).)
- **(iii) Inversione ℘/σ** `kerr_genus2_klein_find_char.py`: **BLOCCATA.**
  Cercando la caratteristica half-integer, `℘22(u)−(r1+r2)` NON è costante
  (spread ~0.74, dovrebbe essere `−η22` costante). Motivo fondamentale:
  `℘ij=−∂i∂j log θ[δ] − ηij`, e la matrice **η** sono i **periodi di 2ª specie**
  (Baker), che **Sage non fornisce** (dà solo la 1ª specie). Nessuna
  caratteristica assorbe η (è additivo). Serve integrare i differenziali
  meromorfi di 2ª specie sui cicli di omologia + costanti di Riemann esatte.

### Conclusione
Le **φ(r) pure** (§1,§2) sono complete e nel paper. L'inversione r(φ) è 2/3
costruita (θ + Abel verificati); l'ultimo passo dipende da η (2ª specie) —
research-grade, direzione futura (Enolski–Baker). Non chiusa in modo affidabile.

Script mattoni: `kerr_genus2_theta.py`, `kerr_genus2_abel_klein.py`,
`kerr_genus2_klein_find_char.py`, `kerr_tau_periods_export.sage.py`.

---

## 8. WKB adiabatico φ(r,A) — forma IBRIDA (FATTO ✓)
A(η) lento lungo orbita (universo espande, E_eff=Ê/A scorre). Ordine dominante
WKB = famiglia CONGELATA (forme chiuse §1,§2 con E_eff istantaneo). Correzione
1° ordine O(A'/A). Forma finale:

  φ(r,A) = φ_0(r;Ê/A)  +  (A'/A)[ Closed(r) + ψ(r) ]  +  O((A'/A)²)

  Closed(r) = −½ Ê · ∂_E φ_0 · η(r)            [CHIUSO]
  ψ(r)      = ½ Ê (ρ − ρ̃)                       [NUMERICO 1D]
              ρ=∫∂_Eφ_0·h dr,  ρ̃=∫η·∂_E F dr,  h=dη/dr,  η(r)=t(r) flusso congelato

### Pezzi CHIUSI (analitici)
- φ_0: forma chiusa Kleinian (§1 t / §2 τ).
- ∂_E φ_0 = A(r)/√R + Σ_{k=0}^4 c_k ∫ r^k/√R  (riduzione 2ª specie completa:
  β_±=0, δ_±=0 → NIENTE terza specie, err=0 verificato ENTRAMBI i rami; gauge c5).
- Closed = −½ Ê ∂_E φ_0 η: prodotto di pezzi chiusi.

### Pezzo NUMERICO (unico irriducibile)
- ψ = ½Ê(ρ−ρ̃): integrale iterato iperellittico lunghezza-2 (polilog iperell.).
  Antisimmetrico ½(ρ−ρ̃) = parte IRRIDUCIBILE (shuffle: ρ+ρ̃=∂_Eφ_0·η chiuso
  simmetrico; ½(ρ−ρ̃) primitivo/irriducibile). ~26% correzione L2 (cresce 9%→28%
  verso turning). Integrato NUMERICAMENTE (trapezio 1D).

### Verifica ENTRAMBI i rami (parità)
Ramo **t** `kerr_adiabatic_phi_hybrid.py` (J=6, clock η=t coordinato dal flusso):
  A'/A=0.005: max|φ_hybrid−φ_full| = 1.39e-06
  A'/A=0.02 : 5.57e-06     A'/A=0.06 : 1.67e-05
Ramo **τ** `kerr_adiabatic_phi_hybrid_tau.py` (J=2.5 scattering, clock η=tempo
proprio η=∫L_τ dr, L_τ=√(Q/w) invariante doranTau.md §2, h=dη/dr=−L_τ analitico):
  A'/A=0.005: max|φ_hybrid−φ_full| = 2.5e-08
  A'/A=0.02 : 1.0e-07      A'/A=0.06 : 3.0e-07
Errore ~lineare in A'/A (residuo trapezio, non struttura). Ibrido = piena.
Forma di φ IDENTICA nei due rami (stesso E_eff=Ê/A); τ ha in più solo il
prefattore A⁻² sul TIMING, non sulla forma. NB τ: backoff dal turning
(F_τ~1/√(Δ−J²w) diverge a r_min) e h analitico → 1e-8.

### Approccio (concordato con utente)
Pezzi analitici trattati come chiusi; UNICO pezzo irriducibile (ψ) integrato
numericamente. Onesto: forma chiusa dove esiste, numerico solo dove dimostrato
irriducibile.

### Nel paper (FATTO ✓)
Sottosezione `sec:adiabatic` "Semi-analytic first-order adiabatic orbit shape"
in Sec V, ENTRAMBI i file (main.tex 43pp, main_prd_revtex.tex 23pp):
- eq:phi-adiab (forma ibrida), eq:dEphi0 (riduzione 2ª specie, no 3ª specie),
  eq:psi-irr (parte irriducibile), tab:adiab-valid (validazione t+τ).
- Framing onesto: ψ NON riducibile alla classe depth-1 (shuffle) MA è polilog
  iperellittico genus-2 (funzione speciale depth-2), valutato numericamente.
  NON scritto "non integrabile".

Script: `ThakurtaMetric/kerr_adiabatic_phi_hybrid.py` (ramo t),
`kerr_adiabatic_phi_hybrid_tau.py` (ramo τ, clock proprio),
`fig_adiabatic_pieces.py` (chiuso 74% vs irriducibile 26%),
`fig_adiabatic_curve.py`, `fig_breathing_families.py` (t+τ), `fig_breathing_wkb.py`
(non-autonomo + errore O(A'/A)).

---

## 9. Livello 3 (polilog iperellittico genus-2) — motore installato (FATTO ✓ primo mattone)
`abelfunctions` 0.2.0 compilato in SageMath 10.9 (vedi memoria
`abelfunctions-sage-install.md`): forzato build x86_64 (Sage in Rosetta) +
patch API deprecate `is_LaurentSeries` ecc. → shim isinstance.
Espone: `RiemannTheta` (con **derivate** via `derivs=[...]`), `RiemannSurface`,
`AbelMap`, `RiemannConstantVector`, `Jacobian`, `differentials`.
Verificato sul reticolo τ (§2): θ + gradiente + Hessiana simmetrica (dà
℘ij=−∂i∂j log θ − ηij), quasi-periodicità 5.8e-16. Sage nativo integro.
NB: eseguire dentro `sage` (non `sage -python`) da cwd neutra.
Resta il muro η (2ª specie, §7) per ℘ij pieno; ma il motore θ-derivate è ora
libreria robusta (non più somma manuale). Acknowledgements software aggiunti a
entrambi i paper (Sage Developers, abelfunctions/Swierczewski, SymPy/NumPy/SciPy/
Matplotlib).

### Mattone ψ (Abel map in avanti) — FATTO ✓
Pivot corretto: obiettivo = chiudere ψ(r) come polilog genus-2 (funzione IN
AVANTI di r), NON r(φ). Script `KerrMetric/kerr_psi_forward_abel.sage`.
Bug abelfunctions: `holomorphic_differentials()`/RS pipeline crasha su Singular
(`integralbasis`) per leading non-monico/coeff grandi (x^6-1 e monic-ish interi
OK; nostra sestica no). AGGIRATO: pipeline RS di **Sage** per τ, A|B (OK 4.9s);
Abel map in avanti u(r) per integrazione diretta di (1,x)/(2√S) conv. f_y=2y,
u=A⁻¹I; θ+derivate da abelfunctions (numpy, ogni τ).
Verifiche (params razionali M=1,a=9/10,E=7/5,J=5/2): A⁻¹B−τ=1e-16; u(r0)=0;
∇θ(0)=1.6e-16 (θ pari); Hessiana≠0. θ,grad,Hess valutati su tutta l'orbita
r∈[4.3,12] → ingredienti kernel Kronecker-Eisenstein genus-2 pronti.
PROSSIMO mattone: assemblare g^(n)(u,τ) e l'integrale iterato lunghezza-2
½(ρ−ρ̃) nelle 1-forme (dEF 2ª specie, L), verificare che riproduca ψ(r) numerica.

### Struttura polilog di ψ — dicotomia t/τ (FATTO ✓, validato)
Prerequisito assemblaggio: le due 1-forme di ψ (ω_a=∂_E F dr, ω_b=clock dη)
vivono sulla curva del ramo? Test simbolico+numerico (`/tmp/oncurve_check.py`,
`tclock_oncurve.py`, validazioni):
- **τ**: β cancella in Q; (dτ/dr)²·S_τ=[r²(r−2M)]² QUADRATO PERFETTO →
  `dτ/dr = r²(r−2M)/√S_τ` RAZIONALE su y²=S_τ (validato 8.9e-16).
  ω_a=∂_E F_τ dr 2ª specie su curva. ⇒ **ψ_τ = polilog iperellittico genus-2 PURO.**
- **t**: β cancella in Q_t; (Q_t/w)·R6=E²r⁶=(Er³)² → √(Q_t/w)=Er³/√R6 on-curve.
  MA il termine B/f porta β=√(2Mr/(r²+a²)) (frame-dragging Doran), NON in √R6:
  `dt/dr = ρ_t/√R6 + c_β√(2Mr/(r²+a²))`, c_β=(1−2Ma²/(rΔ))/f (validato 2.5e-14).
  ⇒ **ψ_t = [polilog genus-2 su R6] + [resto su cover frame-dragging β].**
Fisica: τ frame-independent → curva pulita; t trascina β (velocità fiume Doran)
→ rivestimento √ extra. Dicotomia polilog = dicotomia frame-dragging.
STATO closed form: differenziali ora ESPLICITI e on-curve (τ) / on-curve+cover (t),
validati. Ultimo passo (kernel Kronecker-Eisenstein genus-2 named / q-serie) =
research-grade, non ancora fatto. NON è forma tabulata finale, è la STRUTTURA
polilog rigorosa provata.

### ψ chiude in Kleinian ζ,σ (NON polilog) — η calcolata ✓ (muro §7-iii rotto)
CORREZIONE: ψ NON è polilog. Le sue 1-forme sono 1ª/2ª specie (no 3ª: β_±=δ_±=0).
Iterato antisimmetrico lunghezza-2 di 2ª specie CHIUDE in Kleinian ζ,σ (peso 1,
classe di φ₀). Es. ellittico: ψ=z ζ(z)−2 log σ(z). Mancava solo η (2ª specie).
**η CALCOLATA** (`KerrMetric/kerr_quasiperiods_bel.sage`): modello dispari
(quintica x=1/s, 1 punto ∞) + 2ª-specie canoniche Baker-Enolski-Leykin
(dr_1=(λ3 s+2λ4 s²+3λ5 s³)/4y, dr_2=λ5 s²/4y) + Sage matrix_of_integral_values
(interi). VALIDATA: κ=η ω⁻¹ simmetrica 1.4e-12, Legendre ω'ηᵀ−ωη'ᵀ=−iπ·I.
Il modello pari (deg 6, 2 punti ∞) sbaglia (x²dx/y 3ª specie). Ingredienti ψ
chiuso tutti pronti (u(r), θ+deriv, κ). Resta assemblaggio σ,ζ→ψ + validazione.
PAPER: sec:adiabatic dice "polylog" → correggere in "Kleinian ζ (2ª specie)".

### Natura di ψ — FINALE (decomposizione residui)
Correzione della correzione (onesto). Residui: res_∞(ω_a=∂_E F)=0 (pura 1ª+2ª,
niente 3ª); res_∞(ω_b=dη)=1.063≠0 (dipolo 3ª specie all'∞, perché dη~dr/r, η~log r).
⇒ ψ = [Kleinian ζ,σ peso 1, termini (2ª)×(1ª/2ª), chiudibile con η calcolata]
     + [dilogaritmo iperellittico peso 2, termini (3ª)×(1ª/2ª)=∫(int.abel.)·dlog,
        GENUINO, sourced dal residuo 3ª specie del clock all'∞].
Il "polylog" del paper era giusto nello spirito (c'è un dilog genuino); la mia
"ζ puro" era sbagliata. Ora preciso: ζ-chiudibile + dilog-irriducibile con origine
fisica (crescita log tempo proprio). η (2ª specie) resta utile per la parte ζ.
PAPER: sec:adiabatic "polylog" OK; volendo raffinare in "splits into a Kleinian-ζ
part (weight 1) and a genuine hyperelliptic dilogarithm (weight 2) sourced by the
clock's third-kind residue at infinity".

### Validazione decomposizione ψ vs ODE (vari innesco) — FATTO ✓
`ThakurtaMetric/psi_decomp_launch.py`. ψ = ψ_a (2ª specie, r³/√S → ζ) + ψ_b
(3ª specie, −2Mr²/√S → dilog). Closed+ψ_a+ψ_b = δφ_direct (−Ê∫dEF·η dr, verità
ODE) a **1e-6** per r0 ∈ {14,12,10,8,6.5}. Peso dilog ψ_b: 24.5%→47.4%
(cresce verso turning). VALIDA la struttura per-specie contro l'ODE, ogni innesco.
Livello mancante: valutare ψ_a via ζ(σ) e ψ_b via dilog iperellittico
INDIPENDENTEMENTE (serve caratteristica di Riemann + σ da θ+κ) e rimatchare.

### Livello 2: forma chiusa θ del 3ª-specie — VALIDATA ✓ (cerchio chiuso)
`KerrMetric/kerr_thirdkind_theta_closed.sage`. Il pezzo trascendente di ψ_Li è
l'integrale di 3ª specie η_b=∫h_b, h_b=−2Mr²/√S. Forma chiusa:
  η_b(r) = ρ₀ log[θ[δ](w(r)−e₊)/θ[δ](w(r)−e₋)] + holo,  δ = caratteristica ODD.
VALIDAZIONE (livello differenziale, scan 6 odd δ): d/dr della log-ratio θ,
fittata su base {r^k/√S, k=0..3} (1ª+2ª+3ª specie). δ **#1** unico con residuo
**7.3e-5** (prossimo 7.8e-4, 10×), coeff 3ª specie REALE (−0.2575). Ingredienti:
θ+deriv (abelfunctions), w(r) Abel normalizzato, e₊=w(∞), e₋=−e₊ (base branch pt).
Floor ~1e-4 (troncamento θ + e₋=−e₊ + 40 pt). ⇒ integrale 3ª specie CHIUSO in θ.
Quindi ψ_Li = ½Êρ₀[∂_Eφ₀·L − 2𝓛₂], L=log(θ-ratio) CHIUSO, 𝓛₂=∫L dA il dilog
(peso 2, endpoint). Combinato con L1 (decomp vs ODE 1e-6) e η (κ sym 1e-12):
forma chiusa di ψ validata end-to-end (struttura + pezzi speciali).

### Round-trip primitiva + PAPER aggiornato ✓
Round-trip (`/tmp/roundtrip.sage`): derivando la primitiva chiusa L(r) black-box
(differenze finite) si RIOTTIENE l'integrando algebrico di 3ª specie: dL_fd vs
dL_analitico = 4.5e-8; ricostruzione Σc_k r^k/√S residuo 5.6e-5; c2 reale. Primitiva
CORRETTA. PAPER (main.tex 44pp + PRD 23pp, compilano puliti, 0 undefined):
sec:adiabatic esteso con forme chiuse φ_t/φ_τ(r,A):
- eq:clock-tau (dτ/dr=r²(r−2M)/√S), eq:clock-t (dt/dr=ρ_t/√R6+c_β√(2Mr/(r²+a²)))
- eq:psi-split (ψ=ψ_ζ+ψ_Li, ψ_Li=½Êρ₀[∂_Eφ₀ L−2𝓛₂]), eq:thirdkind-theta (L=log θ-ratio)
- ψ_ζ Kleinian ζ,σ (Legendre 1e-12); L 3ª specie chiuso in θ[δ] (round-trip 5e-8);
  𝓛₂ dilog iperellittico peso-2 (endpoint, serie Kronecker-Eisenstein).

### Decomposizione ANALITICA (no fit) → ψ ha TRE pezzi (non due)
Riduzione 2ª specie di ∂_E F: c_k esatti (razionali in E), @E=7/5
[-0.531,1.979,-0.812,-0.360,0.189] (sympy, identità polinomiale, no fit).
Residuo 3ª specie clock: ρ0=M/(E²-1)^(3/2) ESATTO (analitico), num 1.0631.
Componenti olomorfe via a-periodi esatti (`KerrMetric/kerr_holo_component_check.sage`):
b^A=holo(∂_Eφ0), b^B=holo(clock). det = b^A_0 b^B_1 - b^A_1 b^B_0 = -4.27-7.94i,
|det|/(|bA||bB|)=0.80 ≠ 0.
⇒ **ψ = ψ_ζ (peso1, Kleinian) + ψ_ab (peso2, olomorfo×olomorfo ∫(u1 du2-u2 du1),
regolatore Beilinson) + ψ_Li (peso2, 3ª specie/dilog).** TRE pezzi.
ψ_ab è novità genus≥2 (genus1 ha 1 sola olomorfa → assente); l'ellittico zζ-2logσ
non ce l'ha. Il FIT least-squares nascondeva ψ_ab (assorbito nei coeff liberi);
l'algebra esatta lo rivela. 3ª specie: L=Fay (teoria), ρ0 analitico — NON fit.
PAPER DA CORREGGERE: sec:adiabatic dice ψ=ψ_ζ+ψ_Li (due pezzi) → sono TRE
(aggiungere il termine olomorfo×olomorfo peso-2).

### FORMA ESPLICITA ANALITICA di ψ — VERIFICATA (fix metodologico) ✓✓
`KerrMetric/kerr_psi_explicit_verified.py`. Dopo errori ripetuti (fit, gestione Δ),
fix metodologico = verificare OGNI passo prima di costruirci sopra.
- Riduzione CORRETTA: ∂_E F = N/S^(3/2), N=EJ r⁴(r−2M)²Emu (il Δ di K si CANCELLA).
  Poi 2N=2S𝒜'−𝒜S'+2SM, 𝒜 deg5, M=Σc_k x^k (k=0..4). Verificata dE F(diretto)=
  d(𝒜/√S)+M/√S a 1e-15. (I c_k erano giusti; bug era in 𝒜/costanti del bookkeeping.)
- FORMA ESPLICITA (identità, NO fit): ψ = ½Ê Σ_{k<j} Q_kj W_kj + ½Ê(peso≤1).
  Q_kj = c_k b_j − c_j b_k (ALGEBRICI), b=(0,0,−2M,1,0) (clock), W_kj=∫(U_k dU_j−
  U_j dU_k), U_k=∫x^k/√S (polilog genus-2). Verificata ρ−ρ̃=decomposizione a 4.8e-14.
- Q_01=0 → NIENTE olomorfo×olomorfo (ψ_ab=0, l'artefatto Hodge è risolto).
- Q_02=−2M c_0=1.063, Q_03=c_0, Q_12=−2M c_1, Q_13=c_1, Q_23=c_2+2M c_3,
  Q_24=2M c_4, Q_34=−c_4. c_k funzioni razionali esatte di E → Q_kj SIMBOLICI.
- W_kj = polilog iperellittici genus-2 (peso 2, NON riducibili a peso-1: teorema
  divisore theta Θ=W_{g-1}). ρ_0=M/(E²−1)^{3/2} per il pezzo 3ª specie (L=Fay).
CONCLUSIONE: coefficienti ALGEBRICI (simbolici, no fit, no periodi) × funzioni
polilog genus-2 (transcendenti, endpoint). Questa è la forma chiusa di ψ.

### Soglia di penetrazione ergosfera (ramo t) — diagramma di fase ✓
`KerrMetric/kerr_penetration_threshold.py`, `fig_penetration_threshold`.
Piano (A,J), 4 regimi: plunge(orizzonte)/penetra+rimbalza/scattering/forbidden.
DUE soglie ANALITICHE verificate:
- muro congelamento: A_c^wall = Ehat/sqrt(1-2M/r0) = 1.534 (E_eff<1 -> r_w=2M/(1-E_eff²);
  lancio proibito se r_w<r0). Verticale, indip. da J.
- ergosfera: J_c^+(A)=2M²/a+a+aA²/(2Ehat²) (confine penetra/scatter, diff 1e-15).
Finestra penetrante (J_c^-,J_c^+)~(-8,3.35) matcha dicotomia t di progress §3.
Fisica: A cresce -> muro avanza -> a A_c^wall espelle l'orbita (transizione osservata
nella fig penetranti). Risultato pubblicabile a sé.

### Figure penetrazione nel paper + setup adiabatico Vaidya ✓
Paper (main 48pp, PRD 25pp): aggiunte fig:penetration-phase (diagramma di fase A,J
con A_c^wall e J_c^+(A)) e fig:bounce (orbita penetra-rimbalza J=3.2), paragrafo
"Penetration phase diagram" con eq:Awall.
VAIDYA adiabatico impostato (`VaidyaMetric/vaidya_adiabatic_setup.md`,
`vaidya_dMF_reduction.py`): frozen=Schwarzschild (a=0 di Thakurta-Kerr), genus-2
(NON ellittica: 6 radici distinte; ellittiche sono le geodetiche, non la
brachistocrona vincolata). Parametro lento M(v), Ṁ=dM/dv. Riduzione ∂_M F=N_M/S^(3/2)
VERIFICATA 1e-15, c_k^M dati. Stessa pipeline polilog genus-2; clock v(r) (tempo
avanzato) da esplicitare. No teorema conforme (M(v) non è fattore conforme).

---

## 10. Chiusura W_ij (J generico, genus-2) in funzioni NOMINATE — MATTONE 1 (TK-τ) ✓
Obiettivo (utente): chiudere i W_ij (peso-2, ψ=½Ê Σ Q_ij W_ij) in funzioni speciali
nominate con coeff simbolici, come Brown-Levin Γ̃ sulla separatrice. Caso TK-τ.
Script `KerrMetric/kerr_tau_Wij_reduction.sage` (sympy). Params M=1,a=9/10,E=7/5,J=5/2.

SCHELETRO VERIFICATO (diff 0.00e0): i 5 integrali abeliani U_k=∫r^k dr/√S (k=0..4)
sono TUTTI indipendenti (k=4 provato irriducibile: grado minimo riducibile via forma
esatta = 5 per modello deg-6). Classificazione a r=∞ (a6=E²−1):
- U_0,U_1: 1ª specie (olomorfe) = coordinate Abel u.
- **U_2: UNICO generatore 3ª specie**, residuo 1/√(E²−1) ai due punti r=∞ → log σ-ratio.
  UNICA sorgente di peso-2 genuino.
- U_3,U_4: 2ª specie → Kleinian ζ_i(u).

⇒ Chiusura = ESATTO parallelo separatrice (Weierstrass+Γ̃):
- peso-1 (coppie senza U_2): chiude in Kleinian σ,ζ_i (Legendre/Baker), coeff simbolici
  da c_k, b=(0,0,−2M,1,0), κ. κ validata (`kerr_quasiperiods_bel`, Legendre 1e-12).
- peso-2 irriducibile: UN dilog genus-2 da U_2×(2ª specie), coeff ∝ ρ_0=M/(E²−1)^{3/2}
  (analitico). U_2=log[σ(u−e₊)/σ(u−e₋)] già validato (`kerr_thirdkind_theta_closed`, 5e-8).

Q_ij nonzero (7): Q_02=−2Mc_0, Q_03=c_0, Q_12=−2Mc_1, Q_13=c_1, Q_23=c_2+2Mc_3,
Q_24=2Mc_4, Q_34=−c_4 (c_k razionali in E, simbolici).
PROSSIMO: (a) Sage — U_3,U_4→ζ_i con coeff simbolici, verifica vs U_k diretto;
(b) forma nominata q-serie (Kronecker-Eisenstein/Fay genus-2) del singolo dilog = frontiera.

### Mattone 2a: riduzione peso-1 di U_k a integrali abeliani canonici ✓ (VERIFICATO 1e-14)
`KerrMetric/kerr_tau_Wij_oddmodel_reduce.py`. Modello dispari (quintica Y²=q6=s^6 S(1/s),
s=1/r, 1 punto ∞=r=0). ω_k=r^k dr/√S=-s^{1-k}/Y ds. Riduzione (forme esatte d(s^m Y),
m≤0, cancellano poli ordine≥2 a s=0; residuo 3ª specie n_-1/s tenuto esplicito):
  U_k = [P_k Y]_{r0}^r + c1_k R1 + c2_k R2 + g1_k u1 + g2_k u2 + n_-1,k L
con R_i=∫dr_i (2ª specie BEL), u_i=∫du_i (1ª), L=∫ds/(sY) (3ª specie). COEFF SIMBOLICI:
- ω_2: PURO 3ª specie, n_-1=-1, resto 0 → U_2=-L (letter canonico del dilog).
- ω_3: n_-1=(3-2E²)/(E²-1), c1=-2/(E²-1), P=1/(s(E²-1)).
- ω_4: n_-1=(-625E⁶+1156E⁴-37E²-794)/(200(E²-1)²), c1=3(3-2E²)/(E²-1)²,
       c2=-1/(E²-1), g1=(625E⁴-2581E²+2437)/(200(E²-1)), P=(6E²s+E²-9s-1)/(2s²(E²-1)²).
Tutti i pezzi 3ª specie ∝ stesso L → collassano in UN dilog nel montaggio ψ.
VERIFICHE: identità differenziale 0/7e-15/3.6e-12; identità INTEGRALE U_k(dir) vs
ricostruzione = 1e-14…1e-16 (r=10,8,6). ⇒ peso-1 CHIUSO, coeff simbolici razionali in E.
PROSSIMO (Sage): naming Kleiniano R_i→ζ_i(u), L→log[σ(u-e+)/σ(u-e-)]; poi il dilog.

### Mattone 2b: montaggio ψ nella base canonica ✓ (end-to-end VERIFICATO 1e-15)
`KerrMetric/kerr_tau_Wij_assembly.py`. A=∂_Eφ0=[𝒜/√S]+Σc_k U_k (c_k razionali E),
η=clock=Σb_k U_k, b=(0,0,-2M,1,0). Usando U_k=boundary+Σ_α M_kα V_α (mattone 2a),
V∈{u1,u2,R1,R2,L}: A=A_alg+Σa_α V_α, η=η_alg+Σh_α V_α. Decomposizione A,η verificata
vs diretto 1e-15. ⇒ ψ=½Ê Σ_{α<β} P_αβ w_αβ + T_alg, P_αβ=a_α h_β-a_β h_α SIMBOLICO.
STRUTTURA (clock ha solo 2 letters: R1 con h=-2/(E²-1), L con h=1/(E²-1); source a_L=0
= niente 3ª specie nel source):
- (u1,u2)=0  → niente olo×olo (conferma Q_01=0/ψ_ab=0).
- peso-1 Kleinian (3): (u1,R1),(u2,R1),(R1,R2)  [2ª×1ª, 2ª×2ª].
- DILOG (4, tutte con L): (u1,L),(u2,L),(R1,L),(R2,L) → condividono L = UN dilog genus-2,
  sorgente = residuo 3ª specie del CLOCK (h_L=1/(E²-1)).
- T_alg = termini algebrici (boundary elementari), verificato indipendentemente.
END-TO-END: ψ_dir = ½Ê Σ P_ab w_ab + T_alg a 1e-15 (r=10,8,6.5). Peso dilog ~3% del
trascendente qui (cresce verso turning). ⇒ ψ CHIUSA: [Kleinian peso-1 σ,ζ] + [1 dilog
genus-2] + [algebrico], coeff TUTTI simbolici razionali in E.
PROSSIMO: (Sage) naming w_peso1→ζ_i,σ e L→log[σ(u-e+)/σ(u-e-)]; poi q-serie del dilog.

### Naming Kleiniano: muro divisore-θ per ζ_i nudo (ONESTO, non riuscito con ζ nudo)
`KerrMetric/kerr_tau_Wij_naming{,2,3}.sage` (5 run Sage). Tentato: nominare R_i (2ª specie)
come ζ_i(u) Kleiniano nudo → weight-1 in funzioni Baker tabulate. FALLITO, strutturale.
- Diagnostica OK: ordine cicli coerente (2e-33), τ 1e-16, κ sym 1e-12, Legendre −I/2.
  RiemannConstantVector(abelfunctions) fallisce API ("must be a Place").
- CAUSA: ζ_i(u(r)) sull'orbita REALE cade sul divisore-θ (θ→0 → ζ esplode; dζ/dr O(1)
  complesso vs dR/dr O(0.03) reale; least_squares "residuals not finite"). = muro §7-iii.
- L (3ª specie) invece funziona (5e-8) perché è un RAPPORTO log[θ(u−e+)/θ(u−e−)], δ dispari:
  gli zeri di θ si cancellano num/denom. ζ nudo non ha questa cancellazione.
⇒ La 2ª specie NON è ζ_i nuda. Naming robusto = forma-DIFFERENZA ζ_i(u−e+)−ζ_i(u−e−)/σ-ratios
  agli e_± (i due r=∞, dove ω_3,ω_4 hanno i poli). Il muro §7-iii si AGGIRA (forma-differenza),
  non si sfonda (costante di Riemann). PROSSIMO: implementare forma-differenza, o passare a q-serie.

### Naming forma-DIFFERENZA: FUNZIONA ✓ (muro divisore-θ aggirato)
`KerrMetric/kerr_tau_Wij_diffform.sage`. Riusa convenzione validata di thirdkind (δ dispari,
e_±=±w(r=∞), base branch point, misura 1/√q). Nomina i differenziali canonici agli e_±:
- 3ª specie (U_2): D3=[ζ_δ(w-e+)-ζ_δ(w-e-)]·dw/dr = d/dr log[θδ(w-e+)/θδ(w-e-)] (log-ratio).
- 2ª specie (U_3,U_4): G±=[∇ζ_δ(w-e±)]·dw/dr = ζ_δ shiftata (= ∂ log θδ, NON log nudo).
VERIFICA: r^k/√S generato dalla base {du1,du2,D3,G+,G-} — δ#1 residui: k=2 2.0e-6, k=3 1.5e-6,
k=4 2.7e-4 (floor troncamento θ ~1e-4). ⇒ U_2,U_3,U_4 CHIUSI in forma-differenza θ[δ] agli e_±,
ROBUSTO al divisore-θ (argomenti shiftati). Il muro §7-iii AGGIRATO (non sfondato).
NB per l'utente: il log-ratio (θδ-ratio) basta per la 3ª specie/dilog; la 2ª specie serve la
sua DERIVATA ζ_δ=∂log θδ. Ma tutto da UN oggetto tabulato: θ[δ] dispari ai due punti e_±.

### Naming peso-1 a livello INTEGRALE: VERIFICATO ✓ (coeff da raffinare via residui)
`KerrMetric/kerr_tau_Wij_diffform_integral.sage`. Primitive nominate valutate DIRETTAMENTE
(log θ[δ]-ratio, ζ_δ come valori θ, non integrando derivate). δ#1. Verifica U_k(named) vs
U_k(direct):
- U_2: 5.8e-7..2.3e-6 ; U_3: 1.4e-6..4.1e-6 ; clock η=U3-2M U2: 2.5e-7..1.4e-6  ✓✓
- U_4: 8.6e-5..1.7e-3 (floor troncamento θ + polo ordine alto).
⇒ φ0 e CLOCK chiusi analiticamente in θ[δ] tabulata agli e_± (3ª specie=log θ-ratio;
2ª specie=ζ_δ). Coeff Lrat del clock ≈1.017 (residuo 3ª specie, pulito).
LIMITE ONESTO: coeff globali dal lstsq ENORMI/complessi (u1,u2~1e5, cancellazione) — base
sovracompleta (e_-=-e_+ rende i 4 ζ_δ dipendenti). Naming FUNZIONA ma i coeff puliti simbolici
servono dai residui/parti principali agli e_± (non fit). Raffinamento = prossimo passo.

### Coefficienti simbolici via PARTI PRINCIPALI ✓ (base canonica, cross-check esatto)
`KerrMetric/kerr_tau_Wij_principalparts.py`. Nella base canonica θ[δ] (oggetti a polo unitario
agli e_±: Ω=3ª specie, Z=2ª doppio polo, P=2ª triplo polo), i coeff di U_k sono i coeff di polo
di ω_k=-s^{1-k}/√q6 ds a s=0 = -g_{k-2}, g_i=Taylor di q6^{-1/2}:
  g_0=1/√(E²-1) ; g_1=(2E²-3)/(E²-1)^{3/2} ; g_2=(625E⁶-1156E⁴+37E²+794)/(200√(E²-1)(E²-1)²)
  U_2=-g0 Ω +olo ; U_3=-g0 Z -g1 Ω +olo ; U_4=-g0 P -g1 Z -g2 Ω +olo   (coeff SIMBOLICI).
CROSS-CHECK indipendente vs mattone-2a (residui BEL): n_2a·g0-(-g_k)=0 ESATTO (U_3,U_4).
⇒ risolve il mal-condizionamento del fit globale: coeff FISSI simbolici, gli oggetti Ω,Z,P
(θ[δ] agli e_±) portano la normalizzazione geometrica (come σ,ζ,℘ sulla separatrice).
CLOCK nominato: η=U3-2M U2 = -g0 Z + (2M g0 - g1) Ω + olo, coeff simbolici:
  Z: -1/√(E²-1) ; Ω: 2M/√(E²-1) - (2E²-3)/(E²-1)^{3/2}.
⇒ φ0 e clock CHIUSI in θ[δ] agli e_± con COEFFICIENTI SIMBOLICI. Naming peso-1 COMPLETO.

### q-serie dilog genus-2 — TAPPA 1 (fondamenta) + ostruzione Siegel
`KerrMetric/kerr_tau_dilog_qseries1.sage`. (1) θ[δ] genus-2 somma-nome vs RiemannTheta = 2.4e-9.
(2) ground-truth dilog Λ(r)=∫(A dΩ-Ω dA), Ω=log θ-ratio: Λ(10)=-0.283, Λ(8)=-0.876, Λ(6.5)=-1.679.
(3) NOMI: |q11|=0.045, |q22|=0.009 (piccoli), ma |q12|=6.71>1 (Im τ12=-0.606<0). τ di Sage NON
Siegel-ridotta → serie naïve in q12 diverge (somma reticolare converge, Im τ≻0). 
⇒ q-serie pulita richiede riduzione Siegel Sp(4,Z) (Tappa 2), poi log θ nome (T3), poi kernel
Kronecker-Eisenstein genus-2 (T4=frontiera). NB: il dilog È già nominato (∫log θ-ratio × 2ª specie
= dilog ellittico genus-2 via somma reticolare θ convergente); la q-serie in nomi è raffinamento.

### q-serie tappa 2: riduzione Siegel — serve la trasformazione S
`KerrMetric/kerr_tau_siegel_reduce.py`. Minkowski(Im τ)+shift(Re τ) da soli NON riducono
(q12 resta 6.71). La trasformazione S=-τ^{-1} (inversione modulare/Fricke) dà nomi ridotti
|q11'|=0.054, |q22'|=0.176, |q12'|=0.587 (tutti <1) → q-serie converge in τ'=-τ^{-1}.
NB: la somma reticolare θ converge SEMPRE (Im τ≻0); la riduzione serve solo per la serie di
POTENZE nei nomi. RESTA (frontiera): T3 log θ-ratio come nome-serie in τ'; T4 kernel
Kronecker-Eisenstein genus-2 per Λ (con θ trasformata sotto S). = ricerca aperta.

### AUDIT coefficienti simbolici (genus-2, risposta onesta)
SIMBOLICI e verificati (razionali in E): c_k (source), b=(0,0,-2M,1,0) (clock),
Q_ij=c_k b_j-c_j b_k, coeff mattone-2a (c1_k,c2_k,g1_k,g2_k,n_k in base BEL/odd),
P_ab (assembly), h_L=1/(E²-1), a_L=0, g_i (parti principali, cross-check ESATTO 0).
⇒ TUTTI i coeff di residuo/polo (la fisica) sono simbolici.
NUMERICI (come sulla separatrice): dati delle funzioni speciali — punti e_±, periodi τ,
normalizzazioni di Ω,Z,P; caratteristica δ discreta/esatta. Come σ,ζ,℘ e z_d,z_∞ sulla sep.
CAVO APERTO onesto: la parte OLOMORFA (coeff u1,u2) nel cambio base BEL→canonica θ[δ] mescola
fattori di periodo numerici (Z,P def mod olomorfe). Simbolica in base BEL (mattone-2a), non
ancora ripulita in base canonica. Analogo agli additivi C0/Ce della separatrice (che pure
portavano dati di punto marcato). Estraibile via a-periodi (impor: U_k senza periodo spurio).

### q-serie tappa 3: Ω=log θ-ratio in serie di NOME convergente ✓
`KerrMetric/kerr_tau_dilog_qseries3.sage`. Via S=Fricke (τ'=-τ^{-1}, nomi 0.026/0.088/0.228 <1)
+ caratteristica-zero (δ assorbita nell'argomento, θ0 si trasforma pulita). Formula:
  Ω = [log θ0(t^{-1}ζ1;τ')-log θ0(t^{-1}ζ2;τ')] - iπ(ζ1 t^{-1}ζ1 - ζ2 t^{-1}ζ2) + 2πi a·(e_--e_+)
  ζ_i=(w-e_pm)+τa+b, t=τ. VERIFICA vs diretto (RiemannTheta a τ): 3e-7..1e-10 (mod 2πi ramo log).
⇒ la 3ª specie del dilog HA serie di nome convergente. RESTA T4: kernel Kronecker-Eisenstein
genus-2 per l'integrale iterato Λ=∫(2ª specie)×Ω (frontiera vera).

### Chiusura parte olomorfa: NON estraibile dall'orbita (period-level, come C0/Ce)
`KerrMetric/kerr_tau_Wij_holomorphic.sage`. Tentato: fissare coeff polo ai g_i, risolvere
olomorfa (u1,u2) dal residuo. RISULTATO: ancora mal-condizionato (cond 2e5, coeff ~370 con
cancellazione) ANCHE per U_2 (puro 3ª+olo). CAUSA INTRINSECA: sull'orbita reale du_1=ds/√q e
du_2=s ds/√q sono QUASI PARALLELI (rapporto s ~ costante su arco corto) → α u1+β u2 non
separabile dai dati d'orbita. ⇒ i coeff olomorfi sono PERIOD-LEVEL: determinati dal vincolo
degli a-PERIODI (∮_{a_j}ω_k = Σ pole·∮pole + α∮du1 + β∮du2), NON dall'orbita. STESSO status di
C0/Ce sulla separatrice (principio a-periodi/punto marcato, no fit, ma dipendono dai periodi).
NON sono razionali puri in E (intrinseco a genere≥1). Chiusura esplicita = calcolo a-periodi
(residui simbolici × periodi), step Sage dedicato.

### (A) Parte olomorfa CHIUSA via a-periodi ✓ (sanity esatto)
`KerrMetric/kerr_tau_Wij_aperiods.sage`. (α_k,β_k)=∮_{a_j}ω_k · ω^{-1}, ω_k=x^k dx/y modello
PARI (Sage matrix_of_integral_values, differenziali polinomiali). SANITY ESATTO: k=0→(1,0),
k=1→(0,1). Coeff olomorfi period-level: k=2 (-0.215+1.105j, 0.451-0.102j); k=3 (1.344+3.523j,
4.212+2.326j); k=4 (-0.356+6.721j, 10.127+4.091j). Determinati dagli a-periodi (PRINCIPIO, non
fit; ben condizionato), NON razionali in E — dipendono dai periodi come C0/Ce separatrice.
⇒ DECOMPOSIZIONE CANONICA COMPLETA: U_k = (poli, coeff simbolici g_i) + α_k u1 + β_k u2 + const.
Tutti i coeff determinati da principi (residui simbolici + a-periodi), ZERO fit. Status = separatrice.

### (B) q-serie tappa 4: il dilog È alimentato dalla serie di nome ✓ (milestone)
`KerrMetric/kerr_tau_dilog_qseries4.sage`. Lambda=int Omega dA calcolato con Omega_nome (tappa 3,
serie di nome in tau') vs Omega_diretto (RiemannTheta): diff 1e-9 (r=10,8,6.5). ⇒ la
rappresentazione in nomi alimenta il dilog end-to-end. STRUTTURA derivata:
  Lambda = Sum_{n,m} c_n d_m e^{2pi i(n+m)w}/(2pi i(n+m))  (Kronecker-Eisenstein genus-2)
  c_n=Fourier(Omega=log theta-ratio), d_m=Fourier(2a specie). FRONTIERA APERTA: la resummazione
analitica di questa doppia somma in forma NOMINATA chiusa (kernel Kronecker-Eisenstein genus-2)
= research-grade (letteratura Enriquez/Schlotterer). Il dilog E' nominato+calcolabile via nomi;
la forma-serie-chiusa esplicita resta il pezzo di ricerca.

### (B) q-serie CHIUSA come serie di nome KE genus-2 ✓ (tappa 5, convergenza geometrica)
`KerrMetric/kerr_tau_dilog_qseries5.sage`. Split: Λ = ∫Q dA [ELEMENTARE, quadratica in coord
Abel, ~80%] + ∫L dA [dilog puro, ~20%]. Il pezzo L (log θ0-ratio) HA serie di nome ESPLICITA
GEOMETRICAMENTE convergente: N=1→1e-2, N=2→2e-4, N=3→2e-8, N=4→9e-14, N=6→0. ⇒ q-serie CHIUSA.
FATTO ONESTO: NON riducibile a Li2 classici — genus-2 non ha formula prodotto di Jacobi (triple
product e' genus-1), quindi log θ0 != Σ log(1-x_k) -> polilog genus-2 GENUINO (Enriquez), non Li2.
STATUS = separatrice: là Γ̃ era serie di nome g^(1)=π cot+4πΣ q^{2n}/(1-q^{2n})sin (NON elementare);
qui KE genus-2 (serie di nome 2D, NON classica). Stesso standard: [coeff simbolico P_ab/h_L] ×
[serie di nome, coeff = dati di periodo q'^{Q(n)}]. La q-serie del dilog e' CHIUSA (convergente
esplicita); i suoi coeff interni sono dati di periodo (come le q-potenze di g^(1)), non razionali E.

### T_alg/G_alg ESPLICITO in forma chiusa elementare ✓ (verificato 0 esatto)
`KerrMetric/kerr_tau_Talg_explicit.py`. Prima solo verificato numericamente come resto; ora
CALCOLATO esplicito. G_alg=2 I_el + boundary, I_el=∫A5(r³-2Mr²)/S dr (i due 1/√S -> 1/S ->
RAZIONALE). Forma chiusa: I_el=P(r)+Σ res_i log(r-r_i), P cubica razionale, res_i=A5(r_i)
(r_i³-2M r_i²)/S'(r_i) sui 6 zeri di S (3 nulli: r=0,2M, una radice). Boundary=η·(alg). VERIF
forma chiusa vs diretto = 0 esatto (r=10,8,6.5). ⇒ T_alg ELEMENTARE (polinomio+log)+boundary,
NON irriducibile. Era l'ultimo pezzo asserito-ma-non-calcolato: ora esplicito.

### T_alg/I_el coefficienti SIMBOLICI in E ✓ (Mathematica, cross-check)
`paper/crosscheck_Ialg_symbolic.wl` + `crosscheck_genus2.wl`. I coeff di I_el sono razionali/
algebrici (non period-level come α,β) → resi SIMBOLICI: P(r)=[p_k(E)r^k]/D(E) razionale in E;
log-part = RootSum sul cubico C(r)=rΔ-J²DE con res(x)=A5(x)x²(x-2M)/S'(x) razionale. Residui NULLI
provati simbolici a x=0,2M,2M/(1-E²). Coeff diversi E=7/5 vs 13/10 (non universali, giusto renderli
simbolici). Cross-check Mathematica INDIPENDENTE: g0=1/√(E²-1), g1=(2E²-3)/(E²-1)^{3/2} esatti;
c_k razionali in E; Q_kj pattern esatto; T_alg integrando razionale→elementare. Ora T_alg
INTERAMENTE simbolico (nessun coeff numerico residuo).

### Precisazione: P(r), c_k, Q_kj, g_2 razionali in (M,a,E,J) — NON solo E
`paper/crosscheck_P_params.wl`. Il coeff r^3 di P(r) con a,E,J simbolici dipende da a,E,J tutti
(∂_a,∂_J,∂_E ≠0). Verificato: a=9/10 → -14504578125/443030960699 (= valore prima); a=1/2 diverso.
⇒ P(r),c_k,Q_kj,g_2,res(x) sono razionali in (M,a,E,J), formule simboliche universali (4 variabili).
ECCEZIONE: g_0=1/√(E²-1), g_1=(2E²-3)/(E²-1)^{3/2} dipendono SOLO da E (leading r→∞, solo DE conta).
Corretto GENUS2_CLOSED_FORM.md ("razionale in E" -> "razionale in M,a,E,J").

### Coefficienti SEPARATRICE pienamente SIMBOLICI in (M,a,E,r_d,Jc) ✓ (tutti i rami)
`SEP_COEFF_SYMBOLIC.py`, `SEP_COEFF_SYMBOLIC.md`, `VaidyaMetric/sep_coeff_symbolic.py`.
Formula universale b_i (residui R polo triplo) via DERIVATE DI S a r_d: Q4(rd)=S''(rd)/2,
Q4'(rd)=S'''(rd)/6, Q4''(rd)=S''''(rd)/12, a4=[r^6]curva. b3=h0/s^3, b2=(h1-3a1h0)/s^3,
b1=(h2-3a1h1+(6a1^2-3a2)h0)/s^3, F=N/Q4. r_d,Jc=doppia radice (S(rd)=S'(rd)=0).
Prima erano valutati a params fissi; ora SIMBOLICI in (M,a,E,r_d,Jc). Verifiche: Vaidya tau
match contorno 1e-7 (b1=0.2704,b2=0.0326,b3=0.0099); TK tau Jc=20.328,r_d=-7.130; TK t+
Jc=19.089,TK t- Jc=-18.671 (match noti). a4=E^2-1, e2_zi=1/(E^2-1) (solo E). Vaidya v = stessi
b_i di tau (clock diverso). Riutilizzabile: (M,a,E)->doppia radice->plug formule.

### Cross-check Mathematica indipendente dei b_i separatrice ✓ (1e-16)
`paper/crosscheck_sep_bi.wl`. Via INDIPENDENTE dalla formula h0/s^3: r(t) via InverseSeries
dell'ODE dr/dt=sqrt(Q4), estrazione Laurent della sorgente al polo triplo. Match a precisione
macchina: Vaidya tau b1 7e-16,b2 2e-17,b3 0; TK tau b1 1.5e-14,b2 7e-17,b3 7e-18.
b_i confermati per TRE vie: sympy formula, contorno Python (1e-7), Laurent Mathematica (1e-16).
BUG DEBUG (rigoroso): Q4f=Cancel[S/(r-rd)^2] con r_d NUMERICO -> Cancel non elimina (r-rd)^2 da
poly decimale -> 0/0 -> catena a zero. FIX: Taylor locale Q4^(m)(rd)=m! S^(m+2)(rd)/(m+2)!
(solo derivate di S a rd, niente divisione). 

### Ramo Vaidya v: residui clock SIMBOLICI ✓ (verificati)
`VaidyaMetric/sep_v_clock_residui.py`. Sorgente b_i = Vaidya tau (identica). Clock v_z=E r^3/(r-r_d)
+r sqrt(Q4)/(r-2M): residuo z_d = E r_d^3/s (diff contorno 3.6e-7); residuo orizzonte z=i w_im = 4M
(ESATTO, diff 0, indip dai e_i). Ramo v completo. RESTA: residui clock ramo t (rho_t, poli r±).

### Ramo TK t: residui clock SIMBOLICI ✓ (verificati) - SEPARATRICI COMPLETE
`ThakurtaMetric/sep_t_clock_residui.py`. Clock rho_t=P3+R_Delta/Delta, in z: etpz=rho_t/(r-r_d).
Residui: z_d = rho_t(r_d)/s ; orizzonti z(r±) = sigma R_Delta(r±)/((r±-r∓)(r±-r_d)sqrtQ4(r±)),
sigma=-1 (foglio sqrt Q4). Q4(r±)=R6(r±)/(r±-r_d)^2. INVARIANTE res(r+)+res(r-)=2M (verificato).
Match contorno: z_d,r+,r- tutti 1e-6. ⇒ SEPARATRICI COMPLETE: sorgente b_i (tutti i rami) +
residui clock (tau, v, t tutti) SIMBOLICI e verificati. Restano solo additivi Ce,C0 (period-level).

### Ce,C0 espliciti + natura period-level (risposta concettuale) ✓
`vaidya_sep_C0Ce_closed.py` (Ce,C0 chiuse, verif 1e-8), `sep_periodlevel_test.py` (test 2 params).
Ce=η'(0)+2e1_zd ζ(z_d)-2e2_zi ℘(z_∞)+2e1_zi ζ(z_∞); C0=-Σ[b1 ζ+b2 ℘-b3/2 ℘'](z_∞-a). Coeff
(e_i,b_i) SIMBOLICI; valori ζ,℘ ai punti = period-level. TEST E=7/5 vs 13/10: ω1,z_d,z_∞,ζ(z_d),
℘(z_∞) TUTTI diversi -> period-level NON universali, NON razionali-simbolici, valutati per-curva
(come K(m)). Gerarchia: residui=razionali; e_i,g2,g3=algebrici; z_d,z_∞,ζ,℘,Ce,C0=trascendenti.

### Vaidya tau GENERICO (genus-2, J qualsiasi) - coeff simbolici ✓ (tutti i mattoni)
`VaidyaMetric/vaidya_generic_coeff.py`. Parametro m (massa). Curva S=r(r-2m)DE(r^2(r-2m)-J^2 DE),
sorgente dm F (N_m=S dm K-1/2 K dm S, K=J DE). Mattoni: (1) dm F=N/S^{3/2} =0; (2) c_k^m razionali
in (m,E,J); Q_kj=c_k b_j-c_j b_k (b=(0,0,-2m,1,0)); (3) g0=1/sqrt(E^2-1), g1=m(2E^2-3)/(E^2-1)^{3/2},
g2 razionale; (4) P(r) razionale. Verifica riduzione dm F 1e-15. STESSO schema di TK-tau, coeff
pienamente simbolici in (m,E,J). [g1 mostra dipendenza da m, TK aveva M=1]
