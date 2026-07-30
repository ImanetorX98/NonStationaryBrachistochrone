# Technical Revision Memorandum main14 — change-list (30 lug 2026)

Fonte: `Referee_Memorandum_main14_for_Claude.pdf` (GPT, 10pp). Verdetto: core pubblicabile,
NON splittare ancora; correggi 5 fondamentali, armonizza terminologia/status, rigenera figure,
poi dividi in due paper (I: rail/Kodama+Vaidya sferico; II: Thakurta-Kerr conforme).
Legenda: [ ] aperto · [~] parziale · [x] fatto.

## P0 — Fondamentali — ✅ TUTTI FATTI
- [x] **4.1 Problema terminale definito** (commit 2ae215f). Apertura "between two events" → "from a given
      launch event to a fixed spatial target at a free arrival clock"; nuovo paragrafo "The optimal-control
      problem" (stato, clock, controllo, evento iniziale+target, s_f libero, costo, trasversalità free-clock);
      "fixed endpoints" riservato al two-point BVP (Sec bvp). main+PRD.
- [x] **4.2 HJB non-autonoma** (commit 22f994a) [mio blocco]. Eikonale = solo caso frozen; rail
      non-autonomo obbedisce `∂_s V+min_u{1+∇_x V·f}=0` (HJB estesa su stato+s). Framing come TEOREMA DI
      VERIFICA SUFFICIENTE NON costruito per traiettorie generiche; "candidate minimizer"/PMP extremal
      altrove. main+PRD.
- [x] **4.3 Costato ≠ meccanico** (commit 5c07823) [mia contraddizione]. Tolto "coincidono"; Noether del
      controllo dà il costato, non il meccanico; `dL_mech/dτ=a·∂_φ≠0` (rail forzato). Formula esatta
      VERIFICATA simbolicamente `L_mech=r²u^φ=((Ê²−f)/Ê)J` (eq:costate-mech, `costate_vs_mechanical.py`,
      diff=0); `J_eff=J/A` = normalizzazione del costato, non momento meccanico di Kerr. main+PRD.
- [x] **4.4 App. C canonica** (commit 3102013) [mia derivazione]. `P_r=p_r/A` NON canonica
      (`dr∧dp_r=A dr∧dP_r`): normalizzazione conf. simplettica. Derivazione in (r,p_r) originali:
      canonica `dp_r/ds=-∂_r H_br`, poi P_r + regola catena ⟹ ESATTAMENTE `-H̄_r-αP_r` (eq:normflow);
      `-αP_r` è conseguenza, non trasformazione canonica; sistema normalizzato non-canonico ma equivalente.
      Identità VERIFICATA (`canonical_dilation_check.py`, diff=0). main+PRD.

## P1 — Major — ✅ TUTTI FATTI (commit 6ac905f)
- [x] **4.5 R_*(E,a)**: definizione esplicita `R_*(E,a)=sup{r_0:V(r_0)≤4V(r_min)}` (eq:Rstar), dipende
      da (E,a) non solo E; tutte `R_*(E)`→`R_*(E,a)`. main+PRD.
- [x] **4.6 tracking**: rinominato "Adiabatic derivative of the degeneration family"; NON estremale
      torque-free (J conservato; serve torque esterno per seguire il locus). main+PRD.
- [x] **4.7 J_deg**: App. retitolata "genus-degeneration loci"; Vaidya J=Jc = J_deg algebrico (r_d<0), non
      separatrice; TK Jc^± = separatrice fisica esterna + J_deg. main (PRD già aveva J_deg).
- [x] **4.8 evaporazione**: m'<0 ingoing = "negative-rate continuation" (evaporazione fisica = outgoing).
      main+PRD.
- [x] **4.9 genus-2 peso-due**: App C corretta — U_k peso-1 (Kleiniane), W_jk iterati lunghezza-2 peso-2
      (non riducibili a peso-1, congetturalmente irriducibili), dilog genus-2 congetturale. main.

### (originale P1, per riferimento)
- [x] **4.5 R_*(E) + no-inversion.** R_*(E) mai definito. O definizione esplicita (formula/dominio/
      ipotesi/box) O declassamento coerente (grazing/quarter provato + asintotico statico + interval
      domini elencati + numerico altrove), abstract=conclusione stesso status tiered. `n_t/n_τ=E/f` NON
      è prova di ordinamento. [Sec 6.4, Eq(54-57), Table 5, Abstract, Concl]
- [ ] **4.6 "Separatrix tracking".** J conservato ⟹ una traiettoria non segue la separatrice mobile
      cambiando J senza torque. Rinominare "derivata adiabatica della FAMIGLIA di separatrici lungo
      J=J_sep(λ)"; non è evoluzione di un'estremale torque-free. [fine Sec 6.5; App B.1]
- [ ] **4.7 J_pen / J_sep / J_deg distinti ovunque.** App B chiama ancora "separatrice Vaidya" una
      degenerazione a r_d<0 (= J_deg). J_pen(v_0) soglia dinamica; J_sep separatrice fisica accessibile;
      J_deg doppio root algebrico (anche r_d<0, inaccessibile). "Separatrice" solo per confine reale.
      [Vaidya Sec 4; Sec 6.5; App B; caption; nomi script]
- [ ] **4.8 Evaporazione fisica vs sign-flip ingoing.** Restano posti che chiamano "evaporazione" curve
      ṁ<0 ingoing. ṁ<0 ingoing = "continuazione formale a rate negativo"; rietichettare Fig A11-A12;
      v/u = decomposizione frozen-clock/formale; niente claim quantitativi outgoing finché BVP non
      risolta; tenere outgoing come future direction. [Sec 4.2-4.4; Fig 6-7, A11-A12]
- [ ] **4.9 Status genus-2 peso-due.** App C dice W_jk "chiudono in ζ,σ Kleiniane classe peso-UNO di φ_0"
      — contraddice peso-due. Gerarchia: U_k Abeliani peso-1 (riducibili base Kleiniana); W_jk iterati
      Abeliani lunghezza-2, coeff algebrici; irriducibilità congetturale; completamento single-valued
      aperto; base 5-dim numerica non teorema. "length-two iterated Abelian integral" nei teoremi.
      [Abstract, Table 1, Sec 4.4/5.2, App B, App C]

## P2 — Editoriale / struttura — ✅ TESTO FATTO; immagini+struttura → fase split
FATTO (commit 8b3e421, 15052b2): tabella notazione (tab:notation); conclusione in 3 blocchi
(Proved/Conditional-numerical/Open); fix "App. Appendix X" duplicati (12×); "global attractor"→
"moving kinematic boundary"; J_c^+(A) equazione numerata (eq:Jcplus); Fig 9 forward-ref a Fig 10.
RESTA (richiede rigenerare le IMMAGINI figure via script, o è strutturale → fase split):
- [ ] Fig 10 (label "Eq(40)"→(43); slope legenda 2.12 vs caption 2.05/2.07/2.01 — una sola sorgente).
- [ ] Fig 9 (caption dice full-flow deferred ma Fig 10 lo fornisce → forward ref a Fig 10).
- [ ] Fig 14 (residuo Weierstrass-vs-quadratura E ODE-vs-closed vicino instabilità: due residui distinti).
- [ ] Fig 20 (cross-ref Sec 6.3 → 6.4).
- [ ] Eq per J_c^+(A): equazione numerata propria (non "Eq(48) caption").
- [ ] Cross-ref "App. Appendix A.1"/"App. Appendix C" duplicati; assi "particle ang. mom."→"axial costate J".
- [ ] Residui: specificare norma abs/rel, intervallo, tolleranza solver, precisione, indipendenza.
- [ ] Scala visiva: pannelli appendice troppo piccoli (split/ingrandire).
- [ ] "global attractor"/"freezing wall" → "moving kinematic boundary" (salvo teorema di attrattore).
- [ ] Spostare thrust/fuel in remark/nota; **rimuovere §5.1** quasi-costanti off-equatoriali (o svilupparlo
      a risultato 3D completo — ora sembra un terzo paper) → tenere per paper futuro.
- [ ] Tabella di notazione (Ê, E_eff, J, J_eff, a, A, m, s, η, λ, ε, clock).
- [ ] Non ripetere status in ogni paragrafo; status in theorem labels + Table 1 + conclusione.
- [ ] Conclusione in 3 blocchi: proved / conditional-numerical / open.

## Split (SOLO dopo P0+P1+P2 nel master)
- Paper I "Controlled-rail brachistochrones in non-stationary spacetimes: conformal symmetry, Kodama
  energy, Vaidya dynamics" (~25-35pp): problema controllo+terminali; esistenza/normalità/HJB non-auton.;
  gerarchia W; FLRW degenere; Vaidya (Kodama, costate memory, bounce, adiabatica, sorgente universale
  mass-function); limite: BVP outgoing aperta; App degenerazione algebrica solo se supporta senza dirsi
  separatrice.
- Paper II "Brachistochrones in conformal Kerr spacetimes: ergosphere trichotomy, separatrices, adiabatic
  response" (~30-40pp): indicatrice TK, t≡η; branch arrivo vs τ + dizionario costato/meccanico; Weierstrass/
  Kleiniane frozen + Doran; cusp/trichotomy/penetrazione; same-launch vs fixed-endpoint tiered; first-order
  on/off-shell canonicamente ri-derivato + true-flow; iterati Abeliani lunghezza-2 con terminologia cauta.
- Tenere per dopo: quasi-costanti off-equatoriali/3D; teoria turning-point mobili/Neishtadt crossing;
  framework higher-genus rigoroso + minimalità base.

## Note
- GPT nota 3 lacune reali in pezzi aggiunti da me: 4.2 (HJB autonoma), 4.3 ("coincidono" costato),
  4.4 (canonicità App C). Rigore, non invalidano i numerici.
- 5 P1 = armonizzazione (sostanza fatta, non uniforme su tutto il documento).
- Change-log per sezione richiesto (punto 8.1 del memo): tenere in questo file man mano.
