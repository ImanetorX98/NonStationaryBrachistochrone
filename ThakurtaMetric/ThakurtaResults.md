# Brachistocrone in Thakurta (Schwarzschild conforme) — risultati

**Script**: `thakurta_brachistochrone_sympy.py` (tutte le verifiche passano).
Contesto: `../FLRWmetric/FLRWresults.md` (mossa conforme),
`../SdSMetric/SdSresults.md` (limite statico con Λ),
`../VaidyaMetric/VaidyaResults.md` (Hamiltoniane PMP).

```
g = a(η)² · ĝ_Schw :   ds² = a(η)² [ −f dη² + dr²/f + r²dΩ² ] ,  f = 1−2M/r
```

Il matrimonio esatto dei due semilavorati: mossa conforme FLRW su base
statica Schwarzschild. Come modello fisico di BH cosmologico è discusso
(flusso di calore radiale, dibattito PBH); come laboratorio geometrico è
il caso minimale con **geometria ottica non banale E orologio non banale**.

## R1. Simmetria
`£_{∂η} g = (2a′/a) g`: `∂_η` è Killing **conforme per ogni a(η)**
(verificato); genuino ⇔ `a′ = 0`. Niente stazionarietà: Perlick non
si applica direttamente.

## R2. Cinematica della rotaia conforme `Ê = −u_η`
```
dτ/dη = a²f/Ê ,    γ = Ê/(a√f) ,    v² = 1 − a²f/Ê²
```
**Superficie di congelamento** `a(η)²f(r) = Ê²` — una curva in (η, r) che
unifica le due barriere note:
- limite FLRW (f=1): istante cosmico `a = Ê`;
- limite statico (a=1): raggio `f = E²` (barriera w=0 di SdS/Schwarzschild).

Con `a` crescente la regione accessibile `f < Ê²/a²` si **restringe verso
l'orizzonte**: l'espansione spinge le particelle di rotaia sul buco nero.

## R3. Indicatrice e Hamiltoniane (esatte)
Ellisse senza vento (base statica, `′ = d/dη`):
```
r′ = f v cosθ ,    φ′ = (√f·v/r) sinθ
H_η = v√f·√( f p_r² + J²/r² ) − 1
H_τ = v√f·√( f p_r² + J²/r² ) − a²f/Ê
```
`p_φ = J` conservato; `dH/dη ∝ a′` (forzante conforme); trasversalità
`H = 0` all'arrivo; statico ⇒ H conservata ≡ 0.

## R4. Riduzione tipo Fermat — dove si va OLTRE Perlick
Lungo ogni worldline di rotaia: `dη = dl_opt / v` con `dl_opt` l'arco
della **metrica ottica di Schwarzschild** (`dr²/f² + r²dφ²/f`). Il ramo
η è Fermat con indice
```
n(η, r) = 1/v = Ê / √(Ê² − a(η)² f(r))
```
- `a = 1`: `n = E/√(E²−f)` — **Perlick 1991 recuperato** esattamente;
- `a(η)` generico: `n` dipende **congiuntamente** da η e r ⇒ niente
  orologio globale (il trucco λ(η) di FLRW richiede f=1, omogeneità):
  Fermat non separabile — oltre Perlick in modo essenziale;
- ramo τ: costo `a²f/Ê` pesato in spazio E tempo: mai riducibile.

## R5. Validazione statica (a=1)
Flusso di Hamilton (coordinate di Schwarzschild, non EF) vs forme chiuse:
`dφ/dr` a 10 cifre su r = 4, 6, 9 per **entrambi** i rami
(`𝒦_τ = J`, `𝒦_η = fJ/E`).

## R6. Cattura da espansione (de Sitter conforme, a = −1/(Hη), H=0.02)
Orbita lanciata **in fuga** da r=4 (J=1.3, Ê=1.2, a(η₀)=1):
- sale fino a r = 8.74 mentre la superficie di congelamento scende
  (∞ → 8.74): **contatto** a η = −36.59 con `v → 0`;
- dopo il contatto l'indicatrice degenera a un punto: la particella
  **cavalca la superficie**, che la trascina a `r_freeze → 2M`
  (3.21 → 2.004 entro η → 0).

**Cattura da espansione**: fenomeno assente sia in FLRW (niente
orizzonte) sia in Schwarzschild (niente espansione) — il primo risultato
genuinamente "BH cosmologico" del programma.

## R7. Thakurta ORIGINALE: Kerr conforme (`thakurta_kerr_sympy.py`)
Il paper di Thakurta (1981) è "Kerr in un universo in espansione":
`g = A(η)²·g_Kerr` — il BH cosmologico **rotante** minimale. Equatoriale,
rotaia `Ê = −u_η` (`v̄² = 1−A²f/Ê²`, `P̄ = P + A²(2Ms/r)²/Ê²`):

**Indicatrice = ellisse con vento angolare** (esatta):
```
centro: φ′₀ = (2Ms/r)v̄²/P̄  (frame dragging) ,  r′₀ = 0
R² = f v̄² + P̄φ′₀²  =  Δ·v̄²/P̄        [identità: fP + (2Ms/r)² = Δ]
```

**Hamiltoniane in forma chiusa**:
```
H_η = p_φφ′₀ + R√( (Δ/r²)p_r² + p_φ²/P̄ ) − 1
H_τ = p̃_φφ′₀ + R√( (Δ/r²)p_r² + p̃_φ²/P̄ ) − A²f/Ê ,
      p̃_φ = p_φ − 2MsA²/(rÊ)     (shift gravitomagnetico CONFORME)
```

Validazioni: s→0 = Thakurta-Schwarzschild (numerico, 1e-12); A=1 =
forme chiuse Kerr equatoriali doranTau/doranT (`𝒦_τ = J`,
`𝒦_t = (fJ+2Ms/r)/E`) a 10 cifre su r=4,6,9 (spin 0.9, E=1.2, J=1.3).

**R7a — le due superfici critiche si separano nettamente:**
- `Δ = 0` (orizzonte `r_+ = M+√(M²−s²)`): degenerazione **conformemente
  invariante** — non dipende né da A(η) né da Ê (verificato: radice di
  R² identica a 1e-13 per ogni A). L'indicatrice resta regolare
  attraverso l'ergosfera (a f=0: R² = v̄²Δ/P̄ > 0), come il worldline
  vincolato in Doran;
- `v̄ = 0` (congelamento conforme `A²f = Ê²`): **respira** — scende
  dall'infinito quando A supera Ê (A=1.25→r=25.5; A=3→r=2.38) e
  comprime la regione accessibile nella shell orizzonte–congelamento:
  versione rotante della cattura da espansione.

## R8. Metrica ottica dipendente dal tempo conforme per Thakurta-Kerr
(`thakurta_kerr_optical_sympy.py`)

Risolvendo l'indicatrice per `dη` a spostamento spaziale dato — **Randers
non autonomo in forma chiusa**:

```
dη = n(η,r) · α_K + β_K

α_K² = r²dr²/(fΔ) + Δdφ²/f²          (Riemanniana)
β_K  = −(2Ms/(rf)) dφ                 (1-forma gravitomagnetica)
n(η,r) = 1/v̄ = Ê/√(Ê² − A(η)²f)     (indice)
```

**Teorema (verificato, residui ~1e-14):** `α_K` e `β_K` sono *esattamente*
i dati di Randers del Fermat **nullo** di Kerr (Gibbons–Werner). Tutta la
fisica di massa, rotaia ed espansione sta nell'indice scalare `n(η,r)`,
che moltiplica solo la parte Riemanniana; la 1-forma di trascinamento è
**rigida**: invariante conforme e indipendente da Ê.

Limiti (verificati): `Ê→∞` ⇒ n=1 (ottica nulla di Kerr); `A=1` ⇒
`n = E/√w` (Perlick stazionario su Kerr; EL di F riproduce doranT a 10
cifre); `s=0` ⇒ Thakurta-Schwarzschild; `M=s=0` ⇒ `n = 1/v` FLRW.

Singolarità: `f=0` (ergosfera) — α, β divergono: la *riduzione ottica*
muore all'ergosfera (1/F intrinseco, come in Kerr statico) mentre
l'Hamiltoniana/indicatrice resta regolare fino all'orizzonte (R7a);
`v̄=0` — n→∞ al congelamento conforme (superficie mobile).

È il "Finsler non autonomo" previsto da `nonStationaryBrachi.md` §1.2,
ora esplicito: **Zermelo con paesaggio dipendente dal tempo, in cui
respira solo il termine metrico, non il vento**.

### R8a. I tre rami dalla stessa ellisse (costi diversi, stessa indicatrice)

| ramo | costo L₀ | Hamiltoniana | metrica ottica |
|---|---|---|---|
| η | 1 | H_η (R7) | `dη = n·α_K + β_K` (Randers) |
| t cosmico | A(η) | `H_t = H_η-forma − A(η)` | **stesse curve di η**: `t(η)=∫A dη` monotona ⇒ orologio ottico (coppia η,t degenere come in FLRW) |
| τ | `A²𝔉/Ê` (dipende dal controllo) | H_τ con `p̃_φ` (R7) | `dτ = (A²f/(Êv̄))·α_K` — **Riemanniana PURA** |

Cancellazione verificata (1e-12 su punti casuali dell'ellisse):
sostituendo `dη = nα+β` in `dτ = (A²/Ê)(f dη + (2Ms/r)dφ)`, il termine
`fβ` elide esattamente il gravitomagnetico ⇒ il ramo τ non ha vento —
versione conforme del "ramo τ Riemanniano puro" di doranTau, che quindi
sopravvive all'espansione. Il vento è fenomeno esclusivo dei rami di
tempo d'arrivo.

## Figure
`genera_figure_thakurta.py` → `Thakurtafigures/`: cattura da espansione
(orbita vs superficie + regione accessibile), superfici critiche
Thakurta-Kerr (orizzonte rigido vs congelamento che respira), residui
flusso vs forme chiuse Kerr A=1 (1e-13). In `../PaperFigures/`:
`fig_indicatrici` (concettuale: le indicatrici nei tre casi).

`genera_figure_thakurta_rami.py` → `fig_thakurta_kerr_rami` (R8a in
figura): (a) con |J|=1.3, arrivo r=8: il ramo τ RIMBALZA in coppia
speculare (`Δφ = ±2.039`, periassi uguali a **1e-11** — la Riemannianità
pura non è imposta, emerge dal flusso con p_φ → −p_φ); il ramo η SPIRALA
sull'orizzonte per entrambi i segni (attraversamento marginale
auto-sintonizzato, |J| < J_t^c ≈ 3.4) ed **entrambe le spirali vincono
in senso progrado** (`Δφ = +3.45` e `+2.55`): il trascinamento forza la
co-rotazione (meccanismo doranT). (b) i due scalari conformi `n` e `k_τ`
su α_K rigida, con asintoti al congelamento.

## R9. Inversione del ramo τ e tricotomia all'ergosfera
(`thakurta_kerr_plunge_inversion.py`, figura `fig_thakurta_kerr_plunge_map`)

NOTA terminologica: la tricotomia (doranTau) classifica il comportamento
alla superficie di luce/ERGOSFERA `r_e = 2M`, non all'orizzonte. La
cattura (attraversamento + spirale su r_+) avviene SOLO sulla
separatrice `J = J_c`, a misura nulla; sotto J_c c'è riflessione a
cuspide su r_e, non plunge.

Con A costante, `g = A²g_Kerr` è Kerr riscalato: rotaia Ê ⇔ problema di
Kerr con `E_eff = Ê/A`, `J_eff = J/A` (T_τ̄ = A·T_τ ⇒ p_φ → A·p_φ).
**Raggio di inversione analitico** (tricotomia riscalata):

```
J > J_c(A):  periasse liscio, radice di Δ(r) − (J/A)²[(Ê/A)² − f(r)] = 0
J < J_c(A) = s·A²/Ê:  riflessione a CUSPIDE sull'ergosfera, r_min = r_e = 2M
J = J_c(A):  attraversa l'ergosfera e spirala su r_+ (cattura, misura nulla)
```

**Validazione colormap statica** (griglia 25×25 in (J, A), flusso di
Hamilton vs analitica, cuspide inclusa): residuo relativo **max 7e-12**
— contorni analitici indistinguibili dalla colormap numerica; la linea
`J_c(A)` è la separatrice di cattura. Verificato `H_τ(p_r=0, r_inv)=0`.

**Dinamico (a = −1/(Hη), H=0.02, lancio da r₀=8, A₀ ∈ [1, 1.3])**:
lo shift non adiabatico dell'inversione `(r_min^num − r_qs)/r_qs` è
**sempre positivo, +49% ÷ +120% (mediana +73%)**: inversione
**anticipata** rispetto al quasi-statico — il congelamento in
avvicinamento decelera globalmente la particella e la fa svoltare molto
prima; nessuna cella catturata (la separatrice resta a misura ~nulla su
questa griglia). A differenza di Vaidya, l'adiabaticità fallisce al
primo ordine anche per H piccolo: il drift `dH/dη ∝ A′` è cumulativo
lungo l'infall. Domanda aperta: la dinamica allarga la separatrice di
cattura all'ergosfera in finestra finita (analogo di Vaidya R8)? Serve
uno zoom attorno a `J_c(A)` con classificazione cuspide/attraversamento.

## R10. Plunge inversion t vs τ (protocollo Kerr) e il suo controllo conforme
(`thakurta_kerr_plunge_t_tau.py`, figura `fig_thakurta_kerr_plunge_t_tau`;
complementare: `thakurta_kerr_inversione_t_tau.py`)

Protocollo di `KerrMetric/evidence_t_plunges_deeper.txt`: stesse condizioni
iniziali `(r₀, p_r0, J)` nei flussi delle Hamiltoniane ottiche QUADRATICHE
dei due rami, integrazione fino a `p_r = 0`, confronto di `r_min` rispetto
alla singolarità. Le quadratiche di Thakurta-Kerr (da R8/R8a):

```
H_τ = (1/2k²)[(fΔ/r²)p_r² + (f²/Δ)J²]
H_t = (1/2n²)[(fΔ/r²)p_r² + (f²/Δ)X²],   X = J + 2Ms/(rf)
```

a A=1 coincidono ESATTAMENTE con quelle del progetto Kerr
(`c_p = FΔ·DE/(E²Σ³)` ecc., verificato termine a termine).

**Rapporto delle Hamiltoniane** (potenziali centrifughi a p_r=0):

```
V_t/V_τ = ρ² ,   ρ = A²[ f + 2Ms/(rJ) ]/Ê
```

**Il piano dell'inversione è (a, J)** — spin del BH vs momento angolare
z della particella, come nello studio Kerr; A è il parametro che sposta
il luogo.

- **Riproduzione quantitativa dell'evidence** (A=1, a=0.4, r₀=10,
  p_r=−0.5): Δr(J=0.5) = **+0.0616** (evidence: +0.06), Δr(J=4.5) =
  **−0.417** (evidence: −0.42) — le Hamiltoniane conformi riproducono i
  numeri del codice Kerr indipendente.
- **Colormap (J, a) a A=1**: r_min algebrico (H conservata), ODE a
  **7e-11**; luogo esatto di inversione J_inv(a) da 0.32 (a→0) a 1.75
  (a=0.98), crescente con lo spin (dragging).
- **Formula analitica corretta**: la condizione di inversione è

  ```
  ρ(r*)² = H_t0 / H_τ0        [NON ρ = 1]
  ```

  perché a parità di (r₀, p_r0, J) i due rami partono con valori
  hamiltoniani DIVERSI (coefficienti 1/n² vs 1/k², X vs J): il rapporto
  dei potenziali alla svolta va normalizzato al rapporto di lancio.
  In campo lontano `H_t0/H_τ0 ≈ ρ(r₀)²` ⇒ lettura compatta:
  **inversione dove ρ(r*) = ρ(r₀)** (il profilo ρ torna al valore di
  lancio). Verificata contro il luogo esatto: deviazione **8.5e-13**.
  (La versione ingenua ρ=1 assumeva H_t0 = H_τ0 ed era spostata.)
- **Sovrapposizione numerica/analitica**: colormap e contorno Δr=0
  rifatti con ODE PURA (28×28×2 flussi); luogo di inversione numerico
  per bisezione sui flussi vs formula: coincidenza punto per punto,
  **max |ΔJ| = 2.4e-9** su a ∈ [0.24, 0.96] (per a ≲ 0.2 J_inv scende
  sotto il floor di scansione J=0.3: il dragging svanisce e
  l'inversione si sposta a J→0).
- **Drift conforme**: A cresce ⇒ il luogo di inversione si sposta a J
  più alti (J_inv(a=0.9): 1.58 → 1.68 → 1.93 per A = 1 → 1.1 → 1.2):
  la regione "τ più profondo" si espande con l'espansione; in campo
  lontano ρ → A²/Ê: per **A > √Ê** l'inversione diventa generica
  (E_eff = Ê/A sotto ~1: vince la dilatazione temporale, crossover
  tipo TOV, ora pilotato dal fattore di scala).

### R10b. Teorema: l'inversione conformale esiste ed è raggiungibile
(`inversione_conformale_AJ.py`, `fig_thakurta_kerr_inversione_AJ`)

**CORREZIONE** di una nota precedente ERRATA ("l'inversione non entra
nella finestra accessibile, 0/484 celle"). Il conteggio 0/484 era giusto
per QUELLA finestra (`J ∈ [3.9, 5.5]`), ma la conclusione era falsa: la
finestra cadeva **tra** due inversioni. Come funzione di J (a `a=0.9`,
`A` dato) `Δr = r_t − r_τ` cambia segno DUE volte:
- a **J piccolo** (`J ~ 1.6–1.9`): inversione rotazionale (R10, piano
  `(a,J)`, `ρ(r*)²=H_t0/H_τ0`);
- a **J grande** (`J ~ 16–24`, per `A ≳ 1.2`): **inversione conformale**
  (questo risultato).
Lo scan `[3.9, 5.5]` stava nel mezzo (t più fondo) → 0 celle, ma
l'inversione è eccome accessibile su entrambi i lati.

**Condizione di inversione in forma chiusa** (r_min^τ = r_min^t, cioè
entrambi i rami svoltano allo stesso r; dalle svolte `p_r=0, H=0`):
```
cond(r,A) = A² b Q − P̄ (Ê − A² f) = 0 ,   Q = b v̄² + √(v̄² Δ)
J_inv = P̄/Q
```
**Teorema (esistenza + raggiungibilità)**. Ai bordi in A:
- `cond(A=1) < 0` (Kerr, "t più fondo", nessuna inversione);
- `cond(A_freeze) = Ê r Δ (Ê−1)/(r−2M) > 0` per r>2M, Ê>1 (ESATTO:
  a `v̄²=0` si ha Q=0 e il bracket `Ê−A²f → Ê−Ê² = −Ê(Ê−1)`, col
  fattore `1/(r−2M)`).

Per il teorema del valore intermedio esiste **`A_inv(r) ∈ (1, A_freeze(r))`**
per ogni raggio di svolta r>2M: l'inversione conformale **esiste sempre
ed è SOTTO il congelamento (raggiungibile)**. Far-field `A_inv → √Ê`
(coerente con R10: A>√Ê rende l'inversione generica).

**Verifiche**:
- integrazione diretta dei due rami a `(A_inv, J_inv)`: `r_min^τ =
  r_min^t` a ~1e-4 (es. r*=6: A=1.324, J=16.76; r*=9: A=1.236, J=24.27);
- mappa `(J, A)` estesa (J fino a 30): 1029/3271 celle in inversione;
  contorno numerico `Δr=0` ≡ curva analitica `cond=0` punto per punto.

Ampiezza onesta: nel rosso τ è più fondo di poco (`Δr ~ +0.1÷0.3`), a
J grande (orbite debolmente legate, svolta a r*~6–9 fuori ergosfera):
inversione genuina ma "gentile", non plunge profondo.

### R10c. L'inversione dipende dal PROTOCOLLO: sparisce a estremi fissi
(`KerrMetric/colormap_spin_estremi_fissi.py` → `fig_colormap_spin_estremi_fissi`;
`KerrMetric/colormap_conforme_estremi_fissi.py` → `fig_colormap_conforme_estremi_fissi`;
BVP: `KerrMetric/bvp_kerr_simmetrico.py`, `KerrMetric/bvp_conforme_inversione.py`,
`KerrMetric/bvp_estremi_fissi.py`)

R10/R10b confrontano i due rami a **stesso lancio** (stesso r₀, stesso
istante): là l'inversione esiste, sia rotazionale che conforme. Ma il
confronto fisicamente pulito è a **estremi fissi**: A=(r₀,−Φ), B=(r₀,+Φ)
identici per t e τ. In quel protocollo l'inversione **NON avviene**.

**Teorema (conforme, Schwarzschild via transfer E_eff=Ê/A).** Con gli
indici ottici `n_τ=√(f/(E²−f))`, `n_t=E/√(f(E²−f))` si ha
```
n_t/n_τ = E/f > 1   per ogni r>2M, E>1,   e decrescente in r.
```
La geodetica di Randers/Beltrami dà `r_min` come radice esterna di
`N(r)·r = J`; a Φ fissato n_t più grande ovunque ⟹ il ramo t svolta più
in alto ⟹ `Δr = r_min^t − r_min^τ > 0` **sempre**. Nessuna radice.

**Verifiche numeriche (tutte Δr>0, nessun cambio di segno):**
- **spin**, Kerr flusso di Hamilton, griglia (J,a) 18×20, a∈[0.05,0.95]:
  `Δr ∈ [+0.446, +2.064]`; mappa a bande verticali, gradiente in J
  (angolo estremi), **piatta in a** → lo spin quasi non conta;
- **conforme**, (J,A) via E_eff, griglia 44×40: `Δr ∈ [+0.002, +3.711]`,
  contorno Δr=0 assente (conferma il teorema);
- **BVP diretta** (integrazione dei rami, root-shooting su J a Δφ=2Φ):
  simmetrico `r_τ=3.617 < r_t=4.735`; asimmetrico `r_τ=4.00 < r_t=5.29`.
  A estremi fissi è **τ** a scendere più in basso (opposto a stesso-lancio).

**Estremi ASIMMETRICI** (`colormap_spin_asimmetrico.py`,
`colormap_conforme_asimmetrico.py`; A=(r_A=10,0), B=(r_B=6,Δφ)). Il
verdetto non dipende dalla simmetria: colonne piene su griglia (J,a) e
(J,A) danno di nuovo `Δr>0` ovunque —
spin `Δr ∈ [+0.762, +2.636]`, conforme `Δr ∈ [+0.161, +3.244]`. Mappe a
bande verticali (gradiente in J, piatte in a/A), coerenti col teorema
`n_t/n_τ=E/f` che è **pointwise in r** e quindi indipendente dagli estremi.

**Quadro protocol-dependence:**

| protocollo | spin a | conforme A |
|---|---|---|
| **stesso lancio** | inverte (R10) | inverte (R10b) |
| **estremi fissi** | NO (Δr>0) | NO (teorema n_t/n_τ=E/f>1) |

L'inversione di plunge è dunque un effetto del **protocollo di
confronto**, non una proprietà geometrica invariante delle due
brachistocrone: richiede la libertà sul punto d'arrivo. Fissati entrambi
gli estremi, t sta sempre più in alto di τ.

## R11. Quasi-costanti: cattura universale al muro e programma polare
(`thakurta_kerr_quasicostanti.py`, figura `fig_thakurta_kerr_quasicostante`)

Tentativo: invariante adiabatico `I_r = (1/2π)∮p_r dr` per librazioni
radiali equatoriali nella shell (regime legato `E_eff = Ê/A < 1`), col
piccolo parametro naturale `ε = A′/A` per periodo (l'espansione
dell'Hamiltoniana non autonoma).

**Risultato negativo ma esatto — cattura universale**: il potenziale
radiale `∝ w_A/Δ` è monotono ⇒ la svolta esterna è SEMPRE il muro di
congelamento; lì la velocità orbitale si annulla (`∝√v̄²`) mentre il
muro scende a velocità finita: al primo contatto misurato
`|ṙ_orbita| = 8×10⁻⁴` vs `|ṙ_muro| = 3×10⁻²` (**38×**) ⇒ chi tocca il
muro è catturato, la librazione radiale non esiste, `I_r` non è
definibile. (Rafforza R6/R9: la cattura da espansione è senza scampo
sul ramo τ equatoriale.)

**Dove vive la quasi-costante** (programma per la prossima sessione,
analogo diretto delle quasicostanti Carter di Kerr):
1. `I_θ = ∮p_θ dθ` — l'azione POLARE fuori equatore: il moto in θ è una
   librazione genuina (doppia svolta centrifuga), robusta al muro;
   richiede le Hamiltoniane 3D (indicatrice ellissoidale, stessa
   costruzione);
2. **quasi-Carter conforme**: `Q_eff = Q_Carter(E → Ê/A(η))` — il tensore
   di Killing di Kerr è conforme in Thakurta-Kerr (esatto per la luce);
   per la rotaia massiva drift `dQ/dη ∝ A′` calcolabile al primo ordine
   e testabile col workflow drift-colormap del progetto Kerr.

### R11a. Trasferimento delle espansioni K_t, K_τ al sesto ordine
Le espansioni O(a⁶) del progetto Kerr (`K_tau_expansion.txt`,
`risultati_paper_summary.txt` §1–2: `K_τ = Q_std + a²f₂p_θ² + a⁴f₄p_θ²
+ a⁶f₆p_θ²`, con verdetto numerico: NLO ottimale per τ, Q_std nuda
ottimale per t) si trasferiscono a Thakurta-Kerr **quasi-staticamente
per sostituzione** `E → E_eff = Ê/A(η)` nei coefficienti (riscalamento
conforme esatto a A costante). Due osservazioni:

1. **Il muro di congelamento appare come polo dell'espansione**:
   `DE₀(r) = (E²−1)r² + 2Mr = r²·w(r)` esattamente ⇒ i denominatori
   `D₀·DE₀`, `Δ₀ = D₀²DE₀` ecc. si annullano dove `w_eff = 0`, cioè sul
   muro `r_w = 2M/(1−E_eff²)` — che in Kerr (E>1) non esiste (w>0
   ovunque, solo il polo r=2M) ma in Thakurta-Kerr (E_eff<1) entra dal
   campo lontano e AVANZA verso il buco nero con A(η): il raggio di
   convergenza della serie quasi-costante è mangiato dal congelamento.
2. Il drift totale ha due piccoli parametri: `O(a²)` (residuo Kerr,
   invariato) + `O(A′/A)` (drift conforme): gerarchia testabile con le
   drift-colormap. La saturazione a NLO e la divergenza per a ≳ 0.6
   (verdetti Kerr) vanno ri-mappate in funzione di E_eff.

Nota di coerenza: il criterio di inversione del plunge di
`risultati_paper_summary.txt` §4, `(f_τ/E)²(H₀τ/H₀t)((J−A_φ)/J)² > 1`,
coincide esattamente con il nostro R10 `ρ(r*)² = H_t0/H_τ0` ad A=1.

### R11b. Verifica del trasferimento (`thakurta_kerr_K_expansion.py`)
H 3D off-equatoriale costruita; `K_NLO(E→E_eff)` testata:

- **bracket**: `{K_NLO, H}` resta O(a²) anche a E_eff<1 (rapporto di
  scala 4.00 esatto) — coerente con la non-integrabilità del summary
  Kerr §3: l'NLO riduce il drift lungo le traiettorie, non l'ordine
  puntuale del bracket. Il trasferimento per sostituzione è formalmente
  esatto (riscalamento conforme, A costante).
- **drift test 3D** (24 traiettorie interne r₀∈[4,7], stesso protocollo
  nei due regimi): Kerr di controllo (E=1.2, a=0.4): NLO vince 19/24;
  **Thakurta-Kerr (E_eff=0.96, A=1.25): NLO vince solo 9/24, ratio
  0.999 — la correzione NLO PERDE efficacia nel regime conforme
  E_eff<1** (verdetto tipo ramo-t: meglio Q_std nuda). Plausibile:
  `N₂ = (E²−1)r² + 4Mr − 4M²` cambia struttura di segno per E<1.
- **polo al muro verificato**: f₂ cresce di 100× avvicinandosi a
  `r_w = 2M/(1−E_eff²) = 25.51` (E_eff=0.96).

Aperto: portare il numtest ORIGINALE di Kerr (stessi ensemble/metriche
di drift) nel regime E_eff<1 per confronto 1:1 col loro 0.916; capire
se esiste una correzione NLO adattata al regime conforme (ri-ottimizzare
l'ansatz con il polo al muro).

### R11c. Identità strutturale: dove vale e dove si rompe
**A costante — TEOREMA** (2 righe, off-equatoriale): `u_TK = u_K/A` ⇒
vincolo Ê su ∂_η ≡ vincolo E_eff = Ê/A su ∂_t; `dτ_TK = A dτ_K` ⇒
`T_τ^{TK} = A·T_τ^K(E_eff)` come FUNZIONALI ⇒ `H_τ^{TK} = A⁻²H_τ^K(E_eff)`
e (ramo t, senza prefattore) `n²_{TK} = n²_K(E_eff)` identicamente.
La costruzione perturbativa commuta con la sostituzione: le
quasi-costanti NON possono essere strutturalmente diverse.

**A(η) dinamico — si rompe all'ordine misto a²ε (ε = A′/A)**:
`dK/dλ = {K,H} + ∂K/∂η`; `Q₀` non dipende da E ⇒ niente controtermine
a O(ε); il primo termine nuovo:
`∂K_NLO/∂η = ε·a²(2E_eff²cos²θ − E_eff·∂_E f₂·p_θ²) + …`.
Il controtermine `h` risolve l'equazione di trasporto
`{h, H₀} = −∂K/∂η|_{a²ε}` lungo il flusso imperturbato: struttura
generica FUORI dalla famiglia Kerr (integrali di flusso + polo al muro):

```
K^{TK} = K^{Kerr}(E_eff(η)) + (A′/A)·a²·h(r,θ,p) + O(a²ε², a⁴ε)
```

`h` = il primo oggetto del programma quasi-costanti senza analogo in
Kerr; probabile spiegazione del verdetto R11b.

### R11d. I conti per A(η) non costante (`thakurta_kerr_dynamic_K.py`)
**Hamiltoniana PMP 3D dinamica** (nuova, off-equatoriale, Σ generica):

```
H_τ = p̃_φφ′₀ + R·√( (Δ/Σ)p_r² + p_θ²/Σ + p̃_φ²/Ḡ ) − A²f_Σ/Ê
f_Σ = 1−2Mr/Σ,  b = 2Mar sin²θ/Σ,  Ḡ = G + A²b²/Ê²,
φ′₀ = b·v̄_Σ²/Ḡ,  p̃_φ = J − A²b/Ê,  R² = v̄_Σ²·Δsin²θ/Ḡ
```

con l'identità chiave `f_Σ·G + b² = Δsin²θ` (det del blocco t-φ, sympy
esatta). Riduzione equatoriale ≡ R7 (0.0 su punti casuali).

**Verifiche dinamiche** (`a(η) = −1/(H_cη)`):
- a=0: `L² = p_θ² + J²/sin²θ` conservata a **6e-13** lungo il flusso non
  autonomo (Noether spaziale esatto — validazione della macchina 3D);
- gerarchia dei candidati (a=0.4, A: 1→1.25): `K(E_eff(η))` running
  batte nettamente `K(E congelata)` (6.8e-3 vs 9.6e-3): **la
  sostituzione running è il controtermine di ordine zero**; NLO ≈ Q_std
  nel dinamico (coerente con R11b).

**Sorgente del controtermine (forma compatta, verificata sympy)**:
usando `DE₀ − N₂ = −2M·D₀/r` ⇒ `∂_E f₂ = −2MEcos2θ/(r·DE₀²)`:

```
dK/dη|_espl = ε·a²·S₁ ,   ε = A′/A
S₁ = 2E_eff²·[ cos²θ + M·cos2θ·p_θ²/(r·DE₀²) ]
```

`DE₀ = r²w_eff` ⇒ **polo doppio sul muro di congelamento**: il
controtermine `h` (soluzione di `{h, H₀^{Thak-Schw}} = −S₁`) eredita il
polo — conferma strutturale che h è fuori dalla famiglia Kerr.

### R11e. Trasporto risolto: medie chiuse, controtermine, limiti fisici
(`thakurta_kerr_h_counterterm.py`)

**Medie angolari in forma chiusa** (piano orbitale, `cosθ = sin i·sinψ`,
`k = 1−J²/L²`) — verificate per quadratura a **1e-16**:

```
⟨cos²θ⟩ = k/2 ,  ⟨1/sin²θ⟩ = L/|J| ,  ⟨p_θ²⟩ = L(L−|J|) ,
⟨cos2θ·p_θ²⟩ = |J|(L−|J|)
⟹ ⟨S₁⟩ = E_eff²(1−J²/L²) + 2ME_eff²|J|(L−|J|)/(r·DE₀²)
```

**Risultato strutturale (estende R11)**: per `E_eff < 1` il ramo τ non
ha orbite legate NEMMENO ad A costante: `W_eff(r)` è monotona
decrescente verso il muro ⇒ **il muro è attrattore globale** (endpoint
in η finito con `p_r → ∞`, singolarità di coordinate della riduzione).
Le quasi-costanti dinamiche vivono quindi su ARCHI (scattering,
E_eff>1) o richiedono il ramo t.

**Controtermine e verdetto numerico** (archi di scattering, a=0.15,
ε≈0.002): `K_comp = K(E_eff(η)) − a²∫εS₁dη` rimuove il drift esplicito
*esattamente per costruzione*; su archi realistici però il drift
conforme è **subdominante al residuo Kerr end-to-end**
(6.1e-3 vs 6.7e-3 statico: stesso ordine). La legge secolare ⟨S₁⟩ vale
come media di toro: su archi corti (1–3 oscillazioni polari) il
rapporto esatto/previsione è ~0.50 (spread 0.36–0.83) — mixing
incompleto, fattore atteso.

**Drift map (θ, r)** (`genera_drift_map.py` →
`fig_thakurta_kerr_drift_map`): drift puntuale relativo
`|{K,H₃D} + εa²S₁|/|K|` su griglia (θ, r-log), s=0.9, A=1.25,
ε=0.025, (p_r, p_θ, J) = (−0.3, 1.0, 1.2). Struttura: banda rossa sul
polo r=2M (D₀ di f₂) a ridosso dell'ergosfera; bacino verde
quasi-equatoriale a r moderati; asimmetria nord-sud da p_θ>0; filamenti
verdi = luoghi di cancellazione; regione oltre il muro mascherata.
Range: 1e-9 → 2e2, mediana 7e-3. Variante nel piano **(A, r)** a θ=60°
(`genera_drift_map_A.py` → `fig_thakurta_kerr_drift_map_A`): il muro
`r_w(A)` cala da destra (bordo maschera ≡ curva analitica, check
gratuito), polo r=2M in basso, bacino verde centrale compresso tra i
due poli al crescere di A — il ritratto del "raggio di convergenza
mangiato dal congelamento" (R11a).

**Sintesi dei conti per A(η) non costante**: la quasi-costante dinamica
è `K^{TK}(η) = K^{Kerr}(E_eff(η)) − a²∫ε S₁ dη` con `S₁` chiusa (R11d) e
medie chiuse; i suoi limiti di impiego sono fisici, non tecnici: il
regime legato E_eff<1 non esiste (muro attrattore) e sugli archi il
termine conforme è una correzione fine al residuo Kerr.

## R12. Struttura algebrica delle orbite equatoriali (integrali ellittici?)
Radicando razionalizzato di `dφ/dr` (verificato sympy):

| oggetto | grado radicando | genere | forma chiusa |
|---|---|---|---|
| geodetiche equatoriali | 4 | 1 | ellittica (classica) |
| brachistocrona τ (Schw. E Kerr) | 5, squarefree | **2** | theta di genere 2 (Hackmann–Lämmerzahl) |
| brachistocrona t | 7 | **3** | oltre |
| **separatrice J = J_c = a/E** | 4 (fattorizzazione doranTau `Δ−J_c²w = f(r²+a²/E²)`, residuo 0) | 1 | **ellittica (Weierstrass ℘)** |
| limite nullo | 4 | 1 | ellittica (fotoni) |
| soglie (radici doppie) | degenere | 0 | elementare (coerente con doranT) |

Il fattore di rotaia `√(wf)` alza il genere di 1 rispetto alle
geodetiche; il 𝒦(r) razionale del ramo t di un altro. L'unica
brachistocrona classicamente integrabile è la separatrice penetrante.
Conforme (A costante): `E → Ê/A` ⇒ isola ellittica a `J_c(A) = sA²/Ê`,
attraversata dinamicamente quando `A(η) = √(JÊ/s)`. I tempi T_τ, T_t
condividono il radicando ⇒ stessa classificazione.

### R12a. φ(r) della separatrice in forma chiusa (Weierstrass)
(`../KerrMetric/kerr_separatrix_weierstrass.py`, validazione
`kerr_separatrix_validation.py`: riduzione 5.6e-17, identità ℘ 9e-28,
antiderivata WW 20.53 esatta a 30 cifre)

Sorpresa strutturale: `P₂ = r²+c²` DIVIDE Q₄ ⇒ i residui ai poli
complessi sono nulli: restano solo i poli agli orizzonti r±. Con
`c = a/E`, `Q₄ = r((E²−1)r+2M)(r²+c²)`:

```
A = Ma²/(2E²) ,  B = a²(E²−1)/(12E²)
g₂ = a⁴(E²−1)²/(12E⁴) − a²M²/E²
g₃ = a⁴[36E²M²(1−E²) − a²(E²−1)³]/(216E⁶)

z(r) = ℘⁻¹(A/r + B; g₂, g₃)      [= ∫₀^r dr′/√Q₄]

α± = r±((E²−1)r± + 2M)/(r± − r∓)
λ± = α±/√Q₄(r±) ,   v± = z(r±) ,  ℘(v±) = A/r± + B
Λ₀ = (E²−1) − α₊/r₊ − α₋/r₋ + 2λ₊ζ(v₊) + 2λ₋ζ(v₋)

φ(r) = (a/E)[ Λ₀·z(r) + λ₊ ln(σ(z−v₊)/σ(z+v₊))
                       + λ₋ ln(σ(z−v₋)/σ(z+v₋)) ] + cost
```

Lineare in z + due rapporti di sigma di Weierstrass (equivalente
Legendre: 1 F + 2 Π con caratteristiche agli orizzonti). Conforme:
E → Ê/A(η) in tutti i parametri.

### R12b. Dominio di validità e coordinate di Doran per Thakurta-Kerr
- **La formula R12a vale ATTRAVERSO l'ergosfera**: la cancellazione
  doranTau `Δ−J_c²w = f(r²+c²)` elimina ogni traccia di r_e dal
  radicando (`Q₄(2M) = 4M²(4E²M²+a²) > 0`) e dai poli (che stanno agli
  orizzonti r±, dentro l'ergosfera). È esattamente ciò che rende
  speciale la separatrice. Le orbite sub-critiche (J<J_c) invece hanno
  la cuspide a r_e, dove la RIDUZIONE ottica degenera (worldline
  regolare, come in doranTau §5).
- **Il breakdown è all'ORIZZONTE**, ed è l'artefatto BL noto: i termini
  di terza specie divergono log a r→r₊ (avvolgimento azimutale
  infinito), come per le geodetiche in BL.
- **Doran per Thakurta-Kerr: esiste per sollevamento conforme.**
  `g_TK = A(t_D − F(r))²·g_Kerr-Doran` con la trasformazione di Doran
  del Kerr base (t_D = η + F(r), φ_D = φ + G(r)): il CKV `∂_η ≡ ∂_{t_D}`
  è lo STESSO campo vettoriale (invariato dal cambio di carta), il
  vincolo di rotaia è coordinate-free, e la metrica è regolare
  all'orizzonte; l'unico prezzo è A = A(t_D − F(r)) dipendente anche
  da r nella nuova carta. La separatrice in Doran aggiunge lo shift
  `∫aβ/Δ dr` con radicando CUBICO `2Mr(r²+a²)` (verificato): un secondo
  sistema ellittico (genere 1) — φ_D resta in forma chiusa, combinazione
  di due curve ellittiche (Q₄ e la cubica di Doran).

### R12c. Shift di Doran in forma chiusa: φ_D completa
(`../KerrMetric/kerr_doran_shift_weierstrass.py`; V1-V2 simboliche
esatte, V3 riduzione a 1.1e-16, antiderivata σ/ζ già validata a 30
cifre)

`aβ/Δ·dr = 2Mar·dr/(Δ√C₃)` con `C₃ = 2Mr(r²+a²)` (il fattore (r²+a²)
si cancella: poli solo agli orizzonti). Seconda curva ellittica in
forma di Weierstrass PULITA:

```
r(z) = (2/M)·℘(z; g₂, g₃) ,   g₂ = −M²a² ,  g₃ = 0   (quasi-lemniscatico)
z(r) = ℘⁻¹(Mr/2)
agli orizzonti: r±²+a² = 2Mr± ⇒ √C₃(r±) = 2Mr±, ℘′(v±) = ±M²r±

φ_shift(r) = [a/(r₊−r₋)]·( B₊(z) − B₋(z) )
B_k(z) = 2ζ(v_k)z + ln[σ(z−v_k)/σ(z+v_k)] ,   ℘(v±) = Mr±/2
```

**φ_D(r) = φ_BL(r) [R12a, curva Q₄] + φ_shift(r) [cubica C₃]**: la
separatrice di Doran in forma chiusa su DUE curve ellittiche; le
divergenze log dei due pezzi a r₊ si cancellano (contenuto analitico
della regolarità Doran all'orizzonte). Conforme: E → Ê/A solo nella
parte BL (lo shift non contiene E).

### R12d. Verifica end-to-end della forma chiusa e limite PMP all'ergosfera
(`../KerrMetric/kerr_separatrix_trajectories.py`,
`fig_separatrix_3traiettorie`)

Tre traiettorie della separatrice, stesse condizioni iniziali (r₀=9,
a=0.9, E=1.2), fin dentro l'ergosfera:
- **T3 Weierstrass valutata davvero** (σ,ζ,℘ via theta di Jacobi sul
  reticolo VERO: rombico, τ = 1/2 + 0.4639i, ω₁ = 3.1305; residuo
  identità ℘ 8.4e-15, uniformizzazione 1.5e-11) vs **T2 quadratura**:
  **max |Δφ| = 2.8e-10 attraverso l'ergosfera** — formula verificata
  end-to-end, funzioni speciali incluse;
- **T1 ODE di Hamilton**: segue le altre a ~1e-4/1e-5 e sfiora
  `r_min = 2.000000` (svolta marginale ESATTA a r_e) — non può
  attraversare: **risultato nuovo**: dentro l'ergoregione `H(p) > 0`
  per ogni p (il costo `𝔉 = f + bφ′` cambia segno sul ramo spurio
  dell'ellisse): la descrizione PMP-η del ramo τ TERMINA a r_e — la
  versione hamiltoniana esatta della "singolarità 1/F intrinseca alla
  riduzione" di doranTau. Dentro, la forma chiusa segue l'orbita del
  worldline vincolato (ancorata dalla validazione Doran di doranTau,
  max|Δφ| = 4.4e-6).

### R12e. Galleria: la sovrapposizione a 3 metodi è generale, non speciale
(`../KerrMetric/kerr_separatrix_gallery.py`, `fig_separatrix_gallery`)

Il test R12d ripetuto su QUATTRO parametri (a,E), ciascuno col SUO
reticolo (g₂,g₃,τ diversi). Ricerca automatica del tipo di reticolo
(rombico Re τ=1/2 / rettangolare Re τ=0) con validazione via
uniformizzazione ℘(z(r))=A/r+B e minimizzazione golden-section del
residuo ℘.

| a | E | reticolo τ | residuo ℘ | unif. | \|W−quad\| | graze | \|ODE−q\| |
|---|---|---|---|---|---|---|---|
| 0.9 | 1.20 | 0.5+0.4639i | 6.6e-14 | 1.5e-11 | 2.8e-10 | 2.000000 | 1.5e-4 |
| 0.5 | 1.20 | 0.5+0.4795i | 2.7e-15 | 2.5e-11 | 1.2e-10 | 2.000000 | 2.6e-4 |
| 0.9 | 1.05 | 0.5+0.4901i | 7.0e-14 | 2.1e-11 | 2.5e-10 | 2.000000 | 1.5e-4 |
| 0.7 | 1.50 | 0.5+0.4387i | 3.8e-15 | 6.0e-11 | 2.8e-10 | 2.000000 | 2.2e-4 |

- **Weierstrass valutata ≡ quadratura a ~1e-10 in tutti i casi**: la
  forma chiusa σ/ζ (R12a) non è artefatto di un reticolo particolare —
  regge su reticoli diversi. La superposizione ODE/quadratura/℘ NON è
  un caso fortunato di (0.9,1.2).
- **Reticolo sempre rombico** (Re τ=1/2): per E>1 il discriminante è
  <0 in tutti e quattro. La separatrice τ vive nel regime rombico; il
  ramo rettangolare esiste nel codice ma non è selezionato qui.
- **ODE grazia sempre r=2.000000M = r_e** (ergosfera equatoriale): la
  terminazione PMP-η all'ergosfera (R12d) si conferma su ogni caso —
  fenomeno strutturale, non parametrico. Il residuo \|ODE−quad\|~2e-4
  è tutto nella regione di svolta vicino r_e; fuori è a livello
  quadratura.

### R12f. Tricotomia retrograda: J negativo minimo per la cattura
(`tricotomia_equatoriale.py`, `tricotomia_figura.py`,
`fig_thakurta_tricotomia_Jneg`)

Banda di cattura (l'orbita entrante RAGGIUNGE l'ergosfera `r_e=2M`;
per R12d il ramo τ vi termina) al variare di J, incluso J<0, per i due
rami equatoriali. M=1, a=0.9, Ê=1.2, r0=8.

- **Ramo τ**: banda **SIMMETRICA** `J ∈ [−sA²/Ê, +sA²/Ê]`. Bordi:
  A=1 → `±0.7511` (vs separatrice chiusa `J_c=sA²/Ê=0.75`, diff 1e-3);
  A=1.1 → `±0.9088` (vs 0.9075). La simmetria non è accidentale: la
  fattorizzazione doranTau `Δ − J_c²w = f(r²+c²)` dipende da **`J²`**,
  quindi prograde e retrogrado hanno la STESSA soglia nonostante il
  trascinamento. ⇒ **`J_neg^τ = −sA²/Ê`**.
- **Ramo t (η)**: banda **larga e ASIMMETRICA**. Cattura fino a
  `J≈−8` (retrogrado forte); prograde-catturante oltre `+3`. Bordi
  negativi: A=1 → `−8.05`, A=1.1 → `−8.37`. Sotto `J≈−15` l'orbita non
  parte (nessuna radice entrante).

Fisica: il ramo t è enormemente più *capture-prone* (coerente con "t
affonda di più", cfr. Vaidya R13/R15). Un'orbita retrograda con `J=−5`
è **catturata da t ma diffusa da τ**. Il trascinamento rompe la
simmetria solo per il ramo t; la struttura `J²` di doranTau protegge la
simmetria del ramo τ. `A>1` allarga leggermente entrambe le soglie.

### R12g. La riflessione all'ergosfera è una CUSPIDE (dimostrazione)
(`cuspide_ergosfera.py`, `fig_thakurta_cuspide_ergosfera`)

Dimostrazione analitica che una brachistocrona τ equatoriale con
`|J| ≠ J_c` che raggiunge l'ergosfera vi si riflette con una **cuspide**
(non un periasse liscio). Da
```
dφ/dr = J √(w f) r / (Δ √(Δ − J² w)) ,   w = E² − f
```
vicino a `r_e=2M` (dove `f→0`, `f ≈ (r−r_e)/2M`, `w→E²`, `Δ→a²`):
```
dφ/dr ~ K √(r−r_e) ,  K = J E r_e/(a²√(a²−J²E²)√(2M))
⟹ φ − φ_e ~ (2/3) K (r−r_e)^{3/2}    (potenza 3/2 = CUSPIDE)
```
La cuspide viene da `√f→0` all'ergosfera mentre il radicando di svolta
`√(a²−J²E²)` resta finito per `|J|<J_c`. Tangente **radiale** (`dφ/dr→0`
da entrambi i lati) ⟹ cuspide, non punto angoloso.

**Tricotomia all'ergosfera** (segno di `a²−J²E²`):
- `|J|>J_c` (`a²−J²E²<0`): svolta PRIMA di `r_e`, periasse liscio,
  esponente 1/2 classico (`dφ/dr→∞`);
- `|J|=J_c` (`a²−J²E²=0`, con `J_c=a/E`): `√f` e `√(Δ−J²w)` svaniscono
  insieme ⟹ `dφ/dr` FINITO, attraversamento liscio (separatrice);
- `|J|<J_c` (`a²−J²E²>0`): `dφ/dr~√(r−r_e)` ⟹ CUSPIDE.

Verifiche numeriche: coefficiente `K` simbolico ≡ formula chiusa;
esponente di `dφ/dr` a `r_e`: J=0.5→0.4999, J=0.7→0.4997 (cuspide),
J=J_c→−0.0001 (liscio), J=0.9>J_c→svolta a r=2.141>r_e; e
`Δφ/(r−r_e)^{1.5} → 2K/3` (potenza 3/2). Thakurta-Kerr via `E=Ê/A`,
`J→J/A`, `J_c=sA²/Ê` (r_e conforme-invariante). Mette su base rigorosa
la "riflessione a cuspide" della tricotomia (§ doranTau, R12f).

**Attraversamento in coordinate di Doran** (figura, ramo verde J_c). Il
`√(wf)` va immaginario dentro l'ergosfera (`f<0`), ma sulla separatrice
`√f` si **cancella**: `dφ_BL/dr = J_c√w·r/(Δ√(r²+c²))`, `c=a/E` — regolare
attraverso `r_e`. Per regolarizzare l'ORIZZONTE si passa a Doran:
`φ_D = φ_BL − φ_shift`, `dφ_shift/dr = a√(2Mr/(r²+a²))/Δ` (R12c); a `r_+`
le due divergenze `1/Δ` si cancellano (verificato: BL−shift = +0.036
finito, BL+shift diverge). Così l'orbita `J_c` **attraversa `r_e` e
scende all'orizzonte** con winding finito — la verde compare "dall'altra
parte" dell'ergosfera (impossibile in BL, che vi si ferma).

## Aperti
1. Dinamica sulla superficie di congelamento (dopo il contatto il
   formalismo hamiltoniano degenera: serve descrizione vincolata al
   bordo — analogo del cuspide J<J_c di Kerr?).
2. Mappa di fuga: (J, Ê, η₀) → fuga vs cattura da espansione; soglia
   Ê_c(r₀, η₀).
3. Sultana–Dyer (a ∝ η², sorgente esatta) come variante fisica.
4. Confronto quantitativo con McVittie asintotico (Λ>0 → SdS).
