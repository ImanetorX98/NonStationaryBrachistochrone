# Brachistocrone in Vaidya — risultati

**Script sympy** (tutte le verifiche passano):
- `vaidya_metric_sympy.py` — indagine della metrica, sorgente, orizzonti
- `vaidya_brachistochrone_sympy.py` — EOM worldline vincolato (parametro r)
- `vaidya_brachistochrone_vparam.py` — EOM in parametro v, rimbalzo, timing
- `vaidya_penetration_map.py` — mappa di penetrazione J_c(v₀)
- `vaidya_hamiltonians_sympy.py` — hamiltoniane PMP/Zermelo dei due rami
- `vaidya_optical_metric_sympy.py` — metrica ottica via omotetia (m = μv)
- `genera_figure_vaidya.py` — figure di validazione (`Vaidyafigures/`)

Contesto: `../FLRWmetric/nonStationaryBrachi.md` (piano), `FLRWresults.md`
(caso FLRW), `../KerrMetric/doranTau.md`, `doranT.md` (caso stazionario).

---

## 1. La metrica e cosa descrive

**Vaidya entrante** (coordinate avanzate di Eddington–Finkelstein, c=G=1):

```
ds² = −(1 − 2m(v)/r) dv² + 2 dv dr + r² dΩ²
```

Verificato: `R = 0`, unica componente `G_vv = 2m′(v)/r²`, cioè

```
T_μν = [m′(v)/(4π r²)] l_μ l_ν ,   l_μ = −(dv)_μ ,  l·l = 0,  ∇_l l = 0
```

**polvere nulla radiale entrante** (pioggia sferica di radiazione).

| variante | energia | descrive |
|---|---|---|
| entrante, `m′>0` | NEC/WEC ok (`ρ ∝ m′`) | **BH che accresce** radiazione |
| entrante, `m′<0` | viola NEC | **toy model evaporazione di Hawking** |
| uscente (`u`, `m′(u)<0`) | radiazione positiva uscente | **stella che irradia** (Vaidya 1951) |

Proprietà chiave (verificate):
- `m′=0` ⇒ vuoto ⇒ Schwarzschild in EF avanzate; coordinate regolari all'orizzonte.
- `(£_∂v g)_vv = 2m′/r ≠ 0`: **niente Killing tempo** — il regime non stazionario voluto.
- **Orizzonte apparente**: espansione uscente `θ₊ = (r−2m(v))/r²` ⇒ `r_AH = 2m(v)`,
  velocità `ṙ_AH = 2m′(v)`.
- **Orizzonte eventi ≠ apparente** (numerico, `m = 1+0.02v` su `v∈[0,50]` poi costante):
  durante l'accrescimento `r_EH > r_AH` (gap ~0.2–0.27), poi convergono al
  Schwarzschild finale. L'EH è teleologico: si gonfia *prima* che la massa cada.

---

## 2. EOM della brachistocrona (worldline vincolato non autonomo)

**Vincolo rotaia**: `−u_v = E` fissata. `∂_v` non è Killing: la rotaia impone
il vincolo con forza a componente `f_v ≠ 0` (analogo di `f_η = aH` in FLRW).

**Riduzione** (parametro `λ = r`, tratto a r crescente, `′ = d/dr`).
Dal vincolo `Λ = dτ/dr = (f v′ − 1)/E`; eliminando Λ:

```
φ′² = W(v, v′, r) = [ f v′² − 2v′ − (f v′ − 1)²/E² ] / r²

ramo τ:  F_τ = (f v′ − 1)/E − J √W        (T_τ = ∫ Λ dr)
ramo t:  F_t = v′ − J √W                   (T_v = ∫ v′ dr = Δv)
```

con `J` = momento di Fermat (moltiplicatore per Δφ). EOM: EL in `v(r)`,
secondo ordine, non autonoma.

### R1. Cosa si conserva e cosa no
- **`J` conservato**: la sfericità sopravvive (φ ciclica). Si perde `E`
  (che infatti è imposta dalla rotaia), non `J` — situazione speculare a
  quella attesa: il Killing perso è quello temporale.
- **Forzante ∝ m′ (esatto)**: `F` dipende da `v` solo via `m(v)` ⇒
  `∂F/∂v = (∂F/∂m)·m′(v)` (chain rule; `∂F/∂m` in forma chiusa negli script).

### R2. p_v e la teleologia della rotaia
`p_v = ∂F/∂v′`. Arrivo spaziale con istante `v` libero ⇒ trasversalità
`p_v(r_arrivo) = 0`.

- Stazionario (`m′=0`): `p_v` conservato ⇒ `p_v ≡ 0` sull'orbita.
- Non stazionario: `dp_v/dr ∝ m′` ⇒ `p_v(r) ≠ 0` lungo il percorso pur con
  `p_v(r₁)=0`: la brachistocrona **anticipa l'accrescimento futuro**
  (stessa teleologia dell'orizzonte degli eventi, qui nella dinamica).

### R3. Validazione sul limite Schwarzschild
Con `m` costante e `p_v = 0`, `dφ/dr` riproduce le forme chiuse Kerr a=0
(match numerico a 12 cifre, entrambi i rami):

```
dφ/dr = 𝒦 r √(w f) / ( Δ √(Δ − 𝒦² w) ) ,   Δ = r²f ,  w = E² − f
𝒦_τ = J           (ramo τ)
𝒦_t = f J / E     (ramo t)
```

### R4. Demo numerica (Vaidya lineare)
Salita `r: 3.5 → 10`, arrivo a `v = 40`, `E=1.2`, `J=1.3`,
`m(v) = 1 + 0.01v` vs Schwarzschild `m=1` (integrazione all'indietro
dall'arrivo con `p_v(r₁)=0`; il periasse dell'orbita è a `r≈3.01`, dove la
parametrizzazione in `r` degenera — la partenza sta sopra):

| | Schwarzschild | Vaidya μ=0.01 |
|---|---|---|
| v partenza | 16.103 | 13.044 |
| Δφ | 0.4014 | 0.5293 |
| T_v = Δv | 23.897 | 26.956 |
| T_τ | 7.757 | 7.480 |
| p_v partenza | 0 (esatto) | **0.0188** |

Figure di validazione in `Vaidyafigures/` (`genera_figure_vaidya.py`):
orizzonti EH/AH e raggi autosimilari, residui vs forma chiusa Kerr a=0,
orbite e `p_v(r)`, minimo variazionale diretto di `T_τ` a `Δφ` fisso.

---

## 3. Metrica ottica: si può?

**Vaidya generica: NO.** Serve un Killing tempo (Perlick 1991); non c'è, e
non c'è nemmeno un Killing conforme: `(£_ξ g − 2g)_vv = 2(v m′ − m)/r` con
`ξ = v∂_v + r∂_r`.

**Vaidya lineare `m = μv`: SÌ, via omotetia.** Il residuo si annulla ⇔
`m = μv`: `ξ` è **Killing omotetico** (`£_ξ g = 2g`). Coordinate
autosimilari `s = ln v`, `x = r/v`:

```
g = e^{2s} · ĝ(x) ,    ĝ = (2x − f) ds² + 2 ds dx + x² dφ² ,   f = 1 − 2μ/x
```

`∂_s` è Killing **genuino** di `ĝ` — stessa struttura di FLRW (`g = a² ĝ`),
con `e^s` al posto di `a`.

### R5. Metrica ottica di Randers (brachistocrona nulla in v)
Per i raggi nulli su `ĝ` (`A = −ĝ_ss = f − 2x`):

```
ds_arrivo = [ dx + √( dx² + A x² dφ² ) ] / A
⇒  F = α + β ,   α² = dx²/A² + x² dφ²/A ,   β = dx/A
```

(verificata contro la condizione nulla di `ĝ`). Minimizzare
`s_arrivo = ln v_arrivo` = minimizzare il tempo d'arrivo `v`. Poiché
`v = e^s` è monotono, anche il **ramo t massivo** con rotaia omotetica
`Ê = −u·ξ` si riduce a problema **statico** su `ĝ`. Il ramo τ no:
`∫dτ = ∫e^s dτ̂` non è invariante.

### R6. Superfici di luce = orizzonti autosimilari, accrescimento critico
`A = 0 ⇔ 2x² − x + 2μ = 0`:

```
x± = (1 ± √(1 − 16μ)) / 4 ,     reali ⇔ μ ≤ μ_c = 1/16
```

- `x₋` è l'**orizzonte degli eventi** autosimilare (`r_EH = x₋·v`;
  ipersuperficie nulla `x=cost` ⇔ stessa quadratica; bisezione numerica:
  match a 8 cifre per μ=0.02: `x₋ = 0.04384472`).
- `x_AH = 2μ < x₋`: AH dentro EH, come dev'essere in accrescimento.
- `x₊` = superficie esterna di fuga (attrattore dei raggi uscenti).
- `μ_c = 1/16` è il valore critico noto della letteratura autosimilare
  (Papapetrou / singolarità nude): per `μ > 1/16` niente superfici `x±`.

La singolarità `1/F` della riduzione ottica cade sugli **orizzonti
omotetici** `x±` — l'analogo della coppia superficie di luce/ergosfera di
Kerr, ma qui la superficie interna È l'orizzonte degli eventi.

---

## 3b. Traiettorie in parametro v e finestra temporale
(`vaidya_brachistochrone_vparam.py`, figure `fig_vaidya_bounce`,
`fig_vaidya_timing`)

La parametrizzazione in `r` degenera al periasse (`dv/dr → ∞`). Con `v`
come parametro (`′ = d/dv`):

```
Λ = (f − r′)/E ,   φ′² = W_v = [ f − 2r′ − (f−r′)²/E² ] / r²
F_τ = (f − r′)/E − J√W_v
```

Stesso funzionale (cambio di variabile esatto: `(fv′−1)dr = (f−r′)dv`),
stesse curve estremali, ma EL regolare a `r′ = 0`: si integra
**attraverso il rimbalzo**. Trasversalità presa dalla forma r-param
(`p_v(r₁)=0` ⇒ `r′(v₁) = 1/v′(r₁)`).

**Validazioni (E=1.2, J=1.3, r₁=10):**
- cross-check r-param al passaggio r=3.5: `v`, `Δφ`, `T_τ` coincidono a
  ~1e-7 (Schwarzschild e Vaidya μ=0.01);
- statico: `r_min = 2.72713516` = radice esatta di `r²f − J²(E²−f) = 0`
  a **9e-12**;
- Vaidya μ=0.01: periasse `r_min = 3.010492` a `v = 7.29` — esattamente
  dove la parametrizzazione in r moriva (3.0104915).

**R7. Finestra temporale (scan dell'istante d'arrivo v₁ ∈ [25,70]):**

| μ | r_min(v₁=25) | r_min(v₁=70) | r_min/2m(v_peri) |
|---|---|---|---|
| −0.01 (evapora) | 2.68 | **1.89** | 1.27 → **1.63** (si allontana) |
| 0 (statico) | 2.727 (piatto) | 2.727 | 1.36 costante |
| +0.01 (accresce) | 2.76 | **3.51** | 1.48 → **1.30** (si avvicina) |

Il periasse dipende da *quando* si arriva — fenomeno genuinamente non
stazionario, assente in Kerr. **Controintuitivo**: in evaporazione il
periasse assoluto scende ma l'orizzonte si ritira più in fretta ⇒ la
penetrazione *relativa* peggiora; in accrescimento il periasse sale ma
l'orizzonte cresce più in fretta ⇒ arrivi tardivi sfiorano l'orizzonte
relativo. La "corsa" orbita-orizzonte ha vincitori opposti in assoluto
e in relativo.

## 3c. Mappa di penetrazione (J, v₀)
(`vaidya_penetration_map.py`, figura `fig_vaidya_penetration_map`)

Setup: lancio entrante da `r₀=10` all'istante `v₀`, mira localmente ottima
(`p_v=0` a massa congelata, risolta in forma chiusa; riflessione temporale
`u_in = 2/f − u_out` — le radici radiali di W sommano a `2/f`). Esiti:
riflette (`r′=0` fuori) vs cattura (attraversa `r=2m(v)`). Soglia
`J_c = sup{J : cattura}` per bisezione. Integrazione DOP853 rtol 1e-12
(RK45 1e-10 NON basta per le svolte radenti: r_min sbagliato di 4e-4).

**Validazione statica**: r_min riproduce `2m + J²E²/(2m)` a tutti i J
testati (0.05→2.001797 vs 2.001800; 0.005→2.000019 vs 2.000018);
`J_c = 0` — penetrazione a misura nulla, come da tricotomia Kerr a=0.

**R8. Risultato principale — l'accrescimento apre una finestra finita:**

| caso | J_c(v₀=0) | J_c(v₀=60) |
|---|---|---|
| statico | 0 (misura nulla) | 0 |
| evapora μ=−0.01 | 0 | 0 |
| **accresce μ=+0.01** | **0.919** | **1.286** (cresce con m(v₀)) |

- In accrescimento il ramo τ acquista un **intervallo finito di cattura**
  `J ∈ (0, J_c]` — la struttura che in Kerr era esclusiva del ramo t
  (dicotomia con soglia) emerge qui dalla **dinamica**: le orbite radenti
  indugiano vicino al periasse e l'orizzonte in crescita chiude il gap
  (inseguimento, vedi figura: due orbite quasi identiche a cavallo di J_c
  si separano solo nella fase radente).
- In evaporazione l'orizzonte si ritira più in fretta di quanto l'orbita
  si avvicini: `J_c = 0` come nel caso statico (asimmetria
  accrescimento/evaporazione, coerente con R7).
- Scala empirica: `J_c ≈ (8÷10)·√μ` per μ ∈ [0.0025, 0.02] (v₀=0):
  coerente con la stima `δ_p = J²E²/2m ~ crescita dell'orizzonte durante
  la fase radente`.
- Caveat tecnico: per `J ≲ 0.01` la riduzione degenera
  (`∂²F/∂r′² ∝ J`): EL stiff, esiti non affidabili — le sonde partono
  da J=0.01.

## 3d. Hamiltoniane dei rami t e τ (`vaidya_hamiltonians_sympy.py`)

Il vincolo `−u_v = E` rende l'indicatrice (velocità ammesse a (v,r) fisso)
un'**ellisse** — identità esatta, `w = E²−f`:

```
r′ = (f − E²) + E√w cosθ ,     φ′ = (√w/r) sinθ
```

Navigazione di Zermelo: **vento radiale entrante** `f − E²` (< 0 sempre per
E>1) + ellisse con semiassi `E√w` (radiale) e `√w/r` (angolare). Il max di
Pontryagin su θ dà le Hamiltoniane in **forma chiusa** (p_φ = J):

```
H_v  = p_r(f − E²) − 1 + √w · √( E²p_r² + J²/r² )
H_τ  = p_r(f − E²) − E + √w · √( (Ep_r + 1)² + J²/r² )
```

(ramo τ = ramo v con `p_r → p_r + 1/E` e costo E: il costo dτ/dv dipende
dal controllo). Proprietà:
- `p_φ = J` conservato; `dp_r/dv = −∂H/∂r ≠ 0` (inomogeneità radiale);
- `dH/dv = ∂H/∂v ∝ m′(v)` — il forzante non autonomo in veste hamiltoniana;
- trasversalità: `H = 0` all'arrivo; statico ⇒ `H ≡ 0` (H = −p_v).

**Validazioni**: flusso di Hamilton H_τ riproduce il bounce EL v-param
esattamente (r_min, Δφ, T_τ, v_start a 6 cifre, μ=0 e 0.01); statico
`max|H| = 4e-12`; dinamico `H(partenza) = ∫∂H/∂v dv` a 6 cifre (−0.062126).

**R9 — auto-sintonizzazione del ramo v, generalizzata dinamicamente**:
con gli stessi (E=1.2, J=1.3, arrivo v=40) il ramo v penetra fino a
`r_min = 2.000000` (statico: sfiora ESATTAMENTE l'orizzonte — il
meccanismo doranT `𝒦(r_e) = fJ/E → 0` sopravvive ad a=0) e in
accrescimento fino a `r_min = 1.04` = sfiorare `2m(v)` dell'epoca del
periasse (raggiunto a v presto, quando il BH era piccolo). Il ramo v
"si sintonizza" sull'orizzonte **della sua epoca**.

## 3e. Forme ottiche non autonome per m(v) generica
(`vaidya_optical_nonautonomous_sympy.py`)

Invertendo l'indicatrice (Legendre-inversa delle Hamiltoniane §3d), per
**ogni** m(v) — con `f = f(v,r)`, `w = E²−f`:

```
dτ = √( dr² + f r²dφ² ) / √w                      (Riemann–Jacobi)
dv = dr/f + (E/(f√w)) √( dr² + f r²dφ² )          (Randers: deriva+norma)
```

Verificate contro il vincolo a ~1e-14; con `m′=0` sono **esattamente le
due metriche di Perlick 1991** in EF (τ: identica a `L_τ = √((B²+Af)/w)`
di doranTau con B=1; v: EL riproduce doranT a=0 a 10 cifre): Perlick
generalizzato *verbatim* con `f → f(v,r)`.

**Il punto**: le forme esistono sempre, ma senza CKV nulla è rigido —
deriva e norma respirano entrambe. "Ottico" solo in senso **non
autonomo** (duale lagrangiano delle Hamiltoniane), nessuna geometria
statica sottostante (tranne `m = μv`, omotetia, R5-R6). Gerarchia:

| caso | struttura ottica |
|---|---|
| FLRW | indice separabile → orologio, geometria banale |
| Thakurta-Schw | geometria statica × indice n(η,r) |
| Thakurta-Kerr | geometria statica + vento **rigido** × indice |
| Vaidya μv | geometria statica sul quoziente autosimilare |
| Vaidya m(v) | tutto respira: solo senso debole |

## 3f. Il W della rotaia è il vettore di Kodama
(`vaidya_kodama_sympy.py`)

`K^a = ε^{ab}∇_b r` (2-piano normale alle sfere): in EF avanzate
**K = ∂_v esattamente** — il nostro W. Verificato:
- `K·K = −f` (timelike fuori dall'orizzonte apparente, come un Killing);
- `∇·K = 0` esatta per ogni m(v);
- **miracolo di Kodama**: `J^μ = G^{μν}K_ν = (0, 2m′/r², 0, 0)` con
  `∇·J = 0` per m(v) qualsiasi — corrente conservata SENZA Killing;
- carica associata = massa di Misner–Sharp `= m(v)`.

⇒ Il vincolo `−u_v = E` è `−u·K = E`: **energia di Kodama**. In ogni
spaziotempo sfericamente simmetrico dinamico W ha una definizione
canonica (non è una convenzione di foliazione); Kodama → Killing nel
limite statico, quindi la gerarchia di liceità del §4b del doc di avvio
si raffina: Killing (teorema) → CKV (selezione conforme) → **Kodama
(selezione quasi-locale, sferico dinamico)** → convenzione dichiarata
(nessuna struttura).

### R12. Energia di Kodama conservata lungo la brachistocrona
(`kodama_conservation.py`, `fig_kodama_conservazione`)

Conferma numerica + identità esatta di §3f:
- **[A] simbolico**: `K = ∂_v`, `∇·K = 0`, `K_a = (−f,1,0,0)`, carica di
  Misner–Sharp `= m(v)` — tutto esatto per m(v) qualsiasi.
- **[B] identità della rotaia**: dal solo vincolo `W` cade
  `−g(T,T) = (f u^v−1)²/E²` ⟹ **`−u·K = E` identicamente**. L'energia di
  Kodama È l'invariante di rotaia; non è imposta a mano.
- **[C] verifica numerica** (Vaidya accrescente `m=1+0.02v`, `E=1.2`):
  `|−u·K| = 1.200000000000` costante, `|−u·K − E| < 6.7×10⁻¹⁶`.
- **[C'] non-banalità**: per una geodetica (rotaia spenta) `−u·K`
  driftrebbe a tasso `−m′(u^v)²/r ∝ m′` (qui 0.02–0.14, cumulato 0.33):
  è la simmetria mancante che la rotaia compensa. NON è conservazione di
  Noether — è **controllo** che paga il costo `∝ m′`.

**Test discriminante** (perché Kodama e non altro): nel limite statico
solo il Kodama NON normalizzato `K=∂_v` ha costo nullo. Con `K̂=K/√f`
si avrebbe `−u·K̂ = E_t/√f`, che drifta anche in Schwarzschild statico
(dove esiste un vero Killing) — segno che è la scelta sbagliata. Quindi
`W=∂_v` è fissato dal criterio **statico ⟺ costo-zero**.

### R13. r_min(t) − r_min(τ) ha segno fisso (non rotante)
Ai punti di svolta (`p_r=0`, shell `H=0`), a `a=0`:
- τ: `(E²−f)J² = f r²`
- t: `(E²−f)J² = E² r²/f`

Rapporto locale esatto `J²_t(r)/J²_τ(r) = E²/f(r)² > 1` (esterno,
`f<1≤E`). Con `J²_branch(r)` crescenti ⟹ **`r_min^t < r_min^τ` nel
limite STATICO** (`m′=0`): il ramo t affonda di più, segno bloccato.
L'inversione di plunge di Kerr/Thakurta-Kerr è l'effetto del termine di
trascinamento `2Mas/r`, che deforma `E²/f²` fino ad attraversare 1.
Verificato su griglia `(E,J)` (Cardano/Ferrari per le radici esterne di
cubica/quartica). Il segno "t più fondo" è ROBUSTO anche in dinamica:
l'evaporazione fisica restringe il gap ma NON lo inverte (R16). Per
invertire serve la ROTAZIONE (termine `2Mas/r`) — non esiste inversione
evaporativa fisica.

### R14. Validazione diretta del principio di minimo (tre test)
Tre validazioni variazionali indipendenti, tutte per perturbazione a
estremi fissi con minimo esatto sulla soluzione:
- **FLRW** (`fig_flrw_variazionale`, R2): la retta comovente minimizza
  ENTRAMBI i funzionali t e τ (`argmin` a ε=0 per entrambi).
- **Vaidya non autonoma** (`fig_vaidya_variazionale`): la soluzione EL
  non autonoma minimizza `T_τ` a `Δφ` fisso (`T_min` vs `T_EL` a 1e-8).
- **t vs τ statico** (`verifica_minimo_brachi.py`,
  `fig_verifica_minimo_brachi`): ogni ramo, costruito INDIPENDENTEMENTE
  dalla sua ODE di Beltrami, è il minimo del SUO tempo (`argmin` a ε=0)
  e NON dell'altro (l'altro funzionale è monotòno nella perturbazione).
  Conferma variazionale della dicotomia t/τ (stessa che dà R13).

⇒ Le curve trovate sono effettivamente brachistocrone (minimi veri,
secondo ordine positivo), e t e τ sono curve genuinamente distinte.

### R15. Legge del raggio minimo di plunge t vs τ in Vaidya
(`plunge_vaidya_t_tau.py`)

Ramo t = minimizza il tempo avanzato `v` (a estremi radiali fissi
equivale a minimizzare t: `v = t + r_*`, differenza costante); riduce
al ramo t di Schwarzschild nel limite statico. Legge locale della svolta
(evento di periasse `p_r=0`), forma di R13:
```
τ:  (E²−f_*) J² = f_* r_min²
t:  (E²−f_*) J² = E² r_min²/f_* ,   f_* = 1 − 2m(v_peri)/r_min
```
con `J` momento ottico/Fermat, `E` energia di Kodama (R12).

- **Coerenza m′=0 (Schwarzschild)**: a μ=0 l'integrazione dei due rami
  RIPRODUCE le radici statiche R13 (residuo ~1e-9): `r_τ=7.7573`,
  `r_t=4.7011` (E=1.2, J=8, m=1). Limite continuo, nessuna
  discontinuità. Nota: per J piccolo il ramo t NON gira (cattura,
  coerente con "t più fondo") — serve J grande perché entrambi girino.
- **Segno robusto**: `r_min^t < r_min^τ` per ogni μ∈[−0.02,+0.02].
  L'accrescimento ALLARGA il gap (τ esce, t affonda: gap 1.62→3.06→5.07
  da μ=−0.02 a +0.02), l'evaporazione lo restringe — ma **non lo
  inverte** in regime fisico (R16: l'inversione lineare a μ≈−0.058 è
  artefatto di `m<0`; con `m(v)>0` il segno resta bloccato).
- **Correzione dinamica**: la legge NON è la statica valutata a
  `m(v_peri)` — il residuo di quella forma cresce ~lineare in m′
  (`res/|μ|` ≈ 475–495 per τ, ~costante): è la correzione non autonoma
  dal momento `p_v` teleologico (memoria `∝ m′`, §3f/R12). La svolta
  ricorda l'intera storia di accrescimento, non solo `m(v_peri)`.

Sintesi: stessa FORMA di legge, limite Schwarzschild esatto, segno
robusto, più una correzione `∝ m′` che è la firma della
non-stazionarietà.

### R16. NON esiste inversione di plunge evaporativa fisica [RETTIFICATO]
(`inversione_fisica.py`, `fig_vaidya_no_inversione_evaporazione`)
Rimpiazza la versione precedente (`inversione_evaporazione.py`,
`fig_vaidya_inversione_evaporazione`), che era ERRATA.

**Storia**: la correzione dinamica `∝ m′` di R15 restringe il gap per
`m′<0`; una prima analisi (modello lineare `m=1+μv`, r1=15, v1=40)
sembrava mostrare un'INVERSIONE a `μ_inv=−0.0577` (secondo meccanismo,
evaporativo). **È un ARTEFATTO del regime a massa negativa.**

- Con `m=1+μv` la massa va NEGATIVA per `v>1/|μ|`. A `μ_inv=−0.0577` la
  svolta era a `m≈0.67` (positiva), ma l'ancoraggio della trasversalità
  a `v1=40` era a `m=−1.31` (`m=0` a v=17.3): condizione al bordo in
  regione non fisica (singolarità nuda) che contamina la traiettoria.
- La serie perturbativa `gap(μ)=gap₀+gap₁μ+…` era asintotica **con
  raggio di convergenza `R~|μ_inv|`**: la singolarità È proprio
  l'ingresso in `m≤0`. La radice cubica (−0.053) estrapola la risposta
  fisica a piccolo μ dentro il regime `m<0`.

**Controllo (due modelli, tracciando m_min sull'orbita)**:
- **Lineare, scan 2D (v1, μ)**: OGNI caso con `gap<0` ha `m<0`; OGNI
  caso con `m>0` ovunque ha `gap>0`. Nessuna sovrapposizione ⟹ nessuna
  inversione fisica. Ancorando a `v<0` (m>1) il gap resta largamente
  positivo (+1÷+2.7).
- **Esponenziale `m=e^{−λv}` (m>0 SEMPRE)**: il gap si restringe verso
  un plateau (~0.85) ma **NON inverte mai**, nemmeno a `λ=1` (m(40)=4e-18).

**Ragione fisica**: l'inversione richiede geometria fortemente evaporata
(m piccola) alla svolta; la rampa lineare non sa fare "localmente molto
evaporato ma globalmente m>0" (continua a scendere sotto zero), e
l'esponenziale, tenendo m>0, evapora solo asintoticamente mentre la
svolta resta a m~O(1). In entrambi i casi la svolta avviene dove m è
ancora sostanziale ⟹ ordinamento statico preservato.

**Conclusione**: **l'evaporazione fisica RESTRINGE il gap ma non lo
inverte.** L'ordinamento `r_min^t < r_min^τ` (t più fondo, R13) è
robusto per tutta l'evaporazione dinamica fisica. Esiste **un SOLO
meccanismo di inversione, quello ROTAZIONALE** di Kerr/Thakurta-Kerr
(`ρ(r*)²=H_t0/H_τ0`, R12c, chiuso perché stazionario). L'affermazione
"serve rotazione" di R13/R15 è CORRETTA.

## 4. Prossimi passi

1. ~~Fenomenologia di penetrazione: soglie in J / finestre in v~~ — FATTO
   (§3b timing R7, §3c mappa R8). Resta il sotto-caso evaporante "la
   penetrante può sfuggire" (orbite lanciate dentro che escono quando
   l'orizzonte si ritira; il ramo v con r_min ≈ 2m(v_epoca) in R9 è
   l'indizio che il regime esiste).
2. Brachistocrone nel settore autosimilare: usare la riduzione ottica R5
   per il ramo t massivo con rotaia `Ê = −u·ξ`; tricotomia/dicotomia a `x±`
   in analogia con `J_c` e `J_t^c` di Kerr; raccordo con la mappa R8.
3. Limite quasi-statico `μ → 0`: espansione perturbativa attorno a
   Schwarzschild (analogo dell'espansione in `a` per Kerr); scala
   `J_c ≈ (8÷10)√μ` di R8 come primo dato da spiegare analiticamente.
4. Mappa di penetrazione per il ramo v (auto-sintonizzato, R9): definizione
   di soglia diversa (grazing marginale invece di intervallo?).
5. Articolo: FLRW (base esatta) + Vaidya (fenomenologia dinamica), sezione
   hamiltoniana unificata PMP con trasversalità H=0.
