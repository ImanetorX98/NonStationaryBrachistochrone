# Brachistocrone t e τ in FLRW piatta — risultati

**Script sympy** (tutte le verifiche passano):
- `flrw_brachistochrone_sympy.py` — worldline vincolato, EOM, forza-rotaia
- `flrw_optical_metric_sympy.py` — metrica ottica e orologio λ(η) (R8)
- `flrw_hamiltonians_sympy.py` — hamiltoniane PMP dei rami t e τ (R9)
- `genera_figure_flrw.py` — figure di validazione (`FLRWfigures/`)

Formalismo: *worldline vincolato covariante* (da `../KerrMetric/doranTau.md`),
esteso al caso non stazionario secondo la proposta di
`nonStationaryBrachi.md` §2: vincolo sulla quantità conforme.

## Setup

Metrica FLRW piatta (k=0), tempo conforme η (`dt = a dη`), moto nel piano:

```
ds² = a(η)² [ −dη² + dx² + dy² ]
```

`∂_η` è Killing **conforme** della metrica fisica. Vincolo rotaia:
`Ê = −u_η` fissata (analogo conforme di `E = −u_t` di Kerr).

Lagrangiane (worldline, parametro λ, `Λ = dτ/dλ = a√(η′²−x′²−y′²)`):

```
ramo τ:  L_τ = (1 + μÊ) Λ − μ a² η′
ramo t:  L_t = a η′ + μ (Ê Λ − a² η′)        [T_t = ∫ dt = ∫ a dη]
```

## Risultati (verificati simbolicamente)

### R1. Il vincolo fissa la cinematica, la rotaia redshifta
Dal vincolo `−u_η = Ê`:

```
γ = dt/dτ = Ê/a          (E_phys = γ ∝ 1/a: redshift della rotaia)
v(η) = √(1 − a²/Ê²)      (velocità locale fisica rispetto ai comoventi)
```

`v` dipende **solo dal tempo** (omogeneità), non dalla posizione.

### R2. Degenerazione: i rami t, τ, η coincidono come curve
`x, y` ciclici ⇒ `p_x, p_y` conservati ⇒ **retta comovente** per entrambi
i rami. Sulla soluzione:

```
dτ/dη = a²/Ê ,    dt/dη = a
```

entrambi i funzionali sono crescenti nell'η d'arrivo ⇒ la curva ottima è la
stessa (la retta, cioè la geodetica spaziale 3D). I due rami differiscono solo
nei moltiplicatori:

```
μ_τ(η) = C/(Ê² v) − 1/Ê ,    μ_t(η) = C_t/(Ê² v)
```

**Risposta alla domanda 2 del doc (caso k=0)**: t-, η- e τ-brachistocrona sono
la **stessa curva**; differiscono solo nel valore del tempo di percorrenza.
L'argomento usa solo l'omogeneità ⇒ vale per ogni FLRW (anche k=±1, con
"retta" → geodetica della 3-geometria). In FLRW la fisica della brachistocrona
è tutta **temporale**; la dinamica di forma dell'orbita (Kerr) riappare solo
rompendo l'omogeneità (→ Vaidya).

### R3. Forza-rotaia
Residuo geodetico `f_μ = du_μ/dτ − ½∂_μ g_{αβ} u^α u^β` sulla soluzione:

```
f_η = a′/a = aH     (Hubble conforme: la rotaia "sostiene" Ê contro il redshift)
f_x = −a′/(a v)
f·u = 0             (nessun lavoro proprio — esatto)
du_η/dτ = 0         (vincolo mantenuto)
```

L'eq. di η è soddisfatta identicamente sull'ansatz (identità di gauge):
worldline consistente. Confronto con Kerr: lì `f_t = 0` (E conservata dal
Killing); qui `f_η = aH ≠ 0` è la firma quantitativa della non-stazionarietà.

### R4. Congelamento: barriera temporale, non spaziale
Il moto esiste solo per `a(t) < Ê`. A `a = Ê`: `v → 0` liscio, la particella
si **congela** a τ finito. È l'analogo dinamico della superficie di luce di
Kerr (`F=0`), ma la barriera è nel **tempo** anziché nello spazio.

### R5. de Sitter (forma chiusa), `a(η) = −1/(Hη)`, η ∈ (−∞,0)
Partenza a `a=1` (η₀ = −1/H):

```
congelamento:  a = Ê  a  t_f = ln(Ê)/H
Δx_max · H = ( √(Ê²−1) + asin(1/Ê) − π/2 ) / Ê
T_t(η₁)  = ln a(η₁) / H
T_τ(η₁)  = (−Hη₁ − 1)/(Ê H² η₁) ;   T_τ(congelamento) = (Ê−1)/(ÊH)  finito
```

Limite nullo `Ê → ∞`: `Δx_max·H → 1` — la brachistocrona nulla raggiunge
**esattamente** l'orizzonte comovente di Hubble (`r_H = 1/H` a η₀).
Consistenza conforme Kovner–Perlick verificata (domanda 5 del doc: la
brachistocrona massiva NON penetra l'orizzonte; lo satura solo il limite nullo).

### R6. Rotaia conforme Ê vs rotaia fisica γ₀ (domanda 1 del doc)
Vincolo alternativo `−u·n = γ₀` (n = osservatori comoventi):

- `γ₀` fisso ⇒ `v = √(1−1/γ₀²)` **costante**, nessun congelamento;
  `−u_η = a γ₀` cresce con l'espansione ⇒ la rotaia deve **pompare energia**
  contro il redshift (rotaia attiva).
- `Ê` fisso ⇒ rotaia passiva/conservativa (`f·u = 0`, componente conforme
  fissa `f_η = aH`), con redshift `E_phys ∝ 1/a` e congelamento.

Coincidono solo per `a = cost` (Minkowski). La scelta è fisica, non
matematica: `Ê` è la generalizzazione naturale della rotaia ideale di Kerr
(nessun lavoro netto), `γ₀` modella una rotaia motorizzata.

### R7. Limiti di controllo
- Minkowski `a=1`: `v = √(Ê²−1)/Ê`, `γ = Ê` — standard. ✓
- Limite nullo: conformemente invariante, retta di Minkowski. ✓

### R8. Metrica ottica dal Killing conforme (`flrw_optical_metric_sympy.py`)
Su `ĝ` statica (`g = a² ĝ`) la riduzione ottica esiste ed è **esatta**:

- **Luce**: `ĝ_{ηi} = 0` (niente Randers/rotazione) ⇒ metrica ottica =
  **3-geometria comovente** `dl_k² = dχ² + S_k(χ)² dφ²`, indice `n ≡ 1`.
- **Massive (rotaia Ê)**: Fermat con indice solo temporale `n(η) = 1/v(η)`.
  **Orologio ottico** `λ(η) = ∫ v dη` ⇒ `dλ = dl_k`: moto unit-speed in λ.
  Minimizzare `η_arrivo` ⇔ minimizzare la lunghezza comovente ⇒ metrica
  ottica = geometria comovente per **tutti e tre i rami** (t, τ, η —
  monotonia). Tutta la non-stazionarietà sta nell'orologio, niente nella
  geometria.
- de Sitter: `λ(u) = [√(Ê²−1) − √(u²−1) − acos(1/Ê) + acos(1/u)]/(ÊH)`,
  `λ(η_f) = Δx_max` (identità verificata). k=+1 numerico: brachistocrona =
  cerchio massimo (radiation `a=η`, argmin a ε=0).
- **Degenerazione ottica**: `n → ∞` al congelamento `a → Ê`. In Kerr la
  riduzione ottica degenera su una **superficie spaziale** (`F→0`,
  ergosfera); in FLRW degenera a un **istante cosmico**. Stessa patologia
  `1/F`, ruotata dallo spazio al tempo.

### R9. Hamiltoniane dei rami t e τ (`flrw_hamiltonians_sympy.py`)
Formulazione di controllo ottimo (Pontryagin), tempo cosmico t, controllo =
direzione û, dinamica `dx/dt = (v/a)û`:

```
H_t   = (v/a)|p| − 1
H_τ   = (v/a)|p| − a/Ê          v(t) = √(1 − a²/Ê²)
```

- `dp/dt = −∂H/∂x = 0` (omogeneità) ⇒ direzione fissa ⇒ retta comovente:
  le traiettorie escono dalle eq. di Hamilton (verificato vs forme chiuse
  R5 in de Sitter, 12 cifre).
- `H` **non** conservata: `dH/dt = ∂H/∂t ∝ ȧ` — versione hamiltoniana del
  forzante non autonomo (in Vaidya: `dp_v/dr ∝ m′`).
- Trasversalità (arrivo spaziale, istante libero): `H(t₁) = 0` ⇒
  `|p|_t = a₁/v₁`, `|p|_τ = a₁²/(Êv₁)` — analogo esatto di `p_v(r₁)=0`
  in Vaidya; il moltiplicatore `μ(η)` del worldline è il costate riscalato
  (`C = Êv(1+μÊ)` costante).

## Figure di validazione

In `FLRWfigures/` (`genera_figure_flrw.py`): cinematica del vincolo e
congelamento; worldline comoventi vs orizzonte di Hubble (Δx_max analitico vs
quadratura, err ~2e-6); minimo variazionale diretto di T_t e T_τ sulla retta
comovente (validazione R2).

## Conseguenza per il programma di ricerca

FLRW è il caso "tutto tempo, niente forma": utile come laboratorio per il
vincolo non stazionario (redshift della rotaia, congelamento, forza `f_η=aH`),
ma la fenomenologia di penetrazione con orbita non banale richiede
**inomogeneità spaziale + non stazionarietà** ⇒ prossimo passo: **Vaidya**
(orizzonte apparente mobile `r = 2m(v)`).
