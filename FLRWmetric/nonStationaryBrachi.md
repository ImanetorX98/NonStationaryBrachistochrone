# Brachistocrone in spaziotempi NON stazionari — FLRW e Vaidya

**Documento di avvio ricerca (handoff).** Obiettivo: estendere il formalismo del
*worldline vincolato* (sviluppato per Kerr, dove funziona grazie alla
stazionarietà) a spaziotempi **dinamici**, dove non esiste un vettore di Killing
di tipo tempo. Da qui una nuova chat deve **derivare le equazioni del moto e un
principio variazionale nuovo** per la brachistocrona.

Progetto di origine: `../KerrMetric/` (vedi `doranTau.md`, `doranT.md`,
`doranOffEq.md` per il caso stazionario completo).

---

## 0. Cosa portiamo dal caso stazionario (Kerr)

La brachistocrona relativistica (Perlick 1991): tra due punti **spaziali**, la
curva che minimizza il tempo di percorrenza (proprio `T_τ` o coordinato `T_t`)
per una particella a **energia specifica fissa** `E = −u_t` (rotaia ideale che
conserva l'energia). Due formulazioni equivalenti:

- **Riduzione ottica/Randers** (comprime a una geometria di Finsler 3D fissa sullo
  spazio quoziente). Elegante ma **singolare alla superficie di luce** (ergosfera:
  `g_tt→0`), e **esiste solo se lo spaziotempo è stazionario**.
- **Worldline vincolato** (curva 4D `x^μ(τ)` con moltiplicatore per `−u_t=E`).
  Regolare in coordinate regolari (Doran), integrabile attraverso l'ergosfera.
  Lagrangiana:
  ```
  L = (1 + μ E) Λ + μ g_{tν} ẋ^ν ,     Λ = √(−g_{αβ} ẋ^α ẋ^β) = dτ/dλ
  ```
  Forza-rotaia `f^μ` con `f·u=0` e `f_t=0` (⟹ `E` conservata perché `∂_t` è di Killing).

**Ingredienti che dipendono dalla stazionarietà** (e che qui vanno sostituiti):
- `∂_t` (Killing tempo) → energia conservata `E = −u_t` (il vincolo della rotaia);
- `∂_φ` (Killing assisimmetria) → momento di Fermat `J` conservato;
- il **quoziente** per il flusso di `∂_t` → "spazio" e "tempo di percorrenza" canonici;
- la **riduzione ottica** a Finsler 3D fissa.

---

## 1. Il principio variazionale generale (nucleo invariante)

Il funzionale `∫dτ` è invariante; ciò che va riformulato è il **vincolo**. Tre vie
note, da sviluppare:

1. **Fermat di Kovner–Perlick** (Kovner 1990, Perlick 1990, *Living Reviews* Perlick 2004):
   per la luce, i raggi da una sorgente a una **worldline-osservatore** estremizzano
   il **tempo d'arrivo** (tempo proprio dell'osservatore al ricevimento). Vale in
   spaziotempi **arbitrari**. → La brachistocrona **nulla** (limite `E→∞`) generalizza
   pienamente.

2. **Navigazione di Zermelo / Finsler non-autonomo**: prescrivere il campo di
   velocità locale `v(x,t)` (invece di derivarlo dall'energia). Least-time path =
   geodetica di Randers–Finsler, ora **dipendente dal tempo**. È la generalizzazione
   naturale della brachistocrona **massiva**.

3. **Worldline vincolato covariante**: la Lagrangiana con moltiplicatore, ma con
   vincolo riferito a una **congruenza di osservatori scelta** `n^μ` (foliamento
   `t=cost`). Vincolo: `−u·n = γ` prescritto (velocità locale rispetto a `n`).
   Equazioni **non-autonome**; nessuna riduzione ottica statica; si integra il
   worldline 4D direttamente.

**Tabella di cosa sopravvive:**

| ingrediente | stazionario | non stazionario |
|---|---|---|
| `∫dτ` da minimizzare | ✓ | ✓ (invariante) |
| energia conservata `E=−u_t` | ✓ | ✗ (serve Killing) |
| Finsler 3D fissa | ✓ | ✗ (Finsler non-autonomo) |
| Fermat per la luce | ✓ | ✓ (Kovner–Perlick) |
| "punti spaziali" canonici | ✓ | ✗ (serve scelta di foliamento) |

---

## 2. Caso FLRW (universo in espansione) — struttura conforme

**Metrica** (spazialmente piatta k=0, generalizzabile a k=±1):
```
ds² = −dt² + a(t)² [ dr² + r² dΩ² ]
```
con `a(t)` fattore di scala. **No Killing tempo** (`a` dipende da `t`). MA:

**Tempo conforme** `η`, con `dt = a dη`:
```
ds² = a(η)² [ −dη² + dr² + r² dΩ² ] = a(η)² · η_μν^(Mink)
```
FLRW è **conformemente piatta**. `∂_η` è un **Killing conforme** della metrica
fisica (Killing genuino della metrica conforme statica `ĝ = η_Mink`).

**Conseguenze da sfruttare:**
- **Luce**: le geodetiche nulle sono invarianti conformi → la **brachistocrona
  nulla in FLRW = brachistocrona nulla in Minkowski** (banale) via la
  trasformazione conforme. Buon test di consistenza.
- **Massive**: `Ê = −u_η` è conservata per la metrica **conforme** `ĝ`, ma l'energia
  fisica `E_phys` **redshifta** (`E_phys ∝ 1/a`). Il vincolo della rotaia può essere
  posto sulla **quantità conforme** `Ê` (conservata!) invece che su `E_phys`.

**Principio variazionale proposto (FLRW):** riformulare la brachistocrona nella
metrica conforme statica `ĝ` (dove `∂_η` è Killing genuino), risolvere lì col
formalismo stazionario, poi **ri-mappare** in tempo cosmico `t` via `dt=a dη`.
La velocità locale fisica `v(x,t)` diventa dipendente dal tempo attraverso `a(t)`.

**Domande aperte per la nuova chat (FLRW):**
1. Derivare la Lagrangiana della brachistocrona massiva in FLRW, distinguendo il
   vincolo su `Ê` (conforme, conservato) vs `E_phys` (redshiftato). Quale è la
   "rotaia" fisicamente corretta?
2. La brachistocrona che minimizza il **tempo cosmico** `t` vs quella che minimizza
   il **tempo conforme** `η` vs il **tempo proprio** `τ`: sono diverse? Come?
3. Effetto dell'**espansione** (Hubble `H=ȧ/a`): un'orbita che parte "dopo"
   (universo più grande) è più veloce/lenta? Analogo del θ-focusing?
4. Casi concreti di `a(t)`: radiation (`a∝t^{1/2}`), matter (`a∝t^{2/3}`),
   de Sitter (`a∝e^{Ht}`, che HA un Killing conforme aggiuntivo). de Sitter è il
   più simmetrico — buon punto di partenza.
5. Confine/orizzonte: l'**orizzonte di Hubble** `r_H = 1/H` è l'analogo dinamico
   della superficie di luce? La brachistocrona lo penetra?

---

## 3. Caso Vaidya (buco nero che accresce/evapora) — struttura nulla

**Metrica** (Vaidya entrante, coordinate di Eddington–Finkelstein avanzate):
```
ds² = −(1 − 2m(v)/r) dv² + 2 dv dr + r² dΩ²
```
con `v` = tempo avanzato (`v=cost` sono coni di luce entranti) e `m(v)` = **massa
dinamica**: `dm/dv > 0` accrescimento, `dm/dv < 0` evaporazione (radiazione di
Hawking modellata). **No Killing tempo** (se `m` dipende da `v`).

**Struttura utile:**
- La coordinata `v` è **regolare attraverso l'orizzonte** (come Eddington–Finkelstein).
  Niente singolarità di coordinata al raggio di Schwarzschild dinamico.
- **Orizzonte apparente**: `r = 2m(v)` (dinamico, si muove con `m`). Diverso
  dall'orizzonte degli eventi (globale, non locale in Vaidya).
- Le geodetiche nulle radiali entranti: `v = cost`. Le uscenti: `dr/dv = (1−2m(v)/r)/2`.
- `g_vv = −(1−2m/r)` cambia segno all'orizzonte apparente — analogo di `g_tt` ma
  **dipendente da `v`**.

**Principio variazionale proposto (Vaidya):**
- **Luce**: Fermat di Kovner–Perlick — estremizzare il tempo d'arrivo a una
  worldline osservatore. Diretto, ben posto.
- **Massive**: worldline vincolato covariante. Il vincolo `−u_t=E` NON è conservato
  (`m(v)` rompe la simmetria). Sostituire con: velocità locale prescritta rispetto
  agli osservatori entranti (`v=cost` slicing) OPPURE rispetto agli osservatori
  ZAMO/Doran dinamici. Equazioni **non-autonome** in `v`.

**Domande aperte per la nuova chat (Vaidya):**
1. Derivare le equazioni del moto della brachistocrona in Vaidya col worldline
   vincolato non-autonomo. Quale vincolo sostituisce `−u_t=E`?
2. **Orizzonte apparente in movimento**: se il buco nero **accresce** (`dm/dv>0`),
   l'orizzonte cresce verso la particella — può "inghiottirla" anche se stava
   riflettendo? Se **evapora** (`dm/dv<0`), l'orizzonte si ritira — la penetrante
   può "sfuggire"? È un fenomeno **genuinamente non stazionario**, senza analogo in Kerr.
3. La penetrazione dipende dal **timing** dell'arrivo relativo a `m(v)`: c'è una
   finestra temporale di penetrazione? Un `v_c` critico oltre `J_c`?
4. Caso semplice: **Vaidya lineare** `m(v)=m_0 + μ v` (accrescimento costante) o
   `m(v)=m_0(1 − v/v_evap)` (evaporazione). Integrabili numericamente.
5. Limite quasi-statico (`dm/dv → 0`): recupero del caso Schwarzschild? Espansione
   perturbativa in `dm/dv` (analogo dell'espansione in `a` per Kerr?).

---

## 4. Strategia consigliata per la nuova chat

1. **Partire dal worldline vincolato covariante** (§1.3) — è il formalismo giusto,
   già validato su Kerr, e l'unico che si estende senza stazionarietà.
2. **FLRW de Sitter** come primo caso (massima simmetria residua, Killing conforme
   esplicito) — deriva EOM, verifica il limite Minkowski, poi generalizza a `a(t)`.
3. **Vaidya lineare** come secondo caso (struttura nulla semplice, orizzonte mobile
   calcolabile) — deriva EOM non-autonome, studia penetrazione con orizzonte in moto.
4. **Punto di forza atteso**: fenomeni **genuinamente dinamici** senza analogo
   stazionario — orizzonte che insegue/ritira, finestre temporali di penetrazione,
   redshift della "rotaia". Questo è il contributo nuovo.
5. **Strumenti**: sympy per le metriche/EOM (come in KerrMetric), integrazione
   numerica dei worldline non-autonomi, verifica dei limiti statici.

---

## 4b. RISPOSTA AL §1 (luglio 2026): il principio variazionale generale

Trovato e validato su FLRW, Vaidya, SdS, Thakurta(-Kerr). Tre ingredienti,
nessuno dei quali richiede simmetrie:

1. una funzione tempo v (scelta della foliazione = cosa si minimizza);
2. il vincolo rotaia −u·W = E con W = ∂_v: LINEARE nella velocità;
3. la normalizzazione u·u = −1: quadratica.

Lineare ∩ quadratica = **indicatrice ellittica** a ogni evento ⇒ forma di
Randers/Zermelo `dv = deriva + norma` (Legendre-duale delle Hamiltoniane
PMP con trasversalità H=0). Esiste per QUALUNQUE spaziotempo.

La "quantità conservata" è sostituita dal **vincolo attivo**: E è un
invariante controllato (Noether → controllo). Costo misurabile: f·ξ della
forza-rotaia (Kerr: f_t=0 gratis; FLRW: f_η=aH; Vaidya: f_v ∝ m′).

Gerarchia di cosa compra la simmetria:
- ∂_v Killing → F autonoma (Perlick 1991 recuperato), H conservata;
- CKV/omotetia → scheletro rigido (indice scalare n(η,r), vento rigido
  in Thakurta-Kerr, quoziente autosimilare in Vaidya lineare);
- ∂_φ Killing → momento di Fermat J (Noether puro): integrabilità ODE;
- nessuna → F esiste ma respira tutta; senza ∂_φ servirebbe HJ (PDE).

**Liceità del vincolo (status per il paper):**
1. buona posizione: μ esiste, f·u = 0 identico (= "ideale" covariante:
   la preservazione di u·u = −1 equivale a f·u = 0) — dimostrato;
2. esistenza dei minimi anche non autonomi: teorema di FILIPPOV
   (indicatrice compatta e convessa = la nostra ellisse); PMP non
   richiede autonomia;
3. limite nullo Ê→∞ = Fermat di KOVNER-PERLICK, dimostrato per
   spaziotempi arbitrari (ancoraggio di letteratura);
4. equivalenza con la navigazione di ZERMELO (v(x,t) prescritta): il
   vincolo è una derivazione simmetrica del campo di velocità, non
   un postulato ad hoc.
Status della scelta di W: teorema (Killing) → selezione canonica (CKV:
unica che degenera a Perlick + limite nullo + costo f·ξ minimo) →
selezione quasi-locale (KODAMA, sferico dinamico: -u·K = energia di
Kodama, VaidyaResults 3f) → convenzione da dichiarare (generico;
analogo classico: brachistocrona con filo in movimento). Alternativa
fisica mappata: rotaia motorizzata γ₀ (FLRWresults R6).

**Realizzabilità del controllo (non euristica):**
- costruttiva: f := Du/dτ ha f·u = 0 automaticamente (identità di
  normalizzazione) ed è di tipo magnetico (f^μ = F^{μν}u_ν con F
  antisimmetrico): ogni curva ammissibile è realizzabile da forza
  esterna senza lavoro proprio; mantenere -u·W = E è UNA condizione
  lineare su 3 gradi di libertà di forza ⇒ risolvibile fuori dal
  luogo critico (v̄→0, orizzonte) — che coincide con le superfici
  critiche mappate; μ calcolato esplicitamente (finito).
- variazionale: lettura VAKONOMICA/controllo ottimo (si ottimizza la
  rotaia, non si simula una rotaia data): PMP (nessuna autonomia
  richiesta) + esistenza dei minimi via FILIPPOV (indicatrice compatta
  convessa). Nota da referee: per vincoli lineari in velocità la
  procedura di d'Alembert–Chetaev differisce dalla vakonomica in
  generale; qui non rileva (problema di ottimizzazione di forma) e le
  due collassano nel limite stazionario (Perlick).

Dettagli: `../VaidyaMetric/VaidyaResults.md`,
`../ThakurtaMetric/ThakurtaResults.md`, `../SdSMetric/SdSresults.md`,
`FLRWresults.md`.

## 5. Riferimenti chiave

- **Perlick, V.** (1991), *The brachistochrone problem in a stationary space-time*,
  J. Math. Phys. 32, 3148. [base stazionaria]
- **Kovner, I.** (1990), *Fermat principle in arbitrary gravitational fields*,
  ApJ 351, 114. [Fermat non stazionario per la luce]
- **Perlick, V.** (2004), *Gravitational Lensing from a Spacetime Perspective*,
  Living Rev. Relativity 7, 9. [Fermat generale, review]
- **Gibbons, Herdeiro, Warnick, Werner** (2009), *Stationary metrics and optical
  Zermelo–Randers–Finsler geometry*, PRD 79, 044022. [Zermelo/Randers]
- **Vaidya, P. C.** (1951), *The gravitational field of a radiating star*. [metrica di Vaidya]
- Progetto Kerr di origine: `../KerrMetric/doranTau.md`, `doranT.md`, `doranOffEq.md`.

---

## 6. Contesto dal lavoro su Kerr (per continuità)

Risultati stazionari già stabiliti (da riusare come limite/confronto):
- **Ramo τ** (tempo proprio, Riemanniano puro): penetra l'ergosfera **solo** a
  `J=J_c=a/E` (misura nulla, prograde). Off-equatore riflette; θ-focusing verso
  ~65° = manifestazione di `â_θθ=F→0` (metrica ottica "molle" in θ).
- **Ramo t** (tempo coordinata, Randers): penetra un **intervallo** limitato
  `(J_t^{c,−}, J_t^{c,+})` che include retrograde (frame dragging forza co-rotazione
  a `F=0`). Soglie in **forma chiusa** (elementare + saddle-node). Crossover di
  meccanismo a `a* ≈ M/√2`.
- **Worldline vincolato**: regolare; il ramo t raggiunge l'ergosfera off-equatore
  da ogni θ_0 con velocità/accelerazione finite (`dr/dτ→0` a `F=0`).
- La singolarità `1/F` è **intrinseca alla riduzione ottica** (superficie di luce),
  non alle coordinate. Nessun Carter esatto off-equatore.

L'aspettativa: in FLRW/Vaidya la "superficie di luce" diventa **dinamica**
(orizzonte di Hubble / apparente in moto), e la fenomenologia di penetrazione
acquista una **dimensione temporale** assente nel caso stazionario.
