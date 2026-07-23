# Teorema adiabatico unificato — l'indicatrice che respira

Sintesi analitica delle correzioni adiabatiche (FLRW / Vaidya / Thakurta-Kerr),
collegata al filo conduttore del paper: la **forma dell'indicatrice** (ellisse di
Randers/Zermelo).

## 1. Setup unificato
Indicatrice (eq. ellisse del paper): ẋ(θ)=c+R(θ),
- **centro c** = vento di Zermelo (0 FLRW, radiale entrante Vaidya, azimutale
  frame-dragging Thakurta-Kerr);
- **assi R** = velocità locale (si restringono alla superficie di congelamento).

Il congelato (λ fisso) dà la brachistocrona φ_0(r;λ)=∫F(r;λ)dr sulla curva
spettrale y²=S(r;λ) di genere ≤2. Parametro lento λ = fattore di scala a (FLRW),
massa m(v) (Vaidya), conforme A(η) (Thakurta-Kerr).

## 2. TEOREMA (adiabatico universale)
Se λ=λ(clock) varia lentamente (λ̇ piccolo):

  φ(r) = φ_0(r; λ(clock(r))) + λ̇ · δφ(r) + O(λ̇²),
  δφ = ∫ ∂_λ F · (clock) dr
     = Σ_{k<j} q_{kj} W_{kj} + (prodotti) + (elementari),

con:
- **q_{kj} ALGEBRICI** (razionali nei parametri della metrica; via la riduzione
  ∂_λ F = N_λ/S^{3/2}, verificata: ∂_E F 1e-15 TK, ∂_M F 1e-15 Vaidya);
- **W_{kj}=∫(U_j dU_k − U_k dU_j)** POLILOG genus-2 della curva congelata
  (U_k=∫r^k/√S). **UNIVERSALI** — dipendono solo dalla curva, non da quale
  parametro respira.

La **trascendenza (W_kj) è universale**; i tre spaziotempo differiscono solo in
COSA respira dell'ellisse e in QUALE clock.

## 3. Collegamento alle ellissi — il dizionario del respiro
∂_λ F si spacca geometricamente in **respiro-del-vento** + **respiro-degli-assi**:
  ∂_λ F = ∂_c F · (∂_λ c) + ∂_R F · (∂_λ R).

| | respira | clock | lettera peso-2 (sorgente dilog) |
|---|---|---|---|
| FLRW | assi R isotropi (c=0) | η conforme | NESSUNA (cerchio, genere 0, rami degeneri) |
| Vaidya | **centro c_r** (vento radiale, ∂_m√(2m/r)) | v (tempo avanzato) | **dlog all'ORIZZONTE** r=2m (tortoise di v) |
| Thakurta-Kerr | **assi R** (scala conforme, E_eff=Ê/A) | t / τ | **dipolo all'∞** (crescita log del clock a grande r) |

**Enunciato geometrico:** il dilogaritmo iperellittico (peso 2) è UNIVERSALE come
CLASSE; la sua LETTERA (dove sta il dlog: orizzonte vs infinito) è determinata da
QUALE dato dell'ellisse respira e da QUALE clock:
- respiro del **vento** (c) + clock del tempo **avanzato** → dlog all'**orizzonte**
  (Vaidya);
- respiro degli **assi** (R) + clock t/τ → dipolo all'**infinito** (Thakurta-Kerr).

## 4. Degenerazione FLRW (il cerchio)
FLRW: c=0 (no vento), R isotropi (no eccentricità) ⟹ i rami t,τ,η coincidono
(dizionario del paper) e la struttura polilog collassa (curva genus-0, ottica
piatta conformemente). Il cerchio è il caso degenere: nessun dilog, correzione
banale. È l'origine del filo "shape": cerchio → ellisse-vento → ellisse-rotante.

## 5. Perché è un teorema (non solo analogia)
Le due ossa portanti sono verificate per entrambi i casi non banali:
- riduzione ∂_λ F = N_λ/S^{3/2}, coefficienti algebrici (TK 1e-15, Vaidya 1e-15);
- decomposizione δφ = Σ q_kj W_kj + elem, identità verificata (TK 4e-14, Vaidya 9e-14).
La curva S è genus-2 in entrambi (anche Schwarzschild!); genus-1 solo su separatrici.
Quindi la trascendenza (polilog genus-2) è la stessa; cambia solo la lettera del
dilog per geometria dell'ellisse.

Riferimenti: `ThakurtaMetric/phi_adiabatic_closed_form.md`,
`VaidyaMetric/vaidya_adiabatic_setup.md`, `paper/main.tex` §indicatrix.
