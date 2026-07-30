# La brachistocrona in tempo conforme η in Thakurta–Kerr (estensione referee 4/5)

Sviluppo della η-brachistocrona richiesto dal follow-up main13. Thakurta–Kerr è Kerr conforme
`g = A(η)² g_Kerr`; η = tempo conforme, t = tempo coordinato, τ = tempo proprio. Il paper ha già
H_η (eq:Heta), H_t (eq:Ht-tk), H_τ (eq:Htau-tk) e la riduzione di Randers (eq:randers,
`n_η=Ê/√(Ê²−A²f)`). Questo documento CHIUDE lo sviluppo: dimostra che le soluzioni η (frozen,
adiabatica, separatrice) coincidono con le t come CURVE spaziali, e presenta il gauge η come quello
naturale/conforme.

## Risultato centrale: t ≡ η come curve spaziali (frozen + adiabatica + separatrice)
### (a) Frozen — argomento analitico
`dt = A(η) dη` ⟹ le funzioni di Finsler di arrivo soddisfano `F_t = A·F_η` (letteralmente). A A costante
(slice frozen) le geodetiche di una Finsler sono invarianti per riscalamento costante ⟹ le
brachistocrone η e t sono la STESSA curva. La carica conservata è rietichettata dal fattore conforme:
`J_t = A·J_η` (osservato numericamente, sotto).

### (b) Drift (adiabatica) — argomento di monotonia
Per una curva spaziale γ a estremi fissi, l'arrivo `t_f(γ) = ∫_0^{η_f(γ)} A(η)dη` è funzione
STRETTAMENTE crescente di `η_f(γ)` (perché A>0). Una trasformazione monotona preserva l'argmin ⟹
`argmin_γ η_f = argmin_γ t_f` per QUALSIASI drift A(η). Quindi anche le brachistocrone adiabatiche
(e le separatrici, che sono limiti di famiglie di brachistocrone) coincidono t ≡ η.

### (c) Conseguenza — ci sono DUE famiglie ottimali di arrivo, non tre
- **Ramo di arrivo**: t ≡ η (stessa curva; η è la rappresentazione conforme canonica).
- **Ramo proprio**: τ (genuinamente distinto, peso extra `√(1−v²)`, wind gravitomagnetico che si cancella).

## Perché il gauge η è quello NATURALE (conforme)
- H_η ha costo unitario (`−1`), H_t ha costo `−A`; l'Hamiltoniana estesa `H_ext = p_η + H_η = 0` è la
  forma più pulita (η = parametro di evoluzione geometrico).
- Il drift è una funzione pulita del parametro di evoluzione: `ε = A'/A`, sorgente `A(η)`.
- La correzione adiabatica off-shell è esattamente la `S_D = λ + ∫E H_E dλ` (forma universale conforme,
  vedi `UNIVERSAL_SD_source.md`), con λ = η.
- Regolare attraverso l'ergosfera (`f=0 ⟹ R²=v̄²Δ/P̄>0`); degenera solo all'orizzonte `Δ=0`
  (conformemente invariante).

## Verifica numerica (M=1, a=0.9, Ê=1.4, A=1.3 frozen)
- η-orbita a `J_η=6.0` vs t-orbita: il best-match si ha a **`J_t=7.8000 = A·J_η`** (A=1.3), con
  `max|φ_η(r) − φ_t(r)| = 7.2e-9` (precisione ODE) ⟹ STESSA curva, etichetta di carica diversa.
- Il confronto naive a J uguale dà 0.084 (≠0): non è il confronto corretto (orbite IVP a J fisso su
  shell diverse ≠ stessa geodetica).
- Script: `ThakurtaMetric/eta_brachistochrone_check.py`.

## Cosa NON cambia e cosa cambia
- Le forme chiuse frozen/adiabatica/separatrice della η-brachistocrona SONO quelle della t-branch
  (eq:t-K, eq:t-genus2, ecc.) — nessuna nuova funzione speciale; è la stessa curva.
- Ciò che il gauge η aggiunge: formulazione Hamiltoniana più pulita, drift trasparente, e l'unificazione
  t/τ/η in DUE famiglie (chiude estensioni 4 e 5 insieme).

## Da propagare al paper
Paragrafo/consolidamento "The conformal-time (η) brachistochrone": (a) F_t=A F_η + monotonia ⟹ t≡η;
(b) η = gauge conforme naturale; (c) due famiglie di arrivo/proprio; (d) verifica J_t=A J_η a 7e-9.
