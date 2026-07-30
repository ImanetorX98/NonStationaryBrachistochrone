# Forma UNIVERSALE della sorgente off-shell S_D (30 lug 2026) — TROVATA e VERIFICATA

Risponde all'estensione 3 del referee e al punto di GPT ("S_D non ancora sviluppato in forma
universale"). Risultato: **una forma universale ESISTE**, derivata e verificata (simbolico esatto
+ numerico a precisione macchina).

## Setup
Sorgente della correzione adiabatica first-order (extended Hamiltonian):
`S_D(λ) = ∫_0^λ D H dλ'`, con `D = E_eff ∂_E_eff + J_eff ∂_J_eff + P_r ∂_P_r` (Eulero completo di fase),
`DH = E H_E + J H_J + P_r H_Pr`. Qui H è la Hamiltoniana ottica frozen in momento normalizzato
`P_r=p_r/A`, `E_eff=Ê/A`, `J_eff=J/A` (il fattore conforme entra SOLO via questi riscalamenti).

## Fatto strutturale chiave: H è una forma di Randers (arrival-time Finsler)
Per qualsiasi metrica seed stazionaria assisimmetrica, riscalata conformemente,
```
H = β(r,E) · J  +  √( A_rr(r,E) p_r² + A_φφ(r,E) J² )  −  1
```
(termine lineare in J + radice omogenea grado 1 in (p_r,J) − normalizzazione). Nel codice TK–Kerr:
`β=b v²/P̄²`, `A_rr=(Δ v²/P̄²)(Δ/r²)`, `A_φφ=(Δ v²/P̄²)/P̄²`, `v²=1−f/E²`, `P̄²=P²+b²/E²`.

## Identità universale di Finsler (dimostrata per β,A_rr,A_φφ ARBITRARIE)
Poiché `βJ + √(...)` è omogenea di grado 1 in (p_r, J), Eulero dà
```
(J ∂_J + P_r ∂_P_r) H = H + 1        ← ESATTA, universale (sympy: =0 anche con funzioni astratte)
```
Quindi
```
DH = E H_E + (H + 1)
```

## Collasso sullo shell frozen
L'orbita frozen (ε=0) vive su **H = 0** (verificato: max|H|=2e-15 sull'orbita). Perciò l'integrando
di S_D collassa:
```
DH |_frozen = 1 + E H_E
```
e integrando:
```
┌─────────────────────────────────────────────┐
│   S_D(λ) = λ  +  ∫_0^λ  E_eff H_{E_eff}  dλ' │   ← FORMA UNIVERSALE
└─────────────────────────────────────────────┘
```
- Il termine **λ** è COMPLETAMENTE universale (solo struttura Finsler; vale per ogni brachistocrona
  arrival-time, qualsiasi metrica).
- Il termine **∫ E H_E dλ** è l'UNICO pezzo metrica-specifico ("energy-Euler"), determinato dagli
  STESSI potenziali metrici (f, b, P, Δ) via regole universali: `E ∂_E v² = 2(1−v²)`,
  `E ∂_E P̄² = −2b²/E²`. Nessun input nuovo. Per-metrica riduce alla stessa base di funzioni speciali
  (U_k = ∫r^k/√S, ecc.) già usata — cioè `∫E H_E dλ = ∫E H_E dr/H_Pr`.

## Verifiche (TK–Kerr, M=1,a=0.9,Ê=1.4,J=6,r0=12)
- `(J∂_J+P_r∂_Pr)H − (H+1) = 0` simbolico esatto (anche con funzioni astratte).
- `DH − (E H_E + H + 1) = 0` simbolico esatto.
- `E∂_E v² = 2(1−v²)`, `E∂_E P̄² = −2b²/E²` simbolici esatti.
- Numerico sull'orbita frozen: max|H|=2e-15; max|DH−(1+E H_E)|=2e-15;
  max|S_D_diretto − (λ+∫E H_E dλ)| = 2e-13; valore finale S_D=31.893060 identico.

## Interpretazione / rapporto col paper
- Il paper aveva già `∫P_r H_Pr dλ = ∫p_r dr` (l'azione radiale). Il NUOVO risultato è più forte e
  universale: l'INTERA S_D = λ + (un solo termine energy-Euler), grazie a J H_J + P_r H_Pr = H+1 = 1 su shell.
- Chiude l'estensione 3: la sorgente off-shell HA una forma universale (metrica-indipendente nella
  struttura); GPT aveva ragione che la forma *generale* non era scritta, ma ora lo è: `S_D = λ + ∫E H_E`.

## Ambito e onestà
- Vale per il caso **conformemente stazionario** (D con dilatazione P_r; orbita frozen su H=0) = famiglia TK.
- Per **Vaidya** (non conforme, `m(v)` funzione di massa, `D=Θ=m∂_m`, NIENTE dilatazione P_r) la struttura
  è diversa: lì manca il termine P_r∂_Pr, quindi l'identità H+1 non si applica allo stesso modo. Serve
  l'analogo separato (S_D = ∫ m H_m dλ) — da analizzare a parte.
- La riduzione finale di `∫E H_E dλ` a funzioni speciali resta per-metrica (stessa classe genus-2),
  ma la DECOMPOSIZIONE `S_D = λ + ∫E H_E` è universale e chiusa.

## Caso VAIDYA (mass-function, non conforme) — forma universale DIVERSA, anch'essa trovata
Vaidya frozen = Schwarzschild (a=0, β=0): `H = √(f v²)·√(f p_r² + J²/r²) − 1`, `f=1−2m/r`.
Parametro lento = m (funzione di massa), `D = Θ = m∂_m` (NIENTE dilatazione P_r). Quindi `S_D = ∫ m H_m dλ`.

Tre identità (tutte simboliche esatte, =0):
1. **Auto-similarità** di Schwarzschild: sotto `(r,m)→κ(r,m)` con E,p_r invarianti e J→κJ, H è invariante
   (peso 0) ⟹ `r H_r + J H_J + m H_m = 0`, cioè `m H_m = −(r H_r + J H_J)`.
2. Regola universale `m ∂_m f = f − 1`.
3. La Finsler `(J∂_J + p_r∂_pr)H = H+1` vale ancora (a=0).

Combinando su shell (H=0 ⟹ `J H_J = 1 − p_r H_pr`):
`m H_m = −(r H_r + J H_J) = (p_r H_pr − r H_r) − 1`. E per Hamilton (`dr/dλ=H_pr`, `dp_r/dλ=−H_r`):
`p_r H_pr − r H_r = d(r p_r)/dλ`. Perciò
```
┌──────────────────────────────────┐
│   S_D(λ) = [ r p_r ]_0^λ  −  λ    │   ← FORMA UNIVERSALE VAIDYA (mass-function)
└──────────────────────────────────┘
```
Il termine off-shell è un **termine di bordo** (derivata totale) meno λ — ancora più semplice del conforme.
Verifica numerica (Schwarzschild m=1, Ê=1.4, J=6, r0=12): max|H|=3.5e-15;
max|m H_m + (r H_r+J H_J)|=7e-16; max|S_D − ([r p_r]−λ)| a livello di quadratura; S_D end −7.607436 identico.
Script: `ThakurtaMetric/universal_SD_source_check.py` (blocco TK) + `VaidyaMetric/universal_SD_source_vaidya.py`.

## CONCLUSIONE — c'è una forma universale, ma è DIVERSA per tipo di deformazione (non metrica-specifica)
| | Conforme (TK) | Mass-function (Vaidya) |
|---|---|---|
| operatore lento | `D=E∂_E+J∂_J+P_r∂_Pr` | `D=m∂_m` |
| sorgente | `S_D = λ + ∫E H_E dλ` | `S_D = [r p_r] − λ` |
| pezzo extra | energy-Euler (quadratura, per-metrica) | termine di bordo (derivata totale) |
| origine | omogeneità di Finsler | Finsler + auto-similarità |

- **NON** esiste una singola formula che copra entrambi i tipi: il parametro lento accoppia la geometria in
  modo diverso (riscalamento conforme di TUTTI i momenti vs shift di una funzione metrica con auto-similarità).
- **MA** non è "metrica-specifica": dentro OGNI classe di deformazione la formula è **universale**
  (metrica-indipendente). Tutte le metriche conformemente stazionarie condividono `λ+∫E H_E`; tutte quelle
  a funzione-di-massa auto-simili condividono `[r p_r]−λ`.
- **Scheletro comune**: entrambe hanno il termine universale **±λ** dalla stessa identità di Finsler
  `(J∂_J+p_r∂_pr)H=H+1`. È una parziale unificazione, non una coincidenza.

## Prossimo (opzionale)
- Aggiungere entrambi i risultati al paper (paragrafo "Universal form of the off-shell source"): chiude
  esplicitamente l'estensione 3 e la critica di GPT, e unifica parzialmente TK/Vaidya via il termine ±λ.
