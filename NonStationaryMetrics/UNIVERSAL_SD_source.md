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

## Prossimo (opzionale)
- Aggiungere il risultato al paper (paragrafo "Universal form of the off-shell source"), che chiude
  esplicitamente l'estensione 3 e la critica di GPT.
- Verificare l'analogo per Vaidya (`S_D=∫m H_m dλ`, m∂_m) e cercarne la forma universale sferica.
- Script di verifica: rifà i blocchi qui sopra (da consolidare in `universal_SD_source_check.py`).
