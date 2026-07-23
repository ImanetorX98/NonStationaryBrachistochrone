# Asimmetria accrescimento/evaporazione della correzione adiabatica δφ_V

## Idea fisica
A primo ordine δφ ∝ ṁ ⇒ ṁ→−ṁ flippa solo il segno (banale). L'asimmetria vera:
- **Accrescimento** (ṁ>0): Vaidya **ingoing**, clock = tempo **avanzato** v=t+r∗.
- **Evaporazione** (ṁ<0): Vaidya **outgoing**, clock = tempo **ritardato** u=t−r∗
  (ds²=−f du²−2du dr).

La riduzione congelata ∂_m F è IDENTICA (dipende solo da S); cambia solo il clock, e
v↔u flippa il segno del tortoise r∗.

## Decomposizione (verificata 1e-15)
Con v = E U₃ + r∗,  r∗ = (r−r₀) + 2m ln((r−2m)/(r₀−2m)):

    δφ/ṁ|accr = A∞ + B_hor        δφ/ṁ|evap = A∞ − B_hor
    A∞  = ∫ ∂_m F · E U₃ dr   (SIMMETRICA, polilog all'INFINITO)
    B_hor = ∫ ∂_m F · r∗   dr   (ANTISIMMETRICA, tortoise)

**Asimmetria netta** (a |ṁ| uguale): δφ_accr − δφ_evap = 2|ṁ| B_hor,
portata INTERAMENTE dal tortoise; A∞ (settore E U₃) è comune e si CANCELLA.

## Numeri (M=1, E=1.4, J=2.5)
- A∞ = 8.72756968   (polilog ∞)
- B  = 6.11272429 = B_bulk(4.58262818, 2a specie) + B_hor-log(1.53009610, dilog D_k)
- accr = A+B = 14.84029397  (= clock v, match diretto 3.6e-15)
- evap = A−B =  2.61484540  (= clock u, match diretto 8.9e-16)
- asimmetria = 2B = 12.22544857
- (i tre pezzi 8.73/4.58/1.53 coincidono con lo split del clock in vaidya_clock.py)

## Struttura
B_hor per parti: 2m[A_m ln(r−2m)] − 2m ∫A_m/(r−2m) dr, e ∫A_m/(r−2m) = Σ c_k^M D_k
(dilog d'ORIZZONTE). ⇒ la parte trascendente dell'asimmetria È il dilog D_k già
derivato. L'asimmetria accr/evap = faccia FISICA del pezzo peso-2 d'orizzonte, come la
finestra di penetrazione dell'ergosfera (misura finita per un segno di J, nulla per
l'altro) è la faccia FENOMENOLOGICA della stessa geometria d'orizzonte.

Script: `vaidya_asymmetry.py`. Paper: eq:vaidya-asymmetry (main + PRD).
