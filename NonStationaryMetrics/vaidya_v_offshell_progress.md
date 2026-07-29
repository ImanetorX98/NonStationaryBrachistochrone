# Vaidya ramo v — off-shell generico: struttura verificata + chiusura

## Stato
Struttura **corretta e verificata numericamente** (bug fixato). Assemblaggio simbolico completo
in corso. Curva `S_v=r·DE·Q2v` (genus-2, `Q2v`=quartica a=0), diversa dal ramo τ.

## Bug fixato
Avevo derivato `∂_pr G` dalla `G` GIÀ on-shell → sbagliato. Serve la derivata **off-shell** valutata
on-shell. Con la G on-shell semplificata si perdono termini.

## Struttura verificata
Ingoing `p_r<0` ⇒ `y=K1+p_r K2<0`, `sign(y)=−1`. `K1=−DE/r`, `K2=(r−2m)DE/r²`, `w=DE/r`.
```
kernel = d_pr G/H_pr = A_kernel_v/√S_v ,  A_kernel_v = +E²J·DE·r⁴/Q2v   (PURO 2ª specie)   [1e-16]
inner  = m H_m/H_pr = elem_inn(razionale) + A_inn·sign(y)/√D_v                             [1e-17]
wrap   = -∫(A_kernel_v/√S_v)·Σ dr = block1 + block2                                        [3.5e-18]
```
- **block1** = kernel × Σ_2nd = **genus-2 A+B+C** (stessa macchina del τ: Hermite su Q2v, dilog a
  r=2m [orizzonte] e DE=0).
- **block2** = kernel × Σ_elem = **DOMINANTE** (585%, grande cancellazione con block1). Σ_elem è la
  primitiva elementare della parte razionale dell'inner.

Le "correzioni" precedenti (D_v^{3/2}, parte elementare del kernel) erano **artefatti del bug**.

## Σ_elem CHIUSO simbolicamente + classe weight-2 (block2)
`Σ_elem = ∫elem_inn dr` risolto in forma chiusa (verificato `d/dr=elem_inn`, split totale a **5e-17**):
```
Σ_elem = 4m²/(r−2m) − 2m·log(r−2m)         (UN solo log, all'ORIZZONTE; niente log(DE))
```
Quindi:
```
block2 = -∫(A_kernel_v/√S_v)·[4m²/(r−2m)] dr   →  dilog genus-2 a r=2m (STESSA base del block1)
         + 2m·M^{2m},   M^{2m}=∫(A_kernel_v/√S_v)·log(r−2m) dr   (UNA lettera log-Abeliana, coeff 2m)
```
Molto più pulito del temuto: **una sola** lettera nuova `M^{2m}` (weight-2, Abeliano 2ª specie pesato
dal log all'orizzonte). Via IBP `M^{2m}=Φ_A log(r−2m) − ∫Φ_A/(r−2m)dr`, con `∫Uₖ/(r−2m)dr` la classe
weight-2 distinta (Abeliano/fattore-lineare, senza √S; tabulata). Il resto (E_rat=4m²/(r−2m)) dà dilog
genus-2 a r=2m, stessa base del block1. Radice: costo "−1" tempo avanzato → parte additiva in p_r0.

## Rimane
1. `Σ_elem` simbolico (E_rat, c1, c2) via `sp.integrate(elem_inn)`;
2. block1: A+B+C (riuso macchina τ, kernel A_kernel_v, inner P_inner_v=−m·SECp·r²/(K2²w));
3. block2: parte 2ª specie + `c1 M^{2m}+c2 M^{DE}`;
4. verifica block1+block2 == diretto; physics-anchor (attenzione segno: Vaidya `m` cresce).

## Scripts
`VaidyaMetric/vaidya_v_offshell_structure.py` (struttura verificata 3.5e-18),
`vaidya_v_step1_clean.py`, `vaidya_v_kernel_probe.py`.


## CHIUSO (esplicito)
P_inner_v ridotto: `Σ_2nd = Σ poly_k U_k + a2/(r−2m)² Π²₂ₘ + a1/(r−2m) Π₂ₘ + b/DE Π_DE`,
`a2=32E²m⁵`, `a1=2m²(−E²J²+32E²m²−4m²)`, `b=8m⁴/(E²−1)³` (simbolici). Σ_elem=4m²/(r−2m)−2m·log(r−2m).
Wrap esplicito == diretto **5e-17**; PHYSICS-ANCHORED **7.6e-11**, slope 1.96. Espansione pura W_jk =
stessa macchina τ (block1 A+B+C con letter doppio-polo D²₂ₘ; block2 = terza-specie weight-1 + 2m·M^{2m}).
Script: `VaidyaMetric/vaidya_v_explicit_assembly.py`.
