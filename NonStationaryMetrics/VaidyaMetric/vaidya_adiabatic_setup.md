# Adiabatico in Vaidya — setup (frozen = Schwarzschild genus-2, correzione dM/dv)

Programma: forme chiuse delle brachistocrone t/τ in **Vaidya** (Schwarzschild
dinamico, massa M(v) che irradia), per analogia con Thakurta-Kerr.

## 1. Metrica e brachistocrona
Vaidya entrante:  ds² = −(1−2M(v)/r)dv² + 2 dv dr + r²dΩ².
Sferica, NO spin, NO ergosfera; orizzonte DINAMICO. Coordinata nulla v (tempo
avanzato). Rami: **v-branch** (min tempo avanzato v, analogo di t) e **τ-branch**
(tempo proprio).

## 2. Famiglia CONGELATA = Schwarzschild (limite a=0 di Thakurta-Kerr)
Ponendo a=0 nelle forme genus-2 di Thakurta-Kerr:
- τ: F_τ = J√Emu/√S,  Emu=(E²−1)r+2M,
  S = r(r−2M)Emu(r²(r−2M)−J²Emu)  [SESTICA, genus 2].
- **NON ellittica**: 6 radici distinte (discriminante≠0). La brachistocrona
  VINCOLATA è genus-2 anche in Schwarzschild (le ellittiche in letteratura sono le
  GEODETICHE). Ellittica solo sulle separatrici |J|=Jc.
- Struttura semplifica vs Kerr: orizzonti r_+=2M=r_e, r_-=0; niente frame-dragging.

## 3. Espansione adiabatica (M(v) lento)
Parametro lento: Ṁ = dM/dv. Frozen-coefficient WKB:
  φ(r) = φ_0(r; M(v(r))) + Ṁ · δφ + O(Ṁ²),
  δφ = ∫ ∂_M F · v(r) dr · (fattore),  v(r)=tempo avanzato lungo l'orbita (clock).
STRUTTURA IDENTICA a Thakurta-Kerr con: ∂_M F ↔ ∂_E F, v(r) ↔ η(r), Ṁ ↔ A'/A.

## 4. Riduzione ∂_M F — VERIFICATA ✓ (`vaidya_dMF_reduction.py`)
∂_M F = N_M/S^(3/2), N_M = S ∂_M K − K ∂_M S/2 POLINOMIALE (grado 6),
N_M = J r Emu·(poly grado 4). Riduzione:
  ∂_M F = d(𝒜/√S) + Σ_{k=0}^4 c_k^M r^k/√S,
verificata dM F(diretto)=d(𝒜/√S)+M_poly/√S a **1e-15**.
c_k^M @(M=1,E=7/5,J=5/2) = [−0.8715,−0.3597,1.0354,0.2129,−0.1894].
(Check: N_E=EJr⁴(r−2M)²Emu = limite a=0 di Thakurta-Kerr ✓.)

## 5. Cosa segue (stessa pipeline di Thakurta-Kerr, verificata lì)
- clock v(r): tempo avanzato lungo l'orbita congelata (analogo η=t/τ). Da derivare
  esplicito (dv/dr per la brachistocrona Vaidya).
- decomposizione: δφ = ½ · Ṁ · Σ_{k<j} Q_kj^M W_kj + (prodotti clock×abeliani) +
  (elementari log), Q_kj^M = c_k^M b_j − c_j^M b_k, W_kj = polilog genus-2 (STESSI
  della curva S, dipendono solo dalla curva). Coefficienti ALGEBRICI.
- I W_kj sono già quelli di Thakurta-Kerr a=0 (stessa curva genus-2 S).

## 6. Differenze da Thakurta-Kerr
- NO teorema di trasferimento conforme (M(v) non è fattore conforme): il congelato
  è Schwarzschild(M), non Kerr(E_eff). Ma la struttura adiabatica è la stessa.
- Coordinata nulla v → causalità diversa; Vaidya resta BUCO NERO (vincolo utente).
- Frozen più semplice di struttura (no spin) ma STESSO genere (2).

Script: `vaidya_dMF_reduction.py`. Riferimenti: `../ThakurtaMetric/
phi_adiabatic_closed_form.md`, `../KerrMetric/kerr_psi_explicit_verified.py`.

## 7. Clock v(r) — VERIFICATO ✓ (`vaidya_clock.py`) + risultato NUOVO
Tempo avanzato (M=M(v)): v=t+r_*, quindi
  dv/dr = (E/f) dτ/dr + 1/f = E r³/√S + r/(r−2M)   (verifica simbolica: diff=0).
  dτ/dr = r²(r−2M)/√S (limite a=0, τ-branch).
v(r) = E R_3 + r + 2M ln(r−2M),  R_3=∫r³/√S (2ª specie).
Il tortoise r/(r−2M) → **logaritmo d'orizzonte** 2M ln(r−2M).

**δφ_V = Ṁ ∫∂_M F·v dr** si spacca in TRE (verificato 1.8e-15):
1. polilog genus-2 (da E R_3) — come Thakurta-Kerr.
2. elementari (da r) — 2ª specie + polinomiale.
3. **∫∂_M F·ln(r−2M) dr — DILOGARITMO D'ORIZZONTE, NUOVO.**

RISULTATO: il tempo avanzato di Vaidya inietta un log d'orizzonte (v~ln(r−2M) al
r_+) → un dilogaritmo iperellittico con "lettera" dlog AL R_+ (l'orizzonte),
distinto dal dilog di Thakurta-Kerr (che veniva dal dipolo 3ª specie all'∞).
Fisica: la variazione di massa è pesata dalla divergenza log del tempo avanzato
all'orizzonte. Coefficienti algebrici (c_k^M); i W_kj sono gli stessi (curva S a=0).

## 8. FORMA ESPLICITA CHIUSA di δφ_V — VERIFICATA ✓ (riutilizzabile)
`vaidya_explicit.py`. Verificata a 9.4e-14 (coefficienti DATI, no fit):

  δφ_V/Ṁ = A_M(r) v(r) − E ∫A_M r³/√S dr − ∫A_M r/(r−2M) dr

INGREDIENTI (tutti espliciti):
- A_M = ∂_Mφ_0 = 𝒜^M/√S + Σ_{k=0}^4 c_k^M U_k  (riduzione ∂_M F, c_k^M ALGEBRICI;
  @E=7/5: c_k^M=[−0.8715,−0.3597,1.0354,0.2129,−0.1894]).
- v(r) = E U_3 + (r−r0) + 2M ln((r−2M)/(r0−2M))   (clock tempo avanzato).
- U_k = ∫r^k/√S dr  (integrali abeliani 1ª/2ª specie, sulla curva S a=0 Schwarzschild).

DECOMPOSIZIONE in funzioni NOMINATE (espandendo i due integrali):
  ∫A_M r³/√S = [∫𝒜^M r³/S elem] + Σc_k^M(½U_k U_3 + ½W_{3k}) − const·U_3
  ∫A_M r/(r−2M) = [∫𝒜^M r/((r−2M)√S) = L_{2M} 3ª specie orizzonte]
                 + Σc_k^M·2M·D_k + Σc_k^M(rU_k−U_{k+1}) + elem
- W_{jk}=∫(U_j dU_k−U_k dU_j): polilog genus-2 all'∞ (come Thakurta-Kerr).
- **D_k=∫U_k/(r−2M) dr: dilog genus-2 all'ORIZZONTE (NUOVO in Vaidya).**
- L_{2M}=∫𝒜^M r/((r−2M)√S): 3ª specie all'orizzonte (θ-log, peso 1).

RIUTILIZZO: data la curva S (=Schwarzschild a=0), i c_k^M (algebrici) e le funzioni
nominate {U_k, W_jk, D_k, L_{2M}}, si valuta δφ_V per qualunque M(v) lento.
Trascendenza: W_jk (polilog ∞) + D_k (dilog orizzonte); coefficienti algebrici.
Differenza da Thakurta-Kerr: il dilog è all'ORIZZONTE (tortoise del tempo avanzato),
non all'∞ (dipolo 3ª specie).
