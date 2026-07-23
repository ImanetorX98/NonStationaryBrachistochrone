# φ(r,A) adiabatico — forma semi-analitica chiusa e struttura polilog

Documento consolidato del programma: **forma chiusa delle brachistocrone φ(r)
con fattore conforme A(η) lentamente variabile** (Thakurta–Kerr → Kerr via
E_eff=Ê/A). Rami t (tempo coordinata) e τ (tempo proprio). Tutti i risultati
validati numericamente.

Parametri riferimento: M=1, a=9/10 (o 0.9), E=Ê, J. Curva genus-2 del ramo:
y²=R₆ (t), y²=S_τ (τ).

---

## 1. La forma finale (obiettivo raggiunto a livello di struttura)

```
φ(r,A) = φ₀(r; Ê/A)  +  (A'/A) [ Closed(r) + ψ(r) ]  +  O((A'/A)²)

Closed(r) = −½ Ê · ∂_E φ₀(r) · η(r)          [CHIUSO]
ψ(r)      = ½ Ê ( ρ − ρ̃ )                     [polilog genus-2]
            ρ = ∫ ∂_Eφ₀ · h dr,  ρ̃ = ∫ η · ∂_E F dr,  h = dη/dr
```

η(r) è il **clock** accumulato lungo l'orbita congelata: tempo coordinata t(r)
(ramo t) o tempo proprio τ(r)=∫L_τ dr (ramo τ). Il ramo τ ha in più solo il
prefattore A⁻² sul timing, NON sulla forma; la forma di φ è comune ai due rami.

**Come si arriva a φ(r):** ogni pezzo è una funzione speciale sulla stessa curva
genus-2 → valuti e sommi. Vedi tabella §2.

---

## 2. I pezzi di φ(r) — tutti identificati

| pezzo | forma esplicita | classe | stato |
|-------|-----------------|--------|-------|
| φ₀(r) | ∫F dr, F=dφ/dr (K/√R₆ o K_τ/√S) | 1ª+3ª specie (Kleinian σ,℘,ζ) | ✓ |
| ∂_Eφ₀(r) | A(r)/√R + Σ_{k=0}^4 c_k ∫r^k/√R (no 3ª specie) | 2ª specie | ✓ err=0 |
| η(r) | ∫(dη/dr)dr, clock (vedi §3) | 2ª specie [τ] / 2ª+cover [t] | ✓ |
| Closed | −½Ê ∂_Eφ₀ · η (prodotto dei due sopra) | chiuso | ✓ |
| ψ(r) | ψ_ζ + ψ_ab + ψ_Li (TRE pezzi, vedi §5) | peso1 ζ + DUE peso2 | struttura ✓ |

φ(r) chiusa = [pezzi Kleinian] + (A'/A)[Kleinian + un polilog genus-2].

---

## 3. I differenziali clock dη/dr — resi ESPLICITI e validati

Prerequisito per Closed: dη/dr è on-curve? (test: (dη/dr)²·[sestica] quadrato).

### Ramo τ (tempo proprio) — ON-CURVE puro
- β cancella in Q; `(dτ/dr)²·S_τ = [r²(r−2M)]²` QUADRATO PERFETTO.
- **dτ/dr = r²(r−2M)/√S_τ** — differenziale RAZIONALE su y²=S_τ.
- Validazione: `max|√(Q/w) − r²(r−2M)/√S_τ| = 8.9e-16`.
- ω_a=∂_E F_τ dr (2ª specie su curva) + ω_b=dτ (su curva)
  ⇒ **ψ_τ = polilog iperellittico genus-2 PURO.**

### Ramo t (tempo coordinata) — ON-CURVE + cover frame-dragging
- β cancella in Q_t; `(Q_t/w)·R₆ = E²r⁶ = (Er³)²` → √(Q_t/w)=Er³/√R₆ on-curve.
- MA il termine B/f porta β=√(2Mr/(r²+a²)) (velocità fiume Doran), NON in √R₆:
  ```
  dt/dr = ρ_t(r)/√R₆  +  c_β(r)·√(2Mr/(r²+a²)),   c_β=(1−2Ma²/(rΔ))/f
  ```
- Validazione: `max|dt/dr − ρ_t/√R₆ − c_β·β| = 2.5e-14`.
  ⇒ **ψ_t = [polilog genus-2 su R₆] + [resto su cover frame-dragging β].**

### Fisica
Tempo proprio frame-independent → curva pulita genus-2.
Tempo coordinata trascina β (frame-dragging) → rivestimento √ extra.
**Dicotomia polilog t/τ = dicotomia frame-dragging.**

---

## 4. Validazione della forma ibrida (φ_hybrid == φ_full)

Ibrido = pezzi chiusi (analitici) + ψ numerico. Confronto vs adiabatica piena
(coefficiente non-autonomo −Ê∫∂_E F·η dr). Identità per parti; residuo = trapezio.

| A'/A | ramo t (η=t) | ramo τ (η=τ) |
|------|--------------|--------------|
| 0.005 | 1.4e-06 | 2.5e-08 |
| 0.02  | 5.6e-06 | 1.0e-07 |
| 0.06  | 1.7e-05 | 3.0e-07 |

Script: `kerr_adiabatic_phi_hybrid.py` (t), `kerr_adiabatic_phi_hybrid_tau.py` (τ).

---

## 5. La natura di ψ — decomposizione per specie (FINALE, dai residui)

Storia onesta (due correzioni): paper diceva "polilog"; poi ho detto "ζ puro"
(SBAGLIATO); la decomposizione dei residui dà la verità precisa.

**Residui calcolati** (`/tmp` oncurve, params razionali):
- `res_∞(ω_a=∂_E F) = 0`  → ω_a pura **1ª+2ª specie** (niente 3ª; anche β_±=δ_±=0
  agli orizzonti). ✓
- `res_∞(ω_b=dη) = 1.063` ≠ 0 → dη ha un **dipolo di 3ª specie all'∞**
  (fisica: dη~dr/r a grande r, η~log r; il tempo cresce logaritmicamente).

**ψ = ½Ê ∫(A dB − B dA)**, A=∂_Eφ₀, B=η. Decomposizione ANALITICA ESATTA (no fit,
`KerrMetric/kerr_holo_component_check.sage`) dà TRE pezzi:
```
ψ = ψ_ζ  (peso 1)  Kleinian ζ,σ           termini (2ª)×(1ª/2ª)
  + ψ_ab (peso 2)  olomorfo×olomorfo       det·∫(u1 du2 − u2 du1)   [Beilinson]
  + ψ_Li (peso 2)  3ª specie / dilog        ½Êρ₀[∂_Eφ₀ L − 2𝓛₂]
```
- **c_k** (riduz. 2ª specie ∂_E F): ESATTI, razionali in E (@E=7/5:
  [-0.531,1.979,-0.812,-0.360,0.189]). Analitici, sympy.
- **ρ₀ = M/(E²−1)^(3/2)** ESATTO (residuo 3ª specie clock; @params 1.0631).
- **ψ_ab**: coeff = det(componenti olomorfe) = b^A_0 b^B_1 − b^A_1 b^B_0 =
  −4.27−7.94i ≠ 0 (|det|/|bA||bB|=0.80). Componenti olomorfe da a-periodi esatti.
  NOVITÀ genus≥2: in genus1 c'è 1 sola olomorfa → ψ_ab ASSENTE (ellittico zζ−2logσ
  è pulito). Il FIT lo nascondeva; l'algebra esatta lo rivela.
- **ψ_Li**: L=log(θ-ratio) [teorema di Fay, NON fit], ρ₀ analitico, 𝓛₂=∫L ∂_E F dr
  dilog iperellittico esplicito (integrale singolo di θ, come Li₂).
- STORIA correzioni: paper "polylog" → "ζ puro" (SBAGLIATO) → "ζ+dilog" (2 pezzi,
  INCOMPLETO) → **"ζ + olomorfo×olomorfo + dilog" (3 pezzi, dall'algebra esatta)**.
- ψ_ζ (peso-1) resta da ridurre esplicitamente con le identità ℘ di Baker.

## 5bis. Quasi-periodi η — CALCOLATI e VALIDATI (muro §7-iii rotto)

`KerrMetric/kerr_quasiperiods_bel.sage`. La 2ª specie ora si calcola con Sage:
- Modello **DISPARI** (quintica via x=1/s, 1 punto all'∞): 2ª-specie pulite.
  Il modello pari (deg 6, 2 punti all'∞) dà x²dx/y di 3ª specie → η sbagliato.
- 2ª-specie **canoniche Baker-Enolski-Leykin**:
  `dr_1=(λ3 s+2λ4 s²+3λ5 s³)ds/4y`, `dr_2=λ5 s²ds/4y`.
- Periodi via `matrix_of_integral_values` (interi: no crash Singular).
  Pairing du=[1,s]↔dr=[dr1,dr2].
- **VALIDAZIONE**: κ=η ω⁻¹ **simmetrica a 1.4e-12**; Legendre generalizzata
  `ω'ηᵀ−ωη'ᵀ = −iπ·I` (canonica). ⇒ η corretta.

Ingredienti forma chiusa ψ ORA tutti disponibili: u(r) (Abel avanti), θ+derivate
(abelfunctions), κ/η (qui). PROSSIMO (meccanico): σ=exp(−½uᵀκu)θ, ζ_i=∂_i log σ,
ψ(r)=combinazione chiusa di ζ(u(r)); validare vs ψ numerico.

**PAPER DA CORREGGERE:** `sec:adiabatic` dice ψ "genus-two polylogarithm /
irreducible" — impreciso. ψ chiude in Kleinian ζ,σ (2ª specie), non è polilog.

---

## 6. Motore genus-2 (livello 3) — abelfunctions in Sage

- `abelfunctions` 0.2.0 compilato in SageMath 10.9 (build x86_64 forzato, patch
  API `is_LaurentSeries`; vedi memoria `abelfunctions-sage-install.md`).
- Espone RiemannTheta (con derivate), RiemannSurface, AbelMap, RiemannConstantVector.
- Bug: pipeline RS (holomorphic_differentials→Singular integralbasis) crasha su
  leading non-monico/coeff grandi. AGGIRATO: τ,A|B da RiemannSurface di **Sage**;
  Abel map in avanti u(r) per integrazione diretta; θ+derivate da abelfunctions.
- Verifiche: A⁻¹B−τ=1e-16; u(r₀)=0; ∇θ(0)=1.6e-16 (θ pari); Hessiana simmetrica≠0.
- Script: `KerrMetric/kerr_psi_forward_abel.sage`.

---

## 7. Nel paper (main.tex + main_prd_revtex.tex)

- Sottosezione `sec:adiabatic` "Semi-analytic first-order adiabatic orbit shape":
  eq. ibrida, riduzione 2ª specie, ψ irriducibile, tabella validazione t+τ.
- Framing onesto: ψ NON riducibile a depth-1 (shuffle) ma è polilog genus-2.
- Acknowledgements + bibliografia software: SageMath, abelfunctions (con link),
  NumPy/SciPy/SymPy/Matplotlib (con DOI).

**DA AGGIUNGERE (proposto):** la dicotomia polilog t/τ (§3) — risultato pulito:
proper-time = curva pulita, coord-time = cover frame-dragging.

---

## 8. Cosa resta

1. (opzionale) Figura curve φ(r,A) x-y assemblate dai pezzi analitici, t e τ.
2. Separatrice genus-1: ψ chiude in dilogaritmo ellittico tabulato (win chiuso).
3. Kernel Kronecker-Eisenstein genus-2 named per ψ generico (research-grade).
4. Scrivere la dicotomia polilog t/τ nel paper.

File collegati: `../progress.md` (§8,§9), `../KerrMetric/doranTau.md`,
`../KerrMetric/doranT.md`, script in `../KerrMetric/` e `.`.
