# Algoritmo: forma chiusa esplicita di δφ adiabatico (1° ordine) in funzioni speciali

Procedura **deterministica e riutilizzabile** per ottenere la correzione adiabatica al
prim'ordine δφ come **forma chiusa in funzioni Brown-Levin Γ̃** (dilogaritmi ellittici),
con **coefficienti analitici** (residui, mai fit) e **ogni passo verificato numericamente
prima di costruirci sopra**.

Validata su: separatrice Vaidya, ramo τ (Σ 63 termini = diretto a 9e-10).
Da trasportare a: Thakurta-Kerr, e (con genus-2/Schottky) J-generico.

Riferimenti metodo: Hermite–Ostrogradsky (riduzione algebrica); Mittag-Leffler sul toro;
Chen (integrali iterati); Brown-Levin / Broedel-Duhr-Weinzierl (polilog ellittici);
Kronecker-Eisenstein g^(n); Weinzierl arXiv:2602.09956 (GiNaC eMPL).

---

## Regola d'oro (disciplina, non opzionale)

**Ogni passo produce un oggetto che va VERIFICATO numericamente (≥1e-8) contro la sua
definizione diretta PRIMA di usarlo nel passo successivo.** I check catturano i bug di
codice (segni, mappe di indici, precisione). Se un check fallisce: FERMARSI, isolare, NON
proseguire e NON aggiustare a mano i numeri.

---

## Input

- Metrica → polinomio radiale **S(r)** (da y²=S, brachistocrona vincolata), funzionale
  d'onda **F = K/√S** (K = numeratore del ramo: t, τ, v).
- Parametro adiabatico λ (massa m, o massa Vaidya, ecc.) e clock **η = Σ_j b_j U_j**
  (U_j = ∫ r^j/√S dr; b_j fissati dal clock: τ → η=U_3−2mU_2, ecc.).
- Costanti fisiche (M, E, J).

---

## PASSO 0 — Frozen + condizione di separatrice (se serve genus-1)

- Genus generico: y²=S grado 6 (genus 2). Per lavorare in genus-1 serve la **separatrice**.
- **Condizione separatrice**: S ha radice doppia ⇔ **Res_r(S, S') = 0** ⇔ S(r_d)=S'(r_d)=0.
  Risolvere in J → **J_c** (forma chiusa: es. Vaidya J_c=5√(3011/3072+581√249/9216)).
- Fattorizzare **S = (r−r_d)² Q₄(r)**, Q₄ grado 4 (numpy.polydiv, controllare resto ~1e-12).
  √S = (r−r_d)√Q₄. Radici di Q₄ = e₁<e₂<e₃<e₄ (branch points), leading a₄.
- **Check**: Q₄ ricostruisce S; e_i reali.

---

## PASSO 1 — Toolkit Weierstrass dal reticolo (θ₁)

Da e_i, a₄ costruire il reticolo e σ,ζ,℘ ESATTI (mpmath.jtheta):
```
k² = (e3-e2)(e4-e1)/[(e4-e2)(e3-e1)];  pref = 2/√((e4-e2)(e3-e1))/√a4
ω1 = pref·K(k²) ;  w_im = pref·K(1-k²) ;  τ = i·w_im/ω1 ;  q = e^{iπτ}
η1 = -(π²/12ω1)·θ₁'''(0,q)/θ₁'(0,q)
σ(z)= (2ω1/π) e^{η1 z²/2ω1} θ₁(u)/θ₁'(0),  u=πz/2ω1
ζ(z)= η1 z/ω1 + (π/2ω1) θ₁'(u)/θ₁(u)
℘(z)= -η1/ω1 - (π/2ω1)² [θ₁''(u)/θ₁(u) - (θ₁'(u)/θ₁(u))²]
```
- Coord di Abel: **z(r) = ∫ dr/√Q₄** (base e₄ → z=0). Immagini: z_∞=∫_{e4}^∞, z_d=z_∞+∫_{-∞}^{r_d}.
- **r(z) ESPLICITA**: r(z) = c_r − (1/√a4)[ζ(z−z_∞) − ζ(z+z_∞)], c_r = e₄ − (2/√a4)ζ(z_∞).
- **Check**: r(z(r_test)) = r_test (1e-11).

---

## PASSO 2 — Blocchi weight-1 (U_k, primitive) in σ,ζ,℘

- U_0 = ρ[lnσ(z−z_d) − lnσ(z+z_d)] + C₀ z  (3ª specie a r_d), ρ=1/√Q₄(r_d).
- V_1,V_2 espliciti (3ª/2ª specie a z_∞); **V_{≥3} via ricorsione di Hermite**:
  (2k+4)a₄V_{k+3}+(2k+3)b₃V_{k+2}+(2k+2)b₂V_{k+1}+(2k+1)b₁V_k+2k b₀V_{k-1}=2r^k√Q₄
  ([a₄,b₃,b₂,b₁,b₀]=coeff Q₄). U_k = Σ_{i<k} r_d^{k-1-i}V_i + r_d^k U_0.
- **Check**: U_k(esplicito) = ∫r^k/√S dr (quadratura), 1e-9.

---

## PASSO 3 — Ridurre ∂_λF alla curva: la funzione ellittica R(z)

- ∂_λF = N_λ / S^{3/2}, con **N_λ = S ∂_λK − ½ K ∂_λS** (polinomio in r).
- Sul toro: dr = √Q₄ dz ⇒ **∂_λF dr = R(z) dz**, con
  **R(z) = N_λ(r(z)) / [ (r(z)−r_d)³ · Q₄(r(z)) ]**  (funzione ELLITTICA in z).
- **Check**: R(z) = ∂_λF(r(z))·√Q₄(r(z)) (1e-10).
- ⚠️ NON usare lo split polinomiale eq:psi-split (𝒜 polinomiale): è SINGOLARE alla
  separatrice (S'(r_d)=0 ma N_λ(r_d)≠0 → identità inconsistente). Usare le parti principali.

---

## PASSO 4 — Parti principali di R e η' (Mittag-Leffler, ALGORITMICO)

Ogni funzione ellittica f = C + Σ_a [ b₁ᵃ ζ(z−a) + b₂ᵃ ℘(z−a) − (b₃ᵃ/2) ℘'(z−a) + ... ].

**4a. Trovare i poli** = immagini z dei punti dove il denominatore si annulla:
- di R: z = ±z_d (ordine 3, da (r−r_d)³); z = semiperiodi immagine di e_i con Q₄=0 (ordine 2).
  ⚠️ La mappa e_i → semiperiodo NON è ovvia: e₄→z=0; determinare gli altri (es. Vaidya:
  e₃→z=i·w_im) confrontando il residuo analitico col contorno. Dove N_λ(e_i)=0 → niente polo.
- di η'=(numeratore)/(r−r_d): z=±z_d (ordine 1); z=±z_∞ (ordine = grado numeratore all'∞).

**4b. Coefficienti Laurent b_nᵃ — DUE VIE (devono coincidere):**
- *Contorno* (robusto, sempre): b_n = (1/2πi)∮_{|z−a|=ε} f(z)(z−a)^{n-1}dz  (mpmath.quad).
- *Analitico* (forma chiusa, da nominare i coefficienti): espansione locale.
  Polo semplice a r_d di g(r)/(r−r_d): Res = g(r_d)/√Q₄(r_d) (segno ± sui due fogli ±z_d).
  Polo triplo di N/((r−r_d)³Q₄): con r=r_d+s·w+…, s=√Q₄(r_d), a₁=Q₄'(r_d)/4s, a₂=Q₄''(r_d)/12,
  F=N/Q₄, h₀=F(r_d), h₁=F'(r_d)s, h₂=½F''(r_d)s²+¼F'(r_d)Q₄'(r_d):
    **b₃=h₀/s³, b₂=(h₁−3a₁h₀)/s³, b₁=(h₂−3a₁h₁+(6a₁²−3a₂)h₀)/s³.**
  Polo doppio a semiperiodo (Q₄(e_i)=0): **b₂ = 4N(e_i)/[(e_i−r_d)³ Q₄'(e_i)²].**
- **C (costante)**: C = f(z*) − Σ parti principali, a z* regolare.
- **Check**: f(z) = C + Σ parti principali (1e-8, tutti i punti). ⚠️ dps≥40 vicino ai poli.

---

## PASSO 5 — δφ come integrale iterato length-2 (Chen)

Senza IBP / senza bordo:
  **δφ = ∫_{z0}^z R(z') η(z') dz' = ∫_{z0}^z R(z') [∫_{z0}^{z'} η'(s)ds] dz'**
= integrale iterato length-2 delle due 1-forme ω_R=R dz (esterna) e ω_η'=η' dz (interna).

---

## PASSO 6 — Espansione bilineare → somma di Γ̃ (PROGRAMMATICO, no algebra a mano)

Con R = Σ_i R_i f₁ⁱ (f₁∈{1,ζ_a,℘_a,℘'_a}) e η' = Σ_j η'_j f₂ʲ (f₂∈{1,ζ_b,℘_b}):
  **δφ = Σ_{i,j} R_i η'_j · J[f₁ⁱ, f₂ʲ]**,  con
  **J[f₁,f₂] = ∫_{z0}^z f₁(z') [∫_{z0}^{z'} f₂(s) ds] dz'**  (iterato length-2 di forme Weierstrass).
- Costruire i termini con un **doppio loop** sulle liste (coeff,kind,polo). MAI enumerare a mano.
- **Check**: Σ termini = δφ diretto = ∫∂_λF·η dr (1e-9).

---

## PASSO 7 — Nominare ogni J come Brown-Levin Γ̃(m,n) (tabulabile)

Mappa Weierstrass ↔ Kronecker (û=z/2ω1, ESATTA):
  ζ(z) = (η1/ω1)z + (1/2ω1) g⁽¹⁾(û);  ℘,℘' → g⁽²⁾,g⁽³⁾ (derivando);  lnσ = quadr. + logθ₁(πû).
Kernel di Kronecker (**serie-q TABULABILE**):
  **g⁽¹⁾(x) = π cot(πx) + 4π Σ_{n≥1} q^{2n}/(1−q^{2n}) sin(2πn x)** ,  q=e^{iπτ}.
Funzione Brown-Levin length-2:
  **Γ̃(m,n; x₁,x₂; y) = ∫_0^y g⁽ᵐ⁾(t−x₁) [∫_0^t g⁽ⁿ⁾(s−x₂) ds] dt.**
- Ogni J[f₁,f₂] = Γ̃(m,n; â,b̂) + peso-1 (ζ→m=1, ℘→2, ℘'→3, 1→0; â=a/2ω1).
- **Check**: shuffle Γ̃(x₁,x₂)+Γ̃(x₂,x₁) = [∫g⁽¹⁾][∫g⁽¹⁾] (1e-31); g⁽¹⁾ serie-q = θ₁ (1e-31);
  nucleo dilog ∫logθ₁·g⁽¹⁾ = Γ̃(1,1)+peso-1 (1e-32).

**Risultato**: δφ = (peso-1 elementare) + Σ (residui analitici) · Γ̃(m,n; â,b̂), tutto tabulabile.

---

## Trappole viste (warning per la prossima volta)

1. Selezione radice: `sp.solve(...)[0]` può dare J=0 spurio → filtrare reale >soglia.
2. Punti test FUORI dall'orbita fisica (Q₄<0) → z complesso, quad fallisce. Restare in (e₄,r₀].
3. Segno delle basi ℘,℘': fissare fbasis = funzione *pura* e mettere il segno SOLO nel coeff.
4. **Mappa e_i→semiperiodo**: verificarla col contorno, non indovinarla.
5. **Precisione**: dps≥40 per le doppie quadrature vicino ai poli (dps=30 dà ~1e-4 di errore).
6. Split polinomiale eq:psi-split SINGOLARE alla separatrice → usare parti principali (Passo 4).
7. IBP indipendente da costante additiva di G̃: usare G̃ grezzo (se si usa la forma con bordo).

---

## Estensione a Thakurta-Kerr e genus-2

- **Thakurta-Kerr (a≠0), separatrice**: identico Passo 0-7, con S di TK (ramo t on-curve;
  attenzione alla lettera ∞ iperellittica = terza specie non elementare → Fay, non ln).
- **J-generico (genus 2)**: Passi 1,4,5,6,7 si generalizzano con **Kleinian σ,ζ,℘ (genus 2)**
  e **Schottky-Kronecker** (Broedel-Zerbini arXiv:2406.10051) al posto di Weierstrass+g⁽ⁿ⁾.
  NON coperto da GiNaC (che è genus 1). Le curve J-generiche sono comunque già analitiche in
  *rappresentazione integrale* (verificate); questo pipeline le porta a funzioni speciali esplicite.

---

## File prodotti (separatrice Vaidya, ramo τ)

- `VaidyaMetric/vaidya_separatrix_explicit_Uk.py` — U₀..U₅ (Passo 2).
- `VaidyaMetric/vaidya_sep_G_partialfractions.py` — R, parti principali, G (Passo 3-4).
- `VaidyaMetric/vaidya_sep_weight2_closure.py` — chiusura weight-2 programmatica (Passo 6).
- `VaidyaMetric/vaidya_sep_residui_analitici.py` — residui/coeff ANALITICI (Passo 4b).
- `VaidyaMetric/brown_levin_gamma.py` — Γ̃, serie-q, mappa Weierstrass↔Kronecker (Passo 7).
- `VaidyaMetric/vaidya_sep_deltaphi_full_gamma.py` — **formula 100% termine-per-termine** (Passo 5-6).
- Dettagli e valori: `EXPLICIT_FORMS_PROGRESS.md`.
