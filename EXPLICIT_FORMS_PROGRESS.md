# Progresso: forme chiuse esplicite delle brachistocrone adiabatiche

Stato al: sessione in corso. Obiettivo: **forma chiusa esplicita in FUNZIONI SPECIALI**
(anche non elementari) delle correzioni adiabatiche al primo ordine δφ, testabile
contro la numerica. Metodo: **tutto analitico, coefficienti da residui/algebra, mai
fit; ogni passo verificato numericamente prima di costruirci sopra.**

---

## 1. Mappa 2×2 e stato

| | Vaidya (a=0) | Thakurta-Kerr (a≠0) |
|---|---|---|
| **J generico** (genus 2) | δφ assemblato (v,τ) ✓ rappr. integrale, testato 1e-13/1e-15. Forma in funz. speciali iperellittiche = DA FARE | δφ assemblato (t,τ) ✓ testato 5e-14/9e-13. Idem |
| **Separatrice** (genus 1) | **← QUI.** Settore abeliano/2ª specie ESPLICITO in σ,ζ,℘ ✓. Weight-2 = polilog ellittico (frontiera) | non iniziato (stessa macchina) |

**Brachistocrone adiabatiche al 1° ordine — tutte fatte e verificate (rappr. integrale):**
- TK τ: 5e-14 | TK t (on-curve, cover spurio corretto): 9e-13
- Vaidya v (tempo avanzato): 2e-13 | Vaidya τ (tempo proprio): 9.77e-15

---

## 2. Separatrice Vaidya: dati della curva (M=1, E=1.4, J variabile; separatrice a Jc)

- **Jc = 7.02662374** (radice doppia bracket), **r_d = −3.3637111** (radice doppia)
- Curva ellittica **E: w²=Q4(r)**, radici **{−2.0833, 0, 2, 8.7274}**, a4=0.96
- **k²=0.60672, τ=0.9059733802550 i** (Legendre = Sage, 15 cifre)
- Semiperiodi: om1=0.66913 (reale), om_im=0.60621
- z_∞ (immagine ∞) = 0.23663; z_d (immagine r_d) = 0.46104; z_h (2-torsione orizzonte) = 0.60621 i
- c_r = 0.13406; ρ = 1/√Q4(r_d) = 0.06107; C0 = 0.33327
- Turning fisico sulla separatrice = e4 = 8.7274 → orbita r∈(8.73, 12], z reale
- σ,ζ,℘ da θ₁ (mpmath.jtheta), ℘ ESATTO da θ₁''; reticolo auto-consistente

---

## 3. FATTO: settore abeliano / 1ª-2ª-3ª specie (ESPLICITO, verificato)

Tutte funzioni speciali esplicite (Weierstrass σ,ζ,℘ + algebrico), NO integrali in
quadratura. Coefficienti analitici (residui + ricorsione Hermite). Verifica 1e-9…1e-14.

**Formule (LS ≡ ln[σ(z−z_∞)/σ(z+z_∞)]):**
- U₀ = ρ[lnσ(z−z_d) − lnσ(z+z_d)] + C0·z      (3ª specie a r_d) — 1e-14
- V₁ = c_r·z − (1/√a4)·LS                       (3ª specie ∞)
- V₂ = c_r²z − (2c_r/√a4)LS + (1/a4)[−2ζ(2z_∞)LS − ζ(z−z_∞)−ζ(z+z_∞) + C·z]  (2ª specie)
- V₃,V₄,V₅ : **ricorsione di Hermite** da d/dr(r^k√Q4):
    (2k+4)a4 V_{k+3}+(2k+3)b3 V_{k+2}+(2k+2)b2 V_{k+1}+(2k+1)b1 V_k+2k b0 V_{k-1}=2 r^k√Q4
  con [a4,b3,b2,b1,b0] = coeff di Q4 (np.poly(radici)*a4)
- **U_k = Σ_{i=0}^{k-1} r_d^{k-1-i} V_i + r_d^k U₀**   (dallo split r^k/(r−r_d))
- **Π_h = β ζ(z−z_h) + γ z**,  β=−4/Q4'(2m)=0.07584, γ=−0.198751  (2ª specie orizzonte) — 1e-14

**Script:** `VaidyaMetric/vaidya_separatrix_explicit_Uk.py` (U₀..U₅, verif 1e-9),
`VaidyaMetric/vaidya_separatrix_Pih.py` (Π_h, 1e-14),
`VaidyaMetric/vaidya_ell_dilog_match.py` (U₀ σ,ζ + toolkit).

---

## 4. FRONTIERA: settore weight-2 (polilog ellittico) — NON ancora chiuso

Struttura (fatta l'algebra): W_kj = 2∫U_k dU_j − U_k U_j; ogni termine si riduce a
σ,ζ,℘ TRANNE **∫lnσ(z−a) ζ(z−b) dz = DILOGARITMO ELLITTICO** (Brown-Levin Γ̃ / Zagier D^E).
Analogamente D_k = U_k ln(r−2m) − G_k, G_k=∫ln(r−2m) r^k/√S dr = polilog ellittico.

**Cosa ho implementato (framework, verificato parzialmente):**
- kernel di Kronecker g^(1)(ξ)=θ₁'/θ₁·π + termine Im (single-valued), dispari ✓
- Zagier D^E(z)=Σ_{n≥0}D(qⁿζ)−Σ_{n≥1}D(qⁿ/ζ), e Li₂^ell olomorfo
- coordinate normalizzate ξ_h=τ/2=0.453i, ξ_d=0.3445

**Cosa NON chiude:** il **regulator preciso** — G₀ (reale, ~−0.001…−0.006) non ha
relazione lineare semplice con D^E/Li₂^ell (che sono ordine 2-7). La combinazione
esatta (poli a z_h, z_d E ∞, + correzioni single-valued) richiede la **macchina
completa di Brown-Levin**. È la **frontiera research 2024** (Broedel, Zerbini,
Schottky-Kronecker). NON truccato: il check numerico dice che non torna.

**Script tentativi (framework):** codice inline nelle ultime run (non salvato come file
definitivo — il regulator non chiude ancora).

### AGGIORNAMENTO: GiNaC dilog ellittico VALIDATO (1e-18)

GiNaC 1.8.10 (+CLN 1.3.7, Rosetta `-arch x86_64`) HA il toolkit polilog ellittico
(Weinzierl, arXiv:2602.09956) built-in — **eMPL separato NON serve**.
- Toolchain OK: `iterated_integral(lst{kernels}, y[, N_trunc])`, primo kernel = int. ESTERNA.
- **Normalizzazione FISSATA (verif 1.87e-7, limite solo troncamento):**
  `Kronecker_dz_kernel(2, z_j, tau)` = g^(1)(z-z_j) dz, con
  **g^(1)(z,tau) = d/dz log theta1(pi z, q), reticolo (1,tau), q=e^{i pi tau}.**
  (Mappa alla mia coord Weierstrass: z_ginac = z_phys/(2 om1), tau = w_im/om1 = 0.90597338 i.)
- **Length-2 (dilog ellittico) VERIFICATO 4.38e-18** vs quadratura annidata mpmath:
  I2 = int_0^y g^(1)(t1-za) int_0^{t1} g^(1)(t2-zb). Ordering: 1° kernel lst = esterno.
- Vincolo pratico: polo di g^(1) a z=z_j NON deve stare sul cammino [0,y] (z_j immaginario
  ok); y dentro raggio di convergenza serie (y~0.15 ok, y~0.40 diverge).
- **Compile:** `export PKG_CONFIG_PATH=/usr/local/lib/pkgconfig:$PKG_CONFIG_PATH;
  clang++ -std=c++17 -arch x86_64 f.cpp -o out $(pkg-config --cflags --libs ginac);
  DYLD_LIBRARY_PATH=/usr/local/lib ./out`. Script test in /tmp/ginac_ell/.
- **PROSSIMO:** assemblare W_kj = 2 int U_k dU_j - U_k U_j riducendo ogni ζ->g^(1)+lineare,
  lnσ->int g^(1)+quadratico, pezzo irriducibile = Gamma-tilde(1,1;z_a,z_b;z_end) length-2.
  Poi testare W_kj assemblato vs quadratura diretta.

### NOTA STRATEGICA (mappa 2×2, cosa serve davvero)
- Curve adiab. 1° ordine J-GENERICO: **già analitiche+verificate in rappr. integrale**
  (Vaidya v,τ; TK t,τ). NON richiedono questi conti genus-1.
- Questi conti (GiNaC, genus 1) servono SOLO a portare la **separatrice** a forma esplicita
  in funz. speciali. J-generico esplicito richiederebbe **genus 2 = Schottky-Kronecker**
  (arXiv:2406.10051), NON in GiNaC.

---

## 4bis. WEIGHT-2 CHIUSO + assemblaggio separatrice (VERIFICATO)

**Weight-2 W_kj CHIUSO** (`VaidyaMetric/vaidya_sep_weight2_assembly.py`):
- r(z) = c_r - (1/√a4)[ζ(z-z_∞)-ζ(z+z_∞)] esplicito (elliptic, verif 1e-11)
- f_k(z)=U_k'(z)=r(z)^k/(r(z)-r_d) esplicito (1e-10)
- **W_kj = 2∫U_k dU_j - U_k U_j**, con U_k funz. speciali esplicite (σ,ζ+Hermite) e
  ∫U_k dU_j = **polilog ellittico length-2** (= iterated_integral Kronecker g^(1), verif
  GiNaC 1e-18). W_kj asm vs quadratura diretta: **1e-12…1e-15** su tutte le coppie
  (2,3),(0,2),(0,3),(1,3),(1,2). → il §4 "regulator non chiude" È RISOLTO.

**φ_0|_sep (orbita frozen sulla separatrice) ESPLICITA** (verif 1e-12):
  **φ_0|_sep = Jc[(E²-1) z + ((E²-1)r_d + 2m) U_0(z)]**,  z=∫_{r0}^r dr/√Q4,
  U_0=ρ[lnσ(z-z_d)-lnσ(z+z_d)]+C0 z.

**Jc (condizione separatrice) ESATTO:** doppia radice di S ⇔ Res_r(S,S')=0 ⇔
  S(r_d)=S'(r_d)=0.  **Jc = 5√(3011/3072 + 581√249/9216) = 7.026623740389046**,
  r_d=-3.36371118. (M=1,E=1.4.)

**δφ_τ|_sep(r) ESPLICITO E VERIFICATO end-to-end (1e-7…1e-9)** — MICROSTEP FINALE FATTO.
- Lo SPLIT eq:psi-split (𝒜 polinomiale) è SINGOLARE a Jc (S'(r_d)=0 ma N_m(r_d)≠0 →
  identità 2𝒜'S-𝒜S'+2SM=2N_m inconsistente). CONFERMA §5. Ansatz a gradi fissi NON risolve.
- **RISOLTO per ALGORITMO (Hermite-Ostrogradsky via parti principali sul toro, NIENTE ansatz):**
  `VaidyaMetric/vaidya_sep_G_partialfractions.py`.
  In z: ∂_mF dr = R(z)dz, **R(z)=N_m(r(z))/((r(z)-r_d)³ Q(r(z)))** funzione ellittica.
  Parti principali via contorno b_n^a=(1/2πi)∮R(z)(z-a)^{n-1}dz (analitico). Poli:
  ±z_d (ord.3), z=0 e z=i·w_im (ord.2). Integrazione termine-a-termine:
  **G(z)=∫∂_mF = C₀ z + Σ_a[ b₁^a lnσ(z-a) − b₂^a ζ(z-a) − (b₃^a/2)℘(z-a) ]**
  (usa ∫ζ=lnσ, ∫℘=−ζ, ∫℘'=℘). Verif dG/dz=R a 1e-13.
  Coeff (M=1,E=1.4): C₀=−0.084098;  z_d: b1=+0.27045,b2=+0.032564,b3=+0.0098730;
  −z_d: b1=−0.27045,b2=+0.032564,b3=−0.0098730;  z=0: b2=+0.033998;  z=i·w_im: b2=−0.38946.
- **δφ_τ|_sep(r) = G~(z(r))·η(r) − ∫_{r0}^r G~ η' dr**, η=U_3−2U_2,
  η'-forma=(r³−2r²)/(r−r_d). ASSEMBLATO (G esplicito) = DIRETTO a **1e-7…1e-9** a ogni r.
- **WEIGHT-2 CHIUSO in atomi nominati** (`VaidyaMetric/vaidya_sep_weight2_closure.py`):
  ∫G~ η' dz NON è un integrale irrisolto — è somma di **dilog ellittici**. η' scomposto in
  parti principali (contorno, analitico): poli ±z_d(ord1), ±z_∞(ord2). Prodotto G~·η'
  costruito PROGRAMMATICAMENTE (no algebra a mano) → ogni coppia kind_i×kind_j = 1 atomo:
    **D(a,b)=∫lnσ(z−a)ζ(z−b)dz** (dilog ellittico, = iterated_integral GiNaC, verif 1e-18),
    C(a,b)=∫ζζ, + atomi ℘/z (weight-1 o dilog). Coeff = prodotti di parti principali (ANALITICI).
  Atomi dilog irriducibili: D(±z_d, {±z_d,±z_∞}) con coeff ±1.00232, ±0.28752.
  **Somma atomi = ∫diretto a 1e-9**; **δφ COMPLETO = G~η − Σatomi = diretto a 1.6e-8…4.3e-8.**
  Precisione limitata solo da eps contorno/wpp; nessun fit. FORMA CHIUSA ESPLICITA COMPLETA.
- **COEFFICIENTI DILOG NOMINATI ANALITICAMENTE** (`VaidyaMetric/vaidya_sep_residui_analitici.py`):
  c_ab = Res_a(R)·Res_b(η'), residui in forma chiusa dai dati della curva (verif vs contorno 1e-7..1e-16):
  · **Res_{±z_d}(η') = ±(r_d³−2r_d²)/√Q4(r_d)**  (residuo terza specie a r_d)
  · **Res_{±z_∞}(η'): b₁=∓(1/√a4)(2B+r_d−2), b₂=1/a4** (B=c_r+ζ(2z_∞)/√a4) [b₂=1/a4 ESATTO]
  · **Res_{z_d}(R)** (polo triplo, espansione r=r_d+s·w+..., s=√Q4(r_d), a₁=Q4'(r_d)/(4s), a₂=Q4''(r_d)/12,
    F=N_m/Q4, h₀=F(r_d), h₁=F'(r_d)s, h₂=½F''(r_d)s²+¼F'(r_d)Q4'(r_d)):
      b₃=h₀/s³,  b₂=(h₁−3a₁h₀)/s³,  **b₁=(h₂−3a₁h₁+(6a₁²−3a₂)h₀)/s³**  (Res=b₁).
  Coeff dilog: **D(±z_d,±z_d)=∓1.002315, D(±z_d,±z_∞)=±0.287525** (da formule, non fit).
- **DILOG ELLITTICI NOMINATI:** D(a,b)=∫lnσ(z−a)ζ(z−b)dz = **Brown-Levin Γ̃(1,1;a,b)** = length-2
  Kronecker iterated integral (= GiNaC iterated_integral, verif 1e-18). Poli a,b∈{±z_d,±z_∞}.
  Valori (a r=10): D(z_d,z_d)=−0.11236, D(z_d,−z_d)=+0.06308, D(z_d,z_∞)=−0.31209, D(z_d,−z_∞)=+0.11733,
  D(−z_d,z_d)=−0.06489, D(−z_d,−z_d)=+0.03657, D(−z_d,z_∞)=−0.17968, D(−z_d,−z_∞)=+0.06804.
- **CAVEAT (fisico, già noto):** l'angolo TOTALE diverge al turning: ∂_mF~(r−rt)^{-3/2}
  → δφ~(r−rt)^{-1/2}. È il breakdown adiabatico al punto di svolta (vale per OGNI orbita,
  non solo separatrice; = "Domain of validity" caveat nel paper). δφ(r) accumulato fino a
  r>rt è la quantità esplicita ben definita. NON esiste un singolo numero "δφ|_sep"; il
  "3.7302" precedente era artefatto di cutoff (r=8.9) — RITIRATO.

## 5. Findings / correzioni importanti (verificati)

- **L_2m è SECONDA specie**, NON Fay: 𝒜^m(2m)=−19.6≠0 → integrando ~(r−2m)^{−3/2},
  polo doppio senza residuo → Kleinian ζ. Il genuino 3ª specie orizzonte è la lettera
  dr/(r−2m)=tortoise ln (elementare). **Fay serve in TK (lettera ∞, iperellittica
  r²/√S, non elementare), NON in Vaidya (lettera orizzonte razionale → ln elementare).**
  Paper corretto.
- **t-branch TK: clock ON-CURVE**, il "cover" di frame-dragging era SPURIO. Identità
  (r−2M)Q2+DE(J(r−2M)+2Ma)²=E²r³Δ ⇒ n_t α=E²r³/(f√R6). Corretto + verificato 1e-10.
- **δφ|_sep NON si ottiene mettendo Jc nella riduzione genus-2** (SINGOLARE lì, de Rham
  degenera). Va fatto come limite J→Jc o riduzione ellittica diretta. (L_2m diretto a Jc
  dava ~1e10, assurdo.)
- 𝓘_poly reso esplicito: polinomio + Σ_ρ res·log(r−ρ) sui 6 poli di S, verif 1e-15.

---

## 6. Software per chiudere il weight-2 (in corso di installazione)

Il regulator/polilog ellittico esplicito richiede strumenti dedicati:
- **GiNaC** (C++, `brew install ginac` → 1.8.10 + dipendenza CLN) — INSTALLARE
- **eMPL** (elliptic MPL con argomenti arbitrari): arXiv **2602.09956**, in GiNaC.
  Codice negli ancillary files del paper o repo linkato — DA RECUPERARE
- (genus-2 futuro) **Schottky-Kronecker**: arXiv 2406.10051 (Broedel-Zerbini), ancillary
- PolyLogTools (gitlab.com/pltteam/plt) è genus-0, NON serve al dilog ellittico
Ambiente: clang++ ✓, brew ✓, pkg-config ✓. GiNaC/CLN da installare.

---

## 7. DA DOVE RIPRENDERE

1. **Installare GiNaC** (`brew install ginac`) + recuperare il codice **eMPL** (arXiv:2602.09956).
2. Usare eMPL per valutare/chiudere il **dilog ellittico** ∫lnσ ζ (il regulator) — l'unico
   pezzo mancante della separatrice. Verificare vs G₀ diretto.
3. Assemblare **δφ|_sep = [σ,ζ,℘ espliciti (§3)] + [polilog ellittico (§4)]** come formula
   unica, via **limite J→Jc** (non riduzione a Jc, §5). Testare vs numerica.
4. Trasportare a **J generico (genus 2)** con Schottky-Kronecker, e a **TK** (stessa macchina).
5. Scrivere nel paper la sezione separatrice tutta-esplicita con citazioni (Fay, Baker,
   Buchstaber-Enolskii-Leykin, Brown-Levin, Beilinson-Levin, Zagier, Schottky-Kronecker).

**In sintesi:** settore abeliano ESPLICITO e verificato (§3); manca solo il **dilog
ellittico** del weight-2 (§4), che è frontiera e richiede eMPL/GiNaC. δφ|_sep va
assemblato via limite J→Jc.

## 4ter. BROWN-LEVIN Gamma-tilde SCRITTA ESPLICITA + serie-q TABULABILE (VERIFICATO)

Script `VaidyaMetric/brown_levin_gamma.py`. Risponde a "abbiamo scritto Γ̃ e chiuso con essa?".

**g^(1)(x,τ) = ∂_x log θ₁(πx, q), q=e^{iπτ}, reticolo (1,τ).** Serie-q TABULABILE:
  **g^(1)(x) = π cot(πx) + 4π Σ_{n≥1} [q^{2n}/(1−q^{2n})] sin(2πn x)**   (verif vs θ₁: 1e-31)

**Γ̃(1,1; x1,x2; y) = ∫_0^y g^(1)(t−x1)[∫_0^t g^(1)(s−x2)ds]dt**  (length-2 Kronecker).
  Shuffle Γ̃(x1,x2)+Γ̃(x2,x1)=[∫g1(·−x1)][∫g1(·−x2)] verif 1e-31 (⇒ è l'iterato corretto).

**Mappa Weierstrass↔Kronecker (reticolo Vaidya, ESATTA 1e-31):** con û=z/(2ω1), q=e^{iπτ}, τ=w_im/om1:
  ζ_W(z) = (η1/ω1) z + (1/(2ω1)) g^(1)(z/(2ω1))
  lnσ_W(z) = ln(2ω1/π)−lnθ₁'(0) + η1 z²/(2ω1) + log θ₁(π z/(2ω1))

**Nucleo trascendente del dilog D(a,b) = Γ̃(1,1) + peso-1 (ESATTO 4.5e-32, tutti i poli):**
  ∫_{u0}^u lnθ₁(π(û'−â))·g^(1)(û'−b̂)dû' = Γ̃(1,1;â,b̂;u0,u) + lnθ₁(π(u0−â))·[lnθ₁(π(u−b̂))−lnθ₁(π(u0−b̂))]
  (â=a/(2ω1), b̂=b/(2ω1)). [Attenzione: serve dps≥40 per la doppia quadratura vicino al polo z_∞.]

**CHIUSURA:** D(a,b)=∫lnσ(z−a)ζ(z−b)dz = **Γ̃(1,1;â,b̂)** (coeff 1) + [E₂=∫lnθ₁ (dilog ell.
proprio) + poly×g1 + elementari]. Quindi δφ_τ|_sep = G~η − Σ c_ab[Γ̃(1,1;â,b̂)+peso-1] − (atomi C,℘)
è CHIUSA in funzioni Brown-Levin {Γ̃(1,1), E₂} + elementari, TUTTE con serie-q tabulabile.
Atomi C(a,b)=∫ζζ e ℘-prodotti = stessa famiglia (Γ̃(n1,n2) ordine sup.), riduzione identica.

## 4quater. FORMULA 100% TERMINE-PER-TERMINE (VERIFICATA 1e-9)

Script `VaidyaMetric/vaidya_sep_deltaphi_full_gamma.py`. Forma SENZA bordo (no IBP):
  **δφ_τ|_sep = ∫_{z0}^z R(z')η(z')dz' = Σ_{i,j} R_i η'_j · J[f₁^i, f₂^j]**
- R(z)=Σ_i R_i f₁^i con f₁∈{1, ζ(z-a), ℘(z-a), ℘'(z-a)}, a∈{±z_d, 0, i·w_im};
  coeff R_i = residui ANALITICI (C₀=−0.084099; b₁,b₂,b₃ a ±z_d da formule §residui;
  semiperiodi: b₂ = 4N_m(e_i)/((e_i−r_d)³Q4'(e_i)²), con **e₄→z=0, e₃→z=i·w_im**;
  N_m(e₁)=N_m(e₂)=0 ⇒ niente polo a z=ω1, ω1+i·w_im).
- η'(z)=Σ_j η'_j f₂^j con f₂∈{1, ζ(z-b), ℘(z-b)}, b∈{±z_d, ±z_∞}; coeff analitici (§residui).
- **J[f₁,f₂] = ∫_{z0}^z f₁(z')[∫_{z0}^{z'}f₂(s)ds]dz'** = iterato length-2 di forme Weierstrass
  = **Brown-Levin Γ̃(m,n; â,b̂)** con mappa ordine-g: ζ→g⁽¹⁾, ℘→g⁽²⁾, ℘'→g⁽³⁾, 1→g⁽⁰⁾.
  Tipi presenti: Γ̃(m,n), m∈{0,1,2,3}, n∈{0,1,2}. m=0/n=0 → peso-1/E₂; m,n≥1 → dilog ell.
- **63 termini non nulli. Σ termini = δφ diretto a 9e-10 / 2.6e-9 / 6.5e-9** (r=11/10.5/10, dps=40).
- Reconstruction R: 3.3e-8, η': 1.6e-8. NESSUN integrale irrisolto, nessun fit.
- Ogni Γ̃ valutabile via serie-q g⁽¹⁾ (§4ter) → interamente tabulabile.

## 4quinquies. VALIDAZIONE MULTI-PARAMETRO (PASSATA)

Script `VaidyaMetric/vaidya_sep_validate_multiparam.py`. Gira l'INTERO pipeline (coeff
ANALITICI, non contorno) per vari E (M=1), ognuno con Jc/r_d/curva diversi. Σ63 termini
Γ̃ = ∫diretto ∂_mF·η, a meta' orbita:
  E=1.30: Jc=10.907, r_d=-4.606, turn=11.21, diff=4.4e-10
  E=1.40: Jc= 7.027, r_d=-3.364, turn= 8.73, diff=3.4e-9
  E=1.50: Jc= 4.987, r_d=-2.621, turn= 7.24, diff=7.1e-10
  E=1.60: Jc= 3.764, r_d=-2.128, turn= 6.26, diff=7.1e-9
  (E=1.25: turn=13.2 > finestra r0=12 → orbita fuori range, non testato; alzare r0.)
→ formule chiuse dei residui + assemblaggio Γ̃ VALIDI sull'intera famiglia, non un punto solo.

Nota "riduzione Γ̃(m,n) dei 63 termini": è la riscrittura Weierstrass→Kronecker standard;
NON compattifica (semmai espande per gli shift lineari). Compattificazione vera = sfruttare
simmetria ± (involuzione z→−z: coppie z_d/−z_d, z_∞/−z_∞ con residui ±) + relazioni shuffle
→ combina in dilog ellittici reali, riducendo 63 a poche unita' indipendenti (DA FARE se si
vuole la forma minima).

## 4sexies. COMPATTIFICAZIONE 63 -> 24 (VERIFICATA 1e-9)

Script `VaidyaMetric/vaidya_sep_compact.py`. R(z),η'(z) sono funzioni PARI (r(z) pari,
involuzione iperellittica z->-z). Poli ±z_d, ±z_∞ si combinano in BLOCCHI PARI:
  Z_a=ζ(z-a)-ζ(z+a),  P_a=℘(z-a)+℘(z+a),  P'_a=℘'(z-a)-℘'(z+a)   (tutti pari in z)
  R = C0 + b1 Z_{z_d} + b2 P_{z_d} - (b3/2) P'_{z_d} + b2^(0) ℘(z) + b2^(iw) ℘(z-i w_im)  [6 blk]
  η'= Ce + e1^d Z_{z_d} + e2^∞ P_{z_∞} + e1^∞ Z_{z_∞}                                    [4 blk]
δφ = Σ_{blkR,blkη'} coeff · J[blkR,blkη'],  J=iterato length-2.  **24 prodotti (era 63)**:
  15 dilog ellittici genuini (5×3 non costanti) + 9 con costanti C0/Ce (peso inf/E₂).
Check: R,η' ricostruiti coi blocchi pari (3.3e-8, 1.6e-8); δφ_compatto = diretto a 9e-10…3e-8
(r=11/10/9.2). Stessi coefficienti analitici, solo regruppati per parità → vale a ogni E.
Ulteriore riduzione 24->~15 possibile via relazioni shuffle tra Γ̃, ma rendimenti calanti.

## 4septies. COMPATTIFICAZIONE (a): forma ANTISIMMETRICA / shuffle (VERIFICATA 1e-9)

Script `VaidyaMetric/vaidya_sep_compact_antisym.py`. Decomposizione shuffle dell'iterato:
  J[g,h] = ½(∫g)(∫h) + ½ A[g,h],   A[g,h]=∫(g·∫h − h·∫g)dz = −A[h,g]  (dilog ell. genuino).
Sommando: la parte simmetrica si ricombina in ½(∫R)(∫η')=½ G~ η (peso-1 ELEMENTARE), e:
  **δφ_τ|_sep = ½ G~(z) η(r) + ½ Σ_{i<j} (c_i d_j − c_j d_i) A[e_i,e_j]**
con e = {1, Z_zd, P_zd, P'_zd, ℘_0, ℘_iw, P_z∞, Z_z∞} (8 funzioni PARI), c_i=coeff in R,
d_i=coeff in η'. Le coppie R×R e η'×η' danno (c_i d_j−c_j d_i)=0 → spariscono; diagonali A[g,g]=0.
**21 coppie A sopravvissute** (era 24 prodotti J): 14 dilog ellittici genuini + 7 tipo E₂ (A[1,g]).
Coeff (E=1.4): A[1,Z_zd]=+0.719, A[Z_zd,wp_iw]=−1.443, A[Z_zd,Z_z∞]=+0.288, … (tutti analitici).
δφ_antisym = diretto a 9e-10…3e-8 (r=11/10/9.2). 
Ulteriore riduzione (21→ meno) richiede identità di FAY/Kronecker (oltre lo shuffle), non fatta.

## 4octies. RIDUZIONE DI FAY: rango 5 (16 relazioni) — RISULTATO POSITIVO

Script `VaidyaMetric/vaidya_sep_fay_fast.py`. Metodo rigoroso e VELOCE: le relazioni tra i
21 A[e_i,e_j] modulo peso-1 si trovano dalle DERIVATE A'_k = e_i Pe_j − e_j Pe_i (peso-1 in
forma CHIUSA: Pe = ζ/lnσ/z, nessuna quadratura). Una relazione tra le A' si solleva ESATTA a
una relazione tra gli A (d/dz(ΣλA−ΣμW)=0 ⟹ =cost). Matrice derivate 120×75 (ben posta).
℘' ESATTO da θ₁ (θ₁'''/θ₁ − 3θ₁'θ₁''/θ₁² + 2(θ₁'/θ₁)³), non differenze finite.
- **RANGO di A modulo peso-1 = 5  ⟹  16 relazioni di Fay indipendenti.**
- Gap netto valori singolari A_perp: [1, 1.1e-4, 6.3e-6, 2.3e-6, 1.2e-6 | 9.1e-9, …] (fattore ~130).
- ⟹ la parte trascendente di δφ_τ|_sep ha solo **≈5 dilog ellittici indipendenti** (non 21/63).
- CAVEAT: i 5 SV "non nulli" spaziano 1…1e-6 → conteggio robusto ≈5 (3-5 secondo soglia).
- FATTO: il RANGO. DA FARE: estrarre la FORMA esplicita a 5 termini (5 dilog-base + coeff
  analitici delle 16 relazioni via null-space simbolico) e mapparla su Γ̃(m,n).

## 4nonies. FORMA ESPLICITA A 5 TERMINI (estrazione Fay completata, VERIFICATA)

Script `VaidyaMetric/vaidya_sep_5term.py` (con tqdm + logging). Metodo: algebra lineare
ESATTA sulle DERIVATE (forma chiusa, residuo→0, non fit). QR-pivot sceglie 5 dilog-base;
risolvo δφ_w2' = Σβ_b A_base' + Σγ V'_{peso-1} su 160 punti.
  **δφ_τ|_sep = ½ G~ η + Σ_{b=1}^5 β_b A_base^(b) + [peso-1]**
- 5 dilog-base: A[℘_0,P_z∞], A[℘_0,Z_z∞], A[Z_zd,℘_0], A[1,℘_0], A[P'_zd,P_z∞]
- β = [+0.019246, −0.001957, +0.002138, +0.000115, −0.003883]
- Residuo sistema = 1.4e-10 (condizionamento ℘'/range dinamico; relazione genuina, non fit).
- **VERIFICA end-to-end (5 A_base annidati veri + peso-1) = diretto a 9e-10…3e-8.** ✓
CATENA COMPATTIFICAZIONE: 63 →(parità) 24 →(shuffle) 21 →(Fay) **5 dilog ellittici**.
CAVEAT: β NUMERICI (da algebra lineare esatta ~1e-10), non ancora forma chiusa simbolica;
nominarli richiede le relazioni di Fay simboliche (passo ulteriore, non fatto).

## 4decies. beta SIMBOLICI (Fay analitico): TENTATO, BLOCCATO da condizionamento

Obiettivo: portare i 5 β da numerici (~1e-10) a forma chiusa simbolica.
- Via NUMERICA ad alta precisione (mpmath dps=45 + PSLQ): BLOCCATA. Vicino al turning (z→0,
  r→e4) ℘(z)~1/z² esplode → range dinamico enorme → base peso-1 (46 fn, rango effettivo ~12)
  quasi singolare; ranghi numerici instabili; qr_solve/lu_solve falliscono ("matrix numerically
  singular"). Non si raggiunge la precisione (~40 cifre) per PSLQ. Script `vaidya_sep_beta_hp.py`.
- Via ANALITICA (identità di Fay/addizione Weierstrass per ℘(z-a)ζ(z-b)): derivazione simbolica
  pesante (16 relazioni), non affrontata — error-prone, richiederebbe molto lavoro.
- ESITO ONESTO: β restano NUMERICI (~1e-10). La forma a 5 termini è comunque verificata
  end-to-end a 1e-9 (§4nonies) — solida. La chiusura simbolica dei β è irrisolta.

## 4undecies. MATHEMATICA: rango 5 CONFERMATO ESATTO (residuo 1e-79); beta a 80 cifre

Script `VaidyaMetric/vaidya_sep_fay_symbolic.wl` (wolframscript). Sfrutta precisione arbitraria
(70-80 cifre) + WeierstrassP/Zeta/Sigma built-in + teoremi di addizione.
- Residui IDENTICI a Python (b1zd=0.27044665, b2zd=0.0325642, b3zd=0.0098730), C0=−0.084098. ✓
- **Residuo del sistema = 1.4e-79 ⟹ rango ESATTAMENTE 5** (il 1.4e-10 di Python era solo
  condizionamento double; a 70 cifre svanisce). La compattificazione a 5 dilog è ESATTA.
- β estratti a 80 cifre: [0.047969, 0.012025, −0.271152, −0.201227, −0.0048875] (base come Python;
  nota: differiscono dai β Python perché lo split peso-1 è non-canonico con base W rank-deficient).
- RICONOSCIMENTO SIMBOLICO FALLITO: RootApproximant → razionali giganti (β non algebrici semplici);
  FindIntegerNullVector → solo relazione spuria 24/a4=25 (a4=0.96), nessuna relazione β↔residui.
- ⟹ β NON in forma chiusa nei residui/periodi ovvi. Serve o (a) i valori dei dilog nella base,
  o (b) la riduzione di Fay ANALITICA vera (skeleton in sez.4a del .wl: formule di addizione
  Weierstrass pronte, riduzione algebrica da sviluppare). Sage non aggiungerebbe nulla qui.
- WIN comunque: rango 5 ESATTO (era il dubbio 5-vs-6), e infrastruttura Mathematica pronta.

## 4duodecies. FAY ANALITICA - fondazione simbolica (VERIFICATA), matching finale scoped

Script `VaidyaMetric/vaidya_sep_fay_analytic.wl`. Riduzione simbolica dei dilog via addizione
di Weierstrass, costanti ai poli SIMBOLICHE. FONDAZIONE COMPLETA E VERIFICATA:
- **INSIGHT chiave**: Z_a = zeta(z-a)-zeta(z+a) = -2 zeta(a) + P'(a)/(P(z)-P(a)) -> RAZIONALE in p!
- **Tutti gli 8 e_i sono razionali in p=P(z), pp=P'(z)** (verificato vs diretto 1e-62).
  P(z-+c) = -p-P(c)+1/4((pp-+ppc)/(p-P(c)))^2; Pp_zd via dz[]; ecc.
- **Le 8 primitive Pe_i si riducono a 6 LETTERE** {z, zeta(z), lnsigma(z-+zd), lnsigma(z-+zi)}
  + parte razionale in p (verificato: differenza costante in z, 1e-71).
  (usa zeta(z-+c)=zeta(z)-+zeta(c)+1/2(pp+-ppc)/(p-P(c)); P'(iw)=0.)
- ⟹ A'_k = e_i Pe_j - e_j Pe_i = (razionale in p) x (lettera): struttura ESATTA, niente sampling.
MATCHING FINALE (scoped, NON completato): la relazione Sum lambda_k A'_k = W' (W peso-1) diventa,
per ciascuna lettera, la condizione che il coeff (razionale in p) sia dz-ESATTO (Ostrogradsky).
Le condizioni di non-esattezza -> sistema lineare su lambda -> 16 relazioni -> beta simbolici.
E' ben definito ma e' un'ulteriore implementazione simbolica sostanziosa (test di dz-esattezza
su funzioni razionali di p per ogni lettera). Fondazione pronta; il Solve finale resta da fare.

## 4terdecies. Matching di Ostrogradsky per beta: TENTATO, NON RIUSCITO (residuo 0.72)

Aggiunto in `vaidya_sep_fay_analytic.wl` sez.4 (marcato INCOMPLETO/BUGGATO):
- Idea: mod dz-esatto ogni coeff-lettera (razionale in p) e' caratterizzato da invarianti ai poli
  (residuo t^-1 + doppio-polo t^-2). Match delta_phi vs 5 base -> sistema per beta.
- ESITO: residuo matching = 0.72 (non ~1e-50); beta[4]=beta[5]=0. NON funziona.
- Cause identificate: (1) manca il contenuto COSTANTE del coeff-lettera; (2) polo a z=0
  (da e_5=P(z)) trattato con hack P(10^-30), sbagliato; (3) lo spazio invarianti scelto
  non copre tutto il mod-dz-esatto (serve anche la condizione di PERIODO, non solo residui).
- Il matching corretto = implementare la caratterizzazione completa mod-dz-esatto sulla curva
  ellittica (Ostrogradsky/seconda specie: residui + coeff ℘ + periodi), che e' la riduzione
  degli integrali ellittici in piena regola. Error-prone, non chiuso.

ESITO FINALE onesto sui beta SIMBOLICI: NON RAGGIUNTI. Solidi restano: rango 5 ESATTO (1e-79),
beta NUMERICI a 80 cifre, fondazione analitica verificata (e_i razionali, Pe_i->lettere).
La chiusura simbolica dei beta e' un problema aperto (riduzione integrali ellittici completa).

## 4quaterdecies. C0, Ce in FORMA CHIUSA -> forma a blocchi 100% simbolica nei coefficienti

Script `VaidyaMetric/vaidya_sep_C0Ce_closed.py`. Le uniche due costanti additive numeriche
della forma a blocchi ora sono chiuse (valutazione in un punto regolare + sottrazione parti
principali):
- eta' regolare a z=0 (r=e4):  **Ce = (e4^3-2e4^2)/(e4-r_d) + 2 e1_zd zeta(z_d)
  - 2 e2_zi P(z_inf) + 2 e1_zi zeta(z_inf)**.
- R regolare a z=z_inf E R(z_inf)=0 (decade ~1/r, verif 1.1e-28):
  **C0 = - sum_a [ b1^a zeta(z_inf-a) + b2^a P(z_inf-a) - (b3^a/2) P'(z_inf-a) ]**.
Verifica vs numerico: Ce diff 4.9e-8, C0 diff 2.0e-8 (residuo dal passo finito di P'; esatte).

=> CONCLUSIONE: la forma a blocchi
   delta phi_tau|_sep = 1/2 G~ eta + 1/2 sum (c_i d_j - c_j d_i) A[e_i,e_j]
   ha ORA TUTTI i coefficienti in FORMA CHIUSA (algebrici in m,E + valori Weierstrass ai poli):
   residui b_n^a, e-residui, C0, Ce. E' la forma RIUTILIZZABILE, in funzioni speciali, con
   coefficienti SIMBOLICI. (L'unico pezzo numerico restante e' la COMPRESSIONE opzionale a 5 dilog.)

## 5bis. SIGNIFICATO FISICO adiabatico: Jc(m) mobile + separatrix crossing (chiarimento)

- **Jc NON e' costante**: Jc=m*j(E), j adimensionale (verif Jc/m=7.026624 esatto per m=1,2,3).
  Quindi durante l'evoluzione m(v) la separatrice si SPOSTA.
- delta phi|_sep = risposta al 1o ordine a **J FISSO** (dm F|_J), valutata all'istante J=Jc(m0):
  oggetto ISTANTANEO/locale, NON "orbita che segue la separatrice".
- **Separatrix crossing**: J conservato, Jc(m) mobile -> la particella attraversa la separatrice.
  Accrescimento (m^, Jc^ oltre J) -> CATTURA; evaporazione (m v) -> FUGA. Stesso divario fisico
  dell'asimmetria accr./evap. (eq:vaidya-asymmetry). Regime non-perturbativo, oltre il 1o ordine.
- **Doppia fragilita' sulla separatrice**: (1) periodo radiale DIVERGE (orbita instabile) ->
  l'ipotesi adiabatica stessa si degrada; (2) Jc mobile -> crossing. Matematica = limite controllato
  (genus drop, ellittico); fisica = delicata. delta phi|_sep vive nell'INTORNO, non sul punto singolare.
- Paper: esteso il paragrafo "Domain of validity" (sec:adiabatic) con questo chiarimento.

## 5ter. SEPARATRICE-TRACKING (J=Jc(m)): forma chiusa, STESSA macchina

Script `VaidyaMetric/vaidya_sep_track.py`. Correzione con Jc mobile (orbita che segue la
separatrice istantanea): sorgente = dm F + (Jc/m) dJ F = N_tot/S^{3/2},
  **N_tot = N_m + (Jc/m) N_J**, N_J polinomiale (grado 7, verif), Jc/m = j(E) = 7.026624.
Tutto il resto IDENTICO alla tau: pole-adapted G_track (parti principali), residui b_n^a /
C0 / C_e dalle STESSE formule chiuse (N_m -> N_tot), dilog A[e_i,e_j] invariati (dipendono
dal clock, non da N). Coefficienti tutti algebrici in m,E.
Verifica: dG_track/dz=R_track 2.2e-8; delta phi_track assemblato (IBP, G esplicito) = diretto
int[dmF+(Jc/m)dJF] eta dr a 1e-8..4e-7 (r=11/10/9.2).
=> la correzione separatrice-following e' in forma chiusa esplicita (funz. speciali) con
coefficienti simbolici, esattamente come la J-fissa. Paper: nota in "Domain of validity".

## 5quater. RAMO v (tempo avanzato) sulla separatrice: ASSEMBLATO (verif 1e-8)

Script `VaidyaMetric/vaidya_sep_vbranch_assembly.py`. Stessa curva/R/G~ della tau; clock diverso:
  v = E U_3 + r + 2m ln(r-2m),  v_z = dv/dz = E r^3/(r-r_d) + r sqrt(Q4)/(r-2m).
  - 1o termine (pari, E r^3/(r-r_d)): riproduce i dilog della tau (poli z_d, z_inf).
  - 2o termine (DISPARI sulla curva, r sqrt(Q4)/(r-2m)): polo SEMPLICE all'immagine
    dell'orizzonte 2-torsione z=i w_im (r=2m=e3) -> DILOGARITMO DI ORIZZONTE.
  ln(r-2m) come funzione di z = funzione ellittica (zero doppio a i w_im, poli semplici +-z_inf)
  -> combinazione lnσ (parte reale verificata costante; c'e' termine lineare lambda z di quasi-periodicita').
delta phi_v|_sep = G~(z) v(r) - int_{z0}^z G~ v_z dz. ASSEMBLATO (G~,v_z espliciti) = DIRETTO
  int dm F v dr a 1e-8..6e-8 (r=11/10.5/10). Tutto esplicito in sigma,zeta,P + dilog ellittici,
  coefficienti chiusi (residui, stessa macchina + settore orizzonte).
=> ABBIAMO ENTRAMBI I RAMI (tau e v) sulla separatrice in forma chiusa esplicita con coeff simbolici.
Paper: nota aggiunta in "Domain of validity".

## 6bis. THAKURTA-KERR (a≠0) separatrice, ramo tau: STESSO PIPELINE, verificato 1e-9

Script `ThakurtaMetric/tk_sep_baseline.py`. Scoperta chiave: la correzione adiabatica TK
(parametro = fattore conforme A, via E_eff=Ehat/A) COLLASSA a
   delta phi_TK = -Ehat (A'/A) int eta * dE F dr,
struttura IDENTICA a Vaidya (int clock * sorgente). Quindi riuso il pipeline con:
- sorgente dE F = N_tau/S^{3/2}, **N_tau = E J r^4 (r-2M)^2 DE** (E-derivata, da reproduce_reductions);
- clock **eta = U_3 - 2M U_2** (stesso di Vaidya tau, dtau/dr=r^2(r-2M)/sqrtS);
- curva **S = r(r-2M)DE(rDelta - J^2 DE)**, Delta=r^2-2Mr+a^2 (a=0.9); DE=(E^2-1)r+2M.
Separatrice: radice doppia del cubico (rDelta-J^2 DE). M=1,a=0.9,E=1.2:
  **Jc=20.327866**, r_d=-7.12951, Q4 radici {-4.545,0,2,16.259}, a4=0.44, turning e4=16.259.
  (altra separatrice near-horizon a Jc=0.210, r_d=0.241.)
Verifica: r(z) 8.3e-9; dG/dz=R_tau 3.9e-9 (pole-adapted G, parti principali);
  delta phi assemblato (G esplicito) = diretto int eta dEF dr a 1e-9..4e-8 (r=19/18/17).
=> il rotante chiude con l'IDENTICA macchina genus-1. Residui b_n^a/C0/Ce dalle stesse formule
   (N_m->N_tau). Restano: block assembly weight-2 (dilog), ramo t (curva R6=r Q2 DE).

## 6ter. THAKURTA-KERR tracking (J=Jc(A)): forma chiusa, verificato 1e-8

Script `ThakurtaMetric/tk_sep_track.py`. La separatrice mobile Jc(A) (via E_eff=Ehat/A):
sorgente = dE F + (dJc/dEeff) dJ F = N_tot/S^{3/2}, N_tot = N_tau + (dJc/dE) N_J.
Pezzi in FORMA CHIUSA (tutti verificati):
- **N_tau = E J r^4 (r-2M)^2 DE**  (sorgente E-deriv);
- **N_J = S dJ K - 1/2 K dJ S = r^3 (r-2M)^2 DE^2**  (J-deriv; Delta si CANCELLA, verif N_J - r^3(r-2M)^2DE^2=0);
- **dJc/dEeff = -Ehat Jc r_d / DE(r_d)**  (teorema funzione implicita su g=rDelta-J^2 DE:
  dJc/dE=-g_E/g_J, g_E=-2E J^2 r, g_J=-2J DE). Verif: analitico=-152.960 = finite-diff.
Stessa macchina genus-1 (pole-adapted G). Verifica: dG/dz=R_track 8.9e-10;
delta phi_track assemblato = diretto int[dEF+(dJc/dE)dJF] eta dr a 3e-9..1e-7.
=> anche TK ha ENTRAMBE le versioni (J fisso, J-tracking) in forma chiusa esplicita, coeff simbolici.
Analogo esatto di Vaidya: dJc/dE (TK) al posto di Jc/m (Vaidya). a!=0 rompe lo scaling lineare.
Jc(E_eff) verificato mobile: E=1.15/1.20/1.25/1.30 -> Jc=31.34/20.33/14.53/11.05.

## 6quater. TK tau BLOCK ASSEMBLY: int G~ eta' decomposto in dilog (verif 1e-9)

Script `ThakurtaMetric/tk_sep_blockassembly.py`. Come Vaidya, decompone il weight-2 in atomi:
  delta phi_tau^TK ~ int R_tau eta dz = 1/2 G~ eta + 1/2 sum_{i<j}(c_i d_j - c_j d_i) A[e_i,e_j].
- c_i = residui di R_tau (N_tau) dalle formule chiuse (b1zd,b2zd,b3zd, b2h semiperiodi, C0);
- d_j = residui di eta'=(r^3-2M r^2)/(r-r_d) (e1_zd, e2_zi, e1_zi, Ce);
- A[e_i,e_j] = dilog ellittici antisimmetrici (Brown-Levin), 8 blocchi pari.
Verifica: R,eta' ricostruiti (4e-7,1e-9); **17 atomi A**; block = diretto a 1.3e-9..2e-8 (r=19/18/17).
=> TK tau ora COMPLETO come Vaidya: forma a blocchi esplicita, coeff simbolici (residui chiusi).
NB: log letti da file (unbuffered) per progresso live; tolto finite-diff ridondante altrove.
RESTA: ramo t (curva R6=r Q2 DE, clock 2a+3a specie, separatrice distinta).

## 6quinquies. TK tau TRACKING block assembly (verif 1e-8) -> ramo tau COMPLETO

Script `ThakurtaMetric/tk_sep_track_blockassembly.py` (N_tau -> N_tot=N_tau+(dJc/dE)N_J).
dJc/dE=-152.960 (formula chiusa, no finite-diff -> gira 40s). 13 atomi A (meno dei 17 J-fisso:
con N_tot alcuni coeff si annullano). block = diretto a 3.4e-9..5.3e-8.
Log letti LIVE da file (python3 -u > file, no grep bufferizzato). tqdm ok.
=> TK ramo tau COMPLETO: J-fisso (17 dilog) + J-tracking (13 dilog), forma a blocchi esplicita,
   coeff chiusi (residui). RESTA: ramo t (curva R6, clock 2a+3a specie, separatrice distinta).

## 5quinquies. VAIDYA v BLOCK ASSEMBLY: int G~ v_z decomposto in dilog (verif 1e-7)

Script `VaidyaMetric/vaidya_sep_v_blockassembly.py`. Prima era solo baseline IBP; ora
int G~ v_z esplicito in atomi. delta phi_v = G~(z) v - sum_{i,j} c_i c_j int gf_i gf_j dz.
- v_z = E r^3/(r-r_d) + r sqrt(Q4)/(r-2m); sqrt(Q4)(z) = dr/dz = (1/sqrt a4)[P(z-z_inf)-P(z+z_inf)]
  (meromorfa, corretta per z COMPLESSO: niente hack real-part). Parti principali di v_z via contorno:
  poli +-z_d (res), +-z_inf (ord.2), iw (orizzonte, dal termine DISPARI). Ricostruito 6.7e-8.
- G~ (13 termini principali) x v_z (8 termini) = 104 prodotti; 11 dilog (lnσ×ζ) + C/P + peso-1.
- NB: e' un integrale PRODOTTO int gf_i gf_j (non iterato); bordo con G~ RAW: delta phi_v=Gt(z)v-sum.
Verifica: block = diretto a 5e-8..1e-7 (r=11/10.5/10). Coeff = prodotti di residui chiusi.
=> Vaidya v ora COMPLETO come tau: weight-2 esplicito in dilog nominati (incl. dilog orizzonte).

## 5sexies. VAIDYA v TRACKING block assembly (verif 1e-7) -> VAIDYA COMPLETO

Script `VaidyaMetric/vaidya_sep_v_track_blockassembly.py` (N_m -> N_tot=N_m+(Jc/m)N_J).
N_J = S dJK - 1/2 K dJS = S*DE + J^2 r(r-2m) DE^3; dJc/dm=Jc/m=7.02662 (scaling lineare Vaidya).
Clock v e v_z invariati; cambia solo R via N_tot -> G~_track. block = diretto int[N_tot/S^{3/2}] v dr
a 2.7e-7..7.2e-7 (r=11/10.5/10). 11 dilog + orizzonte, coeff residui chiusi.
=> VAIDYA COMPLETO: tau (fisso+tracking) + v (fisso+tracking), tutti block-assembled, coeff simbolici.
RESTA (per l'intero programma): TK ramo t (curva R6=r Q2 DE, clock 2a+3a specie).

## 6sexies. THAKURTA-KERR ramo t (tempo coordinata): DUE separatrici, baseline verif 1e-7

Script `ThakurtaMetric/tk_t_sep_baseline.py` (prograda) + `tk_t_sep_baseline_retro.py`.
Curva R6=r Q2 DE; clock t = eta2+eta3: eta2'=P3/sqrtR6 (2a specie, b=(8E^2M^3-2(E^2-1)JMa,4E^2M^2,2E^2M,E^2,0)),
eta3'=R_Delta/(Delta sqrtR6) (3a specie, poli agli ORIZZONTI r_pm=1.436,0.564). rho_t=P3+R_Delta/Delta.
Separatrice = radice doppia di Q2 (quartico). ASIMMETRICA (frame-dragging, Q2 lineare in J): DUE separatrici:
- **prograda Jc+=+19.0894427**, r_d=-6.6207, Q4 radici {-4.545,0,1.839,11.403}, turning 11.40;
- **retrograda Jc-=-18.6710563**, r_d=-6.5884, Q4 radici {-4.545,0,2.241,10.935}, turning 10.94.
  |Jc+| != |Jc-| (19.089 vs 18.671). Entrambe r_d<0 (come Vaidya), orbite scattering a r>10.
(Jc,r_d) raffinati con mpmath.findroot su {Q2=0,dQ2/dr=0} (np.roots splittava la doppia in coppia
complessa ±0.037j). R6=(r-r_d)^2 Q4 (resto 3.5e-11 / 5.4e-11).
sorgente N_t=dE(Kt/sqrtR6)R6^{3/2}. G~_t pole-adapted. delta phi_t = G~_t eta_t - int G~_t eta_t':
  prograda = diretto a 3e-8..5e-7; retrograda a 2e-8..4e-7. check dG/dz=R_t ~ 4-5e-9.
=> ramo t chiude con lo stesso pipeline (clock 3a-specie orizzonte gestito). RESTA: block assembly
   (dilog espliciti incl. orizzonte r_pm) + tracking, per entrambe le separatrici.

## 6septies. TK ramo t BLOCK ASSEMBLY (2 separatrici): eta3 + weight-2 espliciti (verif 1e-6)

Script `ThakurtaMetric/tk_t_sep_blockassembly.py` (prograda) + `_retro.py`. Completati (a)+(b):
- (a) eta3 (clock 3a-specie orizzonte) ESPLICITO: eta_t'(dz)=rho_t/(r-r_d)=[P3+R_Delta/Delta]/(r-r_d)
  e' RAZIONALE in r (EVEN, niente sqrtQ4, a differenza di v). Poli: +-z_d, +-z_inf, +-z(r_pm).
  Immagini orizzonti z(r_pm) COMPLESSE (Im=w_im, semiperiodo): r_pm in (e2,e3) mappano li'.
  z(r_pm) trovate con mp.findroot(r_of_z(z)=r_pm). eta_t' ricostruito da parti principali 1-2e-7.
- (b) BLOCK: delta phi_t = G~_t(z) eta_t - sum c_i c_j int gf_i gf_j (prodotto), coeff residui chiusi.
  prograda: 20 dilog (incl. orizzonte z(r_pm)), block=diretto 2e-6..7e-6;
  retrograda: 22 dilog, block=diretto 2e-6..6.5e-6.
  (prec 1e-6, limitata da eps contorno vicino ai poli orizzonte complessi.)
=> TK ramo t COMPLETO (fisso, prograda+retrograda): weight-2 esplicito in dilog, eta3 in lnsigma
   agli orizzonti. RESTA solo il tracking del ramo t (N_t->N_tot con dJc_pm/dE, meccanico).

## 7. CROSS-CHECK MATHEMATICA: identita' esatte + BUG TROVATO E CORRETTO nel ramo t

Script `paper/crosscheck_identities.wl`, `crosscheck_tk_numeric.wl`, `crosscheck_vt_numeric.wl`,
`crosscheck_tkt_flat.wl` (wolframscript, WeierstrassP nativo, indipendente da theta1/mpmath).
- IDENTITA' SIMBOLICHE (tutti i casi) confermate ESATTE: Jc (Vaidya Root=radicale, TK tau, TK t
  prograda +19.089/retrograda -18.671), N_m/N_tau/N_J (=0), N_J(TK)=r^3(r-2M)^2 DE^2, dJc/dE=-EJr/DE,
  scaling Vaidya, asimmetria ramo t (Q2 lineare in J).
- NUMERICO: Vaidya tau (residui, rango 5); TK tau (Jc,r_d,Q4,r(z) 1e-29, delta phi diretto 8 cifre);
  Vaidya v (sqrtQ4=dr/dz 1e-28, delta phi_v diretto); z(r_pm) orizzonti 1e-31.
- **BUG TROVATO (ramo t)**: il Python estraeva P3 con sp.div(numer, Delta) ma il denom di cancel(rho_t)
  e' **2500*Delta** -> P3 era 2500x, P3+RD/Delta != rho_t, e delta phi_t 100x troppo grande.
  FIX: dividere per il denominatore VERO (sp.denom(rho_t)); RD=(rho_t-P3)Delta. Ora P3 leading=E^2=1.44
  (= P3 del paper, b-vector), P3+RD/Delta-rho_t=0. Corretti tk_t_sep_baseline.py/_retro.py +
  blockassembly.py/_retro.py. VALORI CORRETTI: prograda delta phi_t(19)=0.0871960819 (era 8.72),
  retrograda -0.0795252758. block=diretto=Mathematica, 16 dilog (erano 20/22). Verif 2e-8..2e-10.
=> Il cross-check indipendente ha CATTURATO un errore reale (l'assemblato=diretto Python era
   self-consistente col clock buggato, non l'aveva visto). Ora Python==Mathematica.

## 8. TK ramo t TRACKING (separatrice mobile Jc_pm(E)) -- COMPLETO, entrambe le separatrici

Script `ThakurtaMetric/tk_t_sep_track_blockassembly.py`. Segue la separatrice t mentre A varia
(E_eff=Ehat/A): J=Jc_pm(E). Solo la SORGENTE cambia, il clock eta_t (2a+3a specie) e' invariato.
- N_tot = N_t + (dJc_pm/dE) N_J^t,  con N_t=dE(K_t/sqrtR6)R6^{3/2}, N_J^t=dJ(K_t/sqrtR6)R6^{3/2}.
- dJc_pm/dE = -Q2_E/Q2_J valutato a (r_d,Jc): dalla doppia radice di Q2 (Q2_r=0 =>
  Q2_J dJc + Q2_E dE = 0). CONFERMATO con finite-diff (h=1e-6):
    prograda:   dJc/dE = -148.541751  (Jc=+19.0894, r_d=-6.6207)
    retrograda: dJc/dE = +150.328105  (Jc=-18.6711, r_d=-6.5884)
  |dJc/dE| quasi uguale ma NON identico (asimmetria frame-dragging, come |Jc+|!=|Jc-|).
- block assembly (G~ eta_t - sum c_i c_j int gf_i gf_j) vs diretto (N_tot):
    prograda   delta phi_t^track(19)=-4.03767e-01, (18)=-1.875e+00, (17)=-5.003e+00 ; diff ~1e-7
    retrograda delta phi_t^track(19)=+3.889e-01,   (18)=+1.792e+00, (17)=+4.734e+00 ; diff ~2-7e-7
  16 dilog (prograda) / 18 dilog (retrograda), 143 prodotti; coeff residui chiusi su N_tot.
=> PARITA' COMPLETA su tutti i casi: Vaidya (tau,v: fisso+tracking) e TK (tau,t: fisso+tracking,
   entrambe le separatrici). Nessun pezzo meccanico rimasto. (Aperto solo: beta simbolico dei 5-dilog
   Fay, frontiera di ricerca, non richiesto.)

## 9. J GENERICO (off-separatrice, GENUS-2): chiusura W_ij in θ[δ] — TK-τ

Dettaglio completo in `progress.md` §10 + `ThakurtaMetric`/`KerrMetric` script. Off-separatrice
la sestica S resta genus-2 (no doppia radice). Correzione adiabatica ψ = ½Ê Σ Q_ij W_ij.

**STRUTTURA VERIFICATA (mattoni 1-2b, 1e-15):**
- U_0..U_4 = ∫r^k/√S indipendenti; U_2 unico generatore 3ª specie (residuo ∞), U_3,U_4 2ª specie.
- ψ = ½Ê[ 3 iterati Kleinian peso-1 + 1 dilog genus-2 (coppie con L) ] + T_alg(algebrico).
  Clock ha 2 letters (R1 2ª, L 3ª); source senza 3ª (a_L=0). (u1,u2)=0 ⟹ Q_01=0.

**NAMING (forma-differenza θ[δ] agli e_±, aggira muro divisore-θ):**
- Ω=log[θ[δ](w-e+)/θ[δ](w-e-)] (3ª specie); Z=ζ_δ(w-e±), P=℘_δ(w-e±) (2ª specie). δ dispari #1.
- Verifica integrale: U_2,U_3,clock η a 1e-6 (`kerr_tau_Wij_diffform_integral.sage`).
- COEFFICIENTI SIMBOLICI (parti principali, base canonica): U_k coeff=-g_{k-2},
  g0=1/√(E²-1), g1=(2E²-3)/(E²-1)^{3/2}, g2=(625E⁶-1156E⁴+37E²+794)/(200√(E²-1)(E²-1)²).
  Cross-check vs residui BEL mattone-2a = 0 esatto. Clock: η=-g0 Z+(2M g0-g1)Ω+olo.
  ⟹ φ0 e clock CHIUSI in θ[δ] con coeff simbolici (come σ,ζ,℘ separatrice).

**MURI ONESTI:** ζ_i nudo cade sul divisore-θ (§7-iii), aggirato con forma-differenza agli e_±.
**APERTO (frontiera):** q-serie nominata del singolo dilog genus-2 (Kronecker-Eisenstein/Fay),
analogo genus-2 di Γ̃ Brown-Levin. = prossimo lavoro.
