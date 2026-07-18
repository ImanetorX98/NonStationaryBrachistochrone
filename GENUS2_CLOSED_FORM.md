# Forma chiusa — brachistocrona adiabatica genus-2 (Thakurta–Kerr τ, J generico)

Versione **pulita**: coefficienti SIMBOLICI universali × funzioni speciali TABULATE
(valutate per-curva, come K(m) nel pendolo). Nessun coefficiente period-level superfluo.

## Espansione adiabatica
    φ(r,A) = φ₀(r; E_eff) + (Ȧ/A)·δφ(r) + O((Ȧ/A)²),   E_eff = Ê/A

## Curva (genus 2)
    S(r) = r(r−2M)·DE·(rΔ − J²·DE),   DE = (E²−1)r+2M,   Δ = r²−2Mr+a²
    (6 radici semplici per J generico ⇒ genere 2; E ≡ E_eff)

═══════════════════════════════════════════════════════════════════
## FUNZIONI SPECIALI TABULATE (valutate per-curva)
═══════════════════════════════════════════════════════════════════
    U_k(r)  = ∫_{r₀}^r  r'^k dr'/√S            integrali abeliani genus-2 (1ª/2ª/3ª specie)
    Π±(r)   = ∫_{r₀}^r  dr'/[(r'−r±)√S]        3ª specie agli orizzonti r± = M±√(M²−a²)
    W_kj(r) = ∫_{r₀}^r  (U_k dU_j − U_j dU_k)  polilog genus-2 peso-2 (iterato antisimmetrico)

═══════════════════════════════════════════════════════════════════
## ORDINE 0 — orbita congelata
═══════════════════════════════════════════════════════════════════
    ┌──────────────────────────────────────────────┐
    │  φ₀ = C₁·U₁ + C₀·U₀ + α₊·Π₊ + α₋·Π₋           │
    └──────────────────────────────────────────────┘
Coefficienti SIMBOLICI (esatti):
    C₁ = J(E²−1)
    C₀ = 2MJ
    α± = ∓ J a²[ M(E²+1) ± (E²−1)√(M²−a²) ] / [ 2√(M²−a²) ]

═══════════════════════════════════════════════════════════════════
## ORDINE 1 — correzione adiabatica
═══════════════════════════════════════════════════════════════════
    ┌────────────────────────────────────────────────────────┐
    │  δφ = ½ Ê · Σ_{k<j} Q_kj · W_kj  +  ½ Ê · G_alg         │
    └────────────────────────────────────────────────────────┘
Coefficienti SIMBOLICI:
    Q_kj = c_k b_j − c_j b_k          (antisimmetrico)
    b    = (0, 0, −2M, 1, 0)          (vettore clock, dτ/dr = (r³−2Mr²)/√S)
    c_k  = riduzione 2ª specie di ∂_E F: unica soluzione razionale-in-E di
           2N = 2S·𝒜₅' − 𝒜₅·S' + 2S·Σ_k c_k r^k,   N = E J r⁴(r−2M)² DE,  𝒜₅ = poly deg 5

I 7 Q_kj non nulli (in forma simbolica esplicita):
    Q₀₂ = −2M c₀     Q₀₃ = c₀      Q₁₂ = −2M c₁     Q₁₃ = c₁
    Q₂₃ = c₂ + 2M c₃     Q₂₄ = 2M c₄     Q₃₄ = −c₄
    (Q₀₁ = Q₀₄ = Q₁₄ = 0)

G_alg = 2·I_el(r) + [η·(cost − 𝒜₅/√S)] boundary   (ELEMENTARE, coeff SIMBOLICI in E)
        I_el(r) = P(r) + Σ_{x:C(x)=0} res(x)·log(r − x)
          P(r) = [p₃r³+p₂r²+p₁r]/D   coeff RAZIONALI in (M,a,E,J) — NON solo E
          C(r) = rΔ−J²DE = r³−2Mr²+(a²−J²(E²−1))r−2MJ²   (cubico; le sue 3 radici)
          res(x) = 𝒜₅(x)·x²(x−2M)/S'(x)   RAZIONALE in (x,M,a,E,J) (RootSum sul cubico)
          [res(x)=0 provato simbolico a x=0,2M,2M/(1−E²) → solo 3 radici del cubico contribuiscono]
        boundary: η=U₃−2M U₂ (tabulato) × 𝒜₅/√S (𝒜₅ razionale in E).
        VERIFICATO: forma chiusa vs diretto = 0; coeff simbolici cross-check Mathematica
        (crosscheck_Ialg_symbolic.wl); non universali (E=7/5 != E=13/10).

═══════════════════════════════════════════════════════════════════
## STRUTTURA FINE del peso-2 (opzionale, espandendo W_kj in θ)
═══════════════════════════════════════════════════════════════════
Σ Q_kj W_kj  =  [3 iterati Kleiniani peso-1]  +  h_L·Λ  (dilog genus-2)
    h_L = 1/(E²−1)
I 3 iterati Kleiniani chiudono in θ[δ], ζ_δ, σ (funzioni Kleiniane, peso 1).
Λ = polilog ellittico genus-2 (Enriquez): serie di nome Kronecker-Eisenstein,
    geometricamente convergente (N=4 → 10⁻¹³); NON riducibile a Li₂ (genus-2 senza
    formula prodotto di Jacobi). Espansione: Λ = ∫Q dA (elem, ~80%) + ∫L dA (dilog),
    L = log[θ[δ](u−e₊;τ)/θ[δ](u−e₋;τ)] (via S: τ'=−τ⁻¹, char-zero).

═══════════════════════════════════════════════════════════════════
## 📖 LEGENDA — funzioni speciali
═══════════════════════════════════════════════════════════════════
    U_k       integrale abeliano iperellittico genus-2 (k=0,1: 1ª specie;
              k=2: 3ª specie (residuo a r=∞); k=3,4: 2ª specie)
    Π±        integrale abeliano di 3ª specie (poli agli orizzonti r±)
    W_kj      polilogaritmo ellittico genus-2 peso-2 (iterato di 2 abeliani)
    θ[δ](z;τ) theta di Riemann genus-2, caratteristica δ dispari, matrice periodi τ (2×2)
    ζ_δ = ∇log θ[δ]   zeta Kleiniana genus-2 (2ª specie)
    Λ         polilog ellittico genus-2 (Kronecker-Eisenstein, classe Enriquez)
    K(m)      [analogia] integrale ellittico completo — valutato per-caso, universale

═══════════════════════════════════════════════════════════════════
## STATUS ONESTO di ogni pezzo
═══════════════════════════════════════════════════════════════════
🟢 SIMBOLICI UNIVERSALI (formule razionali in M,a,E,J — scrivi una volta):
   C₁, C₀, α±, b, c_k, Q_kj, h_L, g_i, P(r), res(x) — razionali in (M,a,E,J).
   ECCEZIONE: g₀=1/√(E²−1) e g₁=(2E²−3)/(E²−1)^{3/2} dipendono SOLO da E (leading r→∞);
   g₂, c_k, Q_kj, P, res dipendono da tutti (M,a,E,J).
🟡 FUNZIONI TABULATE (valutate per-curva via periodi/θ, come K(m)):
   U_k, Π±, W_kj, θ[δ], ζ_δ, Λ, e i dati di curva τ, e±, δ.

⇒ Universale ciò che deve esserlo (COEFFICIENTI + struttura); per-curva ciò che è
  intrinsecamente speciale (le FUNZIONI). I coefficienti period-level α_k,β_k della
  versione precedente erano una SOVRA-decomposizione di U_k in θ e NON servono: U_k è
  già una funzione tabulata. Come nel pendolo T=4√(L/g)·K(sin²(θ₀/2)): coeff simbolico
  × funzione tabulata valutata per-caso.

═══════════════════════════════════════════════════════════════════
## GERARCHIA period-level (parallela alla separatrice, vedi SEP_COEFF_SYMBOLIC.md)
═══════════════════════════════════════════════════════════════════
| Livello | Oggetti (genus-2) | Natura |
|---|---|---|
| 🟢 RAZIONALE-simbolico in (M,a,E,J) | c_k, Q_kj, g_i, P(r), res(x), h_L, C₁,C₀,α± | formule universali |
| 🟡 ALGEBRICO in (M,a,E,J) | radici di S; radici del cubico C=rΔ−J²DE (per res(x)) | radicali |
| 🔴 TRASCENDENTE (period-level) | matrice τ (2×2); e± (immagini Abel di r=∞); α_k,β_k (a-periodi ∮ω_k·ω⁻¹); i coeff interni q'^{n^T τ' n} della serie del dilog Λ | periodi |

RISPOSTA (identica alla separatrice): i period-level NON sono razionali-simbolici (sono
periodi/valori di θ genus-2), NON universali (cambiano coi parametri — verificato α_k,β_k
diversi J=5/2 vs J=2), si VALUTANO per-curva via procedura universale (S → periodi τ →
Abel e± → θ,ζ_δ). Come K(m): formula universale, valore per-modulo.
CORRISPONDENZA separatrice ↔ genus-2:
  ζ,℘ ai punti (Ce,C0)  ↔  θ[δ],ζ_δ agli e± (U_k, Λ)
  z_d,z_∞ (∫dr/√Q4)     ↔  e± (Abel map ∫du)
  b_i,e_i (residui)      ↔  g_i, c_k, Q_kj (residui/riduzioni)
  α,β genus-2 (a-periodi) sono l'analogo di C0,Ce (additivi period-level).

## Verifiche (in progress.md, tutte committate)
- struttura ψ = Σ Q_kj W_kj + G_alg: end-to-end 1e-15
- Q_kj, c_k, g_i simbolici: cross-check esatto
- W_kj / U_k naming in θ[δ]: integrale 1e-6
- dilog Λ serie di nome: convergenza geometrica N=4→1e-13
- α_k,β_k NON universali (verificato: cambiano J=5/2 vs J=2) ⇒ scartati come superflui
