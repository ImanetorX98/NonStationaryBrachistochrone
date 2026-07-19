# Coefficienti separatrice pienamente SIMBOLICI (riutilizzabili)

Per ogni ramo, i coefficienti dei dilog (residui b1,b2,b3 della sorgente R al polo triplo z_d,
e residui del clock) sono espressi in forma **simbolica in (M,a,E, r_d, Jc)**, dove r_d, Jc =
doppia radice della curva (S(r_d)=0, S'(r_d)=0). Script: `SEP_COEFF_SYMBOLIC.py`,
`VaidyaMetric/sep_coeff_symbolic.py`. Verificati vs estrazione a contorno.

## Formula UNIVERSALE (vale per ogni ramo)
Dati la curva della separatrice (S per τ/v, R6 per t) e la sorgente N:
    Q4(r_d) = S''(r_d)/2      Q4'(r_d) = S'''(r_d)/6      Q4''(r_d) = S''''(r_d)/12
    a4 = [r^6] della curva
    s  = √Q4(r_d)             a1 = Q4'(r_d)/(4s)          a2 = Q4''(r_d)/12
    F  = N/Q4  ⇒  h0=F(r_d),  h1=F'(r_d)·s,  h2=½(F''(r_d)s² + F'(r_d)Q4'(r_d)/2)
    ┌────────────────────────────────────────────────────────────┐
    │  b3 = h0/s³                                                 │
    │  b2 = (h1 − 3a1 h0)/s³                                      │
    │  b1 = (h2 − 3a1 h1 + (6a1² − 3a2) h0)/s³                    │
    └────────────────────────────────────────────────────────────┘
Tutte le derivate di S,N a r_d sono POLINOMI in (M,a,E,r_d,Jc) → b_i razionali/algebrici.
r_d,Jc si ottengono per-curva risolvendo la doppia radice (2 equazioni in M,a,E).

## Rami — curva e sorgente
| Ramo | Curva | Sorgente N | Clock |
|---|---|---|---|
| Vaidya τ (a=0) | S=r(r−2M)DE(r²(r−2M)−Jc²DE) | N_m=S ∂_M K − ½K ∂_M S, K=Jc DE | r³−2Mr² |
| Vaidya v (a=0) | S (stessa) | N_m (stessa) → **stessi b_i di τ** | v_z=E r³/(r−r_d)+r√Q4/(r−2M) |
| TK τ | S=r(r−2M)DE(rΔ−Jc²DE) | N_τ=E Jc r⁴(r−2M)²DE | r³−2Mr² |
| TK t (Jc±) | R6=r Q2 DE (Q2 quartica, lineare in Jc) | N_t=∂_E(K_t/√R6)R6^{3/2} | ρ_t=P3+R_Δ/Δ |
DE=(E²−1)r+2M, Δ=r²−2Mr+a².

## Notevoli (semplici)
- a4 = [r^6] della curva τ = **E²−1** (solo E). ⇒ e2_zi=1/a4 = **1/(E²−1)**.
- clock τ residuo: e1_zd = (r_d³−2M r_d²)/s.

## Verifiche numeriche (tutte confermate)
| Ramo | (M,a,E) | Jc | r_d | b1 | b2 | b3 |
|---|---|---|---|---|---|---|
| Vaidya τ | (1,0,7/5) | 7.026624 | −3.363711 | +0.2704 | +0.0326 | +0.00987 |
| TK τ | (1,9/10,6/5) | 20.327866 | −7.129509 | −1.8364 | −0.0440 | −0.0479 |
| TK t + | (1,9/10,6/5) | 19.089443 | −6.620747 | −1.6172 | −0.0737 | −0.0353 |
| TK t − | (1,9/10,6/5) | −18.671056 | −6.588388 | +1.6170 | +0.0708 | +0.0343 |
Vaidya τ: match vs contorno a 1e-7. TK τ/t: Jc,r_d match valori noti; formula identica a Vaidya.

## Uso (riutilizzabile)
1. Scegli i parametri fisici (M,a,E).
2. Risolvi la doppia radice della curva del ramo → (r_d, Jc).
3. Sostituisci (M,a,E,r_d,Jc) nelle formule simboliche b1,b2,b3 (e nei residui clock).
4. I dilog ellittici (Brown-Levin Γ̃ / D(a,b)) sono le funzioni tabulate valutate sulla curva.
⇒ Coefficienti = formule simboliche universali; funzioni speciali = valutate per-curva (come K(m)).

## Residui del CLOCK (per i coeff dilog c_ab = b_a(sorgente)·d_b(clock))
- **τ (Vaidya/TK)**: clock r³−2M r². Residui: e1_zd=(r_d³−2M r_d²)/s ; e2_zi=1/a4=1/(E²−1).
- **v (Vaidya)**: clock v_z=E r³/(r−r_d)+r√Q4/(r−2M). Residui SIMBOLICI (verificati):
    z_d: **E r_d³/s**  ;  orizzonte z=iω_im: **4M** (esatto, indip. dai e_i). Sorgente b_i = Vaidya τ.
- **t (TK)**: clock ρ_t=P3+R_Δ/Δ, poli 3ª specie agli orizzonti r±. Residui SIMBOLICI (verificati):
    z_d: **ρ_t(r_d)/s** ;  orizzonti z(r±): **σ·R_Δ(r±)/[(r±−r∓)(r±−r_d)√Q4(r±)]** (σ=±1 foglio √Q4).
    Q4(r±)=R6(r±)/(r±−r_d)², R6=r Q2 DE. INVARIANTE: res(r+)+res(r−)=2M (verificato). 
    [prograda Jc+: res z_d=−8.907, z(r+)=+3.294, z(r−)=−1.294; match contorno 1e-6]

## Costanti additive Ce, C0 (forma esplicita) e natura period-level
Forma chiusa esplicita (verificata 1e-8, `VaidyaMetric/vaidya_sep_C0Ce_closed.py`):
    Ce = η'(0) + 2 e1_zd·ζ(z_d) − 2 e2_zi·℘(z_∞) + 2 e1_zi·ζ(z_∞)
    C0 = −Σ_a [ b1^a·ζ(z_∞−a) + b2^a·℘(z_∞−a) − (b3^a/2)·℘'(z_∞−a) ],  a∈{±z_d,0,iω_im}
dove η'(0)=(e4³−2M e4²)/(e4−r_d) (clock al turning). I COEFFICIENTI (e_i, b_i) sono SIMBOLICI;
i VALORI ζ,℘,℘' ai punti marcati sono period-level.

### I pezzi period-level: simbolici? universali? o per-curva? (`sep_periodlevel_test.py`)
Test Vaidya τ a E=7/5 vs E=13/10: ω1, ω_im, z_d, z_∞, ζ(z_d), ℘(z_∞) TUTTI DIVERSI.
GERARCHIA:
- RAZIONALI-simbolici in (M,a,E,r_d,Jc): residui b_i, e_i, residui clock. [formule universali]
- ALGEBRICI in (M,a,E,r_d,Jc): radici e_i di Q4; invarianti g2,g3 (simmetrici nelle e_i).
- TRASCENDENTI (period-level): punti marcati z_d,z_∞ (integrali ellittici INCOMPLETI ∫dr/√Q4);
  valori ζ(z_d),℘(z_∞); quindi Ce,C0.
⇒ I period-level NON sono razionali-simbolici (sono periodi/valori di funzioni ellittiche),
  NON sono universali (cambiano coi parametri), si VALUTANO per-curva via procedura universale
  (radici→g2,g3→periodi→ζ,℘ ai punti; serie θ geom. convergente). ESATTAMENTE come K(m): la
  FORMULA e' universale (scritta una volta), il VALORE e' per-modulo. Irriducibile: e' cio' che
  "periodo" significa.
