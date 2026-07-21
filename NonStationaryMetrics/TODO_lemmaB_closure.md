# TODO — Chiusura del Lemma B (monotonia dell'angolo, no-inversion a estremi fissi)

## Contesto

Il **teorema no-inversion a estremi fissi** (referee #11): a endpoint fissi, la
brachistocrona del ramo `t` resta piu' superficiale di quella `tau`,
`r_min^t > r_min^tau`. Riduzione esatta (parametrizzando per raggio di turning
`r_min`, che fissa `J`; `Phi_br(r_min) = int_{r_min}^{r0} F_br(r;J_br(r_min)) dr`):

- **Lemma A** — `F_t(r;J_t) > F_tau(r;J_tau)` a `r_min` uguale, per ogni `r>r_min`
  ⟹ `Phi_t(r_min) > Phi_tau(r_min)`.  **DIMOSTRATO in forma chiusa** (algebrico).
- **Lemma B** — `Phi_tau(r_min)` strettamente decrescente ⟹ (con A) la tesi.
  **APERTO nel regime fisico.**

Vedi `no_inversion_reduction.py` (funzioni `symbolic_full_proof`,
`lemma_B_partial`, `lemma_B_regime`).

## Stato di Lemma B (`dPhi_tau/dr_min < 0`)

Forma potenziale efficace: `Phi_tau = sqrt(Vmin) int K/sqrt(V - Vmin) dr`,
`V = r Delta / DE`, `K = sqrt(r(r-2M))/Delta`, `W := K/V'`. Con `V = Vmin + sigma^2`
+ una IBP:

```
dPhi/dVmin = (V0-2Vmin) W(Vmin)/(Sigma sqrt(Vmin))
           + (2/sqrt(Vmin)) int_0^Sigma W'(Vmin+sigma^2) g(sigma) dsigma
   g(sigma) = (Sigma - sigma)(sigma + Vmin/Sigma) > 0 su [0,Sigma)   (radici Sigma, -Vmin/Sigma)
   W>0, W'<0, W''>0 (W decrescente E convessa: verificato, r>=2.35 < r_pk)
```

Il secondo termine e' **manifestamente < 0**. Il primo e' `<= 0` sse `V(r0) <= 2 V(r_min)`.

- **PROVATO in forma chiusa**: regime "grazing" `V(r0) <= 2 V(r_min)` (r_min vicino a r0).
- **Esteso (Chebyshev integrale + convessita' di W)**: fino a circa `V(r0) < 4 V(r_min)`.
- **APERTO**: i turning point fisici hanno `V0 ~ 8 Vmin` (es. r_min=4.7, r0=10), fuori
  dai regimi sopra. Il margine e' stretto (|integrale| ~ 2x boundary), quindi servono
  **stime sharp**, non bound grezzi. `Phi_tau` ha un massimo a `r_pk ~ 3M`
  (soglia TRASCENDENTE): nessun certificato polinomiale elementare puo' raggiungerla.

## Da fare — chiusura via BRACHISTOCRONE IN FORMA CHIUSA (periodi)

L'idea: `Phi_tau(r_min)` **e'** un periodo della curva spettrale (le forme chiuse
Weierstrass/theta della brachistocrona sono esattamente `Phi`). Chiudere Lemma B in
forma NON elementare ma rigorosa:

1. **Via periodi / Picard-Fuchs** (self-contained):
   - Esprimere `Phi_tau(r_min)` in forma chiusa: separatrice = Weierstrass `sigma,zeta,wp`
     (genus 1); orbite generiche (scattering, i turning point fisici) = **theta genus-2**.
   - `dPhi_tau/dr_min` = derivata del periodo rispetto al modulo ⟹ soddisfa una
     **ODE lineare di Picard-Fuchs**. Dedurre il segno da struttura ODE + dati asintotici
     (valore e segno a un punto noto, es. r_min->r0 dove Phi->0).
   - NB: i turning fisici sono GENUS-2 ⟹ la monotonia e' proprieta' theta genus-2
     (non-elementare, non piu' facile della stima diretta, ma e' l'ambiente naturale).

2. **Via stime sharp** (analitiche, non numeriche): sfruttare
   - convessita' di W (gia' provata, `W''>0`), PIU'
   - il **decadimento esatto** `W(V) ~ 1/(2V)`, `W'(V) ~ -1/(2V^2)` a V grande
     (forma specifica, non solo convessita'), per controllare la coda `Q` dell'integrale
     nel regime `V0 >> Vmin`. Combinare i regimi (Chebyshev per V0 moderato, decadimento
     per V0 grande) in un'unica disuguaglianza `|integrale| > boundary` per ogni r_min>r_pk.

3. **Via teorema pubblicato** (adattamento): monotonia dell'angolo di deflessione sopra
   la sfera fotonica (Tsukamoto 2016 PRD 94 124001; Bozza). CAVEAT: i risultati sono per
   geodetiche NULL da infinito; il nostro e' brachistocrona MASSIVA, metrica ottica
   Randers/Finsler del ramo tau, r0 FINITO. Serve adattamento, non pura citazione.

## Cosa e' gia' pronto (riusabile)

- `no_inversion_reduction.py`: riduzione + Lemma A (prova completa) + Lemma B (parziale + verifica).
- Convessita' di W provata (`/tmp/Wconvex.py`, da spostare nel repo se serve).
- Forme chiuse delle brachistocrone: separatrice Weierstrass e generiche genus-2 gia'
  nel paper e negli script (`ThakurtaMetric/`, `VaidyaMetric/`, appendici del main).

## Priorita' consigliata

Via 2 (stime sharp) o Via 1 (Picard-Fuchs) — entrambe rigorose non-elementari, chiudono
tutto `r_min > r_pk`. Via 3 la piu' rapida se l'adattamento del teorema regge.
Effort: research-level ma delimitato (NON frontiera come la dimensione del dilog #14).
