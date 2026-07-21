# Lemma B — chiusura della monotonia dell'angolo (no-inversion a estremi fissi)

Record completo di tutte le considerazioni (per compattare il contesto). Il
**teorema no-inversion a estremi fissi** (referee #11): a endpoint fissi
`(r0, +/-Phi)`, la brachistocrona del ramo `t` resta piu' superficiale di quella
`tau`: `r_min^t > r_min^tau`.

## Riduzione esatta (Lemma A + Lemma B => tesi)

Parametrizza ogni ramo per il raggio di turning `r_min` (che fissa `J` via la
condizione di turning: `r Delta = J^2 DE` per tau, `Q2=0` per t). Sia
`Phi_br(r_min) = int_{r_min}^{r0} F_br(r; J_br(r_min)) dr` il semi-angolo spazzato.
- **Lemma A**: `F_t(r;J_t) > F_tau(r;J_tau)` a `r_min` uguale, per `r in (r_min,r0)`
  => `Phi_t(r_min) > Phi_tau(r_min)`.
- **Lemma B**: `Phi_tau(r_min)` strettamente decrescente.
Da A+B: `Phi_tau(r_min^tau) = Phi_t(r_min^t) > Phi_tau(r_min^t)` => (B decrescente)
`r_min^tau < r_min^t`. QED.
NB: a **stesso J** l'ordine e' OPPOSTO (t piu' profondo); il fixed-endpoint e' un
ribaltamento via la mappa `J -> Phi`. Il bound puntuale `n_t/n_tau=E/f>1` NON basta.

Script: `no_inversion_reduction.py`.

## Lemma A — DIMOSTRATO in forma chiusa (algebrico)

`sign(F_t^2-F_tau^2) = sign(N_G)`,
`N_G = (J_t(r-2M)+2Ma)^2 P_tau - J_tau^2(r-2M) Q2`, `P_tau=rDelta-J_tau^2 DE`.
Ridotto mod turning (`J_tau^2=r_min Delta_min/DE_min`, `Q2(r_min,J_t)=0`):
`N_G = (r-r_min) [4r Delta/(r_min-2M)] W_G(r)` con `W_G` LINEARE in r.
Prefattore >0 (r>r_min>2M) => (A) <=> `W_G>0`. `W_G` lineare, `W_G(r_min)>0` e
pendenza `w1>0` perche' la radice prograde `J_t` di `Q2(r_min,.)` sta sotto due soglie:
```
A = coeff J_t^2 di Q2(r_min) = -(r_min-2M)DE_min < 0     (parabola giu', vertice a J<0)
J* = [E^2 r_min(a^2+r_min^2)+2a^2]/(a DE_min) > 0
J**= [E^2 r_min^2(r_min-1)+E^2 a^2 r_min + a^2]/(a DE_min) > 0
Q2(r_min,J*)  = -E^2 r_min^3(E^2 r_min^2+a^2) Delta_min/[a^2 DE_min]                       < 0
Q2(r_min,J**) = -(r_min-2M)(E^2 r_min^2+a^2)(E^2 r_min^2[(r_min-1)^2+a^2]+a^2)/[a^2 DE_min] < 0
=> J* , J** oltre la radice grande => J_t<J*, J_t<J** => W_G(r_min)>0, w1>0 => W_G>0.
```
Valido per ogni `r_min>2M, E>1, a>0`. Funzione `symbolic_full_proof()` in
`no_inversion_reduction.py` (assert simbolici passano). **QED Lemma A.**

## Lemma B — APERTO. Struttura, forme, ingredienti provati

Forma potenziale efficace:
`Phi_tau = sqrt(Vmin) int_{r_min}^{r0} K/sqrt(V-Vmin) dr`,
`V=r Delta/DE` (potenziale, turning `V(r_min)=J^2=Vmin`), `K=sqrt(r(r-2M))/Delta`,
`W := K/V'`. Sostituzione `V=Vmin+sigma^2` + una IBP =>
```
sqrt(Vmin) dPhi/dVmin = (V0-2Vmin) W(Vmin)/Sigma + 2 int_0^Sigma W'(Vmin+sigma^2) g(sigma) dsigma
   g(sigma) = (Sigma-sigma)(sigma + Vmin/Sigma) > 0 su [0,Sigma)   (radici Sigma, -Vmin/Sigma)
   V0=V(r0),  Sigma=sqrt(V0-Vmin)
```
INGREDIENTI PROVATI:
- `W>0`, `W'<0` (verificato), e **`W''>0` (W CONVESSA)** su tutto lo scattering
  (`N_W''` cambia segno a r~2.35 < r_pk). `V'>0` ovunque (r>r_+).
- Il secondo termine (integrale) e' **manifestamente <0** (g>0, W'<0).
- Se `V(r0) <= 2 V(r_min)` (grazing, r_min vicino a r0): primo termine <=0 =>
  `dPhi/dVmin<0` **PROVATO in forma chiusa**. Esteso a ~`V0<4Vmin` con convessita'
  (Chebyshev integrale), ma non oltre in modo pulito.
- `Phi_tau` ha un MASSIMO a `r_pk~3M` (SOGLIA TRASCENDENTE, radice di dPhi=0);
  decrescente per `r_min>r_pk`. I turning point fisici (r_min~4.7-6.4, r0=10) stanno
  sopra r_pk MA a `V0~8 Vmin` >> 4Vmin: fuori dai regimi chiusi.

## SCORCIATOIE ESCLUSE (esplorate e morte)

- **Route 3 (lensing / monotonia deflessione sopra sfera fotonica): MORTA.**
  Il ramo tau NON ha sfera fotonica: `N_V = num(V') = 2((E^2-1)r^3-(E^2-4)r^2-4r+a^2)`
  non ha radici `r>r_+` per nessun (a,E) => `V'>0` ovunque, `V` monotona, nessuna
  orbita circolare instabile => `Phi_tau` NON diverge, non ci sono ipotesi di lensing.
  (La "separatrice tau" a r_d=2.0 di #7/#8 e' il fattore (r-2M) su un turning, non
  un'orbita circolare: `P_tau'(2)!=0`.)
- **Route 1 (Picard-Fuchs / periodi): MORTA per il caso fisico.**
  Per `r0->inf`, `Phi_inf` e' CRESCENTE in r_min (Phi'/Phi>0, ->pi/2), direzione
  OPPOSTA al caso r0 finito (decrescente). E il caso fisico e' `Phi_finito =
  Phi_inf - Phi_coda`, una DIFFERENZA, non un ciclo/periodo chiuso => PF (che vive
  sui periodi r0=inf) da' l'oggetto sbagliato. La curva `y^2=S_tau` e' genus 2
  (sestica irriducibile) anche per a=0, quindi PF sarebbe 4o ordine comunque.
  L'esplorazione numerica dava un'ODE di ordine basso per `Phi_inf` (possibile
  elementare), ma per il r0=inf (direzione sbagliata).

## STRUTTURA CHIAVE del caso fisico (perche' e' hard)

Il no-inversion fisico e' un effetto a **r0 FINITO**. Il risultato a r0->inf ha
direzione OPPOSTA (t piu' profondo). Il termine di coda (da r0 a inf) **ribalta il
segno** di `dPhi/dVmin`: e' una cancellazione delicata. In Form 2 a V0 grande sia il
boundary sia l'integrale sono ~`sqrt(V0) W(Vmin)` (stesso ordine): il segno viene dal
next-order. Questo e' il muro di Route 2.

## ROUTE 2 (stima analitica sharp) — da fare

Target: provare `2 int_0^Sigma (-W') g dsigma > (V0-2Vmin) W(Vmin)/Sigma` per
`r_min>r_pk`, sfruttando:
- convessita' di W (provata) => `-W'` decrescente (concentra il peso vicino al turning);
- decadimento esatto `W(V)~1/(2V)`, `W'(V)~-1/(2V^2)` a V grande (forma specifica);
- combinare i regimi (Chebyshev per V0 moderato, decadimento per V0 grande) in
  un'unica disuguaglianza. ATTENZIONE: cancellazione delicata a V0 grande (entrambi
  i termini stesso ordine); serve il next-order, non bound grezzi.
Alternativa: bound a due regioni con `-W'` esatto vicino a Vmin (non `-W'(2Vmin)`
uniforme, che perde un fattore ~1/4 asintoticamente).

## Script rilevanti
- `no_inversion_reduction.py`: riduzione + Lemma A (`symbolic_full_proof`) +
  Lemma B parziale (`lemma_B_partial`) + regime (`lemma_B_regime`).
- Convessita' di W: `/tmp/Wconvex.py` (spostare nel repo se serve).
- Sfera fotonica assente: `/tmp/photon.py`. PF r0=inf: `/tmp/pf_order.py`,
  `/tmp/pf_identify.py`.

## Stato
Lemma A: **chiuso**. Lemma B: aperto (grazing chiuso; scorciatoie escluse; Route 2
research-level, cancellazione delicata a r0 finito). Il no-inversion e' un teorema
nel regime scattering **modulo Lemma B**.
