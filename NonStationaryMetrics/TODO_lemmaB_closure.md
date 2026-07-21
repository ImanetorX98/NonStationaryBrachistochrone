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

## CASO CONGELATO (Vaidya mu=0 / TK a=0) — riduzione pulita + muro strutturale

Script: `no_inversion_schwarzschild_frozen.py` (+ `/tmp/schw_ineq.py`, `/tmp/check_grazing.py`).
Verificato a 1e-13 che, con `V0=V(r0)` FISSO (estremo fisso):
```
Phi_tau(Vmin) = sqrt(Vmin) int_{Vmin}^{V0} W(V)/sqrt(V-Vmin) dV
              = 2 sqrt(Vmin(V0-Vmin)) I(Vmin),   I(Vmin)=int_0^1 W(Vmin+(V0-Vmin)u^2)du
```
`W(V)=K/V'` fisso (>0, decrescente, convesso); `r(V)` = radice di cubica
`r^3-2r^2-bVr-2V=0` (b=E^2-1) => W ELEMENTARE-algebrico. Log-derivata:
```
Phi'/Phi = (1/2)(1/Vmin - 1/(V0-Vmin)) + I'/I
dPhi/dVmin<0  <=>  (STAR)  -I' > (V0-2Vmin)/(2 Vmin(V0-Vmin)) I,   -I'=int_0^1(-W'(A))(1-u^2)du
```
RISULTATI:
- **(STAR) e' TIGHT a r_min->r_pk+** (margine->0): r_pk e' DEFINITO da Phi'=0 = (STAR)
  con uguaglianza. => soglia genuinamente TRASCENDENTE; **nessun certificato elementare
  uniforme** puo' provare (STAR) su tutto r_min>r_pk. Bound Chebyshev `(2/3)int(-W')du`
  troppo debole (sotto rhs in piu' righe).
- **Sub-risultato ELEMENTARE**: `Vmin>=V0/2` (grazing) => entrambi i pezzi di Phi'/Phi
  <=0 => Phi decrescente. MA raggio V0/2 (~7.5M per r0=10) sta SOPRA i turning fisici
  (~4.7-6.4M): grazing da solo non li copre.
- Punti fisici: `Vmin << V0/2`, tra r_pk e raggio-grazing; (STAR) vale con margine
  positivo ma trascendente.

CONCLUSIONE: anche congelato (Schwarzschild), la chiusura NON e' elementare — la soglia
r_pk e' un punto trascendente dove (STAR) e' uguaglianza. Il "modulo Lemma B" ha ora una
RAGIONE strutturale, non solo mancanza di idee. Chiusura piena richiede o accettare lo
stato verificato-modulo-(STAR), o disuguaglianza sull'angolo di deflessione tight a r_pk
(CAS pesante / stima trascendente ad hoc).

## CASO CONGELATO — FORMA CHIUSA di Phi' + single-crossing (sessione 2026-07-21/22)

Script: `no_inversion_schwarzschild_closedform.py`. Progresso NUOVO oltre a (STAR):

**Forma chiusa di dPhi/dVmin** (verificata ~1e-9 vs diff. finite). Con
`G(x)=int_x^{V0}W(A)/sqrt(A-x)dA`, IBP che sposta d/dV su W (uccide singolarita' sqrt):
`G=2W(V0)sqrt(V0-x)-2 int_x^{V0}W'(t)sqrt(t-x)dt`, poi derivo (nessun bordo singolare):
```
sqrt(x) Phi'(x) = W(V0)(V0-2x)/sqrt(V0-x) + int_x^{V0} W'(t)(2x-t)/sqrt(t-x) dt
```
CONSEGUENZE rigorose:
- **Grazing x>=V0/2**: (V0-2x)<=0, (2x-t)>=0, W'<0 => entrambi <=0 => Phi'<=0. Elementare.
- **Quarto x>=V0/4**: split a t=2x, P1=int_x^{2x}(-W')(2x-t)/sqrt(t-x), P2=int_{2x}^{V0}(...).
  Integrali esatti: int_x^{2x}(2x-t)/sqrt(t-x)dt=(4/3)x^{3/2};
  int_{2x}^{V0}(t-2x)/sqrt(t-x)dt=2[(V0-x)^{3/2}/3 - x sqrt(V0-x)+(2/3)x^{3/2}].
  -W' decrescente (W convessa) => `P1-P2 >= (-W'(2x))(2/3)sqrt(V0-x)(4x-V0) >=0` per x>=V0/4.
  => monotonia controllata per x>=V0/4 (estende grazing V0/2, prima solo abbozzo ~4Vmin).
  **Finestra elementare piena**: se `V0<=4 Vpk` (cioe' r0<=R*(E), ~4.3 per E=1.2) TUTTI i
  turning r>r_pk hanno x>=V0/4 => no-inversion **pienamente elementare per r0 piccolo**.

**Riduzione single-crossing (copre ogni r0)**: Lemma B <=> Phi_tau ha UN SOLO punto
critico. Fatto standard: se a OGNI critico L=log Phi ha L''<0 => al piu' un critico; con
L'(0+)>0, L'(V0-)<0 => esattamente uno => Phi single-peaked => Lemma B. A critico
(L'=0 => I'/I=-(1/2)(a-c), a=1/x, c=1/(V0-x)):
```
L'' = -(3/4)(a^2+c^2) + (1/2)ac + I''/I         (a critico)
(DAGGER)   I''/I < (3/4)(a^2+c^2) - (1/2)ac
```
I=int_0^1 W(A)du, I''=int_0^1 W''(A)(1-u^2)^2du, A=x+(V0-x)u^2. Verificato 40 cifre,
margine >0 STRETTO per ogni (E,r0) finito; margine->0 per r0->inf (il picco r_pk fugge a
inf: r0=inf e' monotono). r_pk~2.94-2.99 (~3M) per tutti gli E testati.

GAP RESIDUO: regione r_pk<r_min<r_quarter (x/V0 in ~0.02..0.25). Serve stima sharp che
sfrutta la CONCENTRAZIONE di -W' vicino al turning (il bound uniforme -W'(2x) perde ~1/4).
Log-convessita' di W: `W W''-(W')^2` ha numeratore polinomiale (in schw_Wode.py) da
certificare positivo (SOS/Sturm) — ingrediente, non chiusura.

## ROUTE B TENTATA (stima concentrazione -W') — MURO CONFERMATO trascendente
Scripts: `/tmp/schw_cheb.py`, `/tmp/schw_cheb2.py` (antideriv. kernel verificata 1e-14).
Riformulazione pulita via W(V0)=W(x)-int_x^{V0}(-W'):  serve
```
int_x^{V0} (-W'(t)) ker(t) dt  >  W(x)(V0-2x)/sqrt(V0-x),
ker(t)=(2x-t)/sqrt(t-x)+(V0-2x)/sqrt(V0-x)  > 0 su (x,V0), ->0 a t=V0, ->+inf a t=x.
```
`ker` e `-W'` ENTRAMBE positive decrescenti => Chebyshev integrale applicabile. Integrali
esatti: int_x^{V0} ker dt = sqrt(V0-x)(V0+2x)/3. Chebyshev singolo intervallo da':
```
int(-W')ker >= (W(x)-W(V0))(V0+2x)/(3 sqrt(V0-x))  =>  chiude sse  2W(x)(4x-V0) >= W(V0)(V0+2x)
```
cioe' soglia ~V0/4 (in pratica ~x/V0≈0.3-0.35 per il termine W(V0)). Chebyshev e' lo
strumento OTTIMALE per due funzioni monotone concordi => non si fa meglio senza sfruttare
la FORZA della concentrazione, non solo la monotonia.
- Split multi-intervallo: converge al vero (che sta sopra bnd) ma margine->0 a r_pk =>
  nessun numero FINITO di split chiude uniformemente.
- Bound uniforme -W'(2x), convessita' W, Phi''<0 a critico (=d/dx[sqrt(x)Phi']<0, unico
  termine ostile A2=int W''(2x-t)/sqrt(t-x)dt concentrato a t=x): tutti stesso muro.
- **Vpk CRESCE con r0** (r0=10->Vpk~2.46; r0=20->6.08; Vpk/V0->0, Vpk~V0^0.54): la regione
  tight INSEGUE r_pk che si sposta => nessun trucco a regione fissa. r0=inf e' monotono
  (nessun critico): margine->0 in quel limite.
CONCLUSIONE ROUTE B: la chiusura NON e' elementare (confermato, non solo congetturato).
Restano solo: (i) proof computer-assistito (interval arithmetic) su range compatto di
(E,r0) = teorema onesto per range fisico stabilito, non pen-and-paper, non tutti r0;
(ii) stima trascendente ad hoc con idea nuova (le vie standard sono esaurite).

## Stato
Lemma A: **chiuso**. Lemma B (congelato Schwarzschild): **forma chiusa di Phi' provata**;
monotonia **elementare per x>=V0/4** (=> pienamente elementare per r0<=R*(E)); ridotta a
single-crossing = **una** disuguaglianza puntuale (DAGGER) a critico, verificata 40 cifre.
GAP: stima trascendente in (r_pk,r_quarter). Kerr congelato: W non elementare (genus-2),
resta verificato-numerico. Il no-inversion e' teorema modulo (DAGGER), con sotto-regime
x>=V0/4 dimostrato.
