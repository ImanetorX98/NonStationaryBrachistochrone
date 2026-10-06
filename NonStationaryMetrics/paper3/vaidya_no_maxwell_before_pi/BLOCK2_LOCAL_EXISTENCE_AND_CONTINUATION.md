# Blocco 2c — esistenza locale, unicità e continuazione sui compatti regolari

**Aggiornamento successivo:** il [blocco 2d](BLOCK2_EXTERIOR_APRIORI_AND_CAPTURE.md)
deduce il compatto dall’evoluzione esterna, esclude gli angoli radiali dal dato
iniziale e certifica l’alternativa cattura/prolungamento all’estremo finito.
Questo documento conserva il checkpoint 2c e i suoi limiti originari.

6 ottobre 2026. Questo checkpoint aggiunge 12 teoremi e porta il progetto a
142 teoremi pubblici. Sono ora formalizzati esistenza locale del sistema di
direzione, unicità locale, permanenza iniziale nel dominio esterno e
continuazione effettiva oltre un estremo finito, sotto confinamento in un
compatto regolare. La continuazione include la prova di incollamento.

Non sono ancora costruiti gli intervalli massimali, né un `DirectionFlow`
futuro globale a partire da ogni dato geometrico. Il teorema geometrico sui
Maxwell resta una formalizzazione parziale.

## 1. Sistema e ipotesi di regolarità

Con \(q=\widehat E^2>1\), \(w=q-1+2m(v)/r\), lo stato è
\(x=(r,\alpha)\) e il campo è

\[
F_q(m,x)=\left(-w+\sqrt{qw}\cos\alpha,
\frac{\sin\alpha}{r}B(q,w,\cos\alpha)\right),
\]

\[
B(q,w,c)=c(2w-q+1)-\frac12\sqrt{q/w}(3w-q+1).
\]

Il sistema è \(x_v=F_q(m(v),x)\). L'esistenza qui dimostrata richiede
\(m\) di classe \(C^1\) vicino al tempo iniziale, \(r_0>0\) e
\(m(v_0)\geq0\). È una condizione sufficiente esplicita, adatta a una
storia di massa liscia; non è una dimostrazione di esistenza sotto la sola
continuità della massa. L'unicità locale, invece, richiede solo continuità
della massa al tempo considerato.

Questi lemmi non abbassano la soglia \(q\geq3/2\) dei confronti di
curvatura e non trattano la connessione di Levi-Civita.

## 2. Esistenza senza una traiettoria fornita

Si introduce un tempo interno \(s\) nello stato aumentato:

\[
\widetilde F_q(x,s)=\big(F_q(m(s),x),1\big).
\]

Il lemma `liftedDirectionField_contDiffAt` ricava la regolarità \(C^1\)
di questo campo dalla regolarità già provata di \(F_q\) in stato e massa
e dalla regolarità di \(m\). Picard–Lindelöf di mathlib costruisce quindi
una curva \((x(v),s(v))\) con dato iniziale \((x_0,v_0)\).

Poiché \(s_v=1\) e \(s(v_0)=v_0\), il teorema del valor medio dà
\(s(v)=v\). Proiettando sulla prima componente si ottiene

\[
\exists\varepsilon>0,\ \exists x:\quad x(v_0)=x_0,\qquad
x_v=F_q(m(v),x(v))\quad\text{per }|v-v_0|<\varepsilon.
\]

`direction_local_solution_exists` non assume una soluzione già esistente.

## 3. Unicità e dominio fisico locale

La regolarità del campo in \((x,m)\) produce una costante Lipschitz e un
intorno comune del dato. La continuità delle due traiettorie, ricavata dalle
loro ODE, e quella della massa le mantengono in tale intorno per tempi vicini.
Per uno stesso valore temporale,

\[
\operatorname{dist}\big((x,m(v)),(y,m(v))\big)
=\operatorname{dist}(x,y)
\]

nella metrica prodotto usata da Lean. Il limite Lipschitz congiunto dà
quindi il limite spaziale necessario a Grönwall.

`direction_local_solution_unique` conclude che due soluzioni con lo stesso
dato iniziale coincidono in un intorno temporale. La costante Lipschitz e la
permanenza in un intorno non sono premesse separate. Il lemma
`direction_local_initial_value_problem` compone esistenza e unicità del germe
di soluzione. Non afferma unicità di estensioni arbitrarie fuori dal loro dominio.

Per un lancio con

\[
m(v_0)>0,\qquad r_0>2m(v_0),\qquad0<\alpha_0<\pi,
\]

`direction_exterior_local_solution_exists` restringe l'intervallo costruito
e conclude, su tutto il nuovo intervallo,

\[
m(v)>0,\qquad r(v)>2m(v),\qquad0<\alpha(v)<\pi.
\]

Questa permanenza iniziale segue dalla continuità e dalle disuguaglianze
strette del lancio. Non è ancora un risultato per ogni tempo futuro.

## 4. Limite a un estremo finito

Sia \(f\) una soluzione su \((a,b)\), \(a<b<\infty\), con

\[
(f(v),m(v))\in C,
\]

dove \(C\subset\mathbb R^2\times\mathbb R\) è compatto e ogni suo
punto soddisfa \(r>0\), \(m\geq0\). La continuità del campo sul
compatto produce una costante finita \(K\) tale che

\[
\|f_v(v)\|\leq K.
\]

Il teorema del valor medio rende \(f\) Lipschitz su \((a,b)\).
`planar_lipschitz_extension` estende separatamente le due coordinate mediante
il teorema di estensione Lipschitz di mathlib e le ricompone. L'estensione
\(F\) è continua e coincide con \(f\) sull'intervallo, quindi

\[
f(v)\longrightarrow F(b)=z\qquad(v\to b^-).
\]

Non si assume il limite. La continuità della massa e la chiusura del compatto
danno \((z,m(b))\in C\), dunque il limite è ancora regolare.
`compact_direction_endpoint_restart` costruisce una nuova soluzione locale
attraverso questo punto, usando \(m\) di classe \(C^1\) vicino a \(b\).

L'estensione Lipschitz ausiliaria \(F\) non viene dichiarata soluzione
fuori dall'intervallo precedente.

## 5. Incollamento e continuazione effettiva

La continuità del campo implica che la derivata di \(F\) converge a
\(F_q(m(b),F(b))\) da sinistra. Il lemma di mathlib sull'estensione della
derivata al bordo produce

\[
D^-F(b)=F_q(m(b),F(b)).
\]

`direction_endpoint_backward_agreement` applica Grönwall all'indietro alla
curva entrante e alla nuova soluzione, con uguale valore in \(b\). Deduciamo
che coincidono su un intervallo \([c,b]\), con \(a<c<b\).

Si definisce allora \(h(v)=f(v)\) per \(v<b\) e \(h(v)=g(v)\)
per \(v\geq b\). Vicino a \(b\) la funzione coincide con \(g\),
anche da sinistra; soddisfa quindi l'ODE anche nel punto di raccordo.

`compact_direction_continues` conclude precisamente

\[
\exists\varepsilon>0,\ \exists h:\quad
h=f\text{ su }(a,b),\qquad
h_v=F_q(m(v),h(v))\text{ su }(a,b+\varepsilon).
\]

Il tratto uscente e la sua coincidenza con quello entrante sono conclusioni.
Non si presume una soluzione globale né un prolungamento già dato.

## 6. Limite del risultato e lavoro successivo

Il confinamento nel compatto regolare è ancora una premessa. Escludere un
estremo finito del dominio massimale richiede ora ricavare questo confinamento
dalle ipotesi del problema. In un intervallo esterno finito, massa positiva
non decrescente e \(r_v\leq q\) suggeriscono

\[
r\geq2m(v_0)>0,\qquad
r(v)\leq r(v_0)+q(v-v_0).
\]

La continuità della massa dà limiti finiti della massa sul tratto; occorre
formalizzare l'intero incastro con gli intervalli massimali e il controllo
\(0<\alpha<\pi\). Queste stime e il loro impiego non sono conclusioni
nuove di questo checkpoint.

Il campo coordinato è regolare anche a \(r=2m>0\): il teorema può
prolungare l'ODE attraverso quel bordo. **Questo non dimostra la permanenza
nel dominio esterno**. La cattura va distinta da una perdita di esistenza.
Serve inoltre collegare le soluzioni costruite al flusso geometrico e ai
vecchi lemmi che assumono un `DirectionFlow` futuro globale.

## 7. Sorgenti e audit

- `LocalDirectionExistence.lean`: 5 teoremi.
- `FiniteEndpointLimit.lean`: 4 teoremi.
- `RegularEndpointContinuation.lean`: 3 teoremi.

`verify.py` compila il progetto e controlla le dipendenze assiomatiche
transitive di tutti i teoremi pubblici. Ambiente fissato: Lean 4.24.0,
mathlib v4.24.0. Log e hash dei sorgenti sono in `verification/`.
