# Blocco 2a — dalla metrica coordinata all'equazione di direzione

**Aggiornamento successivo:** il [blocco 2b](BLOCK2_LOCAL_TIME_AND_NULL_PROPAGATION.md)
costruisce l’inversa temporale locale e certifica propagazione della nullità e
conservazione di L per il flusso Hamiltoniano dichiarato. Questo documento
conserva la descrizione del checkpoint 2a e dei suoi limiti originari.

6 ottobre 2026. Questo documento registra 29 nuovi teoremi Lean, che portano il
progetto a 117 teoremi pubblici. Il blocco 2 è **parzialmente completato**:
è verificata la derivazione coordinata del sistema di direzione, sotto una
riparametrizzazione differenziabile fornita. Non sono ancora certificati
l'esistenza dell'inversa temporale, la continuazione delle soluzioni e il
collegamento geometrico con Jacobi. Il teorema completo sui Maxwell rimane aperto
nella formalizzazione.

## 1. Oggetti e dominio

Poniamo \(q=\widehat E^2\), \(w=q-1+2m/r\). Il ramo usato richiede
\(q>0\), \(r>0\), \(w>0\), \(L>0\). Le applicazioni ai precedenti
lemmi di fuga richiedono inoltre \(q>1\) e le loro specifiche ipotesi sulla
massa e sull'esistenza futura. Il confronto delle curvature conserva la
propria soglia \(q\geq3/2\): questa derivazione non la rimuove.

In coordinate \((v,r,\varphi)\), la forma quadratica dichiarata è

\[
G(k,k)=-wV^2+\frac{(R+wV)^2}{q}+r^2P^2.
\]

Lean dimostra che gli operatori espliciti di abbassamento e innalzamento degli
indici sono inversi e che la forma inversa è

\[
G^{-1}(p,p)=-\frac{P_v^2}{w}+2P_vP_r+(q-w)P_r^2+\frac{L^2}{r^2}.
\]

Questo certifica l'inversione della metrica coordinata. Non è ancora una
costruzione della metrica lorentziana su una varietà né una prova della
corrispondenza tra geodetiche di Levi-Civita e flusso Hamiltoniano.

## 2. Hamiltoniano e ramo nullo futuro

Per \(H=G^{-1}(p,p)/2\), sono verificate le derivate parziali che danno

\[
v_\lambda=-P_v/w+P_r,\qquad
r_\lambda=P_v+(q-w)P_r,\qquad
\varphi_\lambda=L/r^2,
\]

\[
P_{r,\lambda}=-H_r
=\frac{m}{r^2}\left(\frac{P_v^2}{w^2}-P_r^2\right)+\frac{L^2}{r^3}.
\]

La derivata rispetto a \(r\) tiene fisso \(v\), quindi fissa il valore
istantaneo \(m(v)\). **Non assume che la massa sia costante lungo la curva**.
Nel teorema di riparametrizzazione \(m\) e \(P_v\) rappresentano valori
all'evento, non costanti del moto. In particolare, non viene postulata la
conservazione di \(P_v\) nel caso non stazionario.

Con

\[
S=\sqrt{qP_r^2+L^2/r^2}>0,
\]

nullità e \(v_\lambda>0\) selezionano univocamente

\[
P_v=wP_r-\sqrt w\,S,\qquad v_\lambda=S/\sqrt w>0.
\]

Sono certificati sia la nullità e la direzione futura di questo ramo, sia la
sua unicità.

## 3. Ricostruzione dell'angolo e sua dinamica

Definiamo

\[
c=\frac{\sqrt q\,P_r}{S},\qquad s=\frac{L}{rS},\qquad
\alpha=\arccos c.
\]

Lean verifica \(c^2+s^2=1\), \(s>0\), \(0<\alpha<\pi\),
\(\cos\alpha=c\) e \(\sin\alpha=s\). La scelta \(L>0\) fissa
l'orientazione; il ramo \(L<0\) e la simmetria che lo riconduce a questo
non sono trattati da questi nuovi lemmi. Il caso radiale \(L=0\) resta escluso.

Dividendo le equazioni Hamiltoniane per \(v_\lambda\), si ricavano

\[
r_v=-w+\sqrt{qw}\,c,\qquad
\varphi_v=\frac{\sqrt w}{r}s,
\]

\[
P_{r,v}=\frac{mS}{r^2\sqrt w}-\frac{2mP_r}{r^2}
+\frac{\sqrt w L^2}{r^3S}.
\]

La regola della catena, con \(L\) costante, dà

\[
c_v=\frac{\sqrt q L^2}{r^2S^3}
\left(P_{r,v}+\frac{P_r}{r}r_v\right)
=-\frac{s^2}{r}B(q,w,c),
\]

\[
B(q,w,c)=c(2w-q+1)-\frac12\sqrt{q/w}(3w-q+1).
\]

La derivata di \(\arccos\) produce infine

\[
\boxed{\alpha_v=\frac{\sin\alpha}{r}
B(q,w,\cos\alpha).}
\]

L'equazione angolare è una **conclusione**, non una premessa. Non compare una
divisione per \(P_r\) o per \(\cos\alpha\): la derivazione comprende quindi
\(P_r=0\), dove \(r_v=-w\); comprende anche i turning point radiali
\(r_v=0\). Questi sono due casi distinti. Non compare \(m'\) nell'ODE
angolare, benché \(m(v)\) possa variare. La conservazione di \(L\) e la
propagazione della nullità lungo il flusso completo devono ancora essere
formalizzate: i nuovi lemmi assumono rispettivamente un \(L\) fisso e la
nullità all'evento.

## 4. Teorema composto e limite preciso

`ViaB.affine_hamiltonian_to_direction` prende le tre equazioni Hamiltoniane
affini per \(r,P_r,\varphi\), nullità, direzione futura e una mappa
\(\ell(v)\) con derivata

\[
\ell_v=1/v_\lambda.
\]

Conclude le tre equazioni in \(v\) per \(r,\alpha,\varphi\). Il teorema
è puntuale; può essere applicato a ogni evento di una soluzione regolare.
La premessa su \(\ell\) non include l'identità di inversa di una specifica
funzione temporale: basta la regola differenziale per la conclusione puntuale.
Costruire l'inversa locale effettiva rimane un passaggio separato.

Non si assume il sistema angolare per ricavarlo e non si usa l'esistenza
globale per dimostrarla. I vecchi teoremi su `DirectionFlow` continuano a
richiedere soluzioni future globali; questa premessa non è eliminata dal
presente checkpoint.

## 5. Passaggi prioritari ancora da dimostrare

1. Completare il flusso Hamiltoniano dipendente da \(m(v)\), conservazione
   di \(L\) e nullità; costruire l'inversa locale di \(v(\lambda)\) usando
   \(v_\lambda>0\).
2. Dimostrare esistenza locale, unicità e continuazione su intervalli massimali
   per il campo di direzione, con le ipotesi di regolarità di \(m\) esplicite.
3. Collegare il sistema a Jacobi, allo schermo e alla formula delle curvature.
   Occorre verificare il parametro affine comune del confronto Sturm.

Una strategia analitica di continuazione, **ancora non certificata in Lean**,
è lavorare su un intervallo temporale finito con massa continua non decrescente
e positiva. Nel dominio esterno, \(r>2m(v)\geq2m(a)>0\); il limite
\(r_v\leq q\) fornisce un limite superiore finito. Il campo è regolare
sul compatto risultante, incluso il bordo \(r=2m\) come campo coordinato.
L'eventuale uscita dal dominio esterno va distinta dalla perdita di esistenza
dell'ODE. Per massa differenziabile, al bordo

\[
(r-2m)_v=q(\cos\alpha-1)-2m_v<0
\]

per \(0<\alpha<\pi\) e \(m_v\geq0\). Servono inoltre un controllo
del mancato raggiungimento di \(\alpha=0,\pi\) in tempo finito e un
teorema di continuazione: la sola compattezza non è la prova completa.
L'ingresso permanente nel cono \(r>R>2M\), quando le ipotesi sui limiti
della massa valgono, potrà poi escludere l'uscita dal dominio e sostenere
l'esistenza futura. Non va assunto un `DirectionFlow` globale per provare
questo ultimo punto.

## 6. Verifica riproducibile

I moduli aggiunti sono `MetricHamiltonian` (9 teoremi),
`NullHamiltonianDynamics` (14), `DirectionFromHamiltonian` (5),
`AffineReparametrization` (1). `AxiomAudit.lean` copre ogni teorema pubblico;
`verify.py` compila e controlla tutte le dipendenze assiomatiche transitive.
Gli esiti e gli hash dei sorgenti sono in `verification/`.
