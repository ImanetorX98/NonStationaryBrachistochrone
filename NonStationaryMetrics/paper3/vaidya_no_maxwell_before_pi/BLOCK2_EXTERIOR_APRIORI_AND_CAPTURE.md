# Blocco 2d — stime esterne a priori, angolo non radiale e cattura

6 ottobre 2026. Questo checkpoint aggiunge 18 teoremi e porta il progetto a
160 teoremi pubblici. Nel problema di Vaidya via B, il confinamento necessario
alla continuazione viene ora dedotto dall'evoluzione esterna. La permanenza
di \(\alpha\) in \((0,\pi)\) viene dedotta dal dato iniziale, non
assunta sull'intero tratto. È inoltre verificata l'alternativa a un estremo
finito: cattura oppure prolungamento esterno.

Al checkpoint descritto qui non era ancora costruita la soluzione massimale
esterna. Il successivo [blocco 2e](BLOCK2_MAXIMAL_EXTERIOR_AND_GLOBAL_CONE.md)
la costruisce e ricava flussi futuri globali dopo ingresso stretto nel cono,
per il sistema coordinato. Il collegamento completo ai dati geometrici e il
teorema completo sui Maxwell restano aperti nella formalizzazione.

## 1. Ipotesi effettive

Si considera una soluzione esistente \(f=(r,\alpha)\) su \([a,b)\),
con \(a<b<\infty\), per il sistema di direzione già derivato:

\[
q=\widehat E^2>1,\qquad w=q-1+2m(v)/r(v),
\]

\[
r_v=-w+\sqrt{qw}\cos\alpha,\qquad
\alpha_v=\frac{\sin\alpha}{r}B(q,w,\cos\alpha).
\]

La massa è positiva al lancio, \(m(a)>0\), e non decrescente su
\([a,b]\). È prescritta anche al tempo finale e di classe \(C^1\)
vicino a \(b\), dove si costruisce il nuovo tratto. Si assume che il
tratto precedente sia esterno, \(r(v)>2m(v)\), e che il solo angolo
iniziale sia non radiale, \(0<\alpha(a)<\pi\).

Queste ipotesi non includono un limite finale del raggio, un compatto di
confinamento, una soluzione oltre \(b\), né la permanenza angolare futura.
La soglia q≥3/2 dei lemmi di curvatura non viene modificata da questi risultati.
La soluzione precedente e la sua appartenenza al dominio esterno su
\([a,b)\) rimangono premesse: è un teorema di continuazione, non ancora
la costruzione globale della soluzione massimale.

## 2. Raggio e massa: il compatto si ricava dal problema

Il dominio esterno dà

\[
q-1\leq w\leq q.
\]

Il limite già certificato \(r_v\leq q\), con il teorema del valor medio,
produce

\[
r(v)\leq r(a)+q(v-a)\leq r(a)+q(b-a).
\]

La monotonicità e positività della massa danno invece

\[
r(v)\geq2m(a)>0,\qquad0\leq m(v)\leq m(b).
\]

Quando l'angolo è nella striscia chiusa, il compatto esplicito è

\[
C=\big([2m(a),r(a)+q(b-a)]\times[0,\pi]\big)\times[0,m(b)].
\]

Ogni punto di \(C\) ha raggio strettamente positivo e massa non negativa.
`exterior_finite_regular_compact_tube` certifica compattezza, regolarità e
confinamento: non prende un compatto già fornito. La striscia angolare,
inizialmente premessa di questo lemma intermedio, viene eliminata nei
teoremi composti del paragrafo seguente.

## 3. Perché l'angolo non raggiunge 0 o π

Ponendo \(y=\sin\alpha\), la regola della catena dà

\[
y_v=k(v)y,\qquad
k(v)=\frac{\cos\alpha(v)}{r(v)}B(q,w(v),\cos\alpha(v)).
\]

Il limite già provato per \(|B|\), con \(r\geq\rho>0\), implica

\[
|k(v)|\leq\frac{\operatorname{directionBound}(q)}{\rho}.
\]

Questo limite vale per **ogni angolo reale**: usa solo
\(-1\leq\cos\alpha\leq1\), non la permanenza di \(\alpha\)
in una striscia. Nel problema esterno si prende \(\rho=2m(a)\).

Se \(y\) raggiungesse zero in un tempo finito, l'unicità all'indietro
per l'ODE lineare imporrebbe che coincida con la soluzione nulla fino al
lancio. Ciò contraddice \(y(a)>0\). `linear_ode_initial_nonzero` prova
questo passaggio senza assumere il segno di \(y\) lungo il tratto.

Per uscire da \((0,\pi)\), la funzione continua \(\alpha\) dovrebbe
attraversare 0 o π; il teorema dei valori intermedi produrrebbe quindi uno
zero di \(\sin\alpha\). Si conclude

\[
0<\alpha(v)<\pi\qquad\text{su ogni tratto finito esistente esterno}.
\]

`exterior_finite_angle_invariance` non assume la striscia angolare come
premessa. Il teorema composto di continuazione usa questa conclusione per
costruire il compatto e applicare il blocco 2c. Non c'è circolarità fra
permanenza angolare e stima del coefficiente lineare.

## 4. Controllo dell'estremo finito

Il blocco 2c ricava un limite regolare e incolla una soluzione oltre \(b\).
La continuità del tratto prolungato e della massa dà

\[
r(b)\geq2m(b).
\]

Per includere l'angolo al tempo finale, si ripete il confronto lineare su
un tratto \([c,b]\), con \(a<c<b\). L'ODE ora vale anche in \(b\);
i limiti \(r\geq2m(a)\), \(q-1\leq w\leq q\) valgono anche
all'estremo. Si ricava dunque

\[
\boxed{r(b)\geq2m(b),\qquad0<\alpha(b)<\pi.}
\]

In particolare, non basta sapere che l'angolo è interno prima di \(b\):
la stretta interiorità al tempo finale è provata separatamente.

`exterior_finite_endpoint_capture_or_continues` conclude:

- se \(r(b)=2m(b)\), l'estremo è sulla superficie di cattura;
- se \(r(b)>2m(b)\), la continuità produce un intervallo esterno non radiale
  più lungo, sul quale la curva prolungata soddisfa già l'ODE.

Questa alternativa distingue l'uscita dal dominio fisico dalla perdita di
esistenza del campo coordinato.

## 5. Cattura trasversale

Al bordo \(r=2m>0\), si ha \(w=q\). Per un evento con
\(0<\alpha<\pi\), \(m_v=\mu\geq0\), Lean verifica

\[
\frac{d}{dv}(r-2m)=q(\cos\alpha-1)-2\mu<0.
\]

Il lemma `nonradial_capture_gap_derivative` assume l'equazione di direzione,
la derivata della massa all'evento e la non negatività di tale derivata.
La stretta negatività è una conclusione. Il lemma non costruisce ancora il
primo tempo di cattura né formalizza un teorema completo sui tempi di cattura
di famiglie di lanci.

## 6. Il dato iniziale viene conservato

La prima forma dei teoremi di continuazione afferma coincidenza sul tratto
aperto \((a,b)\). `direction_continuation_preserves_launch` completa il
raccordo al tempo iniziale: mantiene la curva precedente per \(v\leq a\)
e usa il prolungamento per \(v>a\). La coincidenza in \((a,b)\) rende
questa funzione uguale alla curva precedente in un intero intorno di \(a\).

`exterior_nonradial_initial_value_continuation` conclude una curva \(H\)
che coincide con \(f\) su **\([a,b)\)** e soddisfa l'ODE su
**\([a,b+\varepsilon)\)**. Dato e derivata al lancio sono quindi
preservati, non soltanto la parte successiva del tratto.

## 7. Prossimi passaggi ancora aperti

1. Costruire una soluzione massimale esterna e l'unicità sull'intero intervallo
   comune; usare l'alternativa finale per dimostrare che un suo estremo finito
   è un tempo di cattura.
2. Collegare l'esclusione della cattura dopo ingresso in un cono uscente alla
   continuazione, così da costruire una soluzione futura globale e poi il
   `DirectionFlow` necessario ai lemmi di fuga già certificati.
3. Completare il collegamento geometrico: metrica, connessione di Levi-Civita,
   schermo/Jacobi, properness, rivestimento e criterio finale sui Maxwell.

Questi punti non sono conclusioni del presente checkpoint. In particolare,
non si usa un `DirectionFlow` futuro globale per provare le stime finite o
la continuazione qui registrate.

## 8. Verifica

`ExteriorFiniteTube.lean` contiene 7 nuovi teoremi;
`AngularFiniteInvariance.lean` ne contiene 11. Tutti sono inclusi nell'audit
delle dipendenze assiomatiche. La verifica usa Lean 4.24.0 e mathlib v4.24.0;
log e hash dei sorgenti sono in `verification/`.
