# Blocco 2h — variazione geodetica e regolarità della connessione

6 ottobre 2026. Il checkpoint aggiunge **12 teoremi** in due moduli al
blocco 2g. Totale: **234 teoremi pubblici**, 49 sorgenti Lean.

## Risultato e limite della certificazione

È formalizzata la differenziazione dell'equazione geodetica rispetto
al parametro di una famiglia. La derivabilità dei coefficienti della
connessione di Fermat è dedotta da quella della massa campionata, del
suo tasso e del raggio, con denominatori regolari. Una massa fisica C²
fornisce le due prime condizioni lungo una coordinata temporale
differenziabile. È verificata anche una famiglia effettiva di geodetiche:
la traslazione angolare di una curva già costruita.

Il risultato generale è **condizionato alla regolarità della famiglia**.
Non costruisce ancora la famiglia differenziabile dei lanci, non prova
lo scambio delle derivate miste e non identifica il coefficiente dello
schermo con la funzione `screenGap`. Il teorema geometrico sui Maxwell
rimane aperto.

## 1. Differenziazione dell'equazione

A un istante affine fissato, ε è il parametro di famiglia. Indichiamo con
Γ(ε) i coefficienti campionati sulla curva, con u(ε) la sua velocità affine
e con b(ε) l'accelerazione. Poniamo

\[
Q^i(\varepsilon)=\sum_{j,k}\Gamma^i{}_{jk}(\varepsilon)
 u^j(\varepsilon)u^k(\varepsilon),\qquad
b^i(\varepsilon)=-Q^i(\varepsilon).
\]

Se DG=∂εΓ e U=∂εu, Lean deriva la regola completa del prodotto:

\[
\partial_\varepsilon Q^i=
\sum_{j,k}\left[
 DG^i{}_{jk}u^ju^k+\Gamma^i{}_{jk}U^ju^k+
 \Gamma^i{}_{jk}u^jU^k\right].
\]

La simmetria Γⁱⱼₖ=Γⁱₖⱼ combina i due ultimi termini. Quando sono
identificate ∂εu=J′ e ∂εb=J″, segue

\[
J^{i\prime\prime}+2\sum_{j,k}\Gamma^i{}_{jk}u^jJ^{k\prime}
 +\sum_{j,k}DG^i{}_{jk}u^ju^k=0.
\]

`connectionQuadratic_hasDerivAt` deriva la contrazione;
`connectionVariation_torsion_free` combina i termini;
`geodesic_family_acceleration_variation` trasferisce la derivata
all'accelerazione usando l'equazione su un intorno di ε.
`coordinate_jacobi_of_commuted_variation` conclude l'identità sopra
con le identificazioni delle derivate come premesse esplicite.

Non è assunta l'equazione linearizzata per dimostrarla. Tuttavia,
i nomi J′ e J″ in quest'ultimo lemma indicano valori forniti dalle
identificazioni: il lemma non prova da solo che siano le derivate
affini di un campo J=∂εx. Questo richiede un teorema distinto sulla
famiglia a due parametri.

## 2. Collegamento alla connessione effettiva di Fermat

`fermat_connectionQuadratic` identifica Q con la contrazione già
usata nell'equazione geodetica verificata nel blocco 2g.

`fermatInverseMetric_family_differentiable` deriva la regolarità
dell'inversa lungo una famiglia per r≠0, w≠0.
`fermatMetricPartial_family_differentiable` tratta le derivate della
metrica, richiedendo anche q≠0 e la derivabilità della massa e del tasso.
`fermatChristoffel_family_differentiable` combina tali risultati
tramite la formula di Koszul. Non assume che Γ sia derivabile.

`fermatChristoffel_sampled_C2_mass_differentiable` applica il risultato
alle funzioni m(v(ε)), m_v(v(ε)), r(ε). Le premesse sono massa C²
all'evento, v e r differenziabili nel parametro e q,r,w non nulli.
Lean deriva la regolarità di m_v tramite la derivata di una funzione C²
e applica la catena. C² è qui una condizione sufficiente; non è
rivendicata come ipotesi minima o necessaria per ogni caso speciale.

`fermat_coordinate_jacobi_of_family` fornisce l'equazione linearizzata
per la connessione effettiva di Fermat. Il termine DG è la **derivata
reale della funzione dei coefficienti lungo la famiglia**, non un
coefficiente arbitrario assunto uguale alla curvatura.

Resta da esplicitare DGⁱⱼₖ=(∂ₗΓⁱⱼₖ)Jˡ per la famiglia dei lanci e
ricomporre l'identità nella forma covariante, con una convenzione di
curvatura dichiarata. Nessuna identificazione con R(J,u)u o con le
curvature dello schermo è ancora certificata da questi moduli.

## 3. Famiglia concreta e limite del controllo di simmetria

`coordinate_null_geodesic_angle_shift` prova che, per ogni costante c,
(v,r,φ+c) soddisfa la stessa equazione nulla affine futura di (v,r,φ).
La connessione e la metrica non dipendono dalla coordinata φ; la
derivata affine di φ+c coincide con quella di φ.

`angle_shift_family_position_derivative` calcola effettivamente
∂ε(v,r,φ+ε)=(0,0,1). Il vettore non si annulla al lancio e non è
il campo scalare di confronto u=(r₀/L)r sinφ, né la variazione dei
lanci dallo stesso evento. È una verifica della simmetria angolare nel
chart planare; non costruisce la variazione fuori piano della geometria
tridimensionale.

## 4. Regolarità: cosa cambia rispetto al blocco precedente

La massa C¹ bastava per le singole geodetiche coordinate costruite.
La differenziazione classica della connessione introduce la variazione
del tasso m_v: non è corretto dedurla dalla sola ipotesi C¹.
Il teorema con massa C² chiude questo punto per i coefficienti,
ma non fornisce automaticamente regolarità congiunta della famiglia,
scambio delle derivate o regolarità dell'orologio affine rispetto al lancio.

Occorre mantenere un unico parametro affine della famiglia, verificando
la dipendenza dell'orologio e dell'inversa dal lancio. La normalizzazione
L>0 dei precedenti costruttori non sostituisce questa prova.

## 5. Audit e prosecuzione

La verifica completa con `verify.py` è passata su tutti i 234 teoremi:
compilazione, copertura dell’audit e controllo degli assiomi transitivi. I log e gli
hash dei sorgenti sono in `verification/`. Sono ammessi soltanto gli
assiomi standard `propext`, `Classical.choice`, `Quot.sound`.

I prossimi obblighi sono:

1. Costruire la famiglia dei lanci su un intervallo comune e provarne
   la regolarità congiunta; includere parametro affine e inversa.
2. Provare le identificazioni delle derivate miste e applicare i lemmi
   di variazione a J=∂εx, con dati iniziali geometrici verificati.
3. Costruire lo schermo nullo e identificare la curvatura e `screenGap`;
   soltanto allora applicare il confronto di Sturm alla geometria.
4. Collegare il risultato al differenziale della mappa di arrivo;
   poi provare properness e conteggio delle fibre.

Papers I/II e lettera CQG restano invariati. Questo sviluppo è per Paper III.


Aggiornamento successivo: il [blocco 2i](BLOCK2_MIXED_DERIVATIVES_AND_ACTUAL_JACOBI.md)
prova gli scambi delle derivate per famiglie C³ fornite e conclude l’equazione
del campo effettivo. La costruzione regolare della famiglia dei lanci resta aperta.
