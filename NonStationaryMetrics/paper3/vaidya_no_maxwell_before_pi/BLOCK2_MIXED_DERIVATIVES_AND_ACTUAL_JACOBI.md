# Blocco 2i — derivate miste e campo di variazione effettivo

6 ottobre 2026. Nuovi moduli: `MixedGeodesicDerivatives` (12 teoremi)
e `FermatC3JacobiFamily` (2 teoremi). Checkpoint di 248 teoremi pubblici
in 51 sorgenti Lean, compresi modulo principale e audit.

## Risultato preciso

Per una **famiglia fornita con regolarità congiunta C³**, le derivate
miste richieste dal blocco 2h sono ora dimostrate mediante Schwarz.
Il campo J è definito come derivata rispetto al lancio, e J′, J″ sono
le sue derivate affini effettive. Una famiglia di geodetiche nulle
coordinate di Fermat soddisfa quindi l'equazione di variazione, senza
assumere separatamente le identificazioni delle derivate miste.

Non è ancora costruita la famiglia C³ dei lanci del problema. L'ipotesi
C³ è una condizione sufficiente di questo blocco, non una conclusione
dall'esistenza delle singole curve né dalla dipendenza continua da Grönwall.
Il coefficiente di curvatura dello schermo e `screenGap` non sono
identificati; il punto 2 completo rimane aperto.

## 1. Parametri e derivate effettive

Scriviamo X(ε,λ) per una componente coordinata della famiglia; ε è
il lancio, λ un unico parametro affine di famiglia. Per una funzione
F:ℝ²→ℝ definiamo

\[
\mathcal P_\varepsilon F=DF(1,0),\qquad
\mathcal P_\lambda F=DF(0,1).
\]

In Lean sono `launchPartial` e `affinePartial`, definite tramite
`fderiv`, non tramite coefficienti assegnati. I lemmi sulle sezioni
verificano le `HasDerivAt` delle funzioni λ↦F(a,λ) e ε↦F(ε,s).
`familyPartial_contDiffAt` prova la perdita di un ordine di regolarità.

`familyPartial_hasFDerivAt` deriva una parziale dalla seconda derivata
Fréchet. `familyPartial_commute` usa il teorema mathlib
`ContDiffAt.isSymmSndFDerivAt`: per F C²,

\[
\mathcal P_\varepsilon\mathcal P_\lambda F
 =\mathcal P_\lambda\mathcal P_\varepsilon F.
\]

Non si aggiunge un assioma di commutazione.

## 2. Velocità, accelerazione e campo J

Per J(λ)=∂εX(a,λ), il primo passaggio verifica

\[
\partial_\varepsilon(\partial_\lambda X)(a,s)=J'(s)
\]

con regolarità C². Per F C³, si applica Schwarz alla funzione
∂λF, che è C², e si differenzia l'uguaglianza del primo scambio
valida su un intorno. Segue

\[
\partial_\varepsilon(\partial_\lambda^2 X)(a,s)=J''(s).
\]

`acceleration_launch_variation_hasDerivAt` certifica il passaggio.
`launch_field_affine_hasDerivAt` e
`launch_field_affine_second_hasDerivAt` mostrano che questi valori
sono rispettivamente la derivata di J e la derivata della sua derivata
ordinaria. Nel secondo lemma l'uguaglianza delle funzioni su un intorno
è stabilita prima di differenziare, evitando di derivare una mera
uguaglianza di valori in un punto.

`affine_slice_second_hasDerivAt` verifica analogamente la seconda
derivata ordinaria di ogni componente della curva geodetica.

## 3. Equazione ottenuta

`C3_geodesic_family_coordinate_jacobi` combina gli scambi dimostrati
con la linearizzazione del blocco 2h. Per u=∂λX e
DG=∂εΓ(X(ε,s)), conclude

\[
J^{i\prime\prime}
 +2\sum_{j,k}\Gamma^i{}_{jk}u^jJ^{k\prime}
 +\sum_{j,k}DG^i{}_{jk}u^ju^k=0.
\]

Qui J′ e J″ sono espresse nella conclusione mediante `deriv` e
`deriv (deriv ...)` del campo effettivo. Le premesse non includono
più le identificazioni ∂εu=J′ e ∂εb=J″.

`fermat_C3_family_actual_coordinate_jacobi` applica la formula alla
connessione di Fermat verificata. Assume una massa fisica C² all'evento,
q,r,w non nulli e la famiglia congiuntamente C³. La derivata di Γ
viene ricavata dal teorema di massa C² del blocco precedente.

`fermat_C3_null_geodesic_family_coordinate_jacobi` parte direttamente
da `CoordinateNullGeodesicAt` per i membri vicini della famiglia.
Deduce le equazioni in termini delle parziali mediante le seconde
derivate effettive delle sezioni, poi applica il risultato precedente.
Nullità e orientazione futura sono contenute nella premessa delle curve;
la linearizzazione in sé richiede soltanto la loro equazione geodetica.

È l'equazione di variazione in coordinate. Restano da ricomporla nella
forma covariante, dichiarare la convenzione di Riemann, proiettarla sullo
schermo e identificare le curvature. DG non è dichiarato uguale a una
curvatura, né a `screenGap`.

## 4. Dato iniziale di posizione

`fixed_launch_event_variation_zero` parte dall'uguaglianza
X(ε,λ₀)=x₀ per i lanci vicini e dalla differenziabilità. Confrontando
la derivata della sezione con quella di una costante conclude

\[
J(\lambda_0)=0.
\]

Il dato non viene assunto direttamente per il campo. Rimangono da
identificare la variazione delle velocità iniziali, la normalizzazione
affine e il dato iniziale dello schermo usato dal confronto di Sturm.

## 5. Regolarità ancora da costruire

Questo blocco **non** prova che massa C² implichi una famiglia dei lanci
congiuntamente C³. Le due ipotesi sono distinte. Le singole geodetiche
locali precedenti e la dipendenza continua non bastano a verificarle.
Non è rivendicata minimalità dell'ipotesi C³: una prova direttamente
variazionale dell'ODE, con regolarità differenziata nei due parametri,
potrebbe consentire ipotesi più deboli.

Il prossimo obbligo è provare la dipendenza differenziabile del flusso
rispetto ai dati, controllare un intervallo comune e trasferire tale
regolarità all'orologio affine e alla sua inversa. Bisogna quindi scegliere
se costruire una famiglia C³ sotto ipotesi sufficientemente forti, oppure
usare le derivate miste necessarie senza pretendere C³ in ogni direzione.

Il chart (v,r,φ) tratta la famiglia planare. Le variazioni fuori piano
e la loro relazione con il campo scalare di riferimento vanno costruite
nella geometria estesa; la sola traslazione di φ non le fornisce.

## 6. Verifica e stato complessivo

L'audit completo con `verify.py` è passato su tutti i 248 teoremi:
compilazione, copertura dell'audit e controllo degli assiomi transitivi. Log e hash sono conservati
in `verification/`. Le dipendenze ammesse sono soltanto `propext`,
`Classical.choice`, `Quot.sound`; nessuna prova viene lasciata incompleta.

Lo scambio delle derivate è chiuso **per famiglie C³ fornite**.
La regolarità della famiglia costruita, lo schermo, l'identificazione
della curvatura e l'applicazione geometrica di Sturm restano aperti.
Properness e conteggio delle fibre appartengono ai passaggi successivi.
Papers I/II e lettera CQG non sono modificati.
