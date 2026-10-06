# Blocco 2e — soluzione massimale esterna ed esistenza globale nel cono

6 ottobre 2026. Il checkpoint aggiunge **19 teoremi**, portando il progetto a
**179 teoremi pubblici e 40 sorgenti Lean**, incluso l’audit assiomatico.
I risultati riguardano il sistema coordinato di direzione di Vaidya.
Non certificano ancora il teorema geometrico completo di esclusione dei Maxwell.

## 1. Problema e ipotesi

Con q = Ê² > 1, f = (r, α), w = q − 1 + 2m(v)/r,

\[
r_v=-w+\sqrt{qw}\cos\alpha,\qquad
\alpha_v=\frac{\sin\alpha}{r}B(q,w,\cos\alpha).
\]

Il dominio fisico qui trattato è m > 0, r > 2m e 0 < α < π.
I casi radiali sono esclusi. La massa è prescritta, senza una conservazione
di P_v né una sostituzione con massa congelata.

L’unicità richiede massa continua sui domini regolari. L’esistenza locale
richiede C¹ vicino al lancio. Per classificare un estremo finito b occorrono
m C¹ vicino a b e massa non decrescente da a in avanti. Per la fuga nel cono
si aggiunge un limite superiore m(v) ≤ M per v ≥ a.

Il teorema finale impacchettato assume m globalmente C¹: è un’ipotesi
sufficiente esplicita, più forte della regolarità strettamente necessaria.
I lemmi intermedi richiedono solo continuità sul dominio unito e C¹ ai
tempi futuri necessari. La soglia q ≥ 3/2 del confronto di curvatura non
viene modificata né dedotta dal requisito q > 1 di questi lemmi.

## 2. Unicità sull’intero intervallo comune

`direction_solution_unique_on_open_connected` considera due soluzioni
regolari sullo stesso aperto connesso, con uguale stato in un punto.
Sul sottospazio dei tempi ammessi, l’insieme dei punti di coincidenza è:

- chiuso, perché entrambe le soluzioni sono continue;
- aperto, perché l’unicità locale produce coincidenza in un intorno di
  ogni punto di coincidenza;
- non vuoto, perché contiene il lancio.

La connessione implica coincidenza ovunque. Non si assume una costante
Lipschitz uniforme sull’intero intervallo. Il corollario per intervalli
aperti viene usato per incollare soluzioni con estremi futuri diversi.

## 3. Costruzione massimale effettiva

`ExteriorSegment q m a l x` contiene un tempo finale b > a e una curva
sull’intervallo (l,b), con l < a fissato, stato iniziale x, ODE e dominio
fisico. `exterior_segment_exists` costruisce almeno un segmento mediante
Picard–Lindelöf; sceglie l’interno del tratto locale come estremo sinistro.

Si definisce

\[
U=\bigcup_{s\text{ ammissibile}}(l,b_s).
\]

I segmenti concordano su tutto il loro intervallo comune. La funzione
`maximalExteriorCurve` sceglie un segmento in ogni punto di U: il valore
non dipende da questa scelta, per l’unicità appena dimostrata.
Vicino a ciascun punto coincide con un segmento fissato, quindi soddisfa
l’ODE e resta nel dominio fisico. Il dato iniziale è conservato.

L’aperto U è connesso, perché tutti gli intervalli contengono a.
Se l’insieme degli estremi ammissibili è limitato superiormente,

\[
U=(l,b_*),\qquad b_*=\sup\{b_s:s\text{ ammissibile}\}>a.
\]

Se è illimitato, U = (l,+∞). Ogni segmento ammissibile è contenuto in U:
questa è la proprietà di massimalità effettivamente formalizzata.
Non si assume l’esistenza di una curva massimale per dedurne l’esistenza.

## 4. Un estremo finito massimale è cattura

`maximal_exterior_finite_endpoint_is_capture` applica alla curva costruita
l’alternativa finita già certificata nel blocco 2d. Il prolungamento
coordinato h esiste attorno all’estremo b_*.

Se h fosse ancora strettamente esterno, l’alternativa fornirebbe un
δ > 0 con permanenza esterna non radiale fino a b_* + δ. Si incolla h
alla curva precedente conservando anche la parte prima del lancio e
la derivata al lancio. La monotonia mantiene positiva la massa sul nuovo
tratto. Si ottiene un nuovo `ExteriorSegment` con estremo b_* + δ,
in contraddizione con la definizione di b_*.

Pertanto

\[
r(b_*)=2m(b_*).
\]

Questo è un endpoint della soluzione esterna. Il campo coordinato può
continuare oltre il bordo; la sua continuazione non equivale a permanenza
nel dominio esterno. Il lemma di cattura trasversale del blocco 2d resta
disponibile, ma non è necessario alla contraddizione per massimalità.

## 5. Il cono elimina l’estremo finito senza circolarità

Per 0 < β < π/2 e q − 1 < q cos²β, il lemma già verificato
`outgoing_cone_margins` costruisce R > 2M, κ > 0 e δ > 0 tali che,
per r ≥ R e 0 ≤ m ≤ M,

\[
B(q,w,\cos\alpha)\le-\kappa,\qquad
r_v\ge\delta\quad\text{quando }\alpha\le\beta.
\]

La prima stima vale per ogni cos α ≤ 1. Il cono stretto è r > R, α < β.
`outgoing_cone_invariant_on_finite_interval` dimostra l’invarianza sui
soli intervalli finiti già esistenti. Non usa un `DirectionFlow` futuro
globale: le barriere accoppiate sono r − R e β − α.

Supponiamo ingresso stretto al tempo c ≥ a, con c ∈ U. Se U avesse
estremo finito b_*, la barriera darebbe r(v) > R per c ≤ v < b_*.
La continuità del prolungamento finale implica r(b_*) ≥ R.
La cattura richiederebbe invece r(b_*) = 2m(b_*) ≤ 2M < R.
Contraddizione. Gli estremi ammissibili sono dunque illimitati.

L’ordine delle deduzioni è:

1. esistenza locale e unicità;
2. costruzione della curva massimale;
3. cattura obbligatoria se l’estremo è finito;
4. barriera del cono sui tratti finiti;
5. impossibilità della cattura dopo ingresso nel cono;
6. esistenza futura globale e costruzione del `DirectionFlow`;
7. applicazione dei precedenti lemmi di crescita e fuga.

Non si applica un teorema di fuga per flussi futuri globali prima di
aver dimostrato l’esistenza di quel flusso.

## 6. Risultati finali utilizzabili

`directionFlow_exists_and_escapes_after_cone_entry` costruisce un
`DirectionFlow` con il dato iniziale originale, massa prescritta e
permanenza strettamente esterna per tutti i tempi futuri. Inoltre,

\[
r(v)-r(c)\ge\delta(v-c)\quad(v\ge c),\qquad r(v)\longrightarrow+\infty.
\]

`directionFlow_exists_from_strict_cone_launch` elimina anche le premesse
intermedie del segmento locale e dei margini del cono. Da m C¹,
m(a) > 0, monotonia futura e m ≤ M, costruisce una soglia R > 2M:
ogni lancio r(a) > R, 0 < α(a) < β produce un flusso futuro globale
esterno e in fuga. Non presuppone alcuna traiettoria già esistente.

## 7. Controlli e limiti per il seguito

`verify.py` compila tutti i moduli, controlla la copertura dei 179 teoremi
nell’audit assiomatico, esclude dichiarazioni non dimostrate e registra
i rispettivi hash. Sono consentiti solo gli assiomi standard di Lean/mathlib:
`propext`, `Classical.choice`, `Quot.sound`.

La revisione della catena ha verificato anche questi punti:

- il supremo viene usato solo nel caso limitato e non vuoto;
- l’estremo sinistro scelto localmente è conservato dal nuovo segmento
  nella contraddizione per massimalità;
- la prova del cono tratta intervalli finiti, non un’esistenza globale
  nascosta nelle sue premesse;
- i margini del cono sono derivati, non introdotti come nuove ipotesi
  nel teorema finale per lanci stretti;
- il confronto con il bordo usa m(b_*) ≤ M anche all’estremo finale.

Il successivo [blocco 2f](BLOCK2_DIRECTION_TO_AFFINE_HAMILTONIAN.md)
ricostruisce dalle soluzioni coordinate il flusso Hamiltoniano nullo affine
locale. Restano l’identificazione tramite connessione e campi di Jacobi,
la properness della mappa di arrivo, il rivestimento/conteggio delle fibre
e il passaggio conclusivo all’esclusione dei Maxwell. Questi 19 teoremi
non certificano tali collegamenti e non provano che ogni lancio esterno
entri nel cono o sfugga. Papers I/II e lettera CQG non sono modificati.
