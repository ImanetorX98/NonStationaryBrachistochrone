# Via B: formalizzazione effettivamente compilata in Lean

6 ottobre 2026, aggiornata con il blocco 2d: confinamento esterno dedotto a priori, permanenza non radiale dall’angolo iniziale, alternativa cattura/prolungamento e conservazione del dato al lancio. **160 teoremi compilati.** Sorgenti di verifica per Paper III, ancora in sviluppo. **Certificazione parziale: non è ancora una formalizzazione del teorema geometrico di esclusione dei Maxwell prima di π.**

## Ambiente e verifica

Lean **4.24.0**, mathlib **v4.24.0**, commit mathlib `f897ebcf72cd16f89ab4577d0c826cd14afaafc7`. Il manifest conserva anche le revisioni delle dipendenze transitive.

Il progetto è nella directory:

```text
NonStationaryMetrics/paper3/vaidya_no_maxwell_before_pi
```

Da questa directory, con Lean e Lake della versione indicata disponibili:

```sh
python3 verify.py
```

Se Lean e Lake non sono nel PATH, indicare la directory dei loro eseguibili:

```sh
python3 verify.py --lean-bin /path/to/lean-4.24.0/bin
```

Predisporre Lean 4.24.0 (per esempio tramite elan), quindi eseguire `lake exe cache get` per scaricare le dipendenze pubbliche e la cache mathlib. Conservare il manifest versionato per usare le revisioni fissate. `verify.py` non aggiorna le dipendenze. `.lake/` contiene dipendenze e artefatti rigenerabili e non viene versionata.

La verifica esegue `lake build` e `lake env lean AxiomAudit.lean`, controlla la copertura di tutti i teoremi pubblici e ammette soltanto le dipendenze assiomatiche standard `propext`, `Classical.choice`, `Quot.sound`. Salva log e hash SHA-256 nella sottodirectory `verification`. Questi controlli **non eliminano le ipotesi esplicite dei teoremi**: occorre leggerle e verificare l'applicazione geometrica.

## Contenuto certificato

| Modulo | Teoremi | Risultato effettivo |
|---|---:|---|
| `ViaB.Algebra` | 4 | Positività del polinomio P e della formula scalare Δ; segni limite del coefficiente B. |
| `ViaB.Wronskian` | 4 | Derivata esatta di W, derivata di z/u, crescita di W su un intervallo di confronto positivo. |
| `ViaB.FirstZero` | 2 | Segno della derivata al primo zero; contraddizione per un candidato primo zero. |
| `ViaB.Sturm` | 4 | Positività iniziale da z′(a)>0; costruzione del primo zero; confronto Sturm completo scalare; positività all'estremo dove u si annulla. |
| `ViaB.Trace` | 3 | Non degenerazione della tangente; costanza locale del numero di fogli di un rivestimento; iniettività da una fibra singola su base connessa. |
| `ViaB.NegativeControls` | 3 | Controlli esatti sotto soglia e sul caso radiale. |
| `ViaB.FarField` | 7 | Continuità, margine negativo uniforme di B, margine radiale uscente e limite superiore r_v≤q. |
| `ViaB.Dynamics` | 5 | Disuguaglianze differenziali, barriere accoppiate e stima logaritmica delle escursioni. |
| `ViaB.ConeEscape` | 5 | Cono uscente invariante, crescita radiale almeno lineare e fuga dopo ingresso stretto nel cono. |
| `ViaB.LingeringLimits` | 4 | Esclusione del limite radiale entrante e del limite α→0, w→q. |
| `ViaB.Excursion` | 2 | Direzione iniziale di un attraversamento uscente; ingresso angolare nel cono durante una grande escursione. |
| `ViaB.LastExit` | 4 | Ultima uscita, velocità non negativa all'uscita e ingresso stretto nel cono. |
| `ViaB.UnboundedEscape` | 2 | Scelta del cono e fuga permanente da raggio illimitato. |
| `ViaB.Barbalat` | 3 | Convergenza della primitiva monotona limitata e decadimento da derivata limitata. |
| `ViaB.BoundedAngle` | 3 | Stima esplicita di B e decadimento di sin α dall'angolo accumulato limitato. |
| `ViaB.AngleLimits` | 2 | Dicotomia α→0 oppure α→π; esclusione del ramo entrante. |
| `ViaB.RadialConvergence` | 4 | Stime della velocità, convergenza di r mediante r+(q/c)φ e convergenza della massa monotona. |
| `ViaB.BoundedLingering` | 4 | Chiusura della contraddizione per angolo finito e divergenza angolare nel caso esterno limitato. |
| `ViaB.EscapeStability` | 3 | Ingresso stretto in un cono fissato; stabilità e apertura della fuga per famiglie con valutazioni continue. |
| `ViaB.FlowDependence` | 7 | Regolarità e Lipschitz locale del campo, stima di Grönwall e dipendenza finita per traiettorie esistenti in un tubo comune. |
| `ViaB.ExteriorDichotomy` | 2 | Angolo limitato implica fuga; alternativa r→∞ oppure φ→∞ nel caso esterno positivo monotono. |
| `ViaB.TubeBootstrap` | 3 | Primo attraversamento escluso dalla stima di Grönwall; continuità e stabilità della fuga da un limite locale del campo, senza assumere permanenza dei vicini nel tubo. |
| `ViaB.CompactFieldBound` | 8 | Uniformità Lipschitz sul riferimento finito regolare; controllo del tubo, dipendenza finita e stabilità della fuga con costanti derivate, apertura dalla continuità iniziale. |
| `ViaB.MetricHamiltonian` | 9 | Metrica coordinata e inversa; ramo nullo futuro e sua unicità. |
| `ViaB.NullHamiltonianDynamics` | 14 | Derivate Hamiltoniane, riduzione del momento radiale e velocità normalizzate. |
| `ViaB.DirectionFromHamiltonian` | 5 | Ricostruzione di α∈(0,π) e derivazione della sua ODE dai momenti. |
| `ViaB.AffineReparametrization` | 1 | Sistema di direzione dalle equazioni affini, con derivata della riparametrizzazione fornita. |
| `ViaB.LocalTimeInverse` | 4 | Inversa temporale costruita da velocità continua positiva; ODE locali con inversa costruita. |
| `ViaB.TimeDependentHamiltonian` | 9 | Derivata temporale di H, cancellazione completa di dH/dλ, propagazione della nullità e conservazione di L. |
| `ViaB.LocalDirectionExistence` | 5 | Esistenza locale con massa C¹, unicità con massa continua e permanenza iniziale nel dominio esterno. |
| `ViaB.FiniteEndpointLimit` | 4 | Estensione Lipschitz ausiliaria, limite finale da velocità limitata e nuovo tratto al limite regolare. |
| `ViaB.RegularEndpointContinuation` | 3 | Derivata sinistra al bordo, coincidenza all’indietro e prolungamento effettivo dell’ODE sui compatti regolari. |
| `ViaB.ExteriorFiniteTube` | 7 | Limiti finiti di r e m, compatto regolare derivato dall’evoluzione esterna e continuazione. |
| `ViaB.AngularFiniteInvariance` | 11 | Striscia non radiale dal dato iniziale; controllo del bordo finale, cattura trasversale, alternativa esterna e raccordo al lancio. |
| **Totale** | **160** | Tutti compilati; nessun `sorry`, `admit` o assioma geometrico aggiunto. |

### Algebra delle curvature

Ponendo q=Ê²:

\[
P=3(2q-3)(q-1)r^2+(33q-38)mr+40m^2.
\]

Lean dimostra P>0 per q≥3/2, m>0, r>0. Dimostra inoltre la positività della funzione definita da

\[
\Delta=\frac{L^2}{r^4(qr+2m-r)^2}
\left[\frac{mP}{r}+m' r(3(q-1)r+4m)\right]
\]

per m′≥0 e L≠0. **L'identificazione di questa funzione con K⊥−K∥ della metrica non è formalizzata.**

Per

\[
B(q,w,c)=c(2w-q+1)-\frac12\sqrt{q/w}(3w-q+1)
\]

sono certificati B(q,q,1)=1/2 per q>0 e B(q,q−1,c)<0 per q>1, c≤1. Il limite inferiore c≥−1 non è necessario per questa seconda disuguaglianza.

### Confronto Sturm

Le equazioni sono scritte in **un unico parametro affine**:

\[
u''+K_\perp u=0,\qquad z''+K_\parallel z=0,
\qquad W=z'u-zu'.
\]

Il teorema `ViaB.sturm_positive` dimostra:

\[
\begin{gathered}
u(a)=z(a)=0,\quad z'(a)>0,\quad
u>0\text{ in }(a,b],\quad K_\perp-K_\parallel>0\text{ in }(a,b)
\\\Longrightarrow z>0\text{ in }(a,b].
\end{gathered}
\]

Sono ipotesi esplicite le quattro relazioni `HasDerivAt` per u, u′, z, z′ su [a,b]. Non si assume z>0 sull'intervallo per ottenere la conclusione. La prova ricava la positività vicino al lancio dalla derivata iniziale, costruisce un eventuale primo zero mediante compattezza e teorema dei valori intermedi, quindi lo esclude con W′=(K⊥−K∥)uz.

`ViaB.sturm_positive_at_reference_zero` sostituisce u>0 in (a,b] con u>0 in (a,b), u(b)=0, u′(b)<0 e conclude z(b)>0. È il lemma scalare utile all'antipodo.

Per applicarlo al problema occorre ancora formalizzare

\[
u=\frac{r_0}{L}r\sin\varphi,\qquad
u'(\lambda_\pi)=-\frac{r_0}{r(\lambda_\pi)},
\]

le ODE dello schermo, i dati iniziali e la relazione tra zero dello schermo e punto coniugato. Il parametro b di Lean non è un angolo: nell'applicazione sarà l'istante affine del bersaglio o dell'antipodo.

### Traccia e fogli

Con A=rψ, B=rλ, C=φψ, D=φλ, il vettore tangente (B,−A) annulla dr, mentre

\[
d\varphi(B,-A)=CB-DA=-(AD-BC)\ne0.
\]

Questo calcolo è certificato anche quando rλ=0. Non contiene ancora la costruzione della varietà T né l'implicazione geometrica da non coniugazione a determinante non nullo.

Per un vero rivestimento f:E→X, Lean dimostra che `Nat.card` della fibra è localmente costante. Su una base preconnessa, se una fibra ha cardinalità 1 allora f è iniettivo. `Nat.card` assegna 0 alle fibre infinite: questo non crea ambiguità nel teorema di iniettività perché la cardinalità propagata è **1**.

Il teorema assume `IsCoveringMap f` e l'esistenza di una fibra singola; **non assume direttamente l'iniettività**, ma queste due premesse devono ancora essere dimostrate per la traccia di Vaidya.

### Controlli negativi

Per q=5/4, m=1, r=20, Lean verifica esattamente P=−45 e Δ<0 con m′=0, L=1. Per L=0 verifica Δ=0. Questi controlli impediscono di allargare automaticamente i lemmi a ogni q>1 o di pretendere stretta positività sui radiali. **Non sono controesempi al teorema geometrico di Maxwell, né provano che q≥3/2 sia una soglia geometrica ottimale.**

## Cosa rimane per il teorema geometrico

1. Completare il collegamento geometrico alla metrica G₃, geodetiche nulle future, parametro affine e schermatura; derivare ODE di Jacobi e formula Δ. I blocchi 2a–2b certificano l’inversione della metrica coordinata, l’equazione di direzione dal suo Hamiltoniano e l’esistenza dell’inversa temporale locale. Nullità e momento angolare sono propagati dal flusso Hamiltoniano completo dichiarato. Il collegamento alla connessione di Levi-Civita resta da formalizzare.
2. Collegare esistenza, continuazione e intervallo angolare alla geometria. Il blocco 2c costruisce ora soluzioni locali per massa C¹ e dimostra continuazione con incollamento sui compatti regolari. Il blocco 2d deduce ora il confinamento dall’evoluzione esterna e la striscia angolare dal dato iniziale; restano costruzione degli intervalli massimali e permanenza futura nel dominio esterno. Il blocco 1 è chiuso: regolarità del campo e compattezza del riferimento producono la stima locale uniforme, senza assumerla. Grönwall controlla il tubo e trasferisce la continuità iniziale alle valutazioni finite e alla stabilità della fuga. `DirectionFlow` assume ancora soluzioni future globali; il lemma fondamentale di uniformità richiede soltanto un riferimento continuo regolare su un tratto finito, senza ODE o esistenza futura.
3. Formalizzare la chiusura della traccia e il controllo degli arrivi tardivi: tempi v e affini, cattura trasversale, limite radiale critico, separazione del ricevitore dal bordo. Da qui ottenere properness.
4. Costruire la traccia come varietà, dimostrare che la mappa angolare è un diffeomorfismo locale, poi che properness implica rivestimento. Dimostrare il conteggio vicino a δ=0 per rT<r0, rT=r0 e rT>r0, compreso il caso vuoto.
5. Costruire il congelamento futuro liscio e provare la località dell'ODE, per rimuovere le ipotesi di massa limitata e ricevitore eternamente esterno.
6. Collegare l'unicità nel rivestimento al problema tridimensionale di primo arrivo, alla riduzione planare, a V_t e all'esclusione dei Maxwell. Formalizzare separatamente l'affermazione sul primo taglio all'antipodo.

**Prossimo obiettivo consigliato:** collegare il sistema alla geometria e al suo intervallo di esistenza, poi chiudere il controllo degli arrivi tardivi e la properness. Non chiamare il teorema geometrico “verificato in Lean” prima di avere chiuso tutti i collegamenti.

## Nuova chiusura dei lemmi ODE

`unbounded_radius_escapes` assume soltanto q>1, M>0, `DirectionFlow q M a` e l'esistenza di raggi arbitrariamente grandi in tempi futuri; conclude r→∞. La monotonicità della massa non serve per questo risultato.

`bounded_exterior_flow_angle_diverges` assume q>1, `DirectionFlow q M a`, m(a)>0, massa non decrescente, 2m≤r≤R e l'equazione φ′=√w sin α/r, dove w=q−1+2m/r; conclude φ→∞. Non richiede un limite radiale tra le ipotesi.

Il passaggio nuovo alla convergenza radiale usa c=√(q−1)/R>0. Se φ fosse limitato, i lemmi di Barbalat darebbero sin α→0 e quindi α→0. Per α≤π/2 valgono r′≥−q sin α e φ′≥c sin α, perciò Y=r+(q/c)φ è non decrescente e limitata. Convergono Y, φ e dunque r. La massa monotona converge; l'ODE radiale forza w→q, mentre il lemma angolare esclude simultaneamente α→0 e w→q.

## Stabilità e dipendenza dai lanci

Il campo `directionField q m (r,α)` è formalmente C¹ e localmente Lipschitz, congiuntamente nello stato e in m, nei punti r>0, m≥0, q>1. Non viene affermato un limite Lipschitz globale fino a r=0.

`trajectory_stays_in_tube` richiede il limite ‖f′−g′‖≤K‖f−g‖ soltanto dove dist(f,g)≤ε, con K≥0. Se dist(f(a),g(a)) exp(K(b−a))<ε, dimostra dist(f(t),g(t))<ε per ogni t∈[a,b]. Il controllo del tubo è quindi una conclusione, non una premessa.

`cone_escape_stable_of_local_field_bound` è la versione precedente con il limite locale uniforme come premessa. `cone_escape_stable_from_regular_reference` ora deriva quel limite internamente: per famiglie di `DirectionFlow` con q fisso, la stessa massa sul tratto [a,s], massa di riferimento continua e dato iniziale continuo, l'ingresso stretto del riferimento implica fuga dei lanci vicini e una crescita radiale comune δ(t−s). Rimane premessa l'esistenza futura delle soluzioni, non la costante Lipschitz né la permanenza dei vicini nel tubo.

`escaping_set_isOpen` è il corollario topologico per famiglie con valutazioni continue a ogni tempo futuro. Non dimostra automaticamente l'apertura nello spazio di tutti i dati geometrici, comprendendo quelli che vengono catturati o cessano di esistere nel dominio.

`escaping_set_isOpen_from_initial` elimina la premessa di continuità a tutti i tempi finiti e la ricava dalla continuità iniziale, per famiglie future globali con q fisso e una stessa storia di massa continua. L'apertura resta relativa a tale famiglia dichiarata.

## Blocco 1: uniformità sul riferimento compatto

`compact_uniform_local_lipschitz` ricopre il compatto con palle di raggio ρᵢ/3, seleziona un sottoricoprimento finito, prende K come massimo delle costanti locali e uno spessore ε>0 minore o uguale a tutti i raggi ρᵢ/3. La disuguaglianza triangolare mantiene il riferimento e lo stato vicino nella stessa palla di raggio ρᵢ. Si conclude una stima tra ciascun punto del compatto e i punti vicini; non si afferma un limite Lipschitz tra tutte le coppie arbitrarie dell'intero tubo.

`finite_reference_uniform_field_bound` applica questa costruzione al riferimento aumentato ((r,α),m) su [a,b], con q>1, continuità, r>0 e m≥0. Non richiede ODE, massa monotona, limiti già noti sul raggio, altre traiettorie o esistenza futura. La continuità di m è esplicita e non è contenuta nel vecchio `DirectionFlow`.

`directionFlow_finite_time_dependence` deriva K e ε dal riferimento e dimostra che ogni altra soluzione esistente con la stessa massa sul tratto e scarto iniziale d₀ exp(K(b−a))<ε resta nel tubo; vale inoltre la stima d(t)≤d₀ exp(K(t−a)). I corollari trasferiscono la continuità iniziale alla valutazione finale e alla fuga stabile. Il blocco non dimostra esistenza o continuazione delle soluzioni.

## Blocco 2a: derivazione coordinata del sistema di direzione

I 29 nuovi teoremi sono descritti in [BLOCK2_COORDINATE_HAMILTONIAN.md](BLOCK2_COORDINATE_HAMILTONIAN.md). La nullità futura seleziona il ramo dei momenti e consente di ricostruire α tramite arccos. La regola della catena deduce l’ODE angolare senza assumerla e senza dividere per il momento radiale o per cos α. La massa è il valore istantaneo m(v); non si assume la conservazione di P_v.

Il teorema composto del blocco 2a prende le equazioni Hamiltoniane affini, L>0, la nullità all’evento e la derivata di una riparametrizzazione; conclude le ODE in v. Il blocco 2b successivo costruisce l’inversa locale e dimostra la conservazione di L e la propagazione della nullità, per soluzioni Hamiltoniane esistenti. Il blocco 2c successivo costruisce soluzioni locali e dimostra continuazione sui compatti regolari. L’esistenza futura globale e i collegamenti di Jacobi restano aperti: il blocco 2 completo e il teorema geometrico non sono ancora certificati.

## Blocco 2b: inversa locale e vincolo nullo

[BLOCK2_LOCAL_TIME_AND_NULL_PROPAGATION.md](BLOCK2_LOCAL_TIME_AND_NULL_PROPAGATION.md) documenta i 13 nuovi teoremi. Da una velocità temporale continua positiva e dalla sua equazione in un intorno si costruisce una vera inversa locale, con entrambe le identità di inversa. `hamiltonian_local_direction_exists` conclude le ODE all’evento senza assumere la riparametrizzazione.

Il flusso completo include P_{v,λ}=−m_v D_mH. La regola della catena verifica dH/dλ=0 anche per massa variabile; il vincolo nullo si propaga dal dato iniziale su un dominio convesso. L_λ=0 implica la conservazione del momento angolare. Queste conclusioni del blocco 2b sono seguite dalla costruzione di soluzioni locali e dalla continuazione sui compatti nel blocco 2c. Resta da eliminare la premessa di esistenza futura globale dei vecchi `DirectionFlow`.

## Blocco 2c: esistenza locale e continuazione

[BLOCK2_LOCAL_EXISTENCE_AND_CONTINUATION.md](BLOCK2_LOCAL_EXISTENCE_AND_CONTINUATION.md) descrive i 12 nuovi teoremi. Picard–Lindelöf costruisce una soluzione locale per massa C¹; l’unicità locale richiede solo massa continua. Per dati strettamente esterni e non radiali, si costruisce un tratto che resta nel dominio esterno.

`compact_direction_continues` dimostra che una soluzione confinata in un compatto regolare dello spazio stato–massa continua oltre un estremo finito, se la massa è C¹ vicino all’estremo. Il limite finale, il nuovo tratto e la loro coincidenza prima dell’estremo vengono tutti dimostrati. Il blocco 2d successivo deduce il confinamento dalle ipotesi esterne e prova l’alternativa finita cattura/prolungamento esterno; restano gli intervalli massimali e il flusso futuro globale. La continuazione del campo coordinato può attraversare r=2m: non certifica da sola permanenza nel dominio esterno.

## Blocco 2d: stime esterne e alternativa finale

[BLOCK2_EXTERIOR_APRIORI_AND_CAPTURE.md](BLOCK2_EXTERIOR_APRIORI_AND_CAPTURE.md) registra 18 nuovi teoremi. Per massa positiva al lancio e non decrescente sul tratto finito, il dominio esterno e l’ODE danno 2m(a)≤r(v)≤r(a)+q(b−a) e 0≤m(v)≤m(b). Il compatto di continuazione viene quindi costruito. La striscia 0<α<π è dedotta dall’angolo iniziale mediante l’ODE lineare di sin α e unicità all’indietro, senza assumerla lungo il tratto.

All’estremo finito il prolungamento ha r(b)≥2m(b) e angolo strettamente non radiale. O arriva al bordo di cattura, oppure continua come soluzione esterna non radiale. Al bordo, con m_v≥0, la derivata di r−2m è strettamente negativa. `exterior_nonradial_initial_value_continuation` conserva anche dato e derivata al lancio. Questi risultati partono da una soluzione esterna esistente sul tratto precedente; resta da costruire la soluzione massimale e ricavare un `DirectionFlow` futuro globale senza assumerlo.
