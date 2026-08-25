# Audit delle citazioni e delle fondamenta teoriche di Paper I e Paper II

**Data:** 25 agosto 2026  
**Autore dell'audit:** Codex  
**Destinatario operativo:** Claude, per l'ultimo passaggio prima della resubmission CQG  
**Intervento sui manoscritti:** nessuno; questo file contiene soltanto valutazioni e modifiche proposte.

## 1. Ambito e versioni controllate

Ho considerato come versioni correnti:

- Paper I: `paper1/paper1_JMP.tex` e `paper1/paper1_JMP.pdf`;
- Paper II: `paper2/paper2.tex` e `paper2/paper2.pdf`;
- risposta alla major review: `paper2/response_to_referees_CQG.tex`;
- bibliografia condivisa: `paper/refs.bib`;
- lettera del referee e documenti di lavorazione già presenti nel progetto;
- i PDF scientifici locali in `Fonti/`, con lettura puntuale dei lavori che sorreggono le affermazioni centrali.

L'inventario completo dei PDF è in `INVENTARIO_BIBLIOGRAFIA_LOCALE_codex.md`. La presente analisi non pretende che ciascuno degli 83 PDF locali abbia lo stesso peso: ho eseguito una lettura semantica approfondita delle fonti da cui dipendono originalità, formulazione variazionale, controllo ottimo, sommersioni, minimalità, correzione adiabatica e classe di funzioni speciali; sulle voci periferiche ho verificato contesto e metadati. Quando manca l'opera esatta citata, lo dichiaro esplicitamente.

## 2. Verdetto esecutivo

Il lavoro è **edificato in modo sostanzialmente corretto sulla letteratura precedente**. La struttura scientifica centrale regge:

1. Perlick fornisce il precedente stazionario a energia fissata, ma non il rail attivamente mantenuto in una geometria non stazionaria.
2. Giannoni--Piccione e collaboratori studiano una brachistocrona relativistica con un vincolo dinamico diverso; le loro condizioni non equivalgono alla conservazione forzata di $-g(u,W)=\hat E$ fuori dal caso Killing.
3. Caponio--Corona--Giambò--Piccione trattano una carica di Noether conservata di un sistema autonomo; non coprono una carica non-Noether mantenuta mediante accelerazione propria.
4. Kovner tratta un principio di Fermat per geodetiche, anche in campi non stazionari, ma non il problema di controllo forzato qui definito.
5. T2 della fonte indicata dal referee non genera né il vincolo del rail, né la dilatazione ottica concreta, né il PMP, né la riduzione abeliana $11\times11$. L'originalità non è quindi assorbita da T2.
6. Le fonti di teoria adiabatica sostengono correttamente il meccanismo qualitativo della separatrix crossing e gli ordini asintotici citati, purché il manoscritto continui a distinguere tali risultati classici dalla verifica specifica del proprio sistema.
7. Le fonti su integrali abeliani ed elliptic polylogarithms sostengono la classificazione prudente adottata nelle parti migliori del testo; non certificano invece da sole l'irriducibilità o l'identificazione canonica del nuovo oggetto di genere due.

Non consiglio tuttavia una resubmission prima di correggere **quattro usi eccessivi o impropri delle fonti**:

- l'attribuzione a Filippov della continuità del controllo su una sola ovale non convessa come insieme;
- il residuo “Fuglede--Ishihara dichotomy” in Paper I;
- la prova incompleta della non-armonicità della proiezione in Paper II;
- l'estensione dalla non-coniugazione locale alla minimalità globale tridimensionale.

Questi punti non distruggono l'impalcatura teorica. Sono soprattutto sovra-formulazioni e dimostrazioni che devono essere ristrette a ciò che è realmente stabilito.

## 3. Matrice delle fonti fondamentali

| Fonte o gruppo | Uso nei paper | Esito | Azione |
|---|---|---|---|
| Perlick (1991) | due brachistocrone a energia fissata; fattori ottici; limite nullo | **Corretto** | Conservare; precisare sempre ramo $t$ o $\tau$ |
| Giannoni--Piccione--Verderesi (1997), Giannoni--Piccione (2002) | brachistocrona relativistica, vincoli e moltiplicatori | **Corretto e centrale** | Conservare la distinzione tra rail attivo e scivolamento senza attrito |
| Giannoni--Piccione--Tausk e Giannoni--Masiello--Piccione | Morse theory, coniugazione, finitezza | **In gran parte corretto** | Non usare i risultati di finitezza per dedurre automaticamente un no-go nello Schwarzschild esterno |
| Caponio--Corona--Giambò--Piccione (2024) | carica di Noether costante e sottovarietà del path space | **Corretto con una precisazione** | Dire che la proprietà di sottovarietà usa anche linearità e regolarità, non la sola conservazione |
| Caponio--Javaloyes (2026) | rassegna del principio di Fermat | **Citazione troppo aggregata** | Collegarla alle formulazioni moderne di Fermat; usare Giannoni per la brachistocrona controllata |
| Kovner (1990) | Fermat non stazionario per massive/massless geodesics | **Corretto** | Conservare la distinzione fra geodetica mass-shell e worldline forzata |
| Filippov (1962), Cesari | esistenza, rilassamento, continuità del controllo | **Parziale; una citazione impropria** | Separare la geometria diretta dell'ovale dal teorema su velocity sets convessi |
| Liberzon | PMP non autonomo e transversality | **Corretto nel contenuto** | Correggere anno dell'edizione e non citare pagine dalla copia lecture-notes come se fossero quelle del libro |
| Fuglede (1978), Ishihara (1979) | semiconformalità e harmonic morphisms | **Paper I errato; Paper II incompleto** | Correggere definizione/teorema e non dedurre non-armonicità dalla sola non-geodeticità delle fibre |
| Lecian, T2 | obiezione di anteriorità del referee | **Non è un antecedente del nucleo** | Rispondere che T2 riformula la semiconformalità; la dilatazione concreta è calcolata direttamente |
| Meena--Zawadzki e Bishop | Clairaut conformal submersions | **Utile ma da qualificare** | Applicare il teorema moderno con tutte le ipotesi; Bishop originale non è disponibile localmente |
| Timofeev; Cary--Escande--Tennyson; Neishtadt | adiabaticità e separatrix crossing | **Corretto come analogia/meccanismo** | Esplicitare la mappa delle ipotesi al sistema radiale e non importare automaticamente stime uniformi |
| Buchstaber--Enolskii--Leykin; Fay; Baker | integrali abeliani, sigma/zeta, terza specie | **Corretto come apparato generale** | Non presentarli come verifica dei coefficienti specifici, che è fornita dagli script/algebra del paper |
| Zagier; Beilinson--Levin; Brown--Levin | elliptic polylogarithms e iterated integrals | **Corretto se formulato prudentemente** | Scrivere “esprimibile nel framework elliptic-polylogarithmic”, non “collassa al Bloch--Wigner” senza mappa esplicita |
| Kodama; Hayward; Abreu--Visser | vettore di Kodama, corrente e massa quasi-locale | **Corretto** | Conservare la distinzione tra carica della corrente e scalare di particella controllato |
| Vaidya; Lindquist--Schwartz--Misner; Ashtekar--Krishnan | segno del flusso e firma del tubo | **Corretto e ben qualificato** | Conservare la distinzione fra accrescimento ingoing ed evaporazione come continuazione formale |
| Mello--Maciel--Zanchin; Kaloper--Kleban--Martin; Sultana--Dyer | Thakurta/McVittie/conformal Schwarzschild | **Sostanzialmente corretto** | Segnalare che tre attribuzioni storiche sono verificate tramite restatement, non tramite gli originali mancanti |

## 4. Fondazione stazionaria: Perlick

Perlick distingue due problemi a energia fissata: minimizzazione del tempo proprio e del tempo coordinato. Nel caso statico entrambi conducono a metriche riemanniane; nel caso stazionario non statico il ramo coordinato include il termine di Coriolis/una-forma. La separazione dei due clock nei manoscritti è dunque ben motivata.

È corretto anche usare Perlick come limite stazionario della costruzione, non come origine del rail non autonomo. La novità non risiede nel solo fattore ottico, che ha un antecedente chiaro, ma nel mantenimento attivo della carica quando il selettore non è Killing e nella conseguente dinamica di controllo.

Da conservare:

> In the stationary limit the construction reduces to Perlick's fixed-energy brachistochrone metric, separately for the coordinate-time and proper-time branches.

Da evitare:

> Perlick's result proves optimality of every controlled extremal.

Perlick stesso distingue extremals, minima e cut points. La minimalità deve provenire da HJB, assenza di coniugati più controllo dei cut points, o argomento globale separato.

## 5. Il precedente più vicino: Giannoni--Piccione

La fonte del 2002 formula le curve di prova imponendo, tra le altre condizioni,

\[
g(\nabla_{\dot\sigma}\dot\sigma,\dot\sigma)=0,
\qquad
g(\nabla_{\dot\sigma}\dot\sigma,Y)=0.
\]

Questo descrive lo scivolamento relativistico senza attrito nella loro formulazione. Nel caso Killing la condizione si integra in una carica costante. In assenza di stazionarietà, la fonte nota che i moltiplicatori non si eliminano in generale e che può emergere una formulazione integro-differenziale.

Il rail attuale impone invece direttamente

\[
-g(u,W)=\hat E
\]

anche quando la quantità non sarebbe conservata geodeticamente, e ammette l'accelerazione necessaria a mantenerla. Questa è una distinzione reale, non terminologica.

La formulazione consigliata è:

> The two charge laws coincide at the Killing rung, where the Giannoni--Piccione constraint integrates to a conserved charge. Away from that rung the variational problems differ: their frictionless-slide condition lets the charge drift, whereas the present rail actively holds it fixed and accounts for the required proper acceleration.

È preferibile a “the two problems coincide precisely when $W$ is Killing”, che può sembrare un'equivalenza globale anche di endpoint, clock e classe ammissibile.

La citazione ai risultati di Morse theory è appropriata quando riguarda indice, coniugazione e molteplicità nelle classi considerate dalle fonti. Non va trasformata in una garanzia diretta di minimalità globale del rail controllato.

## 6. Caponio--Corona--Giambò--Piccione: perché è un precedente vicino ma non risolutivo

La fonte del 2024 studia un Lagrangiano indefinito autonomo con simmetria infinitesima completa e carica di Noether lineare. Le curve critiche hanno carica costante e il relativo insieme è una sottovarietà $C^1$ chiusa del path space; i risultati di esistenza/moltitudine richiedono ipotesi come la pseudocoercività.

Il manoscritto è corretto nel dire che tale impianto non copre direttamente:

- un selettore non simmetrico;
- una carica che non è conservata dal moto libero;
- il mantenimento attivo della carica tramite controllo;
- il termine non autonomo e la risposta off-shell.

Va soltanto affinata questa frase di Paper I:

> ... makes the set of curves of constant Noether charge a $C^1$ closed submanifold of the path space---but the argument uses conservation ...

Proposta:

> ... obtains a $C^1$ closed submanifold of curves with constant Noether charge, using both the conserved Noether structure and the linearity/regularity assumptions on that charge. These hypotheses are absent for the actively maintained non-Noether rail considered here.

## 7. Kovner: citazione corretta

Kovner estende un principio di Fermat a campi gravitazionali arbitrari, per particelle massive e massless, mantenendo la mass shell e variando geodetiche con endpoint appropriati. La fonte sottolinea inoltre che il tempo di arrivo può avere minimo, massimo o sella.

I paper usano correttamente Kovner come antecedente non stazionario, ma separano:

- geodetiche non forzate ed estremalità in Kovner;
- worldline forzate, carica mantenuta, controllo compatto e certificazione della minimalità nel lavoro attuale.

Questa distinzione è una delle parti bibliograficamente più solide e va lasciata invariata.

## 8. P0: uso improprio di Filippov in Paper I

### 8.1 Problema

In `paper1/paper1_JMP.tex`, nella prova del controlled-rail reduction, si dice che l'ovale di velocità è “exactly the situation of Filippov's §III”. La fonte di Filippov mantiene le condizioni della sezione I, dove $R(t,x)$ è un insieme di velocità convesso; aggiunge poi stretta convessità e corrispondenza one-to-one per ottenere la continuità del controllo.

Nel paper l'indicatrice ammissibile è la **frontiera** di un corpo strettamente convesso. Come sottoinsieme del piano, la sola frontiera non è convessa. La conclusione geometrica resta vera: per l'ellisse regolare la funzione supporto ha un unico massimizzatore e tale massimizzatore è liscio per costato non nullo. Ma non è quel teorema di Filippov a provarla nelle ipotesi presentate.

### 8.2 Correzione proposta

Sostituire il paragrafo che inizia con “This is not merely an observation...” con:

> This conclusion follows directly from the geometry of the regular ellipse: for every nonzero covector its support function has a unique maximizer, and the explicit maximizing direction is smooth away from the degeneracy set. Filippov's continuous-control result is closely related but assumes a convex velocity image; our admissible set is the boundary oval itself and is not a convex subset of the plane. We therefore do not invoke that theorem for uniqueness or smoothness, nor do we convexify the mass-shell indicatrix.

La citazione Filippov può restare nella successiva prova di esistenza, dove viene costruito un epigrafo convesso. Anche lì occorre irrobustire un passaggio.

### 8.3 Saturazione dell'epigrafo

La frase “Since $\chi(1)$ is monotone in $t$, any minimizer saturates $t=F$” è troppo rapida, perché $F$ dipende anche da $\chi$. Serve un argomento di confronto per l'ODE scalare: fissato lo stesso arco $x$, la soluzione dell'uguaglianza con lo stesso dato iniziale non supera una traiettoria con $\dot\chi\ge F(x,\dot x,\chi)$, sotto la monotonicità/one-sided Lipschitz effettivamente disponibile. Se la monotonicità in $\chi$ non è garantita, bisogna formulare e provare il lemma di confronto appropriato o usare un'esistenza diretta diversa.

Testo minimo consigliato:

> For a fixed spatial arc, the scalar comparison lemma stated above implies that replacing $t\ge F(x,\nu,\chi)$ by the equality solution cannot increase the terminal clock. Hence a relaxed minimizer admits a saturated representative with the same spatial endpoints and no larger cost.

Il lemma deve essere davvero enunciato con ipotesi sufficienti.

## 9. P0: Fuglede, Ishihara e T2

### 9.1 Paper I contiene ancora l'errore già riconosciuto altrove

`paper1/paper1_JMP.tex` afferma ancora:

> by the Fuglede--Ishihara dichotomy a semiconformal map is either identically zero or conformal ...

In Fuglede questa è sostanzialmente la definizione di semiconformalità sul luogo regolare, non un teorema che conferisce conformalità alla proiezione. T2 della fonte indicata dal referee ripete lo stesso contenuto. Usarlo per provare ciò che è già assunto/definito è circolare.

Sostituzione pronta per Claude:

> Projecting along the selector orbits and equipping the base with the fixed-energy optical metric makes the projection horizontally conformal by direct calculation, with dilation $\Lambda^2=\hat E^2/[f(\hat E^2-f)]$. In Fuglede's terminology this verifies the defining semiconformality condition on the regular horizontal distribution; it is not a consequence of a dichotomy theorem. The stronger Fuglede--Ishihara result concerns harmonic morphisms and requires harmonicity in addition to semiconformality.

### 9.2 Paper II: la sola accelerazione delle fibre non basta nella forma usata

Paper II migliora correttamente la distinzione definizione/teorema, ma poi scrive che Ishihara riduce l'armonicità “for a submersion” alla minimalità delle fibre. L'equivalenza semplice vale per una **Riemannian submersion**. La proiezione del paper è presentata come horizontally conformal con dilatazione non costante. Per una horizontally conformal submersion con target di dimensione diversa da due, la tensione contiene in generale sia la curvatura media delle fibre sia un termine con il gradiente orizzontale della dilatazione. La non-geodeticità delle fibre, da sola, non esclude una cancellazione.

Si aggiunge una seconda cautela: le fonti primarie citate sono riemanniane, mentre lo spazio totale del paper è lorentziano. Il manoscritto lo riconosce, ma subito prima ha già dichiarato il fallimento del harmonic morphism come se il teorema si trasferisse integralmente.

Due opzioni corrette:

1. **Opzione conservativa, raccomandata:** mantenere soltanto la conformalità orizzontale calcolata direttamente e rimuovere la conclusione “not a harmonic morphism”.
2. **Opzione forte:** introdurre esplicitamente una metrica riemanniana associata/threading, calcolare l'intero tension field della mappa con quella metrica e dimostrare che non si annulla. Il solo valore di $\nabla_WW$ non basta.

Testo conservativo:

> The Riemannian harmonic-morphism theorem is not used here, because the total metric is Lorentzian and the projection has nonconstant dilation. What is established without signature or tension-field ambiguities is the pointwise horizontal-conformality identity and its explicit dilation. Determining harmonicity would require specifying an associated Riemannian total-space metric and evaluating the full tension field; the non-geodesicity of the selector fibres is an obstruction term but, for a general horizontally conformal submersion, is not by itself a complete proof.

### 9.3 La risposta al referee attribuisce ancora troppo a T2

In `paper2/response_to_referees_CQG.tex` compare:

> T2 supplies the existence of a conformal factor; it does not supply its value ...

T2 non costruisce né assicura l'esistenza della dilatazione del problema concreto: assume/riformula la semiconformalità. Proposta:

> T2 does not construct the optical projection or its dilation. It restates the horizontal-conformality content once semiconformality is assumed. In the present problem both the regular horizontal projection and the explicit factor $\Lambda^2=\hat E^2/[f(\hat E^2-f)]$ are obtained by solving the rail constraint jointly with the mass shell. T2 therefore neither yields the controlled dynamics nor the Abelian reduction used later.

### 9.4 Originalità rispetto a T2

L'obiezione “trivial application of T2” non è sostenuta dalla fonte. T2 parla del differenziale di una mappa semiconforme. Non contiene:

- il vincolo dinamico $-g(u,W)=\hat E$;
- il costo fisico di mantenerlo;
- la mass shell e l'indicatrice di controllo;
- il massimo hamiltoniano e la transversality non autonoma;
- il sistema lineare $11\times11$ della riduzione di Hermite;
- la separatrice di Kerr/Thakurta--Kerr;
- i termini on-shell e off-shell della risposta adiabatica;
- gli integrali abeliani iterati di peso due.

Il rapporto genuino con T2 è molto più limitato: la metrica ottica rende la proiezione orizzontalmente conforme sul dominio regolare. Questo è un utile inquadramento geometrico, non una derivazione dei risultati.

## 10. P0: minimalità equatoriale e geometria tridimensionale

Le fonti su Morse theory e geometria globale non autorizzano l'attuale salto logico in Paper II.

Nel piano equatoriale congelato, la catena è solida sotto le ipotesi dichiarate:

\[
\text{completezza della superficie} + K\le0
\Longrightarrow
\text{unicità del minimo tra lift fissati}
\Longrightarrow
\text{un minimo per classe di winding equatoriale}.
\]

In tre dimensioni, però, le winding classes equatoriali non restano classi topologiche distinte: un circuito può contrarsi uscendo dal piano. Inoltre la metrica ottica tridimensionale possiede una curvatura sezionale tangenziale positiva in regioni rilevanti, quindi Cartan--Hadamard non è disponibile globalmente.

Paper I prova $\Phi<\pi$ per archi sferici non radiali a singola escursione. Questo, insieme al Jacobi field di rotazione, può escludere il primo coniugato fuori piano su **quegli archi**. Non dimostra la minimalità globale 3D per winding class, né estende automaticamente la calibrazione HJB equatoriale.

Inoltre l'attuale

\[
\text{no conjugate point at all}
\iff K_t<0\ \text{along the arc and}\ \Phi<\pi
\]

è troppo forte. Le condizioni indicate sono sufficienti nel settore considerato, non necessarie: curvatura positiva o a segno variabile su un arco corto non produce obbligatoriamente un punto coniugato.

Testo sostitutivo:

> For a non-radial single-excursion spherical rail arc, Paper I proves $\Phi<\pi$; combined with the in-plane Jacobi estimate, this excludes conjugate points also against the symmetry-generated out-of-plane perturbation on that arc. This is a sufficient local statement. The global-minimality, HJB-calibration and one-minimizer-per-winding-class results remain equatorial, and no unrestricted three-dimensional global-minimality claim is made.

## 11. Finitezza e funzioni convesse

Giannoni--Masiello--Piccione sostengono risultati di finitezza sotto ipotesi precise: global hyperbolicity, una funzione spazialmente propria, strettamente light-convex e invariante, oltre alle condizioni di non-coniugazione rilevanti.

Quando Paper II suggerisce che l'esterno di Schwarzschild “non ha” una tale funzione, la citazione non basta da sola a provarlo. È più rigoroso scrivere:

> We do not assume or exhibit a globally defined proper strictly light-convex function satisfying the hypotheses of the cited finiteness theorem on the optical exterior considered here.

Se si vuole un vero no-go, occorre dimostrare che la topologia/geometria specifica viola una delle ipotesi esatte del teorema, senza passare per una semplice analogia con photon sphere o immagini relativistiche.

## 12. Teoria adiabatica e attraversamento di separatrice

Cary--Escande--Tennyson descrivono precisamente il breakdown dell'approssimazione adiabatica presso una separatrice: il periodo diverge, occorrono regioni before/near/after e compaiono scarti dipendenti dalla fase dell'ordine tipico $O(\varepsilon\log\varepsilon)$ nel problema classico a un grado di libertà lentamente variabile. Timofeev e Neishtadt forniscono risultati correlati in classi definite.

L'uso nel paper è corretto se resta formulato come:

- identificazione dello stesso meccanismo locale;
- non validità della normale espansione adiabatica uniforme sulla separatrice;
- necessità di una teoria di matching dedicata.

Da non scrivere senza ulteriore prova:

> Our rail system is exactly the classical theorem and therefore inherits its error bound.

Testo consigliato:

> After reduction to the frozen radial Hamiltonian, the local loss of adiabaticity has the same separatrix-crossing mechanism as in the classical one-degree-of-freedom theory. The cited estimates are not imported as a uniform bound for the present relativistic system; establishing their hypotheses and the corresponding matching problem is left open.

La fonte esatta `Neishtadt1986` non è disponibile localmente. Le affermazioni centrali sono corroborate da Cary e Timofeev; per una citazione formula-per-formula a Neishtadt occorre ottenere l'originale o una copia editoriale/autore verificabile.

## 13. Funzioni speciali: cosa è sostenuto e cosa no

### 13.1 Parte sostenuta

Buchstaber--Enolskii--Leykin supportano l'uso di sigma, zeta e inversione di Jacobi per integrali iperellittici. Fay supporta la rappresentazione delle differenziali di terza specie mediante prime forms/theta ratios. Brown--Levin forniscono un framework di iterated integrals su curve ellittiche puntate. Zagier e Beilinson--Levin sostengono la classe delle elliptic polylogarithms.

È rigoroso descrivere i $W_{jk}$ come:

> length-two iterated Abelian integrals on the relevant curve.

È anche corretto dichiarare esplicitamente congetturali:

- irriducibilità alla classe di peso uno;
- dimensione minima della base;
- completamento single-valued canonico;
- identificazione con una funzione iperellittica standard universalmente normalizzata.

### 13.2 Formulazioni da attenuare

“Collapses to the tabulated elliptic (Bloch--Wigner/Zagier) dilogarithm” richiede una mappa esplicita fra lettere, punti marcati, normalizzazione e ramo dell'integrale del paper e l'oggetto della fonte. In sua assenza:

> At the genus-one degeneration the length-two object becomes an elliptic iterated integral expressible within the classical elliptic-polylogarithm framework; an explicit identification with a particular normalized Bloch--Wigner/Zagier function is not claimed here.

La frase “not reducible ... (irreducibility conjectural)” è logicamente contraddittoria. Usare:

> not reduced here and conjecturally irreducible to weight one.

### 13.3 Affermazione falsa sulla ricostruibilità testuale

Paper II dichiara che la serie di nome di genere uno $g^{(1)}$ è “printed in Section closed”, ma nel sorgente corrente trovo soltanto menzioni di $g^{(1)}$, non la serie esplicita. Paper I, al contrario, dice correttamente che il testo non stampa una ricetta autonoma di ricostruzione.

Opzioni:

- stampare davvero definizione e serie, con convenzioni e dominio di convergenza;
- oppure sincronizzare Paper II con Paper I:

> the genus-one Kronecker--Eisenstein kernel $g^{(1)}$ identified in Section ...; its numerical nome-series implementation is contained in the archived reproducibility package and is not reproduced as a standalone reconstruction recipe in the manuscript.

## 14. Kodama, Vaidya e Thakurta: audit positivo con tre cautele

La parte Kodama è ben fondata. Kodama introduce la corrente conservata in simmetria sferica; Hayward la collega alla massa di Misner--Sharp; Abreu--Visser precisano che nel limite statico il vettore è in generale parallelo, non necessariamente identico, al Killing vector. Il paper recepisce correttamente questa normalizzazione.

È anche corretta la distinzione tra:

- la carica della corrente $G^{ab}K_b$, che dà la massa quasi-locale;
- lo scalare di particella $-u\cdot K=\hat E$, che qui è un vincolo controllato.

La firma del tubo $r=2m(v)$, con $g(T,T)=4m'(v)$ nella convenzione ingoing, è coerente: accrescimento $m'>0$ dà tubo spacelike; $m'<0$ è la continuazione timelike e non va descritta come il caso fisico principale del chart ingoing.

Per Thakurta/McVittie, le fonti moderne locali corroborano la distinzione fra conformal Schwarzschild/Thakurta, McVittie no-accretion e Sultana--Dyer. Mancano però le versioni esatte di `Thakurta1981`, `DyerHonig1979` e `SultanaDyer2004`; le attribuzioni storiche sono quindi verificate tramite restatement successivi, non direttamente sugli originali.

Cautele:

1. eliminare il paragrafo McVittie duplicato nell'introduzione di Paper II;
2. qualificare qualunque affermazione universale sulla regolarità di $r=2M$ nel caso rotante;
3. citare Mello--Maciel--Zanchin vicino alle affermazioni concretamente verificate, mantenendo le fonti storiche come attribuzione da ricontrollare appena disponibili.

## 15. Solitone: identità corretta, interpretazione fisica da correggere

Con la convenzione

\[
\tfrac12\mathcal L_Xg+\mathrm{Ric}=\lambda g,
\]

lungo una geodetica unitaria l'identità

\[
\frac{d\hat E}{d\tau}=\lambda+\mathrm{Ric}(u,u)
\]

è corretta. La risposta al referee cita però una fonte che usa una convenzione senza il fattore $1/2$; la traduzione va dichiarata esplicitamente, altrimenti cambiano le normalizzazioni di $X$, $\lambda$ e del caso conformal-Killing.

Inoltre $\mathrm{Ric}(u,u)$ non è semplicemente “the energy density the traveller measures”. Con le equazioni di Einstein è una contrazione **trace-adjusted** dello stress-energy, con eventuale termine cosmologico a seconda delle convenzioni. Sostituzione:

> No derivative of $W$ survives: the drift is the soliton constant plus the timelike Ricci-focusing scalar. Through Einstein's equation the latter is related to a trace-adjusted stress-energy contraction measured along $u$, rather than to the local energy density alone.

## 16. Esistenza e normalità: stato delle fondamenta

La strategia di Paper I è plausibile e utile: stato aumentato, insieme esteso convesso, esistenza Filippov--Cesari, poi normalità tramite positività della funzione supporto. Non ritengo che debba essere rimossa. Prima della resubmission deve però essere resa inattaccabile in tre punti:

1. dichiarare $\Omega$ come compact path-connected regular/Lipschitz domain, o assumere separatamente che la classe ammissibile fra gli endpoint sia non vuota; la sola connessione di un compatto arbitrario non basta a garantire un arco Lipschitz;
2. fornire il lemma di confronto che riporta il minimizzatore dell'epigrafo all'uguaglianza $\dot\chi=F$;
3. separare nettamente la prova diretta di unicità del support maximizer dal teorema di Filippov sui velocity sets convessi.

La frase “that substitution is why PMP applies at all” va attenuata. Una formulazione più precisa è:

> The substitution turns the selected rail formulation into a finite-dimensional non-autonomous control system with compact control on the regular domain, to which the stated PMP theorem applies directly.

## 17. Bibliografia e metadati

### P1 obbligatorio per una bibliografia pulita

`paper2/paper2.tex` contiene ancora `\nocite{*}`. Questo inserisce anche voci non usate nel ragionamento e rende difficile dimostrare al referee che ogni riferimento è pertinente. Rimuoverlo e ricompilare BibTeX.

Correzioni in `paper/refs.bib`:

- `GiannoniPiccioneTausk2002`: `pages = {697--724}`;
- `GiannoniPiccione2002`: `pages = {375--423}`;
- `Liberzon`: verificare e portare l'edizione Princeton pubblicata a 2012; la copia locale è un precursore autorizzato, non l'impaginato finale;
- `SchmidtSmith`: chiave fuorviante; gli autori sono correttamente Frittelli--Newman. Aggiungere DOI `10.1103/PhysRevD.59.124001`;
- `NielsenYeom2009`: chiave fuorviante; il lavoro è Nielsen--Yoon (2008);
- `ClaudelViracVirbhadra2001`: refuso nella chiave; aggiungere DOI `10.1063/1.1308507` ed eprint `gr-qc/0005050`;
- `Kovner1990`: verificare che il nome completo sia Israel Kovner, non Isaac, se il campo BibTeX lo espande diversamente.

Le chiavi interne possono restare per evitare rotture, purché i campi mostrati nella bibliografia siano corretti. Una rinomina va fatta solo con sostituzione globale e ricompilazione completa.

Sei nuove fonti locali non sono citate dalla coppia corrente: non vanno aggiunte per riempimento. Una bibliografia rigorosa include soltanto lavori effettivamente usati.

## 18. Cosa Claude dovrebbe lasciare invariato

1. La distinzione netta tra carica mantenuta e carica conservata da Noether.
2. Il confronto con Giannoni--Piccione basato sulle diverse leggi di accelerazione.
3. Il confronto con Kovner basato su geodetiche contro worldline forzate.
4. La gerarchia Killing $\to$ conformal-Killing $\to$ Kodama, con le cautele di normalizzazione già presenti.
5. La separazione tra ramo del tempo coordinato e ramo del tempo proprio.
6. La distinzione fra teorema, verifica simbolica, evidenza numerica e congettura nella classificazione delle funzioni speciali.
7. La risposta adiabatica completa on-shell + off-shell e il test del residuo $O(\varepsilon^2)$, purché restino chiaramente dichiarati dominio regolare e protocollo numerico.
8. La difesa dell'originalità rispetto a T2, ma con la frase corretta: T2 inquadra al massimo la conformalità orizzontale; non “supplies” il fattore né i risultati dinamici.

## 19. Ordine di intervento consigliato a Claude

### P0 — prima della resubmission

1. Rimuovere da Paper I il falso “Fuglede--Ishihara dichotomy”.
2. Correggere la risposta: T2 non fornisce l'esistenza del fattore conforme.
3. Eliminare o dimostrare mediante tension field completo la pretesa che la proiezione non sia un harmonic morphism.
4. Correggere l'“iff” sui punti coniugati e limitare l'estensione 3D alla non-coniugazione locale sugli archi coperti; lasciare minimalità globale, HJB e winding classes nel piano equatoriale.
5. Togliere l'attribuzione impropria a Filippov per l'ovale di frontiera e aggiungere il lemma di saturazione/confronto nella prova di esistenza.
6. Correggere “Ricci scalar = energy density” e sincronizzare la convenzione del solitone con quella della fonte.
7. Correggere l'affermazione secondo cui la serie $g^{(1)}$ sarebbe stampata, oppure stamparla davvero.

### P1 — nello stesso freeze editoriale

1. Rimuovere `\nocite{*}` e ricompilare la bibliografia.
2. Correggere page ranges, DOI, anno di Liberzon e metadati indicati sopra.
3. Attenuare le affermazioni su light-convex functions e sull'identificazione Bloch--Wigner/Zagier.
4. Eliminare il paragrafo McVittie duplicato.
5. Controllare che DOI/release del codice e PDF aggiornato di Paper I puntino davvero alle versioni descritte nella risposta.

### P2 — manutenzione

1. Rinominare eventualmente le chiavi BibTeX fuorvianti con sostituzione globale.
2. Acquisire legalmente gli originali mancanti tramite biblioteca/prestito interbibliotecario.
3. Aggiungere a una suite fail-closed un controllo che vieti `\nocite{*}`, citazioni mancanti e affermazioni di riproducibilità non presenti nel PDF.

## 20. Verdetto finale sulla solidità dell'edificio teorico

La risposta breve è: **sì, l'edificio teorico di base regge**, e la letteratura citata delimita abbastanza chiaramente uno spazio originale per il lavoro. Non emerge una fonte che contenga già il problema completo, né T2 può essere interpretato come una derivazione del rail o della chiusura abeliana.

Le debolezze attuali non sono un crollo del nucleo, ma quattro ponti logici costruiti troppo lunghi rispetto alle fonti: Filippov, harmonic morphisms, minimalità 3D e identificazione speciale troppo forte. Restringendo questi passaggi, il lavoro diventa più rigoroso senza perdere risultati centrali. Anzi, la difesa dell'originalità risulta più credibile perché concede con precisione ciò che appartiene alla letteratura e rivendica soltanto ciò che il manoscritto calcola e dimostra davvero.

