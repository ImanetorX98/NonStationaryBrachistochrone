# Valutazione incrociata di Paper I, Paper II e dell'obiezione basata su T2

**Data dell'audit:** 25 agosto 2026  
**Autore dell'audit:** Codex  
**Scopo:** consegnare a Claude una verifica indipendente dello stato matematico ed editoriale dei due paper, con particolare attenzione all'affermazione del referee secondo cui il lavoro sarebbe una «trivial application» del Teorema T2.

## 1. Verdetto esecutivo

Il nucleo del lavoro **non è un'applicazione banale del Teorema T2 citato dal referee**. T2 riguarda il differenziale di una mappa semiconforme e la sua restrizione al complemento del kernel; l'equazione (9) del Paper II inviato a CQG è invece una riduzione di Hermite di differenziali abeliani su una curva iperellittica. Nella fonte non compare alcun passaggio che possa produrre l'operatore $E\partial_E+J\partial_J$, il radicale della curva spettrale, il numeratore di grado cinque, i coefficienti $c_k$, la decomposizione modulo un differenziale esatto o il sistema lineare $11\times11$. L'affermazione del referee è quindi, allo stato delle fonti, **non dimostrata e matematicamente incompatibile per tipo di oggetto**.

Claude ha compiuto un lavoro sostanziale e in molti punti molto buono: ha ampliato l'audit bibliografico, corretto la replica a T2 in Paper II, delimitato il dominio del controllo, stabilizzato i livelli di evidenza, corretto la classificazione locale e incorporato la correzione completa on-shell/off-shell. Tuttavia il lavoro **non è ancora uniformemente chiuso**. Le correzioni non sono state propagate integralmente da Paper II a Paper I; inoltre restano contraddizioni nel lessico fisico, una duplicazione editoriale, un'incoerenza di normalizzazione del solitone e un disallineamento del DOI riproducibile di Paper I.

La mia conclusione operativa è:

- **originalità rispetto a T2:** difendibile con argomento forte;
- **correttezza del nucleo algebrico controllato:** sostenuta da verifiche esatte;
- **prontezza dei sorgenti per una consegna finale:** non ancora; occorre un ultimo passaggio P0/P1 elencato nella sezione 8;
- **valutazione di Claude:** lavoro scientificamente utile e in gran parte corretto, ma con errori di propagazione e almeno una convenzione matematica da correggere.

## 2. Sorgenti autorevoli usati per il confronto

### Paper I

Il sorgente aggiornato da considerare è:

`paper1/paper1_JMP.tex`

Non ho usato come fonte principale né `paper1/paper1.tex`, che è più vecchio, né la copia congelata inviata a JMP:

`paper1/submitted_JMP_2026-08/paper1_JMP.tex`

Il confronto testuale tra copia inviata e sorgente aggiornato dà **335 righe aggiunte e 19 eliminate**. Il PDF corrente è `paper1/paper1_JMP.pdf` (58 pagine).

### Paper II

Il sorgente aggiornato è:

`paper2/paper2.tex`

La copia congelata inviata a CQG è:

`paper2/submitted_CQG_2026-08-01/paper2.tex`

Il confronto dà **1834 righe aggiunte e 177 eliminate**. La risposta da valutare insieme al manoscritto è:

`paper2/response_to_referees_CQG.tex`

### Referee e fonte T2

- lettera del referee: `Response/PaperII/CQG response/CQGReviewerstask2026I.pdf`;
- fonte contenente T2: `Fonti/Lecian-2026-pseudo-Riemannian-solitons-umbilicity.pdf`, DOI `10.13140/RG.2.2.31647.21921`;
- fonte invocata dal referee per il presunto antecedente di $W$: `Fonti/Lecian-undated-generalized-Schwarzschild-solitons-geodesics.pdf`;
- fonte primaria standard sulla semiconformità: `Fonti/Fuglede-1978-harmonic-morphisms.pdf`, [articolo ufficiale e DOI 10.5802/aif.691](https://www.numdam.org/item/AIF_1978__28_2_107_0/);
- fonte primaria sul mantenimento delle funzioni armoniche: `Fonti/Ishihara-1979-mapping-preserving-harmonic-functions.pdf`, DOI `10.1215/kjm/1250522428`.

## 3. Confronto scientifico tra Paper I e Paper II aggiornati

La divisione del lavoro tra i due paper è ora concettualmente sensata.

### Paper I: fondazioni e casi sferici

Paper I contiene:

1. la definizione del rail attivamente mantenuto, $-g(u,W)=\hat E$;
2. la riduzione a controllo tempo-ottimo su indicatrice compatta e strettamente convessa;
3. esistenza e normalità sul dominio regolare;
4. il criterio HJB di verifica, dichiarato correttamente come condizionale;
5. la gerarchia Killing $\to$ conformal-Killing $\to$ Kodama;
6. FLRW e Vaidya;
7. la riduzione polinomiale dei differenziali abeliani;
8. la risposta adiabatica completa, con termine on-shell e termine off-shell.

I sei risultati condivisi sono numerati manualmente e risultano coerenti nel PDF:

| Risultato di Paper I | Funzione | Uso in Paper II |
|---|---|---|
| Theorem I.1 | riduzione controlled-rail/PMP | specializzazione a Thakurta--Kerr |
| Theorem I.2 | esistenza e normalità | legittimità degli estremali regolari |
| Proposition I.3 | verifica HJB | criterio, non garanzia automatica |
| Lemma I.4 | riduzione polinomiale $11\times11$ | chiusura delle sorgenti abeliane |
| Theorem I.5 | risposta completa di primo ordine | on-shell + off-shell in TK |
| Lemma I.6 | derivata della famiglia di degenerazione | cancellazione del polo dominante |

Paper II cita questi risultati nel ruolo corretto. Non ho trovato uno scambio di numerazione tra I.4 e I.5 nel PDF compilato, nonostante l'ordine fisico dei blocchi nel `.tex` sia diverso: le macro `\Ithm{k}` fissano esplicitamente il numero.

### Paper II: specializzazione rotante

Paper II aggiunge correttamente contenuto che non è una semplice ripetizione di Paper I:

- il caso conformal-Kerr/Thakurta--Kerr;
- l'indicatrice con vento azimutale e frame dragging;
- la delimitazione (r>2M) del problema a controllo compatto;
- le curve spettrali di genere due;
- le separatrici in forma chiusa e la classificazione locale a quattro casi;
- la separatrice retrograda esterna $r_d\simeq3.5139M$;
- le riduzioni di genere due e le degenerazioni ellittiche;
- la risposta adiabatica conformale completa e la verifica $O(\varepsilon^2)$;
- la distinzione tra teorema, prova assistita, evidenza numerica e controllo di limite.

Questa dipendenza è appropriata: Paper II usa le fondazioni di Paper I e dimostra gli elementi specifici della geometria rotante.

## 4. Cosa è migliorato davvero rispetto alle versioni inviate

### Paper I aggiornato rispetto alla copia JMP congelata

Sono miglioramenti reali:

- l'abstract ora afferma che il problema a controllo compatto termina a $r=2m(v)$;
- (m'(v)<0) nella metrica ingoing Vaidya è indicato come continuazione formale, non evaporazione fisica;
- il rapporto con Giannoni--Piccione è molto più preciso: il loro vincolo (g(a,W)=0) lascia derivare la carica, mentre il rail controllato mantiene la carica e paga la spinta;
- Kovner è separato correttamente: geodetiche, solo mass shell ed estremalità, non rail forzato e minimalità;
- PMP è presentato come condizione necessaria, non come prova automatica di ottimalità;
- sono stati aggiunti esistenza/normalità, HJB condizionale, metrica ottica dipendente dal tempo e audit bibliografico;
- la correzione adiabatica completa ha una formulazione canonica e una stima $O(\varepsilon^2)$ su sottoarchi regolari.

### Paper II aggiornato rispetto alla copia CQG congelata

Sono miglioramenti reali:

- il titolo non annuncia più una tricotomia errata;
- la classificazione locale è ora a quattro casi;
- la separatrice retrograda esterna è distinta dalle degenerazioni algebriche non fisiche;
- il dominio di controllo è esplicitato prima delle continuazioni algebriche;
- la risposta on-shell da sola mostra pendenza circa 1, mentre la risposta completa mostra circa $2.12$, coerente con un residuo $O(\varepsilon^2)$;
- la provenienza dei numeri è molto più robusta e centralizzata;
- la sezione su Fuglede/Ishihara in `paper2/paper2.tex` è ora concettualmente corretta;
- la risposta al referee concede il legame geometrico autentico con le sommersioni, senza accettare la falsa dipendenza da T2.

## 5. Analisi forense del Teorema T2

### 5.1 Cosa dice effettivamente T2

La fonte di Lecian introduce una mappa $\psi:(P^m,g)\to(Q^n,h)$, la chiama semiconforme e poi enuncia:

> o $\psi$ è identicamente nulla, oppure è conforme e suriettiva sul complemento di $\ker d\Psi_x$.

Subito dopo identifica il complemento del kernel con lo spazio orizzontale e introduce una dilatazione $\Lambda$ tramite

\[
h(d\psi X,d\psi Y)=\Lambda^2 g(X,Y),\qquad X,Y\in(\ker d\psi)^{\perp}.
\]

Nella formulazione standard di Fuglede, questo contenuto è precisamente la **definizione di mappa semiconforme**: fuori dai punti critici, $d\psi$ ristretto al complemento ortogonale del kernel è conforme e suriettivo. Fuglede dimostra invece un fatto più forte e diverso:

\[
\text{morfismo armonico}
\quad\Longleftrightarrow\quad
\text{mappa semiconforme e armonica}.
\]

Ishihara, Theorem 5.1, parte dall'ipotesi che la mappa preservi le funzioni armoniche e conclude che essa è una pseudo-sommersione conforme armonica, oppure costante quando la dimensione del dominio è minore di quella del codominio. Non parte dalla sola semiconformità.

Quindi T2, letto contro la letteratura primaria, non costruisce una nuova struttura dinamica: al massimo riformula la proprietà puntuale orizzontale di una mappa già presentata come semiconforme. Nella fonte T2 non è accompagnato da una derivazione dell'equazione (9) del paper, né da un calcolo di brachistocrone.

### 5.2 Perché T2 non produce l'equazione (9)

L'equazione (9) della versione inviata di Paper II è

\[
(E\partial_E+J\partial_J)F
=\frac{d}{dr}\left(\frac{\mathcal A(r)}{\mathcal R}\right)
+\sum_{k=0}^{4}c_k\frac{r^k}{\mathcal R}.
\]

Questa è una decomposizione di un differenziale meromorfo sulla curva iperellittica in:

1. un differenziale esatto;
2. una base di differenziali di prima/seconda specie;
3. coefficienti determinati da identità polinomiali.

Nel formalismo aggiornato la decomposizione richiede un numeratore

\[
\mathcal A(r)=\sum_{i=0}^{5}a_i r^i
\]

e il confronto degli undici coefficienti di $r^0,\ldots,r^{10}$, cioè un sistema $11\times11$. Per la sestica di Vaidya il determinante verificato è

\[
\det M=-32(\hat E^2-1)\operatorname{disc}_r S.
\]

T2 non contiene alcuno degli oggetti necessari a questo calcolo. Il disallineamento è formale:

| T2 | Equazione (9) |
|---|---|
| mappa tra varietà | differenziale su curva algebrica |
| differenziale $d\psi_x$ | derivata parametrica $E\partial_E+J\partial_J$ |
| kernel e spazio orizzontale | radicale $\mathcal R$, poli e coomologia abeliana |
| conformalità puntuale | riduzione modulo differenziali esatti |
| nessun coefficiente costruttivo | sistema per $a_i,c_k$ |

Per sostenere una «trivial application» occorrerebbe almeno una mappa esplicita tra ipotesi e oggetti, seguita da una breve catena che produca l'identità. La relazione non è fornita dal referee né dalla fonte.

### 5.3 Il legame autentico con T2/Fuglede

Un legame vero esiste e conviene conservarlo nella risposta.

Per una metrica stazionaria in forma threading, la proiezione lungo le orbite del selettore è orizzontalmente conforme rispetto alla metrica ottica a energia fissata e ha dilatazione

\[
\Lambda^2=\frac{\hat E^2}{f(\hat E^2-f)}.
\]

Questo si verifica direttamente calcolando la metrica sui vettori orizzontali. La teoria delle mappe semiconformi fornisce il linguaggio corretto per descrivere la forma della proiezione; **non fornisce il valore di $\Lambda$, non dimostra che la specifica proiezione sia semiconforme e non genera il problema di controllo**.

Inoltre la proiezione non è in generale un morfismo armonico. Le fibre sono le orbite di $W$; essendo unidimensionali, sono minimali solo se geodetiche. Già nel caso Schwarzschild statico

\[
\nabla_WW=\frac{M}{r^2}f\,\partial_r\ne0,
\]

che è la spinta dell'osservatore statico. È quindi proprio il costo fisico del rail a impedire il passaggio alla struttura armonica più forte. Va mantenuta anche la cautela di firma: Fuglede e Ishihara lavorano su varietà riemanniane; nel paper lo spazio totale è lorentziano, mentre la distribuzione orizzontale usata nel calcolo è spaziale.

### 5.4 La seconda affermazione del referee: $W$ deriverebbe da un precedente Hamiltoniano

La seconda fonte di Lecian studia solitoni di Ricci e traiettorie geodetiche/fotoniche attraverso un approccio Euler--Lagrange/Hamiltoniano. Nel testo non ho trovato:

- un selettore arbitrario $W$ per un rail;
- il vincolo attivamente mantenuto $-g(u,W)=\hat E$;
- la gerarchia Killing--CKV--Kodama;
- il vettore di Kodama;
- l'indicatrice compatta del controllo;
- il support Hamiltonian o la forza propria necessaria a mantenere il rail.

La relazione corretta è più debole: entrambi i programmi incontrano il tensore di deformazione $\mathcal L_Wg$. Se $W$ è anche un campo di solitone, la deriva della carica può essere riscritta in termini di curvatura. Questo è un **caso speciale di contatto**, non una derivazione della definizione del rail e non include i selettori di Thakurta--Kerr o di Vaidya in generale.

## 6. Dove risiede davvero l'originalità

La difesa più forte non consiste nel dire che ogni strumento usato è nuovo. PMP, HJB, metriche di Randers, riduzione di Hermite e integrali abeliani sono strumenti classici. Anche il solo fatto che una proiezione ottica sia orizzontalmente conforme non è nuovo.

Il contenuto originale difendibile è la combinazione strutturale seguente:

1. **scelta del ramo complementare a Giannoni--Piccione:** invece di imporre $g(a,W)=0$ e lasciare derivare la carica, si mantiene $-g(u,W)=\hat E$ e si paga esattamente il difetto di Killing mediante accelerazione propria;
2. **passaggio da geodetiche a controllo attivo:** gli estremali non sono geodetiche di una metrica pesata; il vincolo supplementare produce un insieme locale di velocità ammissibili e un Hamiltoniano come funzione di supporto;
3. **dominio causale del problema:** il controllo compatto esiste solo dove il selettore è timelike e $\hat E>|W|$, e può cessare di essere posto anche se la geometria o le continuazioni algebriche proseguono;
4. **fondazioni non autonome:** esistenza/normalità, stato aumentato, transversality e criterio HJB per un rail non conservato da simmetria;
5. **selettori dinamici:** uso coerente della gerarchia Killing--conformal-Killing--Kodama;
6. **realizzazioni esplicite:** Vaidya e Thakurta--Kerr, con curve spettrali, separatrici e degenerazioni;
7. **risposta adiabatica completa:** distinzione e chiusura del termine on-shell e del termine off-shell, con recupero della pendenza quadratica del residuo;
8. **connessione controllabile tra difetto di simmetria e sforzo:** la quantità che distrugge la conservazione geodetica è la stessa che il controllo deve compensare.

La frase di originalità più robusta sarebbe quindi:

> Il contributo non è una nuova nozione di semiconformità né una nuova versione del principio di Fermat. È una formulazione tempo-ottima per worldline timelike forzate che mantengono attivamente una carica prescritta quando la simmetria necessaria alla sua conservazione non esiste, insieme alla delimitazione del dominio di controllo e alla sua soluzione perturbativa esplicita nei casi Vaidya e Thakurta--Kerr.

Questa è una rivendicazione più precisa e più credibile di «tutto è nuovo». L'audit esclude la dipendenza da T2, ma non sostituisce una ricerca bibliografica universale di priorità; la difesa migliore resta il confronto puntuale con i precedenti pubblicati più vicini, soprattutto Giannoni--Piccione e Caponio--Corona--Giambò--Piccione.

## 7. Verifiche matematiche eseguite in questo audit

Ho rieseguito:

`python3 paper1/verification/verify_paper1_core.py`

Esito: **tutti i controlli passati con residuo esattamente nullo**. In particolare:

- indicatrice di Vaidya e funzione di supporto;
- determinante del sistema $11\times11$ contro il discriminante in tre punti razionali;
- rango 10 e kernel unidimensionale sul luogo di degenerazione;
- identità di omogeneità e cancellazione di $N_{\rm tot}(r_d)$;
- identità canoniche del Theorem I.5 e forma di bordo off-shell.

Questo non costituisce una nuova prova indipendente di ogni frase dei due manoscritti, ma conferma il nucleo algebrico che T2 sarebbe accusato di rendere banale.

Il tentativo di rieseguire `paper2/verification/verify_submersion_link.wls` non ha prodotto output in questa sessione e il processo Wolfram è stato terminato dopo l'attesa. Non considero quindi quel run una verifica completata. Le identità centrali della sommersione sono comunque verificabili direttamente e la loro lettura bibliografica è stata confrontata con Fuglede e Ishihara.

I log LaTeX correnti non mostrano riferimenti indefiniti o errori fatali. Paper II conserva alcuni `Overfull \hbox` fino a circa 23 pt e un warning di font; sono problemi editoriali, non matematici.

## 8. Problemi residui trovati

### P0 — da correggere prima di considerare i testi chiusi

#### 8.1 Paper I conserva l'argomento circolare Fuglede--Ishihara

In `paper1/paper1_JMP.tex:1055-1058` si legge ancora che, «by the Fuglede--Ishihara dichotomy», una mappa semiconforme sarebbe nulla oppure conforme sul complemento del kernel. Paper II e la risposta riconoscono correttamente che questa è la definizione di Fuglede, non un teorema. La correzione non è stata propagata a Paper I.

**Azione:** sostituire quel paragrafo con il calcolo diretto della conformalità orizzontale, seguito dalla distinzione tra semiconformità e morfismo armonico.

#### 8.2 Incoerenza del fattore 2 nella convenzione di solitone

In `paper2/response_to_referees_CQG.tex:750-757` la fonte viene descritta con

\[
\mathcal L_Xg+\mathrm{Ric}=\lambda g,
\]

ma subito dopo si conclude

\[
\dot{\hat E}=\lambda+\mathrm{Ric}(u,u),
\]

che è la formula usata in `paper2/paper2.tex:910-918` per la convenzione

\[
\tfrac12\mathcal L_Xg+\mathrm{Ric}=\lambda g.
\]

Con la convenzione scritta nella risposta, invece, si ottiene

\[
\dot{\hat E}=\tfrac12\big[\lambda+\mathrm{Ric}(u,u)\big].
\]

Le due convenzioni sono convertibili riscalando il campo del solitone, ma la conversione deve essere dichiarata. Anche la condizione CKV cambia da $\mathrm{Ric}=(\lambda-\psi)g$ a $\mathrm{Ric}=(\lambda-2\psi)g$ se si mantiene la convenzione senza $1/2$.

**Azione:** scegliere una convenzione unica, dire esplicitamente come si traduce quella della fonte e sincronizzare manoscritto, risposta e script `verify_soliton_rail.wls`.

#### 8.3 Il protocollo del dominio è contraddetto dal corpo di Paper II

`paper2/paper2.tex:638` prescrive che sotto $r=2M$ non si usino «capture», «plunge» o «bounce» e non si parli di brachistocrona fisica. Tuttavia ricompaiono, tra gli altri:

- `paper2/paper2.tex:1903-1927`;
- `paper2/paper2.tex:2990-3060`;
- `paper2/paper2.tex:3118-3122`;
- l'abstract, con «plunge inversion»;
- `paper2/paper2.tex:3674`, dove una radice interna all'ergosfera viene detta «on the physical arc».

Non è solo lessico: alcune didascalie chiamano «brachistochrones» curve che proseguono fino a (r_+), dopo che il problema di controllo è dichiarato non posto.

**Azione:** troncare la lettura fisica a $r=2M$; mantenere sotto la superficie soltanto «analytic continuation/root pattern/continued branch». Qualificare senza ambiguità quale separatrice sia esterna e quale degenerazione sia interna/algebrica.

#### 8.4 La stessa contraddizione rimane in Paper I

L'abstract di Paper I è stato corretto, ma il corpo usa ancora «plunge», «capture», «horizon surface» e «ergosphere trichotomy»; si vedano `paper1/paper1_JMP.tex:236`, `:612-616`, `:1398-1401`, `:1439`, `:2032-2043`.

Inoltre `paper1/paper1_JMP.tex:594-597` afferma che la spinta diverge avvicinandosi alla «freezing/horizon surface», mentre `:1035-1043` distingue correttamente la divergenza alla freezing surface dalla perdita di compattezza a costo finito sul bordo del dominio. Le due affermazioni non possono restare insieme.

**Azione:** propagare la tassonomia E/C/A e la distinzione tra freezing e perdita di compattezza in tutto Paper I, incluse introduzione, figure e conclusioni.

#### 8.5 Paper I aggiornato cita ancora l'archivio congelato v1.2

`paper1/paper1_JMP.tex:2947-2953` punta a Zenodo v1.2.0, DOI `10.5281/zenodo.22035415`, mentre il repository corrente e Paper II usano v1.5.0, DOI `10.5281/zenodo.22079449`. Poiché il sorgente Paper I è stato modificato dopo la copia inviata, l'archivio v1.2 non può essere assunto automaticamente come riproduzione esatta del testo aggiornato.

**Azione:** o congelare e dichiarare Paper I come versione v1.2, oppure aggiornare citazione, tabella script e data availability a una release che contenga davvero il sorgente e gli script correnti.

### P1 — importanti, ma successive ai P0

#### 8.6 La risposta dice che «T2 supplies the existence of a conformal factor»

`paper2/response_to_referees_CQG.tex` usa questa formulazione. È ancora troppo generosa e tecnicamente imprecisa: T2 non costruisce il fattore e, se si assume già semiconformità nel senso standard, ne riformula la proprietà definitoria.

**Azione:** scrivere che T2 «encodes the horizontal-conformality property once semiconformality is assumed», mentre $\Lambda$ e la semiconformità della proiezione concreta sono ottenute per calcolo diretto.

#### 8.7 Paragrafo McVittie duplicato

In `paper2/paper2.tex:331-355` la descrizione della metrica di McVittie, inclusa la citazione a Kaloper, è ripetuta quasi integralmente.

**Azione:** fondere i due blocchi in un solo paragrafo.

#### 8.8 Documenti di stato non sincronizzati

`Fonti/README.md` contiene la correzione giusta su Fuglede; `POST_REVIEW_CHANGES.md:123-126` e `:385-390`, invece, dicono ancora che la conformalità è «forced by Fuglede--Ishihara». Anche `WORKING_METHOD.md` registra la formula del solitone senza rendere esplicita la convenzione.

**Azione:** aggiornare i documenti durevoli, altrimenti Claude rischia di reintrodurre l'errore in un passaggio successivo.

#### 8.9 Due sorgenti concorrenti per Paper I

`paper1/paper1_JMP.tex` è il sorgente aggiornato, mentre `paper1/paper1.tex` è più vecchio. Il nome `paper1_JMP.tex` può indurre a usare la copia sbagliata o a credere che sia la versione inviata e congelata.

**Azione:** dichiarare in un README quale sia il sorgente canonico e, se possibile, rinominarlo o sincronizzare esplicitamente le varianti dopo la resubmission.

#### 8.10 Richiesta del referee sull'esperimento

La risposta «no experimental validation is available» è onesta e preferibile a una falsa pretesa osservativa, ma non soddisfa completamente una richiesta formulata come «strict». Si può mantenere la distinzione tra validazione del modello gravitazionale ed emulazione analogica, aggiungendo un protocollo possibile con ottica, fluidi o metamateriali tempo-modulati, presentato come test dell'Hamiltoniana efficace e delle leggi di scala, non come riproduzione della gravità di Kerr.

### P2 — pulizia editoriale

- eliminare la nota iniziale di Paper I secondo cui il contenuto fisico sarebbe invariato: non è più vera;
- aggiornare in Paper I «ergosphere trichotomy» alla classificazione corrente;
- risolvere gli `Overfull \hbox` principali di Paper II e il warning `OMS/cmtt`;
- verificare che cover letter, response, sorgente canonico, PDF CQG e DOI vengano rigenerati dallo stesso commit finale.

## 9. Testo suggerito per correggere la parte T2

### Sostituzione consigliata in Paper I

> Projecting along the selector orbits and equipping the base with the fixed-energy optical metric makes the projection horizontally conformal by direct calculation, with dilation $\Lambda^2=\hat E^2/[f(\hat E^2-f)]$. In Fuglede's terminology this verifies the defining condition of semiconformality on the regular set; it is not a consequence of a dichotomy theorem. The stronger Fuglede--Ishihara result concerns harmonic morphisms, which require both semiconformality and harmonicity. The present fibres are not geodesic and hence not minimal even in the static Schwarzschild limit, so the projection is not a harmonic morphism. Only the pointwise horizontal-conformality statement is used here, with the additional caveat that the total spacetime metric is Lorentzian.

### Correzione consigliata nella risposta al referee

> The quoted T2 does not construct the conformal factor. Once a map is assumed semiconformal, its conclusion restates the horizontal conformality of $d\psi$ off the kernel. In our setting both semiconformality and the value of the dilation are obtained directly from the threading metric and the fixed-energy mass-shell constraint. The cited result therefore supplies appropriate terminology for one geometric layer of the construction, but neither equation (9), the controlled-rail Hamiltonian, nor the selector $W$.

### Frase consigliata sulla normalizzazione del solitone

> We use the standard convention $\mathrm{Ric}+\tfrac12\mathcal L_Xg=\lambda g$. The cited preprint writes $\mathrm{Ric}+\mathcal L_{X_{\rm src}}g=\lambda g$; the two agree after $X=2X_{\rm src}$. With our convention, $d\hat E/d\tau|_{\rm geod}=\lambda+\mathrm{Ric}(u,u)$.

## 10. Ordine di lavoro consigliato a Claude

1. **Non modificare le copie congelate** in `submitted_JMP_2026-08` e `submitted_CQG_2026-08-01`.
2. Correggere T2 in Paper I usando la fonte primaria Fuglede.
3. Uniformare la convenzione del solitone tra fonte, risposta, Paper II e script.
4. Eseguire una ricerca globale di `plunge|capture|bounce|horizon|physical arc|trichotomy` in entrambi i paper e applicare il protocollo E/C/A.
5. Separare in ogni punto freezing, perdita di compattezza e continuazione algebrica.
6. Rimuovere la duplicazione McVittie.
7. Aggiornare il DOI di Paper I oppure dichiarare esplicitamente quale release lo riproduce.
8. Sincronizzare `Fonti/README.md`, `POST_REVIEW_CHANGES.md`, `WORKING_METHOD.md` e il handoff finale.
9. Rigenerare Paper I, Paper II, response e cover letter dallo stesso commit.
10. Rieseguire i test simbolici, la provenance completa e una ricerca di riferimenti indefiniti prima del freeze.

## 11. Conclusione finale sull'originalità

Il referee ha identificato un contatto geometrico reale — la proiezione ottica è naturalmente descritta dal linguaggio delle sommersioni orizzontalmente conformi — ma ha esteso quel contatto oltre ciò che T2 dimostra. T2 non contiene il problema di controllo, non impone il rail, non determina il selettore, non produce l'equazione (9) e non chiude la risposta adiabatica.

Il lavoro non è originale perché «usa funzioni speciali» o perché «applica PMP»: questi sono strumenti. È originale, nella forma difendibile emersa dall'audit, perché trasforma la perdita di una carica di simmetria in un vincolo attivamente mantenuto, identifica il relativo dominio di esistenza e la spinta necessaria, e sviluppa questa struttura non autonoma fino a soluzioni e correzioni esplicite in Vaidya e Thakurta--Kerr.

La tesi centrale è quindi difendibile, ma la difesa sarà molto più forte dopo le correzioni P0: soprattutto la rimozione dell'argomento circolare da Paper I, la convenzione del fattore 2 e la coerenza assoluta del dominio fisico.

## 12. Secondo passaggio di audit: verifica puntuale della major revision

Questa sezione integra e, dove necessario, rende più severa la valutazione precedente. Ho confrontato riga per riga:

- il report principale del referee, con i dodici major comments e i dieci minor comments;
- il documento del secondo referee;
- la checklist formale della revised submission;
- `paper2/response_to_referees_CQG.tex`;
- i sorgenti aggiornati `paper1/paper1_JMP.tex` e `paper2/paper2.tex`;
- il PDF pulito, quello evidenziato, la risposta compilata e il pacchetto `paper2/submission_CQG_R1/`;
- le fonti primarie locali su Perlick, Giannoni--Piccione, Fuglede, Ishihara, Myers e i preprint citati dal secondo referee;
- gli script e la provenance della release riproducibile.

Il verdetto aggiornato è più articolato del semplice «Claude ha lavorato bene/male»:

1. **la ricostruzione scientifica principale di Paper II è sostanziale e in larga parte riuscita**;
2. **la classificazione locale a quattro casi e la separatrice esterna sono corrette nei controlli indipendenti svolti**;
3. **la correzione adiabatica completa è riproducibile e mostra il comportamento quadratico promesso**;
4. **alcune correzioni dichiarate nella response non sono state applicate “throughout”**;
5. **una nuova estensione tridimensionale del risultato di minimalità è formulata troppo fortemente e deve essere ridimensionata prima della resubmission**;
6. **il pacchetto CQG già predisposto non è la versione corrente e non deve essere caricato**.

La resubmission ha quindi buone possibilità dopo un ultimo passaggio mirato, ma **non considero ancora sicuro inviare i file nello stato attuale**.

## 13. Matrice di conformità ai dodici major comments del Referee 1

| Commento | Esito dell'implementazione | Valutazione indipendente |
|---|---|---|
| 1. Dominio exterior/contact/continuation | **Parziale** | Il protocollo E/C/A in Paper II è una buona soluzione, ma il corpo e varie didascalie continuano a usare “plunge”, “capture”, “bounce” e “horizon” sotto $r=2M$. Lo stesso problema rimane in Paper I. La response afferma invece che il protocollo è usato “without exception”. |
| 2. Estremale non implica minimo | **Parziale, con una nuova sovra-estensione** | Paper II corregge bene la logica PMP e ottiene risultati forti nel settore congelato equatoriale. Paper I conserva però una pretesa di “local minimizer” basata su perturbazioni numeriche. Inoltre Paper II estende impropriamente il risultato globale equatoriale alla geometria tridimensionale. |
| 3. Protocollo di endpoint e costati | **Sostanzialmente corretto** | La distinzione $J_t=A J_\eta$, i clock e il protocollo a endpoint libero sono ora espliciti. Va mantenuta questa precisione in tutte le figure e nel pacchetto finale. |
| 4. Contraddizione cusp/corner | **Corretto** | La classificazione a quattro casi distingue cuspide semplice, marginale retrogrado asintotico a pendenza finita e marginale progrado attraversante. I numeri e la fattorizzazione sono coerenti con i controlli indipendenti. |
| 5. Separatrice esterna come risultato centrale | **Corretto** | Ho riprodotto indipendentemente $r_d=3.513905124011658M$ e $J_c^-=-8.053516003877020$. Le due disuguaglianze di dominio danno $g(W,W)=-0.430833$ e $\bar v^2=0.700811$. La radice prograda $r_d\simeq1.512292M$ è correttamente interna a $2M$ e non fisica per il selettore scelto. |
| 6. Stato stratificato del fixed-endpoint result | **Corretto nel sorgente corrente** | Teorema analitico, dominio limitato/asintotico, CAP rappresentativo e congettura risultano distinti. La prova assistita resta da consegnare in un pacchetto finale rigenerato e realmente eseguibile. |
| 7. Sorgente adiabatica completa e provenance | **Corretto nel contenuto; tooling migliorabile** | Il termine off-shell è incluso, la finestra $\varepsilon$ e il sottointervallo regolare sono dichiarati, e il run da clone pulito riproduce lo stesso digest. Tuttavia `make_provenance.py --check` non controlla davvero che i file preesistenti siano invariati: li riscrive sempre. |
| 8. Terminologia higher genus | **Quasi corretto** | Il testo adotta la denominazione neutra e dichiara congetturale l'irriducibilità, ma a `paper2/paper2.tex:1645-1648` e `:3620-3625` compare ancora “do not/not reducible” accanto alla qualifica “conjectural”. Va scritto “no such reduction is known; irreducibility is conjectural”. |
| 9. Linguaggio fisico | **Non completato** | La distinzione tra energia mantenuta, costato $J$ e spinta è buona; il lessico delle continuazioni interne non è stato però ripulito globalmente. |
| 10. Riproducibilità | **Molto migliorata, ma non pronta per l'upload** | Manifest di 18 voci: nessun file mancante e nessun hash discordante. Tutte le 22 figure del sorgente corrente si risolvono. Il pacchetto `submission_CQG_R1`, però, è una copia intermedia di 3130 righe/67 pagine contro 3760 righe/79 pagine del corrente. |
| 11. Dipendenze da Paper I | **Parziale** | Esistenza e normalità sono state ampliate e il dominio Vaidya è meglio delimitato. Restano l'argomento circolare Fuglede--Ishihara, il lessico interno e la pretesa di minimalità numerica. |
| 12. Posizionamento e novità | **Corretto nel nucleo** | La distinzione dalla letteratura stazionaria e da Giannoni--Piccione è ora molto più forte. La replica a T2 è difendibile, ma due frasi tecniche della response devono essere corrette per non indebolirla. |

### Esito dei minor comments

M1--M3 e M7--M9 sono implementati in modo sostanzialmente convincente. M4 e M6 sono soltanto parziali perché alcune didascalie e alcuni paragrafi continuano a chiamare fisica una continuazione. M5 è in gran parte rispettato, ma la response ha un problema tipografico in `Protocol~\msProtNames and`: nel PDF compare “Protocol 3and”. M10 richiede ancora un'ultima riduzione e sincronizzazione della conclusione dopo le correzioni scientifiche qui indicate.

## 14. Nuovi risultati della verifica matematica indipendente

### 14.1 Risultati confermati

Il controllo indipendente in SymPy/mpmath, senza usare le formule di output dei notebook Wolfram, conferma:

- l'identità polinomiale dietro il bound di corda di Paper I;
- la positività del polinomio traslato usato per provare $\Phi<\pi$ su una singola escursione radiale;
- per $\hat E=1.4$, $R=10M$, il massimo numerico $\Phi/\pi\simeq0.74320$ vicino a $r_t\simeq2.98M$;
- il polinomio di soglia della curvatura equatoriale e la soglia $\hat E^2=3/2$;
- i valori della separatrice retrograda esterna e della degenerazione prograda interna;
- la formula chiusa di $J_c^-$ valutata al double root;
- la negatività della curvatura scalare ottica, sebbene con una correzione editoriale ai coefficienti riportati nel testo.

Questi controlli danno ragione alla parte migliore del lavoro di Claude: la separatrice e la struttura di curvatura non sono state inserite frettolosamente o soltanto per assecondare il referee.

### 14.2 Errore P0: il criterio tridimensionale non è un “se e solo se”

In `paper2/paper2.tex:1243-1247` compare

\[
\text{no conjugate point at all}
\iff K_t<0\text{ along the arc and }\Phi<\pi.
\]

Il verso sufficiente è valido nel settore non radiale considerato: $K_t<0$ esclude il campo di Jacobi nel piano e $\Phi<\pi$ esclude lo zero del campo di rotazione fuori piano. Il verso necessario è falso in generale. Una curvatura positiva o a segno variabile su un arco sufficientemente corto non produce automaticamente un punto coniugato. Il testo deve usare

\[
K_t\le0\text{ along the arc and }\Phi<\pi
\quad\Longrightarrow\quad
\text{no conjugate point on that arc},
\]

oppure enunciare separatamente il risultato specifico già dimostrato da Paper I per gli archi monotoni.

Vi è anche una qualificazione dimenticata: il campo di Killing scritto nel testo è il modo corretto di controllare la direzione fuori piano per una geodetica **non radiale**, che individua un piano orbitale. Per $J=0$ il generatore scelto può annullarsi lungo l'intera geodetica radiale e non fornisce il campo di Jacobi non banale desiderato. Il caso radiale va escluso dall'enunciato o verificato separatamente.

### 14.3 Errore P0: il bound di Paper I non rende tridimensionale la minimalità globale per winding class

Il passaggio `paper2/paper2.tex:1248-1252` afferma che, poiché Paper I prova $\Phi<\pi$, le Propositions sulla non-coniugazione, sulla minimalità globale e sull'HJB valgono in piena dimensione tre.

Paper I, Theorem `thm:chord`, prova invece esplicitamente $\Phi<\pi$ **su una singola escursione radiale** e precisa che il meccanismo Maxwell richiede più di un'escursione. Gli archi con winding più elevato raggiungono un azimut totale $\pi$, dove il campo di Jacobi rotazionale si annulla. Inoltre l'argomento di Cartan--Hadamard “un minimizzatore per classe di omotopia” è un argomento sulla superficie equatoriale anulare e sul suo universal cover. Nello spazio tridimensionale $\{r>2M\}\times S^2$ tali winding equatoriali non sono classi di omotopia distinte: possono contrarsi uscendo dal piano.

La correzione sicura è:

> For a non-radial single-excursion spherical rail arc, Paper I proves $\Phi<\pi$; combined with the in-plane Jacobi result, this excludes conjugate points also against out-of-plane perturbations on that arc. The global-minimality and one-per-winding-class statements remain equatorial. We make no unrestricted three-dimensional global-minimality claim.

Questo punto è scientificamente importante perché tocca proprio la risposta al major comment 2. Non conviene lasciare al referee una nuova pretesa globale più forte di quella che aveva chiesto di rimuovere.

### 14.4 Coefficienti della curvatura scalare: fattore comune due

In `paper2/paper2.tex:1032-1035` si afferma che il bracket di $R_{\rm opt}$, dopo $r=2M(1+s)$, abbia coefficienti

\[
24\hat E^4,\qquad48\hat E^2(\hat E^2-1),\qquad
8(3\hat E^4-6\hat E^2+8).
\]

I coefficienti del bracket sono invece

\[
12\hat E^4,\qquad24\hat E^2(\hat E^2-1),\qquad
4(3\hat E^4-6\hat E^2+8).
\]

I numeri stampati sono quelli ottenuti moltiplicando anche per il prefattore esterno $2$. Il segno e quindi la conclusione $R_{\rm opt}<0$ restano corretti; è un errore di normalizzazione nella spiegazione, non nel risultato.

### 14.5 Due “iff” puntuali da qualificare

In `paper2/paper2.tex:870-873` si dice che $\partial_\eta\Lambda^2$ “vanishes if and only if $\psi=0$”. Come identità lungo l'intera fibra, il senso voluto è corretto; puntualmente è falso. Per esempio $M=A=1$, $\hat E=1.2$ e

\[
r=\frac{4MA^2}{2A^2-\hat E^2}=7.142857\ldots M
\]

annullano il fattore tra parentesi pur con $\psi\ne0$, e il punto è ancora nel dominio di controllo. Scrivere “vanishes identically along a fibre iff $\psi=0$”.

In `paper2/paper2.tex:892-895` si dice inoltre che il programma vive nel complemento delle ipotesi “umbilical, Clairaut and conformal-Killing”, mentre poche righe prima si è provato esattamente che $W=\partial_\eta$ è conformal Killing per Thakurta--Kerr. Ciò che fallisce è la proprietà Killing/Clairaut e la discesa della dilatazione sul quoziente, non la proprietà conformal-Killing. La frase va riscritta.

### 14.6 Segno della spinta minima

La norma di un'accelerazione non può essere limitata inferiormente da una quantità negativa. Se non si assume esplicitamente $A'>0$, `paper2/paper2.tex:830-839` deve contenere

\[
|a|\ge \frac{|\varepsilon|}{\hat E\bar v},
\]

e analogamente nel testo e nello script. Se il manoscritto vuole trattare solo espansione, deve dichiarare $\varepsilon\ge0$ nell'ipotesi del risultato.

### 14.7 Higher-genus: la cautela è quasi, non completamente, uniforme

La response afferma correttamente che l'irriducibilità a weight one è congetturale. Il corpo contiene però la costruzione grammaticale “do not reduce ... because ...” e, più avanti, “not reducible ... (their irreducibility is conjectural)”. Le due metà non sono logicamente compatibili. Non è una ragione per rimuovere le funzioni speciali o il termine off-shell; basta sostituire la pretesa negativa con lo stato epistemico corretto.

### 14.8 Pretesa eccessiva sulla convex function

`paper2/paper2.tex:1345-1348` afferma in generale che “a black-hole exterior supplies none” riferendosi a una light-convex function. Questa frase universale non è dimostrata dalle fonti citate. È sufficiente scrivere che **lo specifico esterno Schwarzschild con photon sphere non soddisfa l'ipotesi globale confinante invocata**, oppure dare un teorema preciso. Il nesso qualitativo con le immagini relativistiche può restare come motivazione, non come prova.

## 15. Audit di riproducibilità

### 15.1 Risultato positivo: la provenance numerica si riproduce

Da un clone locale pulito del commit corrente ho eseguito:

`python3 NonStationaryMetrics/paper2/provenance/make_provenance.py --check`

Il run è terminato con codice 0 e ha riprodotto:

- tre configurazioni con slope esatta $2.12\pm0.03$;
- slope del run first-order esatto $2.15$;
- sequenze su finestre decrescenti $2.118\to2.068\to2.043$, $2.123\to2.072\to2.045$, $2.122\to2.070\to2.044$;
- finestra $\varepsilon=0.001,0.002,0.004,0.008,0.016$;
- fit su $r/M\in[8,11]$, a distanza minima $6M$ dai loci esclusi;
- lo stesso digest riproducibile del working tree:
  `5ddfa9130887f96783e7d9e9875598413b315a98e8217e70a38feeb48e925ce8`.

Questa è una conferma importante della correzione lineare adiabatica completa, incluso il termine off-shell. Il fatto che l'esponente a finestra finita sia circa $2.12$ e converga verso $2$ restringendo la finestra è coerente con un residuo $c_2\varepsilon^2(1+O(\varepsilon))$.

### 15.2 Bug nel significato di `--check`

Il docstring di `make_provenance.py` promette che `--check` verifichi se gli output siano aggiornati. L'implementazione, invece:

1. riscrive sempre `MANIFEST.tsv`, JSON e frammento TeX;
2. non confronta i file precedenti con quelli appena generati;
3. usa `--check` soltanto per rendere bloccante il test di convergenza delle sottofinestre.

Nel clone pulito il comando è uscito con codice 0 ma ha lasciato modificati JSON e TeX per il timestamp e i tempi macchina. Il digest scientifico è stabile, quindi il contenuto è buono; il flag non è un vero check idempotente. Correzione consigliata: generare in memoria/file temporanei, confrontare il contenuto riproducibile e non considerare timestamp/timing nel diff bloccante.

### 15.3 Il manifest non esegue gli audit bibliografici

`make_provenance.py` include `audit_refs.py` e `audit_by_title.py` nel manifest, ma ne calcola soltanto l'hash: non li esegue. Quindi la presenza delle due righe nel manifest non certifica che la bibliografia abbia superato un controllo.

Inoltre `audit_refs.py` confronta soltanto il **primo numero** del campo pagine. La voce `GiannoniPiccioneTausk2002` contiene `pages={697}` mentre il record ufficiale è `697--724`; lo script interrogato online restituisce comunque “0 entries with a discrepancy”. Il controllo dichiara di confrontare il page range, ma non lo fa.

`audit_by_title.py` restituisce una lista vuota tanto per “nessun match” quanto per errore di rete e non propaga un exit status informativo. Deve distinguere almeno `match`, `no match`, `registry unavailable` e fallire se il registro è indisponibile per tutte le voci.

Il run effettivo con accesso di rete ha riportato 5 casi “flagged” ma ha comunque restituito exit code 0. Almeno due sono falsi abbinamenti del “best match” Crossref: `Chen1977` viene confrontato con un capitolo del 2001 e `Timofeev1978` con un articolo del 1974. Questo conferma che il risultato va sottoposto a verifica umana e che il matching basato essenzialmente sulla prima parola del titolo è troppo debole.

Il run completo di `audit_refs.py` ha segnalato 20/56 voci, ma la maggior parte sono falsi positivi dovuti a iniziali, accenti, `and others` o DOI Zenodo/ResearchGate non risolti da Crossref. Un errore reale è però emerso: `Kovner1990` usa `author={Kovner, Isaac}`, mentre il PDF originale e il registro indicano **Israel Kovner**. Questo va corretto.

### 15.4 Il CAP completo passa da clone pulito

Nel clone pulito ho prima eseguito il self-test di `cap_full.py`, con esito PASS, e poi il runner completo:

`python3 no_inversion_schwarzschild_CAP_r0_10.py`

Esito finale:

`COMPLETE r0=10.0[2.6,6.05] nc=720 S=103 M=57 (753s)`  
`RESULT: True`

Questa verifica risponde direttamente all'obiezione del primo referee sul modulo assente nella vecchia release: nella versione corrente il runner rappresentativo è dependency-complete ed eseguibile da clone pulito. Resta indispensabile inserire nel pacchetto Zenodo/CQG finale lo stesso commit che è stato verificato.

### 15.5 Wolfram non è stato considerato verificato in questa sessione

`wolframscript` è installato, ma perfino una chiamata minimale al kernel non ha prodotto output entro il timeout. Non ho quindi finto che gli script `.wls` fossero stati rieseguiti. Ho trascritto in SymPy le identità ad alto rischio qui sopra e ho ottenuto risultati coerenti. Prima del freeze finale è opportuno far girare la suite Wolfram su una macchina con kernel/licenza funzionanti e conservarne un log.

## 16. Audit bibliografico, DOI e pacchetto di resubmission

### 16.1 `\nocite{*}` va rimosso

`paper2/paper2.tex:3756` contiene `\nocite{*}`. Il `.bbl` ha 99 voci, mentre 74 sono citate esplicitamente; 25 vengono inserite soltanto dallo star. Questo rende la bibliografia più lunga, meno controllabile e include lavori non usati nel ragionamento. Per una resubmission già criticata sulla bibliografia è preferibile rimuovere `\nocite{*}` e conservare solo le voci realmente richiamate.

### 16.2 Metadati bibliografici da correggere

In `paper/refs.bib`:

- `GiannoniPiccioneTausk2002`: `697` deve essere `697--724` ([record ufficiale AIMS](https://www.aimsciences.org/article/doi/10.3934/dcds.2002.8.697));
- `GiannoniPiccione2002`: `375` deve essere `375--423` (il PDF locale ha 49 pagine e l'ultima è numerata 423).
- `Kovner1990`: il nome `Isaac` deve essere `Israel`.

Il secondo referee attribuisce inoltre erroneamente l'articolo arrival-time a tre autori e a CQG; la response fa bene a non copiare quell'errore, ma la voce BibTeX deve comunque essere completa.

### 16.3 Paper I su Zenodo non è la revisione corrente

Paper II e la cover letter puntano al preprint Paper I [`10.5281/zenodo.21781850`](https://zenodo.org/records/21781850). Il record pubblico corrente contiene un solo `paper1.pdf` da 672441 byte, pubblicato il 3 agosto, con abstract ancora basato su “plunge phenomenology” e “dynamical horizon”. Il PDF locale aggiornato è 843501 byte e 58 pagine. La risposta dice al referee che Paper I è stato rivisto, ma il DOI porta ancora a una versione anteriore.

Prima della resubmission:

1. pubblicare una nuova versione Zenodo del preprint Paper I, mantenendo lo stesso concept DOI se possibile;
2. aggiornare `paper/companionI.bib` al nuovo version DOI o usare consapevolmente il concept DOI;
3. aggiornare abstract e metadati Zenodo con il protocollo E/C/A e senza termini fisici ritirati;
4. sincronizzare la data availability di Paper I con la [release codice v1.5.0](https://zenodo.org/records/22079449), DOI `10.5281/zenodo.22079449`, perché il sorgente aggiornato è presente lì mentre il testo cita ancora v1.2.0.

### 16.4 Il pacchetto `submission_CQG_R1` è obsoleto

Il sorgente corrente ha 3760 righe e il PDF 79 pagine. `paper2/submission_CQG_R1/paper2.tex` ha 3130 righe e il PDF 67 pagine: differenza di 630 righe. In particolare conserva una fase intermedia in cui la metrica di arrivo e la curvatura sono attribuite alla $\tau$-branch. **Non caricare questo pacchetto.**

Il PDF pulito corrente e quello highlighted hanno entrambi 79 pagine, che è un buon segnale, ma il pacchetto finale va rigenerato dal medesimo commit solo dopo le correzioni P0. Dovrà includere:

- source TeX pulito e `.bbl` definitivo;
- tutte le 22 figure effettivamente usate;
- frammento di provenance rigenerato;
- clean PDF e highlighted PDF della stessa revisione;
- response aggiornata;
- file di stile necessari, senza copie concorrenti obsolete.

### 16.5 Layout della response

L'ispezione visiva della pagina 5 conferma l'`Overfull \hbox` da circa 106.6 pt: la riga con lo spin scan invade il margine destro. Spezzare la lista su una displayed equation o su due righe. Paper II non ha errori LaTeX fatali o riferimenti indefiniti, ma conserva overfull fino a circa 23 pt, soprattutto nelle tabelle di script e in formule lunghe.

## 17. Correzioni prioritarie da affidare a Claude

### P0 — obbligatorie prima dell'upload

1. Cambiare l'“iff” di `eq:noconj-full` in una condizione sufficiente e qualificare il caso non radiale.
2. Limitare la salita tridimensionale alla non-coniugazione locale su archi a singola escursione; lasciare minimalità globale, HJB e winding classes nel piano equatoriale.
3. Eliminare da Paper I la conclusione di minimalità basata sulla famiglia numerica testata e sostituirla con “finite-dimensional numerical consistency check”.
4. Applicare davvero il protocollo E/C/A a entrambi i paper, comprese didascalie, abstract, appendici e conclusioni.
5. Rimuovere da Paper I il falso “Fuglede--Ishihara dichotomy” e usare il testo suggerito nella sezione 9.
6. Correggere nella response sia “T2 supplies the existence of a conformal factor” sia il fattore due della convenzione del solitone.
7. Pubblicare/collegare la versione aggiornata di Paper I su Zenodo.
8. Rigenerare da zero il pacchetto CQG finale: non usare `submission_CQG_R1`.

### P1 — fortemente consigliate

1. Correggere il fattore due nei coefficienti del bracket di $R_{\rm opt}$.
2. Scrivere $|\varepsilon|$ nel thrust bound o assumere esplicitamente espansione.
3. Sostituire l'“iff” puntuale della derivata della dilatazione con “identically along the fibre”.
4. Correggere la frase secondo cui il programma vive nel complemento del conformal-Killing framework.
5. Indebolire la frase universale sulla light-convex function.
6. Uniformare la terminologia higher-genus allo stato congetturale dichiarato.
7. Rimuovere il doppione McVittie.
8. Rimuovere `\nocite{*}`, completare i page range e rieseguire un audit bibliografico fail-closed.
9. Correggere gli overfull principali e “Protocol 3and”.

### P2 — manutenzione, subito dopo la resubmission

1. Rendere `make_provenance.py --check` idempotente e semanticamente fedele al nome.
2. Far sì che gli audit bibliografici siano eseguiti, non solo hashati.
3. Dichiarare in un README unico quale sorgente Paper I sia canonico.
4. Sincronizzare `POST_REVIEW_CHANGES.md`, `WORKING_METHOD.md` e gli handoff, per evitare che una sessione successiva reintroduca Fuglede o la convenzione errata del solitone.

## 18. Verdetto aggiornato su Claude e sulla probabilità di una resubmission positiva

Claude **non ha lavorato superficialmente sul cuore matematico**. I risultati più importanti richiesti dal primo referee — dominio, quattro casi, separatrice esterna, gerarchia del fixed-endpoint result, sorgente adiabatica completa e provenance — sono stati affrontati con contenuto reale. In particolare, la separatrice esterna e il residuo $O(\varepsilon^2)$ resistono a controlli indipendenti.

Gli errori residui hanno però un pattern preciso: sono soprattutto **sovra-formulazioni e propagazioni incomplete**. La response spesso descrive una correzione più globale di quella effettivamente presente nei sorgenti; la nuova sezione tridimensionale passa da un corretto argomento locale a una conclusione globale non autorizzata; i documenti e i DOI non sono ancora tutti sincronizzati.

Il cuore di originalità resta difendibile:

- T2 non costruisce il rail, il selettore, il support Hamiltonian o l'equazione (9);
- l'analogia con le sommersioni conformi è reale ma riguarda una sola componente geometrica;
- la parte nuova è il controllo attivo di una carica non conservata, il dominio causale della relativa indicatrice e la chiusura della risposta adiabatica on-shell/off-shell nei casi non stazionari.

Una nuova major review può diventare positiva se la resubmission non offre al referee appigli semplici e verificabili. I più pericolosi, oggi, sono la falsa equivalenza in `eq:noconj-full`, l'estensione globale 3D, la frase di minimalità numerica in Paper I, il lessico interno ancora presente e il DOI Paper I obsoleto. Sono correzioni localizzate: non richiedono di rifare la teoria, ma richiedono un ultimo passaggio metodico e un nuovo freeze completo.
