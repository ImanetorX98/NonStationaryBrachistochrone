# Piano di validazione sperimentale e analogica per Paper I e Paper II

**Data:** 25 agosto 2026  
**Autore del piano:** Codex  
**Destinatario operativo:** Claude, per valutare una risposta rigorosa alla richiesta sperimentale del referee  
**Manoscritti di riferimento:** `paper1/paper1_JMP.tex`, `paper2/paper2.tex`, `paper2/response_to_referees_CQG.tex`  
**Intervento eseguito sui paper:** nessuno.

## 1. Risposta breve

Sì: il lavoro ammette protocolli sperimentali falsificabili. Non esiste però un unico esperimento che riproduca contemporaneamente:

- la geometria esatta di Thakurta--Kerr o di Vaidya;
- una particella massive relativistica;
- il rail attivamente mantenuto $-g(u,W)=\hat E$;
- i due clock $t$ e $\tau$;
- la correzione adiabatica on-shell più off-shell;
- la separatrice e la topologia delle winding classes.

La strategia scientificamente corretta è una **validazione a strati**:

1. **emulazione dinamica esatta** dell'Hamiltoniana e del controllo, su un banco robotico o hardware-in-the-loop;
2. **analogo idrodinamico rotante**, con un vortice d'acqua o di elio, per frame dragging, light rings, separatrici e variazione lenta del background;
3. **analogo non stazionario in un fluido quantistico**, per una metrica efficace con fattore temporale controllato;
4. solo in futuro, **test astrofisico**, se esistesse un modello osservativo realistico di un oggetto compatto Thakurta--Kerr.

I primi tre testano aspetti matematici o universali dell'Hamiltoniana efficace. Non testano direttamente le equazioni di Einstein, l'esistenza fisica della soluzione Thakurta--Kerr o l'interpretazione astrofisica dell'oggetto centrale.

## 2. Perché la risposta attuale al referee è insufficiente ma non sbagliata

La risposta corrente dice che non è disponibile una validazione sperimentale diretta e sostituisce tale livello con:

- identità simboliche;
- confronti numerici indipendenti;
- recupero di limiti noti;
- risultati computer-assisted.

Questa è una posizione onesta. Tuttavia il referee aveva chiesto le “techniques of experimental validation”. È possibile soddisfare meglio la richiesta senza fingere un contatto osservativo: basta presentare **tecniche realizzabili** e indicare esattamente quale parte del lavoro ciascuna potrebbe falsificare.

La frase “no experimental validation is available” dovrebbe diventare:

> No direct gravitational or astrophysical validation is presently available, and none is claimed. The effective controlled Hamiltonian nevertheless admits laboratory emulation protocols, described below, which can test its optimal-control structure, its adiabatic response and its separatrix phenomenology without claiming to reproduce the Einstein dynamics of a Thakurta--Kerr compact object.

## 3. I tre significati di “esperimento” da non confondere

### Livello A — emulazione esatta delle equazioni

Un attuatore fisico viene comandato affinché la velocità misurata obbedisca alla stessa indicatrice

\[
\dot x=c(x,\chi)+R(x,\chi)
\begin{pmatrix}
\cos\theta\\
\sin\theta
\end{pmatrix}.
\]

Questo può testare:

- l'unicità del support maximizer;
- la legge PMP;
- le traiettorie a endpoint fissati;
- la minimalità entro una classe misurabile;
- la risposta a una variazione lenta;
- la necessità del termine off-shell;
- la legge di scala del residuo.

È una validazione fisica dell'algoritmo e del sistema dinamico ridotto, ma non una verifica che la natura realizzi spontaneamente la metrica di partenza.

### Livello B — analogo con la stessa classe di Hamiltoniane

Un fluido, un condensato o un mezzo ottico produce un'equazione di dispersione i cui raggi o quasiparticelle vedono una metrica efficace. Non è necessario che la metrica sia identica in tutto lo spaziotempo: può essere sufficiente che l'Hamiltoniana congelata e la sua derivata rispetto al parametro lento coincidano, entro la precisione sperimentale, sul dominio attraversato.

Questo testa l'universalità del meccanismo: frame dragging, turning points, separatrix crossing, perdita di adiabaticità e risposta di primo ordine.

### Livello C — validazione gravitazionale o astrofisica

Richiederebbe osservabili provenienti da un oggetto fisico descritto dalla metrica di Paper II, insieme a un modello della materia, degli orizzonti e dell'emissione. Oggi non è disponibile. Non va promessa nella resubmission.

## 4. Matrice delle piattaforme

| Piattaforma | Cosa testa bene | Come si rende non stazionaria | Limite principale | Priorità |
|---|---|---|---|---:|
| Robot planare/gantry hardware-in-the-loop | indicatrice, PMP, HJB numerica, endpoint, winding classes, off-shell, slope 2 | parametri $A(\chi)$ o $m(\chi)$ programmati e misurati | implementa le equazioni, non la gravità | **1** |
| Floater o robot in vasca con vortice | navigazione attiva, vento Zermelo, pro/retro, separatrice | rampa di pompa e circolazione | il campo di flusso non è Kerr esatto | **2** |
| Onde superficiali in draining vortex | metrica acustica, light rings, superradianza, turning points | rampa lenta della pompa/circolazione | onde null-like; non testano il rail massive | **2** |
| BEC toroidale/2D | metrica espandente, fattore temporale, wave packets, nonadiabaticità | trap o interazione atomica variabili | apparato specialistico; fononi quasi-null | **3** |
| Vortice gigante in He II | geometria rotante con bassa viscosità, bound modes e ringdown | modulazione lenta del propulsore | criogenia e collaborazione necessarie | **3** |
| Polariton BEC con difetti topologici | candidato per quasiparticelle massive Kerr-like | profilo di pump variabile | piattaforma proposta, non protocollo già validato per questo rail | **3** |
| Mezzo ottico temporalmente modulato | Hamiltoniane di raggio, indice $n(x,t)$ | modulazione elettro-ottica o laser | una pura conformal rescaling è invisibile ai raggi nulli ideali | **4** |

## 5. Esperimento consigliato 1: banco planare del rail controllato

### 5.1 Perché è la prima scelta

È il solo protocollo relativamente economico che possa testare **direttamente** il contenuto distintivo dei paper: non un raggio libero, ma un moto attivamente controllato con indicatrice compatta e costo di arrivo.

Il banco può essere realizzato con una delle seguenti opzioni:

- robot olonomo planare con motion capture dall'alto;
- piccolo carrello su tavola ad aria;
- magnete o floater pilotato da bobine sotto una superficie piana;
- piattaforma XY con end-effector e sensore di forza;
- simulatore real-time con attuatore fisico e chiusura di feedback.

La soluzione più pulita è un robot olonomo o una piattaforma XY, perché può realizzare una matrice $R$ anisotropa e non soltanto una velocità di modulo fisso.

### 5.2 Mappa matematica

Si scelga una regione anulare di laboratorio e una scala $(L_0,T_0)$. La posizione misurata $x_{\rm lab}$ viene trasformata nelle coordinate adimensionali del paper. A ogni ciclo di controllo si comanda

\[
\dot x_{\rm cmd}=s_c\left[c(x,\chi)+R(x,\chi)e(\theta)\right],
\qquad
e(\theta)=(\cos\theta,\sin\theta),
\]

dove $s_c=L_0/T_0$. Il controllo reale non deve usare la velocità comandata come dato: la velocità effettiva va ricostruita dal tracking e inserita nell'analisi.

Il banco deve avere due modalità:

1. **open loop**, per misurare errori del modello e instabilità;
2. **closed loop**, per emulare il rail attivamente mantenuto e misurare lo sforzo del feedback.

### 5.3 Test A: ricostruzione dell'indicatrice

Per una griglia di punti $(x_i,\chi_i)$:

1. applicare 32--64 direzioni di controllo $\theta$;
2. misurare $\dot x_{\rm lab}$ dopo il transitorio;
3. fittare centro, semiassi e orientazione dell'ellisse;
4. confrontare con $c$ e $R$ senza rifittare le traiettorie successive.

Osservabili:

- errore ortogonale medio dall'ellisse prevista;
- eccentricità e orientazione;
- area dell'indicatrice;
- distanza dal freezing parameter $\delta$.

Predizione caratteristica:

\[
\text{semiassi}\propto\sqrt{\delta}
\]

quando la superficie di freezing viene avvicinata dal dominio regolare.

### 5.4 Test B: massimizzatore PMP e Hamiltoniana

Per costati di prova $p$ non nulli, calcolare sperimentalmente

\[
\theta_{\rm exp}^*=\arg\max_{\theta}\,p\cdot\dot x_{\rm lab}(\theta).
\]

Confrontare:

- unicità del massimo;
- continuità di $\theta^*$ al variare di $p$;
- valore misurato della funzione supporto;
- formula hamiltoniana chiusa del paper.

Questo è un test particolarmente utile dopo la correzione della citazione a Filippov: la regolarità del massimizzatore viene verificata direttamente sulla geometria effettivamente realizzata.

### 5.5 Test C: winding classes e minimalità equatoriale

La regione di laboratorio deve essere anulare, con un ostacolo centrale. Si selezionano start e target e si costruiscono traiettorie con winding $k=-2,-1,0,1,2$.

Per ogni classe:

1. inizializzare molti controlli diversi;
2. eseguire il boundary-value optimizer;
3. realizzare fisicamente le traiettorie candidate;
4. misurare il costo di arrivo;
5. confrontare con una soluzione HJB o level-set calcolata indipendentemente.

Questo test è precisamente bidimensionale. Non deve essere presentato come validazione della minimalità globale tridimensionale.

### 5.6 Test D: due clock sulla stessa indicatrice

Il banco realizza un'unica dinamica spaziale, mentre il data logger accumula due costi:

\[
T_t=\int \ell_t\,d\chi,
\qquad
T_\tau=\int \ell_\tau\,d\chi.
\]

Si ottimizzano separatamente i due costi con gli stessi endpoint e la stessa indicatrice. Questo può verificare:

- che i rami condividano il vincolo cinematico;
- che costi diversi selezionino estremali diversi;
- l'eventuale inversione dell'ordinamento sotto il protocollo di endpoint dichiarato.

Il secondo accumulatore non è il tempo proprio fisico del robot. È un **cost functional analogico**. Il paper o la risposta devono dirlo.

### 5.7 Test E: correzione adiabatica completa e termine off-shell

Questo è il test più discriminante dell'intero programma.

Si imponga una famiglia lenta, per esempio

\[
A(\chi)=A_0\exp(\varepsilon\chi)
\]

nel caso conformal-Kerr emulato, oppure una famiglia $m(v)$ nel caso Vaidya emulato. Per ciascun $\varepsilon$ si misuri una quantità finale $Q$, preferibilmente:

- angolo totale $\Delta\varphi$;
- turning radius;
- arrival cost;
- costato finale ricostruito.

Si confrontino tre modelli **senza rifit individuale**:

\[
Q_0,
\qquad
Q_0+\varepsilon Q_{\rm on},
\qquad
Q_0+\varepsilon(Q_{\rm on}+Q_{\rm off}).
\]

La firma prevista è:

\[
R_{\rm on}(\varepsilon)\sim\varepsilon,
\qquad
R_{\rm full}(\varepsilon)\sim\varepsilon^2.
\]

Il confronto fra pendenza uno e pendenza due è più convincente del semplice accordo a un solo valore di $\varepsilon$.

### 5.8 Test F: costo del controllo vicino al freezing

Se il banco misura forza, accelerazione o corrente di attuazione, avvicinarsi alla superficie di freezing da valori regolari e testare

\[
|a|_{\min}\propto\frac{\varepsilon}{\hat E\sqrt\delta}.
\]

Il test deve essere interrotto prima della saturazione dell'attuatore. La saturazione va trattata come censura sperimentale, non come conferma della divergenza.

Anche qui la grandezza è l'analogo meccanico del costo di controllo, non un'accelerazione propria gravitazionale misurata direttamente.

## 6. Esperimento consigliato 2: robot o floater in una vasca rotante

### 6.1 Obiettivo

Combinare una dinamica fisica di advezione con un agente che controlla soltanto la direzione o una piccola anisotropia di propulsione. Il flusso fornisce il centro $c$ dell'indicatrice; la propulsione fornisce la parte $R e(\theta)$.

Una vasca anulare con drain e circolazione consente di studiare:

- differenza prograde/retrograde;
- regioni non raggiungibili;
- turning points;
- separatrice tra famiglie di traiettorie;
- winding classes;
- variazione lenta del campo di background.

La fattibilità teorica e tecnologica della navigazione ottima di agenti attivi in campi fluidi, anche non autonomi, è sostenuta dalla letteratura su microswimmers e Zermelo navigation. Daddi-Moussa-Ider, Löwen e Liebchen formulano esplicitamente un problema non autonomo a endpoint fissati e tempo finale libero trattato con PMP e indicano i colloidi attivi programmabili come possibile banco di prova. Gunnarson e collaboratori studiano computazionalmente la navigazione in flussi vorticosi non stazionari e definiscono un protocollo adatto ad agenti robotici dotati di sensori locali. Nessuno dei due lavori costituisce già una realizzazione sperimentale del rail relativistico.

### 6.2 Calibrazione obbligatoria

Usare particle image velocimetry (PIV) o particle tracking velocimetry per misurare

\[
u_r(r,\varphi,t),\qquad u_\varphi(r,\varphi,t)
\]

invece di assumere il profilo ideale $(-D/r,C/r)$. La mappa verso i parametri del paper deve essere costruita dopo questa misura.

Il comando alla pompa non coincide con il parametro fisico: ritardi, isteresi e ricircolo fanno sì che il vero parametro lento sia una quantità inferita dal campo misurato.

### 6.3 Non stazionarietà

Sono possibili tre protocolli:

- rampa lenta della portata del drain $D(t)$;
- rampa lenta della circolazione $C(t)$;
- modulazione lenta della velocità propria dell'agente.

Per testare la struttura di Paper II, è preferibile mantenere il rapporto spaziale del profilo quanto più invariato possibile e cambiare un solo parametro scalare. Se cambiano simultaneamente forma, vorticità e profondità, non esiste più una singola derivata $\partial_\lambda H$ confrontabile con il termine teorico.

### 6.4 Cosa si può rivendicare

Si può rivendicare:

> an active-particle analogue of the controlled optical Hamiltonian in a measured rotating flow.

Non si può rivendicare senza un inverse-design quantitativo:

> an experimental realization of the Thakurta--Kerr spacetime.

## 7. Esperimento consigliato 3: onde su un draining vortex

### 7.1 Stato dell'arte

Le onde superficiali su acqua in moto realizzano una metrica acustica nel regime di lunghe lunghezze d'onda. Sono già stati osservati:

- conversione presso un orizzonte idrodinamico;
- superradianza rotazionale;
- light-ring/ray phenomenology;
- quasinormal ringing;
- una risposta globale lenta del vortice descrivibile con una metrica dinamica.

Questa base sperimentale è molto più vicina a Paper II di un generico “metamateriale”.

### 7.2 Protocollo stazionario

1. Generare un draining vortex stabile.
2. Misurare il campo di velocità con PIV.
3. Eccitare pacchetti narrow-band prograde e retrograde.
4. Ricostruire la superficie con profilometria o imaging ad alta velocità.
5. Estrarre raggi di gruppo, turning radii, frequenze e numeri azimutali.
6. Confrontare la separazione prograde/retrograde e i light rings con l'Hamiltoniana acustica calibrata.

Questo test sostiene il settore Kerr-like e la fenomenologia delle separatrici, non il ramo massive del rail.

### 7.3 Protocollo non stazionario

Ripetere con una rampa lenta della circolazione o della portata. Il parametro adiabatico sperimentale deve essere definito come

\[
\varepsilon_{\rm lab}
=
\frac{T_{\rm orbit}}{T_{\rm ramp}},
\]

o mediante la derivata logaritmica del parametro calibrato, non semplicemente tramite il comando al motore.

Osservabili:

- spostamento del turning radius;
- deriva della frequenza Doppler;
- spostamento dei light rings;
- breakdown della frozen approximation;
- ritardo/isteresi attraversando la separatrice.

Vicino alla separatrice $T_{\rm orbit}$ diverge: nessuna rampa finita è uniformemente adiabatica. Questo rende sperimentalmente visibile proprio il limite della teoria esterna e la necessità di un inner problem di Neishtadt.

### 7.4 Limite fondamentale

Le onde sono excitations null-like della metrica acustica. Una moltiplicazione puramente conforme della metrica non cambia le null geodesics ideali. Per osservare un effetto attribuito al fattore $A(t)$ occorre almeno una delle seguenti condizioni:

- dispersione controllata a frequenza finita;
- un modo con massa efficace;
- una quasiparticella massive;
- un agente attivo;
- una modulazione che cambi non soltanto il prefattore conforme ma anche sound speed o flow profile misurabili.

Per questo il vortice d'acqua non è il test diretto del termine off-shell del rail massive, salvo una mappa hamiltoniana aggiuntiva dimostrata esplicitamente.

## 8. Esperimento consigliato 4: vortice gigante in elio superfluido

Švančara e collaboratori hanno stabilizzato un vortice gigante in He II, misurato onde micrometriche sulla superficie ed estratto bound states e firme di ringdown in una geometria rotante efficace. La viscosità ridotta e il grande numero di quanta di circolazione rendono questa piattaforma la più vicina oggi a una geometria rotante analogica pulita.

Estensione proposta:

1. mantenere una configurazione stazionaria per calibrare il campo;
2. modulare lentamente la frequenza del propulsore entro il regime stabile;
3. misurare in tempo reale circolazione, forma della superficie e spettro $(\omega,m)$;
4. ricostruire la posizione delle barriere/ring states come funzione del tempo;
5. confrontare frozen prediction e correzione di primo ordine;
6. attraversare deliberatamente un valore critico e misurare il breakdown adiabatico.

Questo richiede una collaborazione con un laboratorio criogenico. Non è un primo esperimento da costruire autonomamente, ma è il miglior obiettivo scientifico esterno per Paper II.

## 9. Esperimento consigliato 5: BEC espandente o fluido quantistico

### 9.1 Perché è rilevante

Esperimenti esistenti hanno già:

- realizzato condensati toroidali in espansione con redshift dei fononi;
- misurato attenuazione/amplificazione di Hubble;
- costruito in un BEC 2D metriche efficaci con curvatura spaziale e temporale controllata;
- osservato particle-pair production durante espansioni programmate.

Questi risultati dimostrano che un background adiabaticamente o rapidamente variabile non è soltanto un'idea teorica.

### 9.2 Cosa testare

Un BEC toroidale o 2D può testare bene il settore FLRW/conformal-time di Paper I:

- dipendenza dalla storia $A(t)$;
- redshift e deriva dell'energia efficace;
- confronto frozen/slowly varying;
- transizione da regime adiabatico a nonadiabatico;
- simmetria fra espansione e contrazione.

Per avvicinarsi a Paper II occorrerebbe aggiungere una circolazione quantizzata o un pump con momento angolare e misurare modi co- e contro-rotanti.

### 9.3 Audit del kymograph BEC attualmente generato

Lo script `Experimental/ring_bec_poc_codex.py` usa attualmente il modello

\[
q_m''+\gamma q_m'+\left(\frac{m c_s}{R(t)}\right)^2q_m=0,
\qquad
\delta n_{\rm toy}(\theta,t)=
\sum_m q_m(t)\cos(m\theta+\phi_m).
\]

Il kymograph risultante è numericamente coerente con questa equazione ed è utile come illustrazione deterministica di interferenza fra modi stazionari e redshift temporale. Non è però ancora corretto chiamarlo previsione quantitativa della densità di un BEC toroidale reale.

Il modello sperimentale di Banik et al. evolve innanzitutto la perturbazione di fase:

\[
\ddot{\delta\phi}_m+
\left[2\gamma(t)+\gamma_H\frac{\dot R}{R}\right]
\dot{\delta\phi}_m+
\omega_m^2(t)\delta\phi_m=0,
\]

con

\[
\omega_m(t)=\frac{m c_\theta(t)}{R(t)},\qquad
c_\theta(t)=c_{\theta,i}
\left(\frac{R(t)}{R_i}\right)^{-\alpha/2}.
\]

Nel limite Thomas--Fermi di anello sottile, il volume efficace scala come $V\propto R^\alpha$ e la previsione ideale è $\gamma_H=\alpha$. Questa uguaglianza non va però imposta ai dati: Banik et al. misurano uno scostamento rispetto alla previsione ideale e trattano $\gamma_H$ come parametro sperimentale.

La densità misurata non coincide con $\delta\phi_m$. Dalla relazione di Gross--Pitaevskii usata nell'analisi sperimentale,

\[
\partial_t\delta\phi_m
=-\frac{g_{\rm eff}}{\hbar}\frac{\delta n_m}{R^\alpha},
\]

segue

\[
\delta n_m(t)=-\frac{\hbar}{g_{\rm eff}}
R^\alpha(t)\dot{\delta\phi}_m(t),
\qquad
\delta n_{1D}(\theta,t)=
\sum_m\delta n_m(t)\sin(m\theta+\delta\theta_m).
\]

Lo script corrente deve quindi essere descritto come **toy-model phonon kymograph** finché non saranno introdotti:

- il termine di Hubble $\gamma_H\dot R/R$;
- la dipendenza di $c_\theta(t)$ dal raggio e dalla densità;
- la ricostruzione della densità da $R^\alpha\dot{\delta\phi}$;
- una scelta esplicita fra densità assoluta $\delta n_{1D}$ e contrasto normalizzato $\delta n_{1D}/n_{0,1D}$;
- parametri dimensionali e una pipeline sintetica compatibile con l'imaging.

Per i parametri illustrativi correnti,

\[
\max_t\left|\frac{\dot R}{R}\right|=0.0522857,
\]

e il parametro di lentezza è

\[
\max_t\frac{|\dot R/R|}{\omega_m}
=
\begin{cases}
0.0100, & m=6,\\
0.00667, & m=9.
\end{cases}
\]

Il kymograph è dunque nel regime lentamente variabile del modello semplificato. Questi numeri vanno ricalcolati dopo aver introdotto $c_\theta(t)$ e la calibrazione in unità fisiche.

I modi $m=6$ e $m=9$ sono stati scelti per rendere leggibile l'interferenza, non perché siano già validati per un apparato. Il primo confronto realistico dovrebbe usare il modo quasi puro $m=1$ dell'esperimento di Banik et al. Modi superiori sono ammissibili soltanto se

\[
k_m\xi=\frac{m\xi}{R}\ll1,
\]

se restano separati dalle eccitazioni radiali e se l'ampiezza rimane nel regime lineare.

Poiché l'ascissa del kymograph è l'angolo comovente $\theta$, il numero d'onda angolare $m$ non cambia durante l'espansione. Il redshift deve apparire come rallentamento dell'oscillazione nel tempo; la lunghezza d'onda fisica $2\pi R(t)/m$ cresce invece con il raggio. L'assenza di un allargamento orizzontale nel grafico in coordinate $\theta$ non è quindi un errore.

Le equazioni e la procedura sperimentale sono documentate nelle fonti locali [Eckel et al. 2018](Fonti/EckelEtAl-2018-expanding-BEC-universe.pdf) e [Banik et al. 2022](Fonti/BanikEtAl-2022-Hubble-attenuation-amplification-BEC.pdf), rispettivamente [DOI 10.1103/PhysRevX.8.021021](https://doi.org/10.1103/PhysRevX.8.021021) e [DOI 10.1103/PhysRevLett.128.090401](https://doi.org/10.1103/PhysRevLett.128.090401).

### 9.4 Protocollo minimo per un kymograph confrontabile con i dati

1. Misurare o prescrivere $R(t)$ e propagare la relativa incertezza.
2. Calibrare $c_\theta(R)$ su anelli stazionari e stimare $\alpha$.
3. Preparare inizialmente un solo modo $m=1$ a piccola ampiezza.
4. Integrare l'equazione della fase includendo $Q(R)$, $\gamma_H$ e la fase iniziale.
5. Ricostruire $\delta n_{1D}$, sottraendo il background senza fonone e integrando radialmente come nell'esperimento.
6. Produrre separatamente kymograph della fase, densità assoluta e contrasto normalizzato.
7. Confrontare espansione, anello congelato e contrazione; solo dopo aggiungere più modi o entrare nel regime nonadiabatico.

Questo protocollo testa redshift, attenuazione/amplificazione di Hubble e breakdown adiabatico del settore analogico. Non misura da solo il termine off-shell del rail massive di Paper I/II.

### 9.5 Quasiparticelle massive

La proposta di Solnyshkov e collaboratori usa difetti topologici in un polariton condensate come particelle massive in un analogo Kerr. Sarebbe la piattaforma concettualmente più vicina al rail timelike, ma non va presentata come una realizzazione sperimentale già pronta del presente protocollo.

Un progetto realistico richiederebbe:

- stabilizzazione e tracking del difetto;
- ricostruzione della metrica efficace;
- feedback per mantenere l'analogo di $\hat E$;
- pump profile lentamente variabile;
- confronto fra moto libero e moto controllato.

## 10. Perché l'ottica temporalmente modulata non è la prima scelta

Mezzi con indice $n(t)$ o $n(x,t)$ possono realizzare metriche efficaci tempo-dipendenti. Sono utili per testare Hamiltoniane di raggio, mixing modale e particle creation. Tuttavia:

- i fotoni sono nulli;
- le null geodesics non vedono un puro fattore conforme ideale;
- dispersione e material response possono dominare il piccolo termine adiabatico cercato;
- un setup puramente ottico non misura naturalmente il costo di un rail attivo.

L'ottica diventa competitiva se si usa un waveguide cutoff o un band structure che dia massa efficace, oppure se si realizza direttamente il propagatore hamiltoniano con una rete fotonica programmabile. In quel caso si tratterebbe ancora di un'emulazione del livello A/B, non della gravità fisica.

## 11. Protocollo statistico per il test del termine off-shell

### 11.1 Disegno preregistrato

Prima dei dati finali fissare:

- osservabile primaria $Q$;
- intervallo di $\varepsilon$;
- dominio regolare escluso da turning point, freezing e separatrice;
- modelli $Q_0$, $Q_{\rm on}$ e $Q_{\rm full}$;
- criteri di esclusione;
- metodo di fit;
- soglia del noise floor.

### 11.2 Campionamento suggerito

Come pilot:

- 7--10 valori di $\varepsilon$ approssimativamente log-spaced su almeno un decennio;
- almeno 15--30 repliche per valore se il rumore meccanico è rilevante;
- 3--5 condizioni iniziali tenute fuori dalla calibrazione;
- rampe sia positive sia negative;
- ripetizione stazionaria $\varepsilon=0$ a inizio, metà e fine sessione.

I numeri finali devono essere scelti da una power analysis basata sul pilot, non copiati meccanicamente da questo documento.

### 11.3 Fit

Per ciascun modello calcolare

\[
R(\varepsilon)=|Q_{\rm obs}-Q_{\rm model}|.
\]

Fittare

\[
\log R=\alpha+\beta\log\varepsilon
\]

includendo:

- incertezza sull'asse $\varepsilon$ dovuta alla calibrazione del ramp;
- correlazione fra repliche che condividono la stessa calibrazione;
- bootstrap per gli intervalli su $\beta$;
- censura dei punti sotto il noise floor;
- confronto con modelli a pendenza fissata $\beta=1$ e $\beta=2$.

Esito atteso:

- on-shell-only incompatibile con pendenza due;
- full correction compatibile con pendenza due;
- stationary/frozen model con errore di ordine inferiore.

Il valore $2.00$ o $2.12$ del software non deve essere imposto al fit sperimentale.

## 12. Controlli negativi e simmetrie

1. **Frozen control:** $A'=0$ o parametro di flusso costante.
2. **Ablation:** rimuovere deliberatamente il termine off-shell dal predittore.
3. **Ramp reversal:** $\varepsilon\to-\varepsilon$.
4. **Rotation reversal:** $a\to-a$, $J\to-J$ o $C\to-C$.
5. **Endpoint swap:** verificare la dipendenza dal protocollo P1/P2 senza confonderli.
6. **Clock swap:** stessi dati cinematici, costi $t$ e $\tau$ accumulati separatamente.
7. **Blind initial conditions:** la previsione viene prodotta prima di rivelare i dati di tracking.
8. **Independent calibration:** il campo di flusso o l'indicatrice vengono misurati in sessioni diverse da quelle usate per la traiettoria.
9. **Noise injection:** latenza, rumore di posizione e errore del parametro lento vengono variati per testare la robustezza della slope.
10. **Boundary control:** cambiare la dimensione della vasca o il riflettore esterno per separare effetti fisici da riflessioni finite.

## 13. Error budget minimo

Ogni pubblicazione sperimentale dovrebbe riportare:

- errore di posizione e velocità;
- latenza e bandwidth del feedback;
- accuratezza della matrice $R$ e del centro $c$;
- deriva termica/meccanica;
- ripetibilità del parametro lento;
- viscosità e dispersione, per i fluidi;
- PIV calibration e spatial smoothing;
- effetto delle pareti e riflessioni;
- differenza fra comando al motore e campo realmente misurato;
- saturazione dell'attuatore vicino al freezing;
- propagazione delle incertezze fino a $\beta$.

## 14. Separatrice: protocollo distinto dal test adiabatico regolare

Non mescolare nello stesso fit:

- punti uniformemente lontani dalla separatrice, dove si testa $O(\varepsilon^2)$;
- punti che attraversano la separatrice, dove l'espansione regolare non è uniforme.

Per il crossing:

1. preparare condizioni iniziali con fase controllata;
2. far muovere lentamente il parametro critico attraverso la traiettoria;
3. misurare lobe/capture outcome, tempo di permanenza e phase jump;
4. ripetere per molte fasi iniziali;
5. confrontare con una teoria inner di separatrix crossing, non con la sola formula esterna.

Un mancato slope-two vicino alla separatrice non falsifica il teorema locale del paper; conferma invece che si è usciti dal suo dominio uniforme.

## 15. Che cosa ciascun risultato potrebbe falsificare

| Misura | Predizione falsificabile | Parte del lavoro |
|---|---|---|
| ellisse locale | centro/assi/orientazione e collasso $\sqrt\delta$ | rail reduction |
| support maximizer | unicità e formula hamiltoniana | PMP/indicatrice |
| costi a endpoint fissati | differenza fra rami e protocol dependence | $t$ vs $\tau$ |
| winding classes | un candidato globale per classe nel banco 2D | minimalità equatoriale |
| control effort | scaling $\delta^{-1/2}$ | thrust bound |
| residual slope | 1 on-shell-only, 2 full | termine off-shell |
| pro/retro threshold | split e separatrice esterna | Paper II rotante |
| ramp sign | cambio di segno della risposta lineare | derivata adiabatica |
| separatrix crossing | breakdown nonuniforme e phase dependence | limite della teoria esterna |

## 16. Cosa non deve essere dichiarato

Non scrivere:

- “the experiment validates the Thakurta--Kerr spacetime”;
- “water waves are massive controlled particles”;
- “the accumulated $\tau$ cost is the robot's proper time”;
- “a draining vortex is exactly Kerr”;
- “slope two proves the full gravitational theory”;
- “an experiment on null phonons tests a pure conformal factor”;
- “a finite set of trajectories proves global minimality in an infinite-dimensional class”.

Scrivere invece:

- “tests the reduced controlled Hamiltonian”;
- “laboratory emulation of the effective dynamics”;
- “analogue of the rotating optical geometry”;
- “falsifiable test of the first-order correction and its residual scaling”;
- “does not test Einstein's field equations or the astrophysical realization”.

## 17. Testo pronto per la risposta al referee

```tex
\reply We agree that the manuscript should state not only what has been
validated numerically, but also what could be tested experimentally.  No direct
astrophysical validation of the Thakurta--Kerr compact object is presently
available, and we claim none.  The reduced dynamics nevertheless admits two
falsifiable laboratory protocols.

First, the controlled rail can be emulated by a planar active agent whose measured
velocity indicatrix is calibrated to
$\dot x=c(x,\chi)+R(x,\chi)(\cos\theta,\sin\theta)$.  A slowly programmed parameter
$A(\chi)$ then permits a direct comparison of the frozen, on-shell-only and complete
first-order predictions.  The discriminating observable is the residual scaling:
the on-shell truncation predicts a remaining $O(\varepsilon)$ discrepancy, whereas
the complete on-shell plus off-shell correction predicts $O(\varepsilon^2)$ on
regular arcs.  The same platform can measure the collapse of the indicatrix and the
control effort near freezing, and can test separate winding classes in an annular
two-dimensional domain.

Second, rotating-fluid experiments provide an analogue test of the Kerr-like
sector.  Surface waves on draining vortices have already displayed rotational
superradiance, light-ring and ringdown phenomenology, while a giant quantum vortex
in superfluid helium has realized a low-viscosity rotating effective geometry.
Slowly ramping the measured circulation or drain rate would test frozen-orbit
drift, prograde/retrograde turning structures and the loss of uniform adiabaticity
at a moving separatrix.  Such wave experiments probe a null/acoustic analogue and
not the massive actively forced rail; this limitation is explicit.

These protocols would validate the effective Hamiltonian and its asymptotic
response, not Einstein's equations or the existence of an astrophysical
Thakurta--Kerr object.  We have added this distinction and the corresponding
experimental references to the outlook.
```

## 18. Testo breve per il corpo/conclusioni

```tex
\paragraph{Experimental outlook.}
No direct gravitational validation of the controlled Thakurta--Kerr rail is
presently available.  Its reduced Hamiltonian is nevertheless experimentally
falsifiable.  A planar active agent can be feedback-controlled so that its
measured velocity indicatrix realizes the affine oval of the rail; a slow ramp of
the calibrated background parameter then tests the distinctive residual scaling,
$O(\varepsilon)$ after the on-shell term alone and $O(\varepsilon^2)$ after the
complete off-shell correction.  Complementarily, rotating water or superfluid
vortices can test the Kerr-like prograde/retrograde and moving-separatrix
phenomenology using surface-wave packets.  The latter are null/acoustic analogues,
not massive forced worldlines, so neither protocol is claimed as a realization of
the Einstein dynamics of a Thakurta--Kerr compact object.
```

## 19. Fonti primarie consigliate

### Disponibilità locale verificata

Aggiornamento del 25 agosto 2026: i ventuno lavori elencati in questa
sezione sono ora consultabili direttamente nella cartella `Fonti/`. I PDF
sono stati controllati con `file` e `pdfinfo`; tutti risultano documenti PDF
validi e leggibili. Quando la versione editoriale non era liberamente
scaricabile, è stata archiviata la versione autore/arXiv open-access. Il tipo
di evidenza è indicato esplicitamente per non confondere un esperimento, una
proposta sperimentale e un lavoro puramente teorico.

| N. | Riferimento breve | File locale in `Fonti/` | Evidenza utile al progetto |
|---:|---|---|---|
| 1 | Unruh 1981 | [Unruh-1981-experimental-black-hole-evaporation.pdf](Fonti/Unruh-1981-experimental-black-hole-evaporation.pdf) | fondazione teorica dell'analogo acustico |
| 2 | Weinfurtner et al. 2011 | [WeinfurtnerEtAl-2011-stimulated-Hawking-emission.pdf](Fonti/WeinfurtnerEtAl-2011-stimulated-Hawking-emission.pdf) | esperimento in acqua |
| 3 | Torres et al. 2017 | [TorresEtAl-2017-rotational-superradiance-vortex.pdf](Fonti/TorresEtAl-2017-rotational-superradiance-vortex.pdf) | esperimento in vortice d'acqua rotante |
| 4 | Torres et al. 2018 | [TorresEtAl-2018-waves-vortex-rays-rings-resonances.pdf](Fonti/TorresEtAl-2018-waves-vortex-rays-rings-resonances.pdf) | teoria, ray tracing e confronto con dati di vasca |
| 5 | Torres et al. 2020 | [TorresEtAl-2020-quasinormal-modes-analogue-black-hole.pdf](Fonti/TorresEtAl-2020-quasinormal-modes-analogue-black-hole.pdf) | esperimento di ringdown in vortice |
| 6 | Patrick et al. 2021 | [PatrickEtAl-2021-backreaction-analogue-black-hole.pdf](Fonti/PatrickEtAl-2021-backreaction-analogue-black-hole.pdf) | esperimento di backreaction e metrica dinamica |
| 7 | Švančara et al. 2024 | [SvancaraEtAl-2024-giant-quantum-vortex-curved-spacetime.pdf](Fonti/SvancaraEtAl-2024-giant-quantum-vortex-curved-spacetime.pdf) | esperimento in elio superfluido |
| 8 | Eckel et al. 2018 | [EckelEtAl-2018-expanding-BEC-universe.pdf](Fonti/EckelEtAl-2018-expanding-BEC-universe.pdf) | esperimento BEC in espansione |
| 9 | Banik et al. 2022 | [BanikEtAl-2022-Hubble-attenuation-amplification-BEC.pdf](Fonti/BanikEtAl-2022-Hubble-attenuation-amplification-BEC.pdf) | esperimento BEC in espansione e contrazione |
| 10 | Viermann et al. 2022 | [ViermannEtAl-2022-quantum-field-simulator-curved-spacetime.pdf](Fonti/ViermannEtAl-2022-quantum-field-simulator-curved-spacetime.pdf) | simulatore quantistico 2D sperimentale |
| 11 | Daddi-Moussa-Ider et al. 2021 | [DaddiMoussaIderEtAl-2021-optimal-microswimmer-navigation.pdf](Fonti/DaddiMoussaIderEtAl-2021-optimal-microswimmer-navigation.pdf) | teoria di navigazione/PMP con proposta di test |
| 12 | Gunnarson et al. 2021 | [GunnarsonEtAl-2021-learning-navigation-vortical-flows.pdf](Fonti/GunnarsonEtAl-2021-learning-navigation-vortical-flows.pdf) | simulazione e benchmark di optimal control |
| 13 | Solnyshkov et al. 2019 | [SolnyshkovEtAl-2019-quantum-Kerr-Penrose-BEC.pdf](Fonti/SolnyshkovEtAl-2019-quantum-Kerr-Penrose-BEC.pdf) | proposta teorica Kerr-like con quasiparticelle massive |
| 14 | Westerberg et al. 2014 | [WesterbergEtAl-2014-time-dependent-optical-media.pdf](Fonti/WesterbergEtAl-2014-time-dependent-optical-media.pdf) | proposta ottica tempo-dipendente |
| 15 | Patrick et al. 2018 | [PatrickEtAl-2018-quasibound-states-draining-vortex.pdf](Fonti/PatrickEtAl-2018-quasibound-states-draining-vortex.pdf) | teoria di modi massivi/trappolati in vortice rotante |
| 16 | Bossard et al. 2023 | [BossardEtAl-2023-hydrodynamic-analogue-horizons-laser-cavities.pdf](Fonti/BossardEtAl-2023-hydrodynamic-analogue-horizons-laser-cavities.pdf) | validazione in canale d'acqua e progetto di cavità a due ostacoli |
| 17 | Jaskula et al. 2012 | [JaskulaEtAl-2012-dynamical-Casimir-BEC.pdf](Fonti/JaskulaEtAl-2012-dynamical-Casimir-BEC.pdf) | esperimento BEC con velocità del suono modulata |
| 18 | Steinhauer et al. 2022 | [SteinhauerEtAl-2022-cosmological-particle-creation-fluid-light.pdf](Fonti/SteinhauerEtAl-2022-cosmological-particle-creation-fluid-light.pdf) | esperimento con quench in un fluido quantistico di luce |
| 19 | Liebchen e Löwen 2019 | [LiebchenLoewen-2019-optimal-active-particle-navigation.pdf](Fonti/LiebchenLoewen-2019-optimal-active-particle-navigation.pdf) | teoria: Fermat generalizzato e soluzioni esatte in flussi vorticosi |
| 20 | Yang et al. 2025 | [YangEtAl-2025-dynamic-flow-control-active-matter.pdf](Fonti/YangEtAl-2025-dynamic-flow-control-active-matter.pdf) | esperimento: campi di flusso spaziotemporali programmati con luce |
| 21 | Haeufle et al. 2016 | [HaeufleEtAl-2016-external-control-self-propelled-particles.pdf](Fonti/HaeufleEtAl-2016-external-control-self-propelled-particles.pdf) | esperimento: colloidi Janus controllati verso un target |

La voce 20 è la nuova precedenza sperimentale più concreta per costruire un
background realmente programmabile nel tempo. Non realizza direttamente la
metrica del paper, ma dimostra il componente hardware che prima mancava:
generazione, composizione e riconfigurazione ottica di flussi vorticosi con
trasporto misurabile di particelle. La voce 21 è il precedente sperimentale
più diretto per il controllo feedback di un agente materiale verso un target.
La voce 19 è invece il confronto teorico più vicino alla brachistocrona: deve
essere usata per chiarire con precisione che cosa il formalismo attuale
aggiunga al principio di Fermat/Zermelo già noto per particelle attive.

### Fondazione delle metriche acustiche

1. W. G. Unruh, “Experimental Black-Hole Evaporation?”, *Phys. Rev. Lett.* **46**, 1351--1353 (1981), [DOI 10.1103/PhysRevLett.46.1351](https://doi.org/10.1103/PhysRevLett.46.1351). Fondazione dell'analogo acustico.

### Esperimenti con acqua e vortici rotanti

2. S. Weinfurtner et al., “Measurement of Stimulated Hawking Emission in an Analogue System”, *Phys. Rev. Lett.* **106**, 021302 (2011), [DOI 10.1103/PhysRevLett.106.021302](https://doi.org/10.1103/PhysRevLett.106.021302). Dimostra l'uso quantitativo di onde superficiali e orizzonti idrodinamici.

3. T. Torres et al., “Rotational superradiant scattering in a vortex flow”, *Nature Physics* **13**, 833--836 (2017), [DOI 10.1038/nphys4151](https://doi.org/10.1038/nphys4151). Evidenza sperimentale della componente rotante.

4. T. Torres et al., “Waves on a vortex: rays, rings and resonances”, *J. Fluid Mech.* **857**, 291--311 (2018), [DOI 10.1017/jfm.2018.752](https://doi.org/10.1017/jfm.2018.752). Deriva le due famiglie co- e counter-rotating di orbite instabili e confronta il ray tracing con dati di vasca.

5. T. Torres et al., “Quasinormal Mode Oscillations in an Analogue Black Hole Experiment”, *Phys. Rev. Lett.* **125**, 011301 (2020), [DOI 10.1103/PhysRevLett.125.011301](https://doi.org/10.1103/PhysRevLett.125.011301). Mostra che ringdown e barriere efficaci sono misurabili in un vortice.

6. S. Patrick et al., “Backreaction in an Analogue Black Hole Experiment”, *Phys. Rev. Lett.* **126**, 041105 (2021), [DOI 10.1103/PhysRevLett.126.041105](https://doi.org/10.1103/PhysRevLett.126.041105). Particolarmente importante per un background effettivamente dinamico.

### Vortice quantistico rotante

7. P. Švančara et al., “Rotating curved spacetime signatures from a giant quantum vortex”, *Nature* **628**, 66--70 (2024), [DOI 10.1038/s41586-024-07176-8](https://doi.org/10.1038/s41586-024-07176-8). Piattaforma più vicina al settore rotante a bassa viscosità.

### Background non stazionari in BEC

8. S. Eckel et al., “A Rapidly Expanding Bose--Einstein Condensate: An Expanding Universe in the Lab”, *Phys. Rev. X* **8**, 021021 (2018), [DOI 10.1103/PhysRevX.8.021021](https://doi.org/10.1103/PhysRevX.8.021021). Realizzazione sperimentale di espansione e redshift dei fononi.

9. S. Banik et al., “Accurate Determination of Hubble Attenuation and Amplification in Expanding and Contracting Cold-Atom Universes”, *Phys. Rev. Lett.* **128**, 090401 (2022), [DOI 10.1103/PhysRevLett.128.090401](https://doi.org/10.1103/PhysRevLett.128.090401). Misura quantitativa della risposta a espansione e contrazione.

10. C. Viermann et al., “Quantum field simulator for dynamics in curved spacetime”, *Nature* **611**, 260--264 (2022), [DOI 10.1038/s41586-022-05313-9](https://doi.org/10.1038/s41586-022-05313-9). Dimostra un simulatore 2D con curvatura e storia temporale configurabili.

### Agenti attivi e quasiparticelle massive

11. A. Daddi-Moussa-Ider, H. Löwen and B. Liebchen, “Hydrodynamics can determine the optimal route for microswimmer navigation”, *Communications Physics* **4**, 15 (2021), [DOI 10.1038/s42005-021-00522-6](https://doi.org/10.1038/s42005-021-00522-6). Formula PMP, endpoint fissati e tempo finale libero anche nel caso non autonomo; propone active colloids programmabili.

12. P. Gunnarson et al., “Learning efficient navigation in vortical flow fields”, *Nature Communications* **12**, 7143 (2021), [DOI 10.1038/s41467-021-27015-y](https://doi.org/10.1038/s41467-021-27015-y). Studio computazionale di agenti a velocità propria fissa in flussi non stazionari con confronto contro optimal control; supporta il design del protocollo, non costituisce una validazione sperimentale del rail.

13. D. D. Solnyshkov et al., “Quantum analogue of a Kerr black hole and the Penrose effect in a Bose--Einstein condensate”, *Phys. Rev. B* **99**, 214511 (2019), [DOI 10.1103/PhysRevB.99.214511](https://doi.org/10.1103/PhysRevB.99.214511). Proposta Kerr-like in cui difetti topologici agiscono come test particles massive; da citare come proposta, non come validazione già eseguita del rail.

### Ottica tempo-dipendente, solo come alternativa

14. N. Westerberg et al., “Experimental quantum cosmology in time-dependent optical media”, arXiv:1403.5910, [arXiv](https://arxiv.org/abs/1403.5910). Proposta di metriche con $n(t)$; utile per motivare la fattibilità della modulazione temporale, con le cautele sulla conformal invariance dei raggi nulli.

### Altri lavori mirati emersi dalla ricerca

15. S. Patrick et al., “Black Hole Quasibound States from a Draining Bathtub Vortex Flow”, *Phys. Rev. Lett.* **121**, 061101 (2018), [DOI 10.1103/PhysRevLett.121.061101](https://doi.org/10.1103/PhysRevLett.121.061101). Mostra che la vorticità del core produce un termine di massa locale e modi quasi-legati longevi. È un precedente teorico, non la misura di una separatrice massiva, ma rende più plausibile cercare osservabili massive/trappolate in una piattaforma a vortice.

16. A. Bossard et al., “How to create analogue black hole or white fountain horizons and LASER cavities in experimental free surface hydrodynamics?”, arXiv:2307.11022 (2023), [arXiv](https://arxiv.org/abs/2307.11022). È una fonte pratica per il dimensionamento di un canale d'acqua: valida sperimentalmente le leggi del flusso transcritico con un ostacolo e propone una cavità a due ostacoli. La parte a due ostacoli resta prospettica e non va descritta come risultato già osservato.

17. J.-C. Jaskula et al., “Acoustic Analog to the Dynamical Casimir Effect in a Bose--Einstein Condensate”, *Phys. Rev. Lett.* **109**, 220401 (2012), [DOI 10.1103/PhysRevLett.109.220401](https://doi.org/10.1103/PhysRevLett.109.220401). Modula sperimentalmente densità e velocità del suono cambiando la rigidezza della trappola e misura eccitazioni correlate. È una dimostrazione diretta che un parametro metrico efficace può essere variato nel tempo, ma misura fononi e non worldline massive.

18. J. Steinhauer et al., “Analogue cosmological particle creation in an ultracold quantum fluid of light”, *Nature Communications* **13**, 2890 (2022), [DOI 10.1038/s41467-022-30603-1](https://doi.org/10.1038/s41467-022-30603-1). Realizza un quench delle interazioni in un fluido quantistico di luce e misura la produzione di particelle analoghe. È una seconda piattaforma non stazionaria, più semplice di una BEC atomica, con la cautela che una coordinata di propagazione svolge il ruolo del tempo efficace.

19. B. Liebchen and H. Löwen, “Optimal Control Strategies for Active Particle Navigation”, arXiv:1901.08382 (2019), pubblicato come “Optimal navigation strategies for active particles”, *EPL* **127**, 34003 (2019), [DOI 10.1209/0295-5075/127/34003](https://doi.org/10.1209/0295-5075/127/34003). Formula un principio di Fermat variazionale per particelle attive e ricava traiettorie esatte anche in flussi vorticosi. È il confronto concettuale obbligatorio per delimitare l'originalità del rail relativistico, del vincolo di energia e della correzione adiabatica off-shell.

20. F. Yang et al., “Dynamic flow control through active matter programming language”, *Nature Materials* **24**, 615--625 (2025), [DOI 10.1038/s41563-024-02090-w](https://doi.org/10.1038/s41563-024-02090-w). Dimostra sperimentalmente campi di flusso micrometrici, vorticosi e riconfigurabili nello spazio e nel tempo mediante reti motore--microtubulo controllate con luce. È la migliore fonte trovata per giustificare un banco attivo con background programmabile, pur restando necessaria una specifica procedura di calibrazione per realizzare l'ovale affine del rail.

21. D. F. B. Haeufle et al., “External control strategies for self-propelled particles: Optimizing navigational efficiency in the presence of limited resources”, *Phys. Rev. E* **94**, 012617 (2016), [DOI 10.1103/PhysRevE.94.012617](https://doi.org/10.1103/PhysRevE.94.012617). È il precedente sperimentale più vicino al demonstrator materiale: colloidi Janus di silice/carbonio in una miscela acqua--2,6-lutidina, imaging CCD e velocità di propulsione regolata con l'intensità luminosa in funzione di posizione e orientazione rispetto al target. Non realizza ancora un flusso Kerr-like né l'indicatrice del rail, ma dimostra tracking, feedback e controllo della motilità necessari al protocollo.

## 20. Sottoinsieme minimo da citare nel paper

Non è necessario aggiungere tutti i ventuno riferimenti. Per un solo paragrafo di outlook sono sufficienti:

- Unruh 1981 — fondazione dell'analogo acustico;
- Torres et al. 2017 — settore rotante sperimentale;
- Patrick et al. 2021 — metrica analogica dinamica;
- Švančara et al. 2024 — vortice quantistico rotante;
- Viermann et al. 2022 — background tempo-dipendente controllato;
- Daddi-Moussa-Ider et al. 2021 — navigazione attiva non autonoma/PMP;
- Haeufle et al. 2016 — controllo feedback sperimentale di un colloide attivo.

Se si vuole proporre un demonstrator davvero costruibile, aggiungerei anche Yang et al. 2025 e Bossard et al. 2023. Se invece il punto è delimitare l'originalità matematica, aggiungerei Liebchen--Löwen 2019 e Patrick et al. 2018. Eckel, Banik, Weinfurtner, Jaskula, Steinhauer e Solnyshkov possono entrare se il testo dedica più di un paragrafo all'argomento.

## 21. Voci BibTeX pronte per Claude

Queste voci sono proposte per `paper/refs.bib` soltanto se le rispettive fonti vengono effettivamente citate nel manoscritto. Non sono state aggiunte automaticamente.

```bibtex
@article{Unruh1981Analogue,
  author = {Unruh, W. G.},
  title = {Experimental Black-Hole Evaporation?},
  journal = {Phys. Rev. Lett.},
  volume = {46},
  pages = {1351--1353},
  year = {1981},
  doi = {10.1103/PhysRevLett.46.1351}
}

@article{TorresEtAl2017Superradiance,
  author = {Torres, Theo and Patrick, Sam and Coutant, Antonin and Richartz, Mauricio and Tedford, Edmund W. and Weinfurtner, Silke},
  title = {Rotational superradiant scattering in a vortex flow},
  journal = {Nature Phys.},
  volume = {13},
  pages = {833--836},
  year = {2017},
  doi = {10.1038/nphys4151}
}

@article{PatrickEtAl2021Backreaction,
  author = {Patrick, Sam and Goodhew, Harry and Gooding, Cisco and Weinfurtner, Silke},
  title = {Backreaction in an Analogue Black Hole Experiment},
  journal = {Phys. Rev. Lett.},
  volume = {126},
  pages = {041105},
  year = {2021},
  doi = {10.1103/PhysRevLett.126.041105}
}

@article{SvancaraEtAl2024Vortex,
  author = {\v{S}van\v{c}ara, Patrik and Smaniotto, Pietro and Solidoro, Leonardo and MacDonald, James F. and Patrick, Sam and Gregory, Ruth and Barenghi, Carlo F. and Weinfurtner, Silke},
  title = {Rotating curved spacetime signatures from a giant quantum vortex},
  journal = {Nature},
  volume = {628},
  pages = {66--70},
  year = {2024},
  doi = {10.1038/s41586-024-07176-8}
}

@article{ViermannEtAl2022CurvedSpacetime,
  author = {Viermann, Celia and Sparn, Marius and Liebster, Nikolas and Hans, Maurus and Kath, Elinor and Parra-L\'opez, \'Alvaro and Tolosa-Sime\'on, Mireia and S\'anchez-Kuntz, Natalia and Haas, Tobias and Strobel, Helmut and Floerchinger, Stefan and Oberthaler, Markus K.},
  title = {Quantum field simulator for dynamics in curved spacetime},
  journal = {Nature},
  volume = {611},
  pages = {260--264},
  year = {2022},
  doi = {10.1038/s41586-022-05313-9}
}

@article{DaddiMoussaIderEtAl2021Navigation,
  author = {Daddi-Moussa-Ider, Abdallah and L\"owen, Hartmut and Liebchen, Benno},
  title = {Hydrodynamics can determine the optimal route for microswimmer navigation},
  journal = {Commun. Phys.},
  volume = {4},
  pages = {15},
  year = {2021},
  doi = {10.1038/s42005-021-00522-6}
}

@article{PatrickEtAl2018Quasibound,
  author = {Patrick, Sam and Coutant, Antonin and Richartz, Mauricio and Weinfurtner, Silke},
  title = {Black Hole Quasibound States from a Draining Bathtub Vortex Flow},
  journal = {Phys. Rev. Lett.},
  volume = {121},
  pages = {061101},
  year = {2018},
  doi = {10.1103/PhysRevLett.121.061101}
}

@misc{BossardEtAl2023HydrodynamicHorizons,
  author = {Bossard, Alexis and James, Nicolas and Aucouturier, Camille and Fourdrinoy, Johan and Robertson, Scott and Rousseaux, Germain},
  title = {How to create analogue black hole or white fountain horizons and {LASER} cavities in experimental free surface hydrodynamics?},
  year = {2023},
  eprint = {2307.11022},
  archivePrefix = {arXiv},
  primaryClass = {physics.flu-dyn}
}

@article{JaskulaEtAl2012DynamicalCasimir,
  author = {Jaskula, Jean-Christophe and Partridge, Guthrie B. and Bonneau, Marie and Lopes, Raphael and Ruaudel, Josselin and Boiron, Denis and Westbrook, Christoph I.},
  title = {Acoustic Analog to the Dynamical Casimir Effect in a Bose--Einstein Condensate},
  journal = {Phys. Rev. Lett.},
  volume = {109},
  pages = {220401},
  year = {2012},
  doi = {10.1103/PhysRevLett.109.220401}
}

@article{SteinhauerEtAl2022ParticleCreation,
  author = {Steinhauer, Jeff and Abuzarli, Murad and Aladjidi, Tangui and Bienaime, Tom and Piekarski, Clara and Liu, Wei and Giacobino, Elisabeth and Bramati, Alberto and Glorieux, Quentin},
  title = {Analogue cosmological particle creation in an ultracold quantum fluid of light},
  journal = {Nature Commun.},
  volume = {13},
  pages = {2890},
  year = {2022},
  doi = {10.1038/s41467-022-30603-1}
}

@article{LiebchenLoewen2019Navigation,
  author = {Liebchen, Benno and L\"owen, Hartmut},
  title = {Optimal navigation strategies for active particles},
  journal = {EPL},
  volume = {127},
  pages = {34003},
  year = {2019},
  doi = {10.1209/0295-5075/127/34003}
}

@article{YangEtAl2025DynamicFlow,
  author = {Yang, Fan and Liu, Shichen and Lee, Heun Jin and Phillips, Rob and Thomson, Matt},
  title = {Dynamic flow control through active matter programming language},
  journal = {Nature Mater.},
  volume = {24},
  pages = {615--625},
  year = {2025},
  doi = {10.1038/s41563-024-02090-w}
}

@article{HaeufleEtAl2016ExternalControl,
  author = {Haeufle, Daniel F. B. and B\"auerle, Tobias and Steiner, Jakob and Bremicker, Lena and Schmitt, Syn and Bechinger, Clemens},
  title = {External control strategies for self-propelled particles: Optimizing navigational efficiency in the presence of limited resources},
  journal = {Phys. Rev. E},
  volume = {94},
  pages = {012617},
  year = {2016},
  doi = {10.1103/PhysRevE.94.012617}
}
```

## 22. Roadmap concreta

### Fase 0 — una settimana: definizione del demonstrator

- scegliere una sola metrica/famiglia di parametri;
- scegliere una sola osservabile primaria, preferibilmente $\Delta\varphi$;
- scegliere il dominio regolare lontano dalla separatrice;
- derivare $c$ e $R$ in coordinate di laboratorio;
- produrre un digital twin con rumore e latenza realistici.

### Fase 1 — due/quattro settimane: hardware-in-the-loop

- motion tracking;
- calibrazione dell'ellisse;
- traiettorie stazionarie;
- prima verifica PMP;
- acquisizione del pilot per scegliere il range di $\varepsilon$.

### Fase 2 — quattro/otto settimane: flagship off-shell

- dataset preregistrato;
- on-shell vs full;
- fit delle slope;
- ablation, ramp reversal e condizioni blind;
- archivio di dati grezzi, firmware e analisi.

### Fase 3 — collaborazione esterna

- inviare un proposal di due pagine a un gruppo di analogue gravity;
- offrire le Hamiltoniane, le separatrici e le predizioni dimensionali;
- chiedere se un parametro del loro vortice può essere modulato lentamente mantenendo il profilo calibrabile;
- definire prima la mappa $H_{\rm lab}\leftrightarrow H_{\rm paper}$ e solo dopo l'esperimento.

## 23. Verdetto finale

Il miglior esperimento immediatamente producibile non è una vasca che “crea Thakurta--Kerr”, ma un **demonstrator planare del controlled rail**. Può verificare proprio ciò che è originale nel lavoro: indicatrice, controllo, endpoint, termine off-shell e slope due.

Il miglior esperimento esterno per Paper II è invece un **vortice rotante misurato**, preferibilmente superfluido o, più accessibilmente, una vasca con draining vortex. Può verificare la fenomenologia Kerr-like e la separatrix crossing, ma non il carattere massive/forzato senza un agente attivo o una quasiparticella massive.

La combinazione dei due protocolli sarebbe scientificamente forte: il primo testa esattamente la matematica del rail; il secondo mostra che la struttura rotante e non stazionaria ha un analogo fisico indipendente. Presentata con questi limiti, la proposta risponde al referee senza sovravendere il risultato.
