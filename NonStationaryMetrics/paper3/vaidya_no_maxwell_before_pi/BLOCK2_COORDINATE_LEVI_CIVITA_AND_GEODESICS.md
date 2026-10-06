# Blocco 2g — connessione coordinata e geodetiche affini

6 ottobre 2026. Checkpoint compilato: **222 teoremi pubblici**, di cui
23 nuovi in tre moduli. Lean 4.24.0, mathlib v4.24.0.

## Risultato e ambito

Il passaggio dalle ODE di direzione all'equazione geodetica della metrica
efficace di Fermat è ora verificato nel chart (v,r,φ). La connessione
è ricavata dalle derivate della metrica mediante Koszul e caratterizzata
come unica connessione coordinata senza torsione e compatibile con la metrica.
Una derivazione separata dalle equazioni Hamiltoniane dà la stessa accelerazione.

Questo risultato non conclude ancora l'esclusione geometrica dei Maxwell:
restano il ponte di Jacobi e curvatura, la mappa di arrivo e la properness.
La metrica efficace di Fermat non va confusa con la metrica fisica di Vaidya.

## 1. Metrica, derivate e unicità della connessione

Poniamo q costante, w=q−1+2m(v)/r e μ=m_v. La matrice è

\[
g=\begin{pmatrix}
-w+w^2/q&w/q&0\\
w/q&1/q&0\\
0&0&r^2
\end{pmatrix},\qquad
g^{-1}=\begin{pmatrix}
-1/w&1&0\\
1&q-w&0\\
0&0&1/r^2
\end{pmatrix}.
\]

Lean verifica entrambi i prodotti inversi per q,r,w non nulli. Le derivate
coordinate non nulle sono

\[
\begin{aligned}
\partial_v g_{vv}&=(2w/q-1)2\mu/r,&
\partial_v g_{vr}&=2\mu/(qr),\\
\partial_r g_{vv}&=(2w/q-1)(-2m/r^2),&
\partial_r g_{vr}&=-2m/(qr^2),\\
\partial_r g_{\varphi\varphi}&=2r.
\end{aligned}
\]

Si includono le componenti simmetriche; le derivate angolari sono zero.
Queste matrici non sono soltanto definizioni algebriche: i teoremi
`fermatMetric_time_derivative`, `fermatMetric_radius_derivative` e
`fermatMetric_angle_derivative` verificano le rispettive `HasDerivAt`.

Con D_k g_ij così verificata,

\[
C_{ijk}=\tfrac12(D_jg_{ik}+D_kg_{ij}-D_ig_{jk}),\qquad
\Gamma^i{}_{jk}=\sum_h g^{ih}C_{hjk}.
\]

`fermatChristoffel_torsion_free` e
`fermatChristoffel_metric_compatible` dimostrano

\[
\Gamma^i{}_{jk}=\Gamma^i{}_{kj},\qquad
D_kg_{ij}=\sum_h(g_{hj}\Gamma^h{}_{ik}+g_{ih}\Gamma^h{}_{jk}).
\]

`coordinateKoszul_unique` prova l'unicità dei coefficienti abbassati;
`fermatChristoffel_unique` la trasferisce alla connessione con indice alto
usando la metrica inversa. Non si assume l'unicità come assioma geometrico.

## 2. Contrazione e confronto indipendente con Hamilton

Per velocità u=(V,R,P), la contrazione abbassata è

\[
\begin{aligned}
K_v&=(2w/q-1)\frac{\mu}{r}V^2
 -\frac{2m(2w/q-1)}{r^2}VR-\frac{2m}{qr^2}R^2,\\
K_r&=\left(\frac{2\mu}{qr}+\frac{m(2w/q-1)}{r^2}\right)V^2-rP^2,\\
K_\varphi&=2rRP.
\end{aligned}
\]

Qui K_i sono coefficienti della contrazione, **non curvature dello schermo**.
Il risultato della formula di Koszul è

\[
\Gamma(u,u)=(-K_v/w+K_r,\ K_v+(q-w)K_r,\ K_\varphi/r^2).
\]

Separatamente, dalle velocità Hamiltoniane
V=−P_v/w+P_r, R=P_v+(q−w)P_r, P=L/r², poniamo
F=−H_r e G=H_m. La catena temporale dà

\[
\dot w=\frac{2\mu}{r}V-\frac{2m}{r^2}R,
\]

quindi

\[
\begin{aligned}
\dot V&=F+\mu G/w+P_v\dot w/w^2,\\
\dot R&=-\mu G+(q-w)F-\dot w P_r,\\
\dot P&=-2LR/r^3.
\end{aligned}
\]

`hamiltonian_raised_velocity_derivative` verifica queste derivate lungo
il flusso. `hamiltonian_acceleration_equals_connection` verifica
algebricamente che tale accelerazione è −Γ(u,u).
Quest'ultima identità non richiede H=0: la nullità è un vincolo separato.
La connessione non è definita a partire dall'accelerazione Hamiltoniana;
le due costruzioni vengono confrontate solo dopo averle derivate.

## 3. Seconde derivate effettive e applicazione alle curve costruite

Da un intorno in cui valgono tutte le equazioni Hamiltoniane e la
regolarità della massa, `hamiltonian_germ_coordinate_geodesic` deduce

\[
\ddot x^i=-\Gamma^i{}_{jk}(x)\dot x^j\dot x^k.
\]

Le seconde derivate sono espresse tramite `HasDerivAt` delle derivate
coordinate effettive. Non si presuppone una curva C² che già soddisfi
l'equazione cercata. La struttura `CoordinateNullGeodesicAt` registra
questa equazione, g(ẋ,ẋ)=0 e v̇>0.

`direction_coordinate_null_geodesic_germ_exists` parte dalle ODE
coordinate, q>0, L>0, r>0, w>0, 0<α<π e una derivata della massa
su un intorno. Costruisce il parametro affine λ, l'inversa v e la
primitiva angolare, e conclude la struttura appena descritta su uno
stesso intorno affine. Non assume i momenti né l'orologio affine.

`maximal_exterior_coordinate_null_geodesic_germ_exists` applica il
risultato a ogni evento della curva massimale esterna già costruita,
con q>1, L>0 e massa C¹ sul dominio. Le ODE e la regolarità dello
stato sono dedotte dalla costruzione massimale.

## 4. Verifica e limiti

`verify.py` ha compilato tutti i 222 teoremi e controllato la copertura
dell'audit. Nessun `sorry`, `admit` o assioma geometrico aggiunto;
solo `propext`, `Classical.choice`, `Quot.sound` tra gli assiomi ammessi.
Log e hash dei sorgenti sono conservati in `verification/`.

È una certificazione nel chart, senza istanziare ancora una connessione
su una varietà astratta. Gli intorni affini sono locali: non si deduce
λ→∞ da v→∞. L>0 seleziona il ramo non radiale e una normalizzazione
dell'orientazione; i radiali L=0 non sono inclusi nel teorema composto.

## 5. Prossimi passaggi, senza circolarità

1. Provare la dipendenza differenziabile dal lancio, controllando la
   regolarità della massa necessaria. La dipendenza continua da Grönwall
   non basta a produrre una variazione di Jacobi.
2. Derivare la variazione dell'equazione geodetica, proiettarla sullo
   schermo nullo e identificare i coefficienti geometrici. In particolare,
   la positività della funzione `screenGap` già certificata non prova
   ancora che essa sia la differenza delle curvature effettive.
3. Collegare Jacobi, dati iniziali e determinante della mappa di arrivo;
   verificare la trasversalità del bersaglio senza dividere automaticamente
   per una velocità radiale che potrebbe annullarsi.
4. Provare la properness controllando fuga, cattura e limiti radiali dei
   lanci, poi ottenere rivestimento e conteggio delle fibre.

La soglia q≥3/2 appartiene alla positività scalare della curvatura;
q>1 appartiene ai risultati dinamici esterni. Non sono intercambiabili.
Papers I/II e la lettera CQG non sono stati modificati da questo blocco.
