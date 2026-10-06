# Blocco 2f — dalla direzione al flusso Hamiltoniano nullo affine

6 ottobre 2026. **20 nuovi teoremi**, **199 complessivi**, **44 sorgenti
Lean** inclusi radice e audit. Questo passaggio costruisce il collegamento
inverso che mancava: dalle soluzioni del sistema di direzione alle equazioni
Hamiltoniane affini della metrica di Fermat coordinata di Vaidya.

Il risultato è locale, su un intero intorno dello stesso parametro affine.
Non è ancora l’identificazione formale tramite Levi-Civita/Jacobi né il
teorema completo di esclusione dei Maxwell prima di π.

## 1. Dati e metrica effettivamente trattati

Si pone q = Ê² > 0, w = q − 1 + 2m(v)/r. La ricostruzione richiede
r > 0, w > 0 e 0 < α < π. Si sceglie una normalizzazione L > 0 costante.
Per la curva massimale esterna del blocco 2e, q > 1 e m > 0 garantiscono
w > 0 e tutti i requisiti sullo stato sono già dedotti dalla costruzione.

La forma quadratica coordinata e il suo Hamiltoniano sono quelli già
formalizzati nel progetto:

\[
g(\dot\gamma,\dot\gamma)
=-w\dot v^2+\frac{(\dot r+w\dot v)^2}{q}+r^2\dot\phi^2,
\]

\[
H=\frac12\left(-\frac{P_v^2}{w}+2P_vP_r+(q-w)P_r^2+
\frac{L^2}{r^2}\right).
\]

Si tratta della metrica di Fermat del problema, non di una dichiarazione
che queste siano le geodetiche nulle della metrica fisica originale di
Vaidya. L’identità fra questa metrica e la sua inversa è già certificata
in `MetricHamiltonian.lean`.

## 2. Ricostruzione algebrica, regolare alle inversioni radiali

Dal sistema di direzione

\[
r_v=-w+\sqrt{qw}\cos\alpha,\qquad
\alpha_v=\frac{\sin\alpha}{r}B(q,w,\cos\alpha)
\]

si definiscono

\[
P_r=\frac{L\cos\alpha}{\sqrt q\,r\sin\alpha},\qquad
P_v=wP_r-\frac{L\sqrt w}{r\sin\alpha},\qquad
V=\frac{L}{r\sqrt w\sin\alpha}>0.
\]

Lean verifica S = √(qP_r² + L²/r²) = L/(r sin α), la ricostruzione
di seno e coseno, H = 0 e −P_v/w + P_r = V > 0.
Nessun passaggio divide per r_v, P_r o cos α. Sono quindi ammessi
sia i punti di inversione radiale sia α = π/2; restano esclusi i limiti
radiali sin α = 0.

L > 0 è una normalizzazione del momento angolare. Il teorema non tratta
direttamente L = 0 o L < 0. La soglia q ≥ 3/2 del confronto di curvatura
resta distinta dalle ipotesi di questa ricostruzione.

## 3. Equazione radiale e forza temporale

La derivata della formula di P_r, usando le ODE di r e α, dà

\[
(P_r)_v=\frac{-H_r}{V}.
\]

L’uguaglianza viene dimostrata tramite la combinazione radiale già
verificata; non si usa una massa congelata. La forza radiale contiene
il valore istantaneo di m, mentre il tasso m_v entra nella forza temporale.

Per quest’ultima, prima si verifica la differenziabilità della formula
esplicita di P_v. Si applica poi la regola della catena ad H lungo la
curva ricostruita, usando H = 0 **algebricamente** in un intorno e le
equazioni radiali già dimostrate. Con μ = m_v si ottiene

\[
0=H_m\mu-H_r^{\rm force}r_v+V(P_v)_v+
H_{P_r}(P_r)_v,
\qquad
(P_v)_v=-\frac{\mu H_m}{V},
\]

dove H_r^{force} = −H_r e

\[
H_m=\frac1r\left(\frac{P_v^2}{w^2}-P_r^2\right).
\]

I due termini radiali si cancellano perché r_v = H_{P_r}/V e
(P_r)_v = H_r^{force}/V. Questo elimina l’incognita della derivata di
P_v senza presupporre la sua equazione.

**Controllo di non circolarità:** qui la nullità non viene ottenuta dal
lemma di propagazione Hamiltoniana, che richiederebbe già la forza
temporale. Viene dalla formula esplicita del covettore in ogni evento.
Non si assume che P_v sia conservato.

## 4. Parametro affine e angolo costruiti

La velocità del parametro affine è positiva:

\[
\frac{d\lambda}{dv}=F(v)=\frac{r\sqrt w\sin\alpha}{L}=\frac1V.
\]

`continuous_interval_primitive_exists` costruisce una primitiva usando
l’integrale di una estensione continua ottenuta con un clamp ai bordi
di un intervallo compatto locale. Nel tratto interno il clamp coincide
con l’identità; non modifica il sistema fisico.

Si costruisce quindi λ con λ(a) = 0 e una vera inversa locale v, con
v(0) = a ed entrambe le identità di inversa in un intorno. La derivata
dell’inversa viene verificata **nello stesso intorno**, non soltanto
all’evento iniziale:

\[
\frac{dv}{d\lambda}=V(v(\lambda)).
\]

Separatamente viene costruita una primitiva angolare con φ(a) = 0 e

\[
\phi_v=\frac{\sqrt w}{r}\sin\alpha.
\]

Non sono premesse né l’esistenza di λ e della sua inversa, né quella
di φ. La scelta φ(a) = 0 fissa l’origine angolare; un’aggiunta costante
non altera le equazioni.

## 5. Teorema composto e applicazione alla curva massimale

`direction_affine_hamiltonian_germ_exists` assume le ODE di direzione
e una massa differenziabile in un intorno dell’evento regolare. Con
le funzioni appena costruite conclude, nello stesso intorno affine,

\[
v_\lambda=H_{P_v},\quad r_\lambda=H_{P_r},\quad
\phi_\lambda=H_L,\quad
(P_v)_\lambda=-m_v H_m,\quad (P_r)_\lambda=-H_r,\quad L_\lambda=0,
\]

oltre a H = 0 e v_λ > 0. La struttura `HamiltonianLiftAt` registra
tutte queste equazioni e i due vincoli, senza omettere il momento temporale.

`maximal_exterior_affine_hamiltonian_germ_exists` applica il risultato
alla `maximalExteriorCurve`. Per m C¹ sul suo dominio aperto e q > 1,
deduce internamente l’ODE, r > 0, w > 0 e la striscia angolare.
La traiettoria Hamiltoniana affine locale è quindi una conclusione
della costruzione massimale coordinata.

## 6. Verifica e stato dei prossimi passi

La verifica completa compila 199 teoremi, controlla la copertura
dell’audit assiomatico e conserva log e hash SHA-256. Nessuna prova
viene lasciata in sospeso; si ammettono solo gli assiomi standard di
Lean/mathlib. Sono stati ricontrollati il segno della forza temporale,
il reciproco delle due velocità temporali, la positività dei denominatori
e la validità delle equazioni in un intorno della stessa inversa.

I prossimi passaggi sono:

1. Identificare formalmente il flusso Hamiltoniano della metrica di
   Fermat con la sua equazione geodetica e con i campi di Jacobi;
   i risultati scalari di Sturm non sostituiscono questa identificazione.
2. Costruire la mappa di arrivo e trasferirvi regolarità e assenza di
   degenerazione locale sotto le ipotesi del teorema via B.
3. Provare la properness controllando separatamente fuga, cattura e
   limiti dei dati di lancio; la sola fuga nel cono non basta.
4. Applicare il rivestimento e il conteggio delle fibre già disponibili
   in forma astratta, quindi chiudere il passaggio ai Maxwell.

Non è affermata una completezza affine futura: un dominio illimitato
in v non implica automaticamente λ → +∞. Non sono ancora formalizzati
una connessione di Levi-Civita, il ponte completo di Jacobi, la properness
o la conclusione geometrica sui Maxwell. Papers I/II e la lettera CQG
restano invariati; questi sviluppi appartengono a Paper III.
