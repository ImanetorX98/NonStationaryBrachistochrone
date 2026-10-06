# Blocco 2b — inversa temporale costruita e propagazione della nullità

6 ottobre 2026. Questo checkpoint aggiunge 13 teoremi ai 117 del blocco 2a.
Il progetto raggiunge 130 teoremi pubblici. Il risultato chiude due premesse
che il checkpoint precedente aveva lasciato esplicite: l'esistenza di una
vera inversa temporale locale e la propagazione della nullità per il flusso
Hamiltoniano completo non stazionario. Non chiude esistenza e continuazione
delle soluzioni, né il collegamento con le geodetiche di Levi-Civita e Jacobi.

Il documento [BLOCK2_COORDINATE_HAMILTONIAN.md](BLOCK2_COORDINATE_HAMILTONIAN.md)
conserva lo stato del blocco 2a. Le sue indicazioni di lavoro ancora aperto
sull'inversa e sulla nullità sono superate dai risultati qui descritti.

## 1. Costruzione dell'inversa: nessuna premessa di invertibilità

Il lemma `local_time_inverse_exists` assume

\[
\forall\lambda\text{ vicino ad }a:\quad
v_\lambda(\lambda)=V(\lambda),\qquad
V\text{ continua in }a,\qquad V(a)>0.
\]

La continuità della derivata rende \(v\) strettamente differenziabile in
\(a\). Il teorema della funzione inversa di mathlib **costruisce** una
mappa \(\ell\) e dimostra

\[
\ell(v(a))=a,\qquad
\ell(v(\lambda))=\lambda\text{ localmente vicino ad }a,
\]

\[
v(\ell(t))=t\text{ localmente vicino a }v(a),\qquad
\ell_v(v(a))=\frac1{V(a)}.
\]

La sola derivata positiva in un punto non è la premessa del lemma: sono
esplicite l'equazione temporale in un intorno e la continuità di \(V\).
Non si presume una mappa inversa né la sua derivata.

Per il sistema Hamiltoniano,

\[
V(\lambda)=-\frac{P_v(\lambda)}{w(\lambda)}+P_r(\lambda),\qquad
w(\lambda)=q-1+\frac{2M(\lambda)}{r(\lambda)}.
\]

`affineTimeSpeed_continuousAt` deriva la continuità di \(V\) dalla
continuità di \(M,r,P_v,P_r\), con \(r(a)\ne0\), \(w(a)\ne0\).
`hamiltonian_time_inverse_exists` applica quindi la costruzione a una
soluzione futura regolare. \(M\) è la massa lungo il parametro affine;
nell'applicazione \(M(\lambda)=m(v(\lambda))\).

## 2. Collegamento con le ODE già derivate

`hamiltonian_local_direction_exists` compone l'inversa appena costruita con
il blocco 2a. Assume \(q>0,r(a)>0,w(a)>0,L>0\), nullità all'evento,
direzione futura, continuità degli stati e le equazioni Hamiltoniane affini.
Conclude l'esistenza dell'inversa locale e le tre ODE in \(v(a)\):

\[
r_v=-w+\sqrt{qw}\cos\alpha,
\]

\[
\alpha_v=\frac{\sin\alpha}{r}\left[
\cos\alpha(2w-q+1)-\frac12\sqrt{q/w}(3w-q+1)\right],
\qquad \varphi_v=\frac{\sqrt w}{r}\sin\alpha.
\]

L'angolo è \(\alpha=\arccos(\sqrt q P_r/S)\),
\(S=\sqrt{qP_r^2+L^2/r^2}\). La conclusione differenziale del teorema
composto è puntuale, mentre le identità di inversa valgono in intorni.
Applicare il risultato a ogni evento non equivale ancora a costruire una
soluzione futura globale o a incollare formalmente tutte le carte temporali.

## 3. Non stazionarietà: il momento temporale deve evolvere

Per

\[
H=\frac12\left[-\frac{P_v^2}{w}+2P_vP_r+(q-w)P_r^2+\frac{L^2}{r^2}\right],
\qquad D_mH=\frac1r\left(\frac{P_v^2}{w^2}-P_r^2\right),
\]

Lean certifica la derivata rispetto alla massa e, mediante composizione,

\[
H_v=m_vD_mH,\qquad P_{v,\lambda}=-m_vD_mH.
\]

La seconda equazione è la corrispondente equazione Hamiltoniana, usata come
premessa del teorema sul flusso. **Non viene imposta la conservazione di
\(P_v\)**. Il parametro \(q=\widehat E^2\) resta quello della metrica;
il valore conservato \(H\) è il vincolo Hamiltoniano, non un'asserzione
di conservazione dell'energia fisica temporale nel caso non stazionario.

`mass_along_affine_derivative` dimostra, per \(M=m\circ v\),

\[
M_\lambda=m_vv_\lambda.
\]

`hamiltonian_along_path_derivative` verifica la regola della catena completa
lungo un cammino arbitrario differenziabile. Inserendo le equazioni
Hamiltoniane e \(L_\lambda=0\), si ottiene

\[
\frac{dH}{d\lambda}
=D_mH\,m_vv_\lambda+H_r r_\lambda
+H_{P_v}P_{v,\lambda}+H_{P_r}P_{r,\lambda}
+H_L L_\lambda=0.
\]

La cancellazione include il termine di massa; il risultato non assume
nullità. Il lemma `hamiltonian_flow_null_propagates` usa poi il teorema
del valor medio su un dominio convesso per concludere

\[
G^{-1}(p(a),p(a))=0\quad\Longrightarrow\quad
G^{-1}(p(b),p(b))=0.
\]

La nullità iniziale è l'unica premessa di nullità. La positività di
\(v_\lambda\) all'evento resta una distinta ipotesi di orientazione futura.

## 4. Momento angolare e limiti dell'applicazione

Il modulo verifica che \(H\) non dipende da \(\varphi\) e che la sua
derivata parziale angolare è zero. Dall'equazione Hamiltoniana
\(L_\lambda=0\), `angular_momentum_constant` conclude
\(L(b)=L(a)\) su un dominio convesso.

Le soluzioni Hamiltoniane sono ancora fornite nelle premesse, non costruite.
Le ipotesi valgono sul dominio indicato; per propagare un livello bastano
\(r\ne0,w\ne0\). Il ramo di direzione richiede inoltre
\(q>0,r>0,w>0,L>0\). I nuovi lemmi non cambiano la soglia
\(q\geq3/2\) dei lemmi sulle curvature.

La propagazione della nullità e la costanza di \(L\) sono composte con le
equazioni Hamiltoniane nei loro teoremi. Rimane da confezionare il flusso
geometrico complessivo e usarlo per costruire le famiglie `DirectionFlow`
su intervalli massimali, senza assumere a priori esistenza futura globale.

## 5. Prossimo passaggio

Dimostrare esistenza locale e unicità del sistema di direzione per una storia
di massa regolare, quindi un teorema di continuazione sui compatti del dominio.
Separare l'uscita dal dominio esterno \(r>2m(v)\) dalla perdita di esistenza
dell'ODE. Soltanto dopo questo passaggio si potranno applicare i lemmi di fuga
senza ricevere una soluzione futura globale come premessa.

Rimangono separati la connessione di Levi-Civita, lo schermo/Jacobi, properness,
rivestimento e il collegamento finale ai Maxwell. Non sono conclusioni di
questo checkpoint.

## 6. Sorgenti e verifica

- `ViaB/LocalTimeInverse.lean`: 4 nuovi teoremi.
- `ViaB/TimeDependentHamiltonian.lean`: 9 nuovi teoremi.
- `AxiomAudit.lean`: copertura di tutti i teoremi pubblici.
- `verification/`: compilazione, dipendenze assiomatiche e hash dei sorgenti.

La verifica usa Lean 4.24.0 e mathlib v4.24.0. I soli assiomi ammessi sono
`propext`, `Classical.choice`, `Quot.sound`. Le ipotesi esplicite dei teoremi
restano parte dell'enunciato anche dopo il controllo del kernel.
