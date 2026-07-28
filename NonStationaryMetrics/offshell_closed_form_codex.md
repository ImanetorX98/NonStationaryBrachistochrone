# Chiusura del termine off-shell di dilatazione

## Esito

Il contributo off-shell generato dalla dilatazione del momento radiale nella branca
\(t/\eta\) di Thakurta--Kerr ammette una forma chiusa generale in
\((M,a,E,J)\), purché per «forma chiusa» si intenda la forma canonica nelle
funzioni abeliane iterate della curva iperellittica di genere due. In generale
non esiste una riduzione a funzioni elementari, né alle sole funzioni
Kleiniane di peso uno.

La formula derivata qui chiude il contributo che proviene dal termine
\(P_r H_{P_r}\) dell'operatore completo

\[
 \mathcal D H=\Theta H+P_rH_{P_r}.
\]

Più precisamente, chiude

\[
 \Phi_{\mathrm{off,dil}}(r)
 =-\int_{r_0}^{r}
   \left.\frac{\partial_{P_r}G_0}{H_{0,P_r}}\right|_{\rm shell}(x)
   \left(\int_{r_0}^{x}p_{r0}(z)\,dz\right)dx .
\]

Il contributo costruito con \(\Theta H\), già trattato nella riduzione
precedente del progetto, deve essere sommato a questa espressione per ottenere
l'intero termine \(S_{\mathcal D}\). Quindi «FULL» nel seguito significa:
tutti i blocchi della nuova parte di dilatazione, non l'intera correzione
adiabatica del paper presa isolatamente.

L'implementazione verificabile è:

[`KerrSessionScripts/offshell_tbranch_closed_form_codex.py`](KerrSessionScripts/offshell_tbranch_closed_form_codex.py)

## 1. Curva e dati algebrici

Poniamo

\[
 D(r)=(E^2-1)r+2M,\qquad
 \Delta(r)=r^2-2Mr+a^2
\]

e

\[
\begin{aligned}
Q_2(r)={}&E^2r^4
 +(J^2-E^2J^2+E^2a^2)r^2\\
&+\bigl(2E^2J^2M-4E^2JMa+2E^2Ma^2
        -4J^2M+4JMa\bigr)r\\
&+4M^2(J-a)^2 .
\end{aligned}
\]

La curva spettrale è

\[
 y^2=S(r),\qquad S(r)=rD(r)Q_2(r).
\]

Si sceglie un punto base \(r_0\), un cammino che non attraversi zeri o poli e
un foglio continuo di \(y=\sqrt S\). Tutte le primitive riportate sotto si
annullano in \(r_0\).

Il momento radiale congelato e il kernel esterno sono

\[
 p_{r0}=\frac{y}{\Delta D},
 \qquad
 K_{\rm off}
 =\left.\frac{\partial_{P_r}G_0}{H_{0,P_r}}\right|_{\rm shell}
 =\frac{A(r)}{y},
\]

\[
 A(r)=\frac{E^2Jr^4D(r)}{Q_2(r)}.
\]

## 2. Chiusura dell'azione radiale

L'identità algebrica essenziale è

\[
 \frac{S}{\Delta D}
 =\sum_{k=0}^{3}a_kr^k+\frac{C_\Delta(r)}{\Delta(r)}.
\]

I coefficienti sono

\[
\begin{aligned}
a_0&=2M\left[4E^2M^2-J^2+2(1-E^2)Ja\right],\\
a_1&=4E^2M^2-(E^2-1)J^2,\\
a_2&=2E^2M,\\
a_3&=E^2.
\end{aligned}
\]

Il numeratore residuo, lineare in \(r\), è

\[
\begin{aligned}
C_\Delta(r)={}&
\Big[
 E^2J^2a^2-8E^2JM^2a+16E^2M^4-4E^2M^2a^2\\
&\hspace{22mm}-J^2a^2+4M^2a^2
\Big]r\\
&+4E^2JMa^3-8E^2M^3a^2+2J^2Ma^2-4JMa^3 .
\end{aligned}
\]

Questo chiarisce un punto lasciato ambiguo in alcuni script precedenti: il
resto della divisione originale è \(D(r)C_\Delta(r)\). Il fattore \(D\) si
cancella contro il denominatore. Non esiste quindi una lettera di terza specie
nel punto

\[
 r_E=-\frac{2M}{E^2-1}.
\]

Gli unici poli di terza specie sono, nel caso non estremo,

\[
 r_\pm=M\pm\sqrt{M^2-a^2},
\qquad
 \rho_\pm=\frac{C_\Delta(r_\pm)}{r_\pm-r_\mp}.
\]

Definiamo

\[
 U_k(r)=\int_{r_0}^{r}\frac{x^k\,dx}{y(x)},\qquad
 \Pi_q(r)=\int_{r_0}^{r}\frac{dx}{(x-q)y(x)}.
\]

Allora l'azione radiale è

\[
 \boxed{
 I(r)=\int_{r_0}^{r}p_{r0}(x)\,dx
 =\sum_{k=0}^{3}a_kU_k(r)
  +\rho_+\Pi_{r_+}(r)+\rho_-\Pi_{r_-}(r).
 }
\]

Indicheremo con

\[
 \Pi(r)=\rho_+\Pi_{r_+}(r)+\rho_-\Pi_{r_-}(r)
\]

la parte di terza specie.

## 3. Riduzione di Hermite del kernel

La divisione del numeratore del kernel dà

\[
 \frac{E^2Jr^4D}{Q_2}=JD+\frac{N}{Q_2},
\]

\[
 N(r)=-JD(r)\bigl[Q_2(r)-E^2r^4\bigr].
\]

Cerchiamo un polinomio

\[
 P(r)=P_0+P_1r+P_2r^2+P_3r^3
\]

che soddisfi

\[
 \boxed{
 P(r)\,T(r)\,Q_2'(r)\equiv-2N(r)\pmod{Q_2(r)},
 \qquad T(r)=rD(r).
 }
\]

Per una curva generica questa congruenza ha soluzione unica. La specificazione
dei quattro coefficienti è completamente esplicita: se

\[
 H_{\ell m}
 =[r^\ell]\left(r^mTQ_2'\bmod Q_2\right),
\qquad
 b_\ell=-2[r^\ell]N,
\qquad 0\leq\ell,m\leq3,
\]

allora

\[
 (P_0,P_1,P_2,P_3)^{\mathsf T}=H^{-1}b.
\]

Equivalentemente, senza alcuna inversione implicita,

\[
 P_m=\frac{\det H_m}{\det H},
\]

dove \(H_m\) si ottiene sostituendo con \(b\) la colonna \(m\)-esima di \(H\).
Questa è una formula algebrica chiusa in \((M,a,E,J)\). Espandere i quattro
rapporti di determinanti produce espressioni molto più lunghe, ma non aggiunge
informazione matematica ed è numericamente meno stabile.

Definiamo poi

\[
 R_H(r)=
 \frac{
 2N-2P'Q_2T-PQ_2T'+PQ_2'T
 }{2Q_2}.
\]

La congruenza precedente garantisce che \(R_H\) sia un polinomio di grado non
superiore a quattro. Scrivendo

\[
 JD+R_H=\sum_{j=0}^{4}g_jr^j
\]

si ottiene l'identità differenziale

\[
 \boxed{
 \frac{A(r)}{y(r)}\,dr
 =d\!\left(\frac{P(r)y(r)}{Q_2(r)}\right)
  +\sum_{j=0}^{4}g_j\frac{r^j\,dr}{y(r)} .
 }
\]

Poniamo quindi

\[
 B(r)=\frac{P(r)y(r)}{Q_2(r)},\qquad
 \Phi_G(r)=\sum_{j=0}^{4}g_jU_j(r).
\]

## 4. Primitive elementare residua

L'integrazione per parti del differenziale esatto \(dB\) genera

\[
 E_P(r)=\int_{r_0}^{r}\frac{xP(x)}{\Delta(x)}\,dx.
\]

Questa parte è veramente elementare. Se

\[
 \frac{rP(r)}{\Delta(r)}
 =q_P(r)+\frac{c_+}{r-r_+}+\frac{c_-}{r-r_-},
\]

allora

\[
 c_\pm=\frac{r_\pm P(r_\pm)}{r_\pm-r_\mp}
\]

e

\[
\boxed{
\begin{aligned}
E_P(r)={}&Q_P(r)-Q_P(r_0)\\
&+c_+\log\frac{r-r_+}{r_0-r_+}
 +c_-\log\frac{r-r_-}{r_0-r_-},
\qquad Q_P'=q_P .
\end{aligned}
}
\]

Per cammini complessi i logaritmi devono essere continuati sullo stesso ramo
del cammino. Nel dominio reale esterno a \(r_+\) la formula è reale con i
logaritmi reali ordinari.

## 5. Formula chiusa completa del wrap di dilatazione

Introduciamo le lettere iterate

\[
 W_{jk}(r)
 =\int_{r_0}^{r}
 \left(U_k\,dU_j-U_j\,dU_k\right),
\]

\[
 D_{j,q}(r)
 =\int_{r_0}^{r}
 U_j(x)\frac{dx}{(x-q)y(x)}.
\]

Le \(W_{jk}\) sono antisimmetriche:
\(W_{jk}=-W_{kj}\). Le \(D_{j,q}\) sono le lettere miste di seconda/terza
specie responsabili del dilogaritmo iperellittico.

La forma chiusa cercata è

\[
\boxed{
\begin{aligned}
\Phi_{\mathrm{off,dil}}(r)
={}&-B(r)I(r)+E_P(r)\\
&-\frac12\sum_{j=0}^{4}\sum_{k=0}^{3}
 g_ja_k\left[U_j(r)U_k(r)+W_{jk}(r)\right]\\
&-\Phi_G(r)\Pi(r)
 +\sum_{q\in\{r_+,r_-\}}\rho_q
   \sum_{j=0}^{4}g_jD_{j,q}(r).
\end{aligned}
}
\]

Non è necessario aggiungere costanti: tutte le lettere e \(E_P\) sono
normalizzate a zero in \(r_0\), mentre \(I(r_0)=0\).

### Dimostrazione sintetica

Dal paragrafo 3,

\[
 K_{\rm off}\,dr=dB+d\Phi_G.
\]

Pertanto

\[
 -\int I\,dB=-BI+\int B\,dI.
\]

Usando \(dI=p_{r0}dr=y\,dr/(\Delta D)\), \(B=Py/Q_2\) e
\(y^2=rDQ_2\), segue

\[
 B\,dI=\frac{Pr}{\Delta}\,dr,
\]

che produce \(E_P\).

Per la parte polinomiale dell'azione,

\[
 \int U_k\,dU_j
 =\frac12\left(U_jU_k+W_{jk}\right).
\]

Per la parte di terza specie,

\[
 -\int \Pi\,d\Phi_G
 =-\Pi\Phi_G+\int\Phi_G\,d\Pi,
\]

e l'ultimo integrale è precisamente

\[
 \sum_q\rho_q\sum_jg_jD_{j,q}.
\]

La somma dei tre passaggi dà la formula nel riquadro.

## 6. Che cosa significa qui «forma chiusa»

Il risultato è chiuso nel seguente senso rigoroso:

1. tutti i coefficienti sono funzioni algebriche esplicite di
   \((M,a,E,J)\);
2. non rimane alcun kernel razionale non classificato;
3. la parte di peso uno è espressa in integrali abeliani standard sulla
   singola curva \(y^2=S\);
4. la parte di peso due è espressa nelle lettere canoniche
   \(W_{jk}\) e \(D_{j,r_\pm}\), cioè in polilogaritmi iperellittici di genere
   due;
5. l'unica primitiva elementare residua è stata risolta esplicitamente in
   polinomi e logaritmi.

Non sarebbe corretto affermare che il risultato si riduca genericamente alle
sole \(\sigma,\zeta,\wp\) Kleiniane di peso uno. I termini \(D_{j,r_\pm}\)
sono genuinamente di peso due. Una realizzazione mediante serie theta/nome,
utile per la valutazione numerica veloce, è una rappresentazione delle stesse
funzioni e non una nuova riduzione algebrica.

## 7. Verifiche eseguite

Lo script `_codex.py` esegue:

- l'identità simbolica generale per la riduzione dell'azione;
- la divisione simbolica generale del kernel;
- la verifica che il sistema \(4\times4\) codifichi esattamente la congruenza
  di Hermite;
- la riduzione di Hermite esatta, con aritmetica razionale, in due punti
  parametrici indipendenti;
- il confronto dell'azione ridotta con l'integrale diretto;
- il confronto della primitiva di Hermite con l'integrale diretto;
- il confronto della parte elementare logaritmica con la quadratura;
- il confronto dell'intera formula chiusa con l'integrale annidato che
  definisce \(\Phi_{\mathrm{off,dil}}\).

Risultati dell'esecuzione:

| Parametri | intervallo | errore azione | errore Hermite | errore elementare | errore assemblaggio |
|---|---:|---:|---:|---:|---:|
| \(M=1,\ a=9/10,\ E=7/5,\ J=6\) | \(11\to6.5\) | \(8.9\times10^{-16}\) | \(5.2\times10^{-15}\) | \(2.8\times10^{-14}\) | \(5.0\times10^{-14}\) |
| \(M=6/5,\ a=4/5,\ E=3/2,\ J=7\) | \(13\to7\) | \(3.6\times10^{-15}\) | \(7.2\times10^{-16}\) | \(7.1\times10^{-15}\) | \(4.4\times10^{-15}\) |

I valori completi del termine sono:

\[
\begin{aligned}
\Phi_{\mathrm{direct}}&=-1.985249985987,
&\Phi_{\mathrm{closed}}&=-1.985249985987,\\
\Phi_{\mathrm{direct}}&=-3.867100883395,
&\Phi_{\mathrm{closed}}&=-3.867100883395.
\end{aligned}
\]

Queste verifiche controllano l'identità matematica dell'assemblaggio. Non
sostituiscono la verifica separata della derivazione fisica
extended-Hamiltonian, già affidata nel progetto al confronto con il flusso
non autonomo.

## 8. Ipotesi e casi non ancora coperti

La formula nel riquadro vale direttamente nel caso generico:

- \(S\) è square-free;
- \(\gcd(TQ_2',Q_2)=1\), equivalently \(\det H\ne0\);
- \(M^2\ne a^2\), così \(r_+\ne r_-\);
- il cammino non attraversa una radice di \(S\) o un polo di terza specie.

Restano da trattare mediante limiti o riduzioni dedicate:

1. **caso estremo \(M^2=a^2\):** i due poli \(r_\pm\) confluiscono e la
   decomposizione semplice va sostituita con una lettera a polo doppio;
2. **separatrice \(J\to J_c\):** la curva diventa singolare e il limite deve
   essere combinato con il tracking della radice mobile. La presente formula
   chiude il termine generico prima del limite, ma non risolve da sola il
   residuo a polo triplo segnalato nel paper;
3. **branca \(\tau\):** ha una curva e coefficienti differenti; la stessa
   strategia è applicabile, ma questa derivazione non costituisce una prova
   della sua formula;
4. **implementazione theta/nome:** le lettere di peso due sono verificate qui
   mediante le loro quadrature canoniche. Un valutatore indipendente basato
   sulle serie theta sarebbe un upgrade numerico, non una correzione della
   formula.

## 9. Conseguenza per il paper

Per la branca generica \(t/\eta\), non è più corretto descrivere il wrap di
dilatazione come «non assemblato»: l'assemblaggio matematico è ora completo
nella base iterata abeliana.

È invece rigoroso scrivere che:

- la chiusura è iperellittica di genere due e di peso due;
- i coefficienti sono generali in \((M,a,E,J)\);
- la valutazione theta/nome esplicita può ancora essere implementata;
- il limite di separatrice con tracking della radice mobile resta un problema
  distinto e aperto;
- il risultato riguarda la parte di dilatazione e deve essere sommato agli
  altri blocchi della correzione adiabatica completa.
