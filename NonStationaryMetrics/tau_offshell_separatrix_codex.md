# Chiusura off-shell della branca tau e audit delle separatrici

## Risultato complessivo

Sono stati completati i quattro passaggi richiesti:

1. il wrap off-shell di dilatazione della branca \(\tau\) è stato derivato dal
   vero Hamiltoniano proprio e chiuso con coefficienti simbolici generali in
   \((M,a,E,J)\);
2. la forma speciale è stata confrontata con la quadratura annidata in due
   configurazioni indipendenti;
3. lo stesso sottoblocco è stato estratto dalla perturbazione canonica e
   collegato alla pendenza fisica \(2.073471\);
4. le degenerazioni di separatrice sono state classificate e ne è stata
   determinata la struttura singolare.

File di verifica:

- `KerrSessionScripts/offshell_taubranch_closed_form_codex.py`;
- `KerrSessionScripts/offshell_taubranch_physics_anchor_codex.py`;
- `ThakurtaMetric/separatrix_offshell_audit_codex.py`.

## 1. Hamiltoniano tau corretto

La derivazione parte dall'Hamiltoniano proprio congelato

\[
H_\tau=
\widetilde J\,\varphi'_0+
\sqrt{\frac{\Delta v}{\bar P}}\,
\sqrt{\frac{\Delta}{r^2}P_r^2+\frac{\widetilde J^2}{\bar P}}
-\frac{f}{E},
\]

\[
f=1-\frac{2M}{r},\qquad
v=1-\frac{f}{E^2},\qquad
\widetilde J=J-\frac{2Ma}{Er}.
\]

Un vecchio controllo in `adiabatic_offshell_reduction.py` impiega invece
\(-1\) come ultimo termine anche nel ramo \(\tau\). Quella configurazione non
è quindi la riduzione del medesimo Hamiltoniano proprio usato dal test fisico.
I nuovi file `_codex` non utilizzano tale scorciatoia.

## 2. Curva spettrale tau

Definiamo

\[
D(r)=(E^2-1)r+2M,\qquad
\Delta(r)=r^2-2Mr+a^2,
\]

\[
Q_3(r)=r\Delta(r)-J^2D(r).
\]

Sulla shell congelata:

\[
\boxed{
y^2=S_\tau(r)=r(r-2M)D(r)Q_3(r)
}
\]

e, per l'arco entrante,

\[
\boxed{
P_{r0}=-\frac{y}{\Delta D}.
}
\]

La scelta \(P_{r0}=+y/(\Delta D)\) definisce invece la convenzione algebrica
positiva usata nell'assemblaggio. I due risultati differiscono per il segno
globale dell'azione radiale.

## 3. Azione radiale tau

L'identità simbolica generale è

\[
\frac{S_\tau}{\Delta D}
=\sum_{k=0}^{3}a_kr^k+\frac{C_\Delta(r)}{\Delta(r)},
\]

con

\[
\boxed{
a_0=-2MJ^2,\qquad
a_1=-(E^2-1)J^2,\qquad
a_2=-2M,\qquad
a_3=1
}
\]

e

\[
C_\Delta(r)=J^2a^2D(r).
\]

Nel caso non estremo

\[
r_\pm=M\pm\sqrt{M^2-a^2},
\]

\[
\rho_\pm=\frac{C_\Delta(r_\pm)}{r_\pm-r_\mp}.
\]

Ponendo

\[
U_k(r)=\int_{r_0}^{r}\frac{x^k\,dx}{y(x)},\qquad
\Pi_q(r)=\int_{r_0}^{r}\frac{dx}{(x-q)y(x)},
\]

l'azione a radicale positivo è

\[
\boxed{
I_+(r)=
\sum_{k=0}^{3}a_kU_k(r)
+\rho_+\Pi_{r_+}(r)+\rho_-\Pi_{r_-}(r).
}
\]

L'azione fisica entrante è \(I_{\rm in}=-I_+\).

## 4. Kernel off-shell tau

La cancellazione dei radicali interni fornisce

\[
\boxed{
\left.
\frac{\partial_{P_r}G_\tau}{H_{\tau,P_r}}
\right|_{\rm shell}
=
\frac{A_\tau(r)}{y(r)}
}
\]

con

\[
\boxed{
A_\tau(r)=
\frac{r^2D(r)\,[EJr-2Ma]}{Q_3(r)}.
}
\]

Questa identità è stata controllata direttamente contro le derivate del
Hamiltoniano a precisione macchina.

La divisione polinomiale è

\[
\frac{r^2D(EJr-2Ma)}{Q_3}
=A_{\rm pol}(r)+\frac{N_\tau(r)}{Q_3(r)},
\]

\[
A_{\rm pol}
=EJ(E^2-1)r
+2E^3JM-2E^2Ma+2Ma.
\]

## 5. Riduzione di Hermite

Si cerca

\[
P(r)=P_0+P_1r+P_2r^2
\]

tale che

\[
\boxed{
P(r)\,T(r)\,Q_3'(r)\equiv-2N_\tau(r)\pmod{Q_3(r)},
}
\]

\[
T(r)=r(r-2M)D(r).
\]

I coefficienti sono determinati dal sistema simbolico \(3\times3\)

\[
H_{\ell m}
=[r^\ell]\left(r^mTQ_3'\bmod Q_3\right),
\qquad
b_\ell=-2[r^\ell]N_\tau,
\]

\[
(P_0,P_1,P_2)^{\mathsf T}=H^{-1}b.
\]

Quindi

\[
P_m=\frac{\det H_m}{\det H}
\]

è una funzione razionale esplicita di \((M,a,E,J)\).

Ponendo

\[
R_H=
\frac{
2N_\tau-2P'Q_3T-PQ_3T'+PQ_3'T
}{2Q_3},
\]

la congruenza garantisce che \(R_H\) sia polinomiale. Definendo

\[
A_{\rm pol}+R_H=\sum_{j=0}^{4}g_jr^j,
\]

si ottiene

\[
\boxed{
\frac{A_\tau}{y}\,dr
=d\!\left(\frac{Py}{Q_3}\right)
+\sum_{j=0}^{4}g_j\frac{r^j\,dr}{y}.
}
\]

## 6. Formula speciale completa tau

Definiamo

\[
B=\frac{Py}{Q_3},\qquad
\Phi_G=\sum_{j=0}^{4}g_jU_j,
\]

\[
\Pi=\rho_+\Pi_{r_+}+\rho_-\Pi_{r_-},
\]

\[
W_{jk}
=\int_{r_0}^{r}(U_k\,dU_j-U_j\,dU_k),
\]

\[
\mathcal D_{j,q}
=\int_{r_0}^{r}
U_j(x)\frac{dx}{(x-q)y(x)}.
\]

La parte elementare è

\[
E_P^\tau(r)
=\int_{r_0}^{r}
\frac{x(x-2M)P(x)}{\Delta(x)}\,dx,
\]

che si riduce a un polinomio più logaritmi in \(r-r_\pm\).

Per la convenzione \(P_r=+y/(\Delta D)\):

\[
\boxed{
\begin{aligned}
\Phi_{\tau,\mathrm{dil}}^{(+)}
={}&-BI_++E_P^\tau\\
&-\frac12\sum_{j=0}^{4}\sum_{k=0}^{3}
g_ja_k\left(U_jU_k+W_{jk}\right)\\
&-\Phi_G\Pi
+\sum_{q=r_\pm}\rho_q
\sum_{j=0}^{4}g_j\mathcal D_{j,q}.
\end{aligned}
}
\]

Per l'arco fisico entrante:

\[
\boxed{
\Phi_{\tau,\mathrm{dil}}^{\rm in}
=-\Phi_{\tau,\mathrm{dil}}^{(+)}.
}
\]

## 7. Verifica contro la quadratura

### Caso 1

\[
M=1,\quad a=\frac9{10},\quad
E=\frac75,\quad J=\frac52,\quad 12\to6.
\]

\[
\Phi_{\rm direct}^{(+)}
=-2.544217613489,
\]

\[
\Phi_{\rm closed}^{(+)}
=-2.544217613489,
\]

con errore

\[
1.02\times10^{-14}.
\]

### Caso 2

\[
M=\frac65,\quad a=\frac45,\quad
E=\frac32,\quad J=3,\quad14\to8.
\]

\[
\Phi_{\rm direct}^{(+)}
=-1.979110119216,
\]

\[
\Phi_{\rm closed}^{(+)}
=-1.979110119216,
\]

con errore

\[
4.00\times10^{-15}.
\]

Sono inoltre verificati separatamente:

| Controllo | caso 1 | caso 2 |
|---|---:|---:|
| azione radiale | \(0\) | \(8.9\times10^{-16}\) |
| Hermite | \(9.4\times10^{-16}\) | \(1.9\times10^{-15}\) |
| primitiva elementare | \(1.4\times10^{-14}\) | \(7.1\times10^{-15}\) |

## 8. Collegamento alla fisica

Sul medesimo arco fisico \(12\to6\):

\[
\Phi_{\tau,\mathrm{dil}}^{\rm closed,in}
=+2.544217613489,
\]

mentre l'estrazione indipendente dal vettore di perturbazione canonica dà

\[
\Phi_{\tau,\mathrm{dil}}^{\rm PT}
=+2.544217621260.
\]

La differenza

\[
7.77\times10^{-9}
\]

è il residuo delle due integrazioni cumulative trapezoidali impiegate per
estrarre il sottoblocco dalla soluzione ODE.

Per la correzione completa:

| \(\varepsilon\) | residuo contro il flusso originale |
|---:|---:|
| \(10^{-3}\) | \(1.70914\times10^{-7}\) |
| \(2\times10^{-3}\) | \(6.96661\times10^{-7}\) |
| \(4\times10^{-3}\) | \(2.91376\times10^{-6}\) |
| \(8\times10^{-3}\) | \(1.27710\times10^{-5}\) |

Il fit log--log dà

\[
\boxed{p_\tau=2.073471}.
\]

Senza dilatazione:

\[
\boxed{p_{\tau,\mathrm{old}}=1.001141}.
\]

Le variabili canoniche originali e quelle normalizzate con il termine
\(-\alpha P_r\) concordano entro

\[
1.96\times10^{-15}.
\]

## 9. Degenerazione tau \(J=a/E\)

Nel caso

\[
J=\frac aE
\]

si ha esattamente

\[
\boxed{
Q_3(r)
=(r-2M)\frac{E^2r^2+a^2}{E^2}.
}
\]

La sestica degenera in

\[
S_\tau
=(r-2M)^2Q_4(r),
\]

\[
\boxed{
Q_4(r)
=\frac{rD(r)(E^2r^2+a^2)}{E^2}.
}
\]

Ponendo \(y=(r-2M)Y\), \(Y^2=Q_4\), il numeratore del kernel contiene
lo stesso fattore \(r-2M\), e

\[
\boxed{
A_{\tau,\mathrm{sep}}
=\frac{E^2ar^2D}{E^2r^2+a^2}.
}
\]

Di conseguenza

\[
K_{\rm off}
=\frac{A_{\tau,\mathrm{sep}}}{(r-2M)Y}
\sim\frac{\text{costante}}{r-2M}.
\]

Non compare il polo cubico.

L'azione ellittica si riduce a

\[
dI_+
=R_I(r)\frac{dr}{Y},
\]

\[
\boxed{
R_I(r)=
\frac{r(r-2M)(E^2r^2+a^2)}{E^2\Delta(r)}.
}
\]

La divisione simbolica è

\[
R_I=
\left(r^2-a^2+\frac{a^2}{E^2}\right)
-\frac{
a^2[2E^2Mr-E^2a^2+a^2]
}{E^2\Delta}.
\]

Definendo l'integrale ellittico iterato

\[
\mathbb E[R_1,R_2](r)
=\int_{r_0}^{r}
\left[
\int_{r_0}^{x}R_1(z)\frac{dz}{Y(z)}
\right]
R_2(x)\frac{dx}{Y(x)},
\]

la chiusura canonica è

\[
\boxed{
\Phi_{\tau,\mathrm{sep,dil}}^{(+)}
=-\mathbb E[R_I,R_K],
}
\]

\[
R_K(r)
=\frac{E^2ar^2D}
{(r-2M)(E^2r^2+a^2)}.
\]

I poli di \(R_I\) sono \(r_\pm\); quelli di \(R_K\) sono \(2M\) e i
due punti complessi \(\pm ia/E\). Pertanto le primitive di peso uno sono
combinazioni standard di \(\wp,\zeta,\sigma\), mentre
\(\mathbb E[R_I,R_K]\) è un dilogaritmo ellittico.

Questa degenerazione è quindi chiusa nella classe ellittica.

### Interpretazione corretta

\(J=a/E\) è una collisione tra il fattore \(r-2M\) della curva e una
radice di \(Q_3\). Nel linguaggio del referaggio è una degenerazione
algebrica/soglia ergosferica, non automaticamente la separatrice esterna
del sistema dinamico.

Per \(M=1,a=0.9,E=1.4\), le vere doppie radici di \(Q_3\) a \(J>0\)
sono:

\[
\begin{array}{c|c|c}
J & r_d & \text{classificazione}\\ \hline
0.204287249715 & 0.233308738193 & \text{interna a }r_-\\
7.180679377301 & -3.416755239862 & \text{raggio negativo}.
\end{array}
\]

Nessuna delle due è una separatrice esterna per questa configurazione.

## 10. Separatrice fisica t/eta

Per

\[
M=1,\qquad a=0.9,\qquad E=1.2
\]

la doppia radice fisica è

\[
\boxed{
J_c=2.9363507418674620791,
\qquad
r_d=1.5122917132047833253.
}
\]

La famiglia di doppie radici soddisfa

\[
\frac{dJ_c}{dE}=-0.1150284688798347261,
\]

\[
\frac{dr_d}{dE}=0.0514609415651983001.
\]

Localmente, con \(x=r-r_d\),

\[
K_{\rm off}
=\frac{\kappa_3}{x^3}+O(x^{-2}),
\]

\[
\kappa_3=0.4140411280683596912.
\]

L'azione fisica accumulata tende a

\[
I_d=34.75044976579203842\neq0.
\]

Quindi

\[
-K_{\rm off}I
\sim-\frac{14.38811542191140141}{x^3},
\]

e la primitiva esterna contiene

\[
\boxed{
\Phi_{\mathrm{off,dil}}
\sim\frac{7.194057710955700706}{x^2}.
}
\]

## 11. Conseguenza per il tracking

La forma congelata di separatrice ha comportamento

\[
\varphi_0\sim c\log|x|.
\]

Una traslazione liscia della radice

\[
r_d\longrightarrow r_d+\varepsilon\,\delta r_d
\]

produce

\[
\delta\varphi_{\rm shift}
\sim-\frac{c\,\delta r_d}{x},
\]

cioè solamente un polo semplice. Non può cancellare da sola il polo doppio
\(x^{-2}\) del wrap di dilatazione.

Una chiusura uniforme è possibile soltanto in uno dei due modi:

1. la trasformazione canonica completa nella coordinata comovente
   \(x=r-r_d(s)\) produce un nuovo integrando cubico con coefficiente
   esattamente opposto a \(14.3881154219\);
2. il coefficiente non cancella, e allora la perturbazione regolare fallisce
   nella regione

   \[
   |r-r_d|=O(\sqrt{\varepsilon}),
   \]

   che deve essere trattata con un problema interno e matching asintotico.

Il solo tracking \(dJ_c/dE\) o la sola sostituzione \(r-r_d(s)\) non
costituiscono una chiusura.

## 12. Stato finale

- La branca generica \(t/\eta\) è chiusa in forma iperellittica.
- La branca generica \(\tau\) è ora chiusa allo stesso livello.
- Entrambe sono collegate a verifiche fisiche con pendenza circa due.
- La degenerazione \(\tau\), \(J=a/E\), è chiusa in forma ellittica, ma va
  chiamata soglia ergosferica/algebrica.
- La separatrice fisica \(t/\eta\) è chiusa come espressione esterna
  ellittica lontano da \(r_d\), ma non possiede una correzione lineare
  uniforme ottenibile mediante una semplice traslazione della radice.
- Il prossimo problema matematico è quindi la derivazione del Hamiltoniano
  canonico nella coordinata comovente e, se il polo non cancella, la
  costruzione esplicita del boundary layer \(O(\sqrt\varepsilon)\).

## 13. Addendum: termine canonico terminale a polo cubico

La derivazione successiva nella coordinata comovente ha separato due problemi
al contorno distinti.

Per il problema canonico con arrivo libero e ancoraggio terminale sulla
separatrice, la costante d'integrazione della variazione di shell impone

\[
I_T(r)=\int_{r_d}^{r}p_{r0}(u)\,du
=I_L(r)-I_L(r_d).
\]

Rispetto alla formula ancorata al lancio ciò aggiunge

\[
\frac{d\Phi_{\rm can}}{dr}
=I_dK_{\rm off}(r),
\qquad
I_d=I_L(r_d).
\]

Il suo coefficiente cubico è \(+\kappa_3I_d\), esattamente opposto a
\(-\kappa_3I_d\) del wrap ancorato al lancio. La somma ha soltanto un polo
semplice e una primitiva logaritmica.

Questa cancellazione non contraddice le conclusioni delle sezioni 10--12:
quelle riguardano il problema IVP a lancio fissato, per il quale non è lecito
cambiare la costante d'integrazione e il boundary layer resta necessario.

La derivazione, i coefficienti simbolici e i controlli su entrambe le branche
fisiche sono in:

- `separatrix_canonical_cubic_codex.md`;
- `ThakurtaMetric/separatrix_canonical_cubic_codex.py`.
