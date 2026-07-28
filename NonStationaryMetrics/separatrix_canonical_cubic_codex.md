# Termine canonico a polo cubico sulla separatrice \(t/\eta\)

## Risultato

Il termine canonico richiesto può essere aggiunto in modo non arbitrario
quando il problema è formulato come problema variazionale con arrivo libero e
ancoraggio terminale sulla separatrice.

Sia

\[
x=r-r_d,
\qquad
Q_2(r)=(r-r_d)^2q(r),
\qquad
Y^2=rDq,
\qquad
D=(E^2-1)r+2M.
\]

Sul foglio algebrico positivo,

\[
\sqrt S=xY,
\]

mentre per il ramo fisico entrante

\[
p_{r0}(r)=-\frac{xY}{\Delta D},
\qquad
\Delta=r^2-2Mr+a^2.
\]

Il kernel off-shell è

\[
K_{\rm off}(r)
=
\frac{E^2Jr^4D}
{x^3q(r)Y(r)}
=\frac{\kappa_3}{x^3}+O(x^{-2}),
\]

con

\[
\boxed{
\kappa_3=
\frac{E^2Jr_d^4D_d}
{q_d\sqrt{r_dD_dq_d}},
\qquad
q_d=\frac12Q_2''(r_d).
}
\]

## Ancoraggio al lancio e ancoraggio terminale

La formula generica verificata sulle sottotratte regolari usa l'azione
ancorata al lancio

\[
I_L(r)=\int_{r_0}^{r}p_{r0}(u)\,du.
\]

Alla doppia radice,

\[
I_d=I_L(r_d)
\]

è finita e generalmente non nulla. Il wrap di dilatazione è

\[
\frac{d\Phi_{\rm dil,L}}{dr}
=-K_{\rm off}(r)I_L(r),
\]

e pertanto contiene

\[
\frac{d\Phi_{\rm dil,L}}{dr}
\sim-\frac{\kappa_3I_d}{x^3}.
\]

Nel problema canonico di separatrice con condizione terminale, la variazione
della shell deve invece essere ancorata al punto terminale:

\[
I_T(r)
=I_L(r)-I_d
=\int_{r_d}^{r}p_{r0}(u)\,du.
\]

Questo cambio dell'integrazione costante aggiunge precisamente

\[
\boxed{
\frac{d\Phi_{\rm can}}{dr}
=I_dK_{\rm off}(r).
}
\]

Il suo coefficiente cubico è

\[
\boxed{
C_3^{\rm can}=+\kappa_3I_d,
}
\]

opposto al coefficiente

\[
C_3^{\rm dil,L}=-\kappa_3I_d.
\]

Quindi

\[
C_3^{\rm dil,L}+C_3^{\rm can}=0
\]

esattamente.

La forma combinata è

\[
\boxed{
\frac{d\Phi_{\rm dil,T}}{dr}
=-K_{\rm off}(r)
\left[I_L(r)-I_d\right]
=-K_{\rm off}(r)
\int_{r_d}^{r}p_{r0}(u)\,du.
}
\]

## Ordine del polo residuo

Poiché

\[
p_{r0}(r)
=p_1x+O(x^2),
\qquad
p_1=-\frac{\sqrt{r_dD_dq_d}}{\Delta_dD_d},
\]

si ha

\[
I_T(r)=\frac{p_1}{2}x^2+O(x^3).
\]

Pertanto

\[
\frac{d\Phi_{\rm dil,T}}{dr}
=\frac{C_{\log}}{x}+O(1),
\]

dove

\[
\boxed{
C_{\log}
=-\frac{\kappa_3p_1}{2}
=\frac{E^2Jr_d^4}{2q_d\Delta_d}.
}
\]

La primitiva non ha più un polo doppio: resta il comportamento

\[
\Phi_{\rm dil,T}
=C_{\log}\log|r-r_d|+O(1),
\]

cioè lo stesso tipo di non uniformità logaritmica della forma congelata.

## Interpretazione canonica

La coordinata comovente è generata da

\[
F_2(r,P_x;s)=[r-r_d(s)]P_x,
\qquad
x=r-r_d(s),
\]

e il nuovo Hamiltoniano contiene

\[
K=H-\dot r_dP_x.
\]

Il termine \(-\dot r_dP_x\) fissa coerentemente il confronto a \(x\)
costante, ma da solo non autorizza una sottrazione locale scelta a piacere.
Il contributo cubico sopra nasce dalla costante d'integrazione selezionata
dalla condizione canonica terminale. Con questa condizione,
\(I_T(r_d)=0\), il termine è univocamente determinato e il suo segno non è
libero.

Equivalentemente,

\[
\Phi_{\rm can}(r)
=I_d\int_{r_0}^{r}K_{\rm off}(u)\,du.
\]

L'integrale di \(K_{\rm off}\) è un integrale abeliano ellittico di seconda
e terza specie; il polo cubico del differenziale produce un termine
\(\wp\) nella primitiva. Anche \(I_d\) è un'azione ellittica. Il
controtermine canonico è quindi chiuso nella stessa base
\(\wp,\zeta,\sigma\) già usata per la separatrice.

## Verifica numerica

La verifica è implementata in
`ThakurtaMetric/separatrix_canonical_cubic_codex.py`.

Per ciascuna delle due doppie radici fisiche a

\[
M=1,\qquad a=0.9,\qquad E=1.2,
\]

lo script controlla:

1. \(x^3K_{\rm off}\to\kappa_3\);
2. \(I_T/x^2\to p_1/2\);
3. \(x^3\,d\Phi_{\rm dil,L}/dr\to-\kappa_3I_d\);
4. \(x^3\,d\Phi_{\rm can}/dr\to+\kappa_3I_d\);
5. il coefficiente cubico della somma tende a zero;
6. \(x\,d\Phi_{\rm dil,T}/dr\to C_{\log}\).

Sono trattate sia la branca prograda ergosferica
\((J_c\simeq2.93635,r_d\simeq1.51229)\), sia la branca retrograda esterna
\((J_c\simeq-8.05352,r_d\simeq3.51391)\).

I coefficienti ottenuti sono:

| branca | \(I_d\) | \(C_3^{\rm dil,L}\) | \(C_3^{\rm can}\) | \(C_{\log}\) |
|---|---:|---:|---:|---:|
| prograda ergosferica | \(34.7504497658\) | \(-14.3881154219\) | \(+14.3881154219\) | \(+8.9097450808\) |
| retrograda esterna | \(25.6576965285\) | \(+64.5218683294\) | \(-64.5218683294\) | \(-1.8184607884\) |

## Limite della conclusione

La cancellazione algebrica è rigorosa una volta imposta la condizione
canonica terminale e rende chiusa la struttura singolare esterna della
separatrice.

Non bisogna aggiungere lo stesso termine al problema a lancio fissato usato
nei test generici: in quel problema la condizione corretta è
\(I_L(r_0)=0\), il polo doppio non viene rimosso e la regione
\(|r-r_d|=O(\sqrt\varepsilon)\) richiede ancora un'analisi interna.

Per trasformare questo risultato in un teorema completo del paper resta
necessario derivare la condizione terminale dal problema variazionale
free-arrival/fixed-endpoint completo, inclusa la correzione \(J_1\) imposta
dall'endpoint angolare.
