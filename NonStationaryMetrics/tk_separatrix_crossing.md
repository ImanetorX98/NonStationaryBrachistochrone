# Separatrice TK t-branch off-shell — natura del problema: SEPARATRIX CROSSING (non un controtermine)

## Risultato definitivo
Il "controtermine di tracking `dr_d/dλ`" che il paper (App C) ipotizzava **NON esiste come tale**.
Test concreto (`ThakurtaMetric/tk_sep_tracking_counterterm.py`): sotto l'evoluzione adiabatica il
raggio di dilatazione `(E_eff,J_eff)=(Ê/A, J/A)` **lascia la separatrice**, quindi il doppio root
non si muove — si **DISSOLVE** (Q2 perde il doppio root fuori dalla separatrice). L'orbita
**ATTRAVERSA** la separatrice.

## Numeri (E=1.2, a=0.9)
- Separatrice: `Jc=2.93635, r_d=1.51229`, `dJc/dE=−0.115`.
- A=1.05: `J_eff=Jc/A=2.797`, ma `Jc(E_eff=Ê/A)=2.943` → `J_eff<Jc` → SOTTO la separatrice (crossed).
- A=1.10: `J_eff=2.669` vs `Jc(E_eff)=2.950` → sempre più sotto.
- Meccanismo: `E_eff=Ê/A` scende con A; `dJc/dE<0` ⟹ `Jc(E_eff)` sale; `J_eff=J/A` scende ⟹ l'orbita
  cade sotto la separatrice mobile.

## Interpretazione (corretta)
La non-uniformità off-shell a `r_d` (polo di potenza `1/(r−r_d)²`, residuo ∝ ΔS(r_d)=34.75, dal KERNEL
`A~1/Q2²` non dalla sorgente) NON è un polo residuo cancellabile da un tracking algebrico. È il
segnale che l'espansione adiabatica al primo ordine è **non-uniforme all'attraversamento della
separatrice** — il problema di **separatrix crossing** (Neishtadt, Timofeev): il cambiamento
dell'invariante adiabatico attraverso la separatrice ha struttura logaritmica/salto, richiede analisi
di strato limite (inner/outer matched asymptotics), NON una forma chiusa a coefficienti simbolici.

## Conseguenza per il paper (App C)
Sostituire "the moving-root separatrix tracking must be extended to cancel the residual triple pole"
con: "*sotto l'evoluzione adiabatica l'orbita attraversa la separatrice (il doppio root si dissolve, non
si muove lungo il raggio di dilatazione); la non-uniformità off-shell a r_d è un fenomeno di
separatrix crossing (Neishtadt), fuori dallo scope della riduzione a coefficienti simbolici — non un
controtermine algebrico mancante.*" Onesto e definitivo: chiude la domanda identificandone la natura.

## Cosa resterebbe (ricerca a sé, opzionale)
L'analisi completa di separatrix-crossing (jump dell'invariante adiabatico, forma matched) è un
problema classico ma non banale; darebbe la correzione uniforme vicino a r_d. Non necessaria per la
correttezza del paper: basta dichiarare che il limite off-shell è separatrix crossing, non chiuso in
forma speciale. Vaidya immune (r_d off-path, nessun crossing).
