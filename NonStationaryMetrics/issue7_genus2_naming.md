# Issue 7 — genus-2 Kronecker–Eisenstein naming: status + fix

## Esito
In gran parte **GIÀ affrontato** nel testo; corretti gli over-claim dell'abstract. L'identificazione
COMPLETA con la classe higher-genus canonica (Baune/D'Hoker) resta **frontiera** (progetto a sé),
onestamente dichiarata "open" nel paper.

## Cosa era già onesto (main.tex §generic t-branch, righe ~1218-1251)
- Le `W_kj` chiamate "length-two iterated Abelian integrals" (termine neutro richiesto dal referee).
- Dichiara esplicitamente che "genus-two dilogarithm" richiederebbe il framework higher-genus
  (Baune/D'Hoker: base points, paths, tangential reg., monodromy, single-valued completion) **prima**
  di essere theorem-level, e NON lo rivendica.
- Brown–Levin/Kronecker–Eisenstein citati come analogo **ellittico g=1**, non come identità genus-2.
- Irriducibilità e dimensione base tenute **congetturali**; rank-5 "strong evidence, NOT a proof".
- Ostruzione theta-divisor `Θ=W_{g-1}` spiegata (perché NON riduce a ζ,σ Kleiniane peso-1).
- Le occorrenze "Kronecker–Eisenstein" nel contesto SEPARATRICE (genus-1) sono CORRETTE (Brown–Levin
  è genuinamente genus-1 lì).

## Fix applicati (abstract, main.tex + PRD)
1. "Every symbolic identity is independently confirmed in Mathematica" → "**Algebraic** coefficient
   identities in Mathematica; the genus-two theta/nome evaluations use **Sage/abelfunctions**"
   (allinea con VERIFICATION_STATUS: θ genus-2 non nativa in Mathematica).
2. "explicit off-shell closed form ... remain open" → "off-shell assemblato in forma chiusa (iterati
   abeliani depth-2 + un pezzo elementare) su TUTTE le 4 branch, verificato vs vero flusso a O(ε²);
   solo l'identificazione con classe higher-genus nominata è deferred" (riflette la chiusura fatta).
Entrambi compilano (main 70pp, PRD 33pp, 0 undefined).

## Resta APERTO (frontiera, dichiarato)
- Teorema di identificazione delle `W_kj`/`D_{j,root}` genus-2 con la classe Baune/D'Hoker
  (arXiv:2306.08644, 2407.11476): rappresentazione theta esplicita + mappa letter-by-letter +
  convergenza/monodromia. Richiede l'env Sage con `abelfunctions` (non installato qui) ed è materiale
  da paper dedicato. Il paper lo dichiara onestamente "open".
- Terminologia U_k: nel contesto separatrice (genus-1) "third-kind elliptic" è corretto; nel generico
  (genus-2) il paper usa già "iterated Abelian integrals". OK.
