# Issue 5 — Lemma B / no-inversion a estremi fissi: status ONESTO + fix di testo

## Status logico (da propagare al paper)
- **Lemma A — DIMOSTRATO** (analitico, chiuso): `sign(F_t²−F_tau²)=sign(N_G)`,
  `N_G=(r−r_min)[4rΔ/(r_min−2M)]W_G`, `W_G` lineare e positivo nel dominio esterno.
  Script `no_inversion_reduction.py`. (referee: "Lemma A is a genuine algebraic result".)
- **Lemma B — PARZIALE**: monotonia di `Φ_tau(r_min)` provata in un sotto-regime (grazing + stima
  elementare `V_min≥V_0/4`), criterio single-crossing `(DAGGER)` derivato e verificato numericamente ad
  alta precisione; l'intervallo tra il picco di deflessione `r_pk` e il regime elementare resta APERTO;
  per Kerr la monotonia resta numerica; la transizione è transcendente (niente certificato polinomiale).

## Tesi
A+B ⟹ fixed-endpoint no-inversion (`r_min^tau < r_min^t`). Con Lemma B solo parziale, la tesi è
**CONDIZIONALE** fuori dal sotto-regime provato.

## Fix di TESTO richiesti (referee Issue 5) — per propagazione
1. NON dire che Lemma B è "un teorema nel regime di scattering": è provato nel sotto-regime,
   numericamente supportato altrove. Usare linguaggio condizionale.
2. Togliere "never" / "mai inversione" salvo prova globale.
3. Togliere dalla tabella (protocollo) il rapporto puntuale `n_t/n_tau=E/f>1` come "prova": il testo
   stesso spiega che il rapporto puntuale NON ordina i minimi di due funzionali diversi (la mappa
   `J→Φ` ribalta l'ordine). Contraddizione tabella-derivazione.
4. Opzioni referee: (A) prova analitica completa; (B) teorema computer-assisted su dominio compatto
   (interval arithmetic); (C) teorema condizionale + evidenza numerica. **Consigliato: (C)** con abstract,
   caption, tabella, conclusioni allineati allo stesso status.

Riferimenti: `TODO_lemmaB_closure.md`, `no_inversion_reduction.py`, `no_inversion_schwarzschild_*.py`.
