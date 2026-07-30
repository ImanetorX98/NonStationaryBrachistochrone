# Referee follow-up su main12 — analisi + azioni (30 lug 2026)

Fonte: `Referee_Followup_main12_Adiabatic_Reconstructibility.pdf` (5 pp).

## Verdetto del referee (POSITIVO)
"Responded successfully to the principal technical referee criticism." La dinamica
first-order Thakurta–Kerr corretta ora include la dilation del momento radiale
(`−(A_s/A)P_r`, Euler completo `D=E∂_E+J∂_J+P_r∂_Pr`) e passa il test O(ε²) vs vero
flusso. Frozen brachistochrones, on-shell adiabatico, e risposta locale alla
separatrice sono ricostruibili. Restano azioni editoriali/di scoping, non errori.

## §2 — CONTRADDIZIONE INTERNA (blocker, azione 1)
Abstract dice off-shell "assembled in closed form on all four branches"; Table 1
(tab:claims righe 156-157) dice "Explicit off-shell special-function closed form —
**open**; assembly deferred". Non possono coesistere.
- Realtà nostra: 4 branch (TK t/τ, Vaidya τ/v) ASSEMBLATE+verificate (1e-14…5e-17)+
  physics-anchored; solo il *naming* canonico theta/polylog resta (ostruzione
  theta-divisor). → **Opzione 1 del referee** (assembly completo; solo naming/
  single-valued completion/irriducibilità/basi aperti). Table 1 va allineata all'abstract.

## §10 — Azioni raccomandate (8)
1. [ ] Eliminare contraddizione abstract/status vs Table 1 sull'assembly off-shell. (§2)
2. [ ] Dire esplicitamente che la formula first-order completa è uniformemente valida
   solo su sub-archi regolari con `H_Pr≠0` (al turning point `H_Pr=dr/ds=0` → non uniforme).
3. [ ] Aggiungere sottosezione con le CONVENZIONI di curva: base point, sheet, cuts,
   orientamento path, costanti additive per ogni `U_k`, `W_kj`, e ogni lettera di 3ª specie.
4. [ ] Chiarire quali ingredienti sono ricostruibili dal solo PDF vs quali richiedono
   il repository (in particolare la q-series theta-nome + normalizzazione, da stampare per intero).
5. [ ] Aggiungere la correzione fixed-endpoint `J_1 = −δΦ_dyn/∂_J Φ_0` (da
   `J=J_0+εJ_1+O(ε²)`), OPPURE dichiarare che le formule attuali sono fixed-charge/same-launch.
6. [ ] Separare "correzione adiabatica ISTANTANEA alla separatrice" da "crossing dinamico/
   long-time following di una separatrice in movimento" (§6.3: non stabilito; serve inner
   problem di Neishtadt + jump/capture law + matched asymptotics). Già abbiamo il linguaggio
   separatrix-crossing — verificare che i due siano distinti chiaramente.
7. [ ] Presentare η (tempo conforme) come la rappresentazione CANONICA del ramo di ARRIVO
   di Thakurta–Kerr (NON una terza famiglia spaziale indipendente: `dt=A(η)dη`, monotono ⟹
   `argmin η_f = argmin t_f`, stesso path). τ = ramo proper-time genuinamente distinto.
8. [ ] NON introdurre un tempo t "alla Schwarzschild" come terzo clock canonico in Vaidya
   (t non canonico quando `m=m(v)`: il tortoise `r_*` dipende da v). Trattare u (retarded)
   come problema outgoing-Vaidya SEPARATO (geometria/Hamiltoniane/costati/endpoint diversi).

## §3-5 — Ricostruibilità (verdetti, per riferimento)
- Frozen brachistochrones: ricostruibili (Hamiltoniana/integrando + radicale spettrale + numeratore).
- On-shell adiabatico: ricostruibile (identità polinomiale `2A'S−AS'+2S Σc_k r^k=2N_λ`,
  sistema lineare 11×11, niente fit; poi `U_k`, `W_kj`, pesi algebrici + resti log elementari).
- Complete first-order: ricostruibile come nested quadrature su archi regolari (`H_Pr≠0`);
  qualificazioni = convenzioni di curva (az.3), theta-nome q-series (az.4), turning point (az.2).

## Note di priorità
- Azioni 1,2,6,7,8: editoriali/scoping (testo), rapide-medie.
- Azione 5: piccola-media (aggiungere formula J_1 o dichiarare fixed-charge).
- Azione 3: la più sostanziale (nuova sottosezione convenzioni di curva).
- Azione 4: media (paragrafo PDF-alone vs repo + stampa q-series).
- Tutto da propagare a ENTRAMBI main.tex (iopart) e main_prd_revtex.tex (PRD).
