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

## §10 — Azioni raccomandate (8) — ✅ TUTTE INDIRIZZATE
1. [x] Contraddizione abstract/status vs Table 1 RISOLTA (opzione 1): off-shell assemblato
   come depth-two iterated Abelian integrals + elementare su 4 branch; solo naming deferito.
   Abstract + intro + riga Table 1 allineati (main+PRD). Commit bc2b85c.
2. [x] Validità uniforme solo su sub-archi regolari `H_Pr≠0` GIÀ presente (main 1405-1406 /
   PRD 1213: "reduction holds on compact sub-arcs bounded away from the turning point").
3. [x] Convenzioni di curva: paragrafo "Curve conventions (fixed once for all letters)" in
   App. (main+PRD): sheet `√S>0` arco esterno, cuts tra radici reali, base point `r_0` +
   orientamento crescente-r, costanti additive pinnate a riferimento regolare (cancellano in δφ).
   Commit b225b9e.
4. [x] PDF-solo vs repository: paragrafo "What is reconstructible from the text alone" (main) /
   "Reconstructibility" (PRD): testo fissa tutto l'algebrico + la FORMA della lettera genus-2 +
   q-series genus-1; i NUMERI genus-2 (τ, Siegel, nome series) via Sage/abelfunctions. Commit b225b9e.
5. [x] Fixed-charge/fixed-endpoint: paragrafo "Fixed-charge versus fixed-endpoint" + eq
   `J_1=−δφ_dyn/∂_J φ_0` (eq:J1-fixedendpoint) in main; versione condensata PRD. Commit 7c57bc8.
6. [x] Istantanea vs crossing GIÀ separate (main 1971-1989): "instantaneous, local first-order
   statement about the neighbourhood of the separatrix" vs "long-time separatrix-following";
   periodo radiale diverge, ipotesi adiabatica degrada, Neishtadt crossing "beyond first order".
7. [x] η = rappresentazione canonica del ramo di ARRIVO (non terza famiglia): `dt=A(η)dη`
   monotono ⟹ `argmin t_f=argmin η_f`; τ distinto; Vaidya nativo = v. main+PRD. Commit 7c57bc8.
8. [x] Vaidya: clock canonici = v (advanced null, Kodama `K=∂_v`) + τ; il t del confronto di
   profondità è il t frozen (Schwarzschild-like), NON canonico sotto `m=m(v)`; u (retarded) =
   problema outgoing-Vaidya SEPARATO. main+PRD (+ 700-712 già presente). Commit b225b9e.

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
