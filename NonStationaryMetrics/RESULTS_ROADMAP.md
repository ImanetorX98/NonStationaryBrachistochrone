# Roadmap risultati analitici — da curare tutti

Priorità e stato dei risultati analitici chiusi da sviluppare.

## In corso
- [ ] **Teorema adiabatico unificato** (FLRW / Vaidya / Thakurta-Kerr) —
  stessa struttura WKB, parametri lenti diversi (scala a, massa m(v), conforme A),
  coefficienti algebrici, trascendenza nei polilog genus-≤2.
  **+ collegamento alle ellissi iniziali (indicatrici di Randers):** il filo "shape"
  FLRW=cerchio (no wind) → Vaidya/Kerr = ellissi (con wind) che "respirano" col
  parametro lento.

## Da fare (tutti da curare)
- [x] **Separatrice → dilog ellittico (Bloch-Wigner).** FATTO: Jc=7.0266, r_d=-3.3637,
  E:w²=Q4 {-2.083,0,2,8.727}, k²=0.6067, τ=0.90597 i (Legendre=Sage, 15 cifre),
  orizzonte=2-torsione. **PESO-1 forma chiusa σ,ζ verificata 1e-14**:
  U_0=ρ[lnσ(z−z_d)−lnσ(z+z_d)]+Cz, ρ=1/√Q4(r_d) e C algebrici (derivati, base al
  punto di Weierstrass ⇒ involuzione z→−z). PESO-2: circolo D_0=U_0 ln(r−2m)−G_0
  verificato 1e-14; G_0=dilog ellittico ∫lnσ ζ; q-serie Zagier D^E implementata.
  UNICO aperto (teoria pura): normalizzazione esatta regulator che lega D_0 fisico a
  D^E. Vedi vaidya_ell_dilog_match.py + vaidya_separatrix_bloch_wigner.md.
- [x] **Asimmetria accrescimento/evaporazione.** FATTO (eq:vaidya-asymmetry, verif
  1e-15): accr=ingoing/avanzato v=t+r∗, evap=outgoing/ritardato u=t−r∗; ∂_m F identico,
  v↔u flippa r∗. δφ/ṁ|accr=A∞+B, δφ/ṁ|evap=A∞−B; asimmetria netta 2B portata SOLO dal
  tortoise (A∞ polilog-∞ si cancella); parte trascendente = dilog D_k d'orizzonte.
  Legata alla finestra di penetrazione (faccia fisica vs fenomenologica dell'orizzonte).
  A∞=8.728, B=6.113, accr=14.840, evap=2.615. Vedi vaidya_asymmetry.md.

## Fatti (riferimento)
- [x] Thakurta-Kerr: φ(r,A) chiuso, δφ = ½Ê Σ Q_kj W_kj + elementari (verif 4-5e-14).
- [x] ψ = polilog genus-2 (NON riducibile a peso-1: divisore theta Θ=W_{g-1}).
- [x] Quasi-periodi η (2ª specie), κ simmetrica 1e-12; L=Fay (3ª specie, round-trip 5e-8).
- [x] Vaidya adiabatico: frozen Schwarzschild genus-2 (a=0), ∂_M F ridotto 1e-15,
  clock v(r) tempo-avanzato, δφ_V esplicito 9.4e-14, dilog d'ORIZZONTE D_k.
- [x] δφ_V FORMA COMPLETA assemblata (eq:vaidya-full, verif 2e-13): tutti i building
  block espliciti U_k,W_3k,D_k,L_2m,I_poly con coeff algebrici c_k^M,a_j^M.
- [x] D_k CHIUSO (eq:horizon-dilog): D_k = U_k ln(r-2m) - G_k, G_k=∫ln(r-2m)r^k/√S
  dilog iperellittico ancorato all'ORIZZONTE (lettera dlog residuo 2 al branch point
  r=2m), compagno di W_kj, irriducibile (divisore theta). Verif IBP 1e-12.
  Separatrice → collasso ellittico Bloch-Wigner (struttura verificata).
- [x] Soglia penetrazione ergosfera: A_c^wall=Ê/√(1-2M/r0), J_c^+(A) (1e-15).
