# Brachistocrone relativistiche in spaziotempi non stazionari — scheletro del paper

Indice dei risultati e delle figure. Dettagli nei md di cartella:
`FLRWmetric/FLRWresults.md`, `VaidyaMetric/VaidyaResults.md`,
`SdSMetric/SdSresults.md`, `ThakurtaMetric/ThakurtaResults.md`;
principio generale e liceità in `FLRWmetric/nonStationaryBrachi.md` §4b.

---

## §1 Introduzione
Perlick 1991: brachistocrone a rotaia in spaziotempi stazionari, riduzione
ottica singolare alle superfici di luce. Obiettivo: recuperarlo e
oltrepassarlo (spaziotempi dinamici, attraversamento delle superfici
critiche). Origine: programma Kerr (`KerrMetric/doranTau.md`, `doranT.md`).

## §2 Il principio: rotaia ideale come invariante controllato
- vincolo `−u·W = E` lineare in velocità + normalizzazione ⇒
  **indicatrice ellittica** a ogni evento (fig_indicatrici,
  `PaperFigures/`);
- liceità: f·u = 0 identico (costruttivo, forza tipo-magnetica
  `F^{μν}u_ν`); esistenza minimi = Filippov (compattezza+convessità);
  PMP con trasversalità H=0; lettura vakonomica (nota d'Alembert);
- selezione canonica di W in gerarchia: Killing (teorema, rotaia
  gratuita) → CKV (conforme) → **Kodama** (sferico dinamico,
  `VaidyaMetric/vaidya_kodama_sympy.py`) → convenzione dichiarata;
- costo di controllo `f·ξ` come misura della simmetria mancante
  (0 in Kerr/SdS; aH in FLRW; ∝m′ in Vaidya; ∝A′ in Thakurta);
- limiti-teorema: nullo = Kovner–Perlick (spaziotempi arbitrari),
  stazionario = Perlick 1991.

## §3 FLRW — caso base esattamente risolubile
R1–R9 in `FLRWresults.md`. Degenerazione (brachistocrona = geodetica
comovente), redshift della rotaia, congelamento a=Ê, forza f_η=aH,
metrica ottica con orologio λ(η), hamiltoniane PMP, forme chiuse de
Sitter, limite nullo = orizzonte di Hubble.
Figure (`FLRWfigures/`): cinematica, worldlines vs orizzonte
(err 2e-6), minimo variazionale diretto (entrambi i rami).

## §4 Vaidya — fenomenologia dinamica (W = Kodama)
§1–3f in `VaidyaResults.md`. EOM non autonome (J sopravvive, forzante
∝m′, memoria teleologica p_v), parametro v attraverso il periasse,
hamiltoniane Zermelo (vento radiale f−E²), metrica ottica via omotetia
(m=μv: x±, μ_c=1/16) e forme ottiche deboli per m(v) generico
(= Perlick con f→f(v,r)).
Risultati principali: **finestra temporale del periasse** (R7) e
**mappa di penetrazione** (R8: accrescimento ⇒ J_c finito ≈ (8÷10)√μ,
evaporazione ⇒ J_c=0); auto-sintonizzazione dinamica del ramo v (R9).
Figure (`Vaidyafigures/`): orizzonti EH/AH + autosimilari, residui Kerr
a=0 (1e-12..17), traiettorie+p_v, variazionale, bounce, timing,
penetration map.

## §5 Buchi neri cosmologici — SdS e Thakurta(-Kerr)
- **SdS** (`SdSresults.md`): Killing genuino solo nella vasca; f_max,
  r_ph=3M esatta, doppie barriere w=0; dicotomia massimale (τ in
  scatola, t attraversa entrambi gli orizzonti); svolte = radici esatte
  (1e-12).
- **Thakurta-Schw** (`ThakurtaResults.md` R1–R6): CKV, superficie di
  congelamento a²f=Ê² che unifica le barriere FLRW/statiche,
  **cattura da espansione** (fig_thakurta_cattura).
- **Thakurta-Kerr** (R7–R8a): indicatrice con vento; identità
  fP+(2Ms/r)²=Δ ⇒ R²=Δv̄²/P̄; orizzonte RIGIDO (conformemente
  invariante) vs congelamento che RESPIRA (fig_thakurta_kerr_superfici);
  hamiltoniane chiuse con shift gravitomagnetico conforme; validazione
  su doranTau/doranT a 10 cifre (fig_thakurta_kerr_residui);
  **plunge inversion t vs τ** (R10): rapporto delle Hamiltoniane
  ρ = A²[f+2Ms/(rJ)]/Ê, inversione controllata dal fattore conforme,
  asintoto A_inv = √Ê (fig_thakurta_kerr_plunge_t_tau).

## §6 Metriche ottiche non stazionarie — la gerarchia
Tabella (VaidyaResults §3e): orologio FLRW → indice n(η,r) su scheletro
statico (Thakurta) → vento rigido + indice (Thakurta-Kerr:
`dη = n·α_K + β_K` con α,β = Randers NULLO di Kerr, R8) → quoziente
autosimilare (Vaidya μv) → senso debole (Vaidya generico).
**R8a**: dalla stessa ellisse, tre costi: rami η,t degeneri (orologio
t(η)=∫A dη); ramo τ **Riemanniano puro** (`dτ = k_τ·α_K`, il
gravitomagnetico cancella esattamente): il vento è esclusivo dei tempi
d'arrivo, con/senza espansione.
Figura: fig_thakurta_kerr_rami (τ speculare 1e-11 vs η spirali
trascinate prograde; scalari conformi n, k_τ).

## §6b Quasi-costanti in spaziotempi conformi (ThakurtaResults R11–R11e)
- Kerr: espansioni K_τ, K_t a O(a⁶) (K_tau_expansion.txt): trasferimento
  a Thakurta-Kerr ESATTO per A costante (teorema di riscalamento,
  R11c) per sostituzione E → E_eff = Ê/A; verdetto numerico però
  ribaltato a E_eff<1 (NLO perde efficacia, R11b); il muro di
  congelamento entra come polo dei coefficienti (DE₀ = r²w_eff).
- A(η) dinamico: rottura all'ordine misto a²ε (ε = A′/A); Hamiltoniana
  PMP 3D dinamica (R11d, identità f_ΣG+b²=Δsin²θ, L² a 6e-13);
  sorgente chiusa S₁ = 2E_eff²[cos²θ + Mcos2θ·p_θ²/(rDE₀²)]; medie
  angolari chiuse (⟨cos2θp_θ²⟩=|J|(L−|J|), 1e-16); controtermine
  K^{TK} = K^{Kerr}(E_eff(η)) − a²∫εS₁dη (R11e).
- Risultato strutturale: per E_eff<1 il muro è ATTRATTORE GLOBALE del
  ramo τ anche ad A costante (niente orbite legate, p_r→∞ in η finito):
  la quasi-costante dinamica vive su archi di scattering, dove il
  termine conforme è correzione fine al residuo Kerr.
- Figura: fig_thakurta_kerr_quasicostante (cattura al muro).

## §6c Struttura algebrica (ThakurtaResults R12)
Geodetiche eq.: genere 1; brachistocrone τ: genere 2 (rotaia √(wf));
t: genere 3; separatrice J_c: ELLITTICA (fattorizzazione doranTau);
soglie: elementari. Conforme: isola ellittica mobile J_c(A)=sA²/Ê.
φ(r) della separatrice ESPLICITA (R12a): lineare in z + due rapporti
di σ di Weierstrass, poli agli orizzonti, g₂ g₃ chiusi; validata a
5.6e-17 / 30 cifre; shift di Doran chiuso (R12c: seconda curva
g₂=−M²a², g₃=0, pesi ±a/(r₊−r₋)): φ_D regolare attraverso ergosfera
e orizzonte (appendice naturale del paper); verifica end-to-end
con σ,ζ,℘ reali (theta, reticolo rombico): 2.8e-10 attraverso
l'ergosfera (fig_separatrix_3traiettorie); la PMP-η del ramo τ
termina a r_e (H(p)>0 dentro): hamiltonianizzazione della
singolarità intrinseca della riduzione. Generalità (R12e,
fig_separatrix_gallery): la sovrapposizione a 3 metodi tiene su
4 parametri (a,E), ciascuno col suo reticolo (|W−quad| ~1e-10,
reticolo sempre rombico per E>1, graze a r_e=2M in ogni caso) —
non è un caso speciale.

## §7 Conclusioni e sviluppi
Fuga dall'evaporazione (Vaidya m′<0 dall'interno), dinamica sulla
superficie di congelamento, settore autosimilare di Vaidya lineare
con rotaia omotetica, Vaidya–de Sitter, scala J_c(μ) analitica;
teoria del primo ordine per lo shift non adiabatico dell'inversione
(R9: +73% mediano); zoom sulla separatrice di cattura all'ergosfera
J_c(A)=sA²/Ê; ramo t dinamico per le quasi-costanti (attraversa invece
di congelarsi: orbite lunghe dove il controtermine conforme domina);
NLO ri-ottimizzato per il regime E_eff<1 (base con polo al muro);
porting del numtest Kerr originale sotto E_eff=1.

## Inventario figure
(tutte rigenerate single-column, label inglesi, PDF vettoriale +PNG,
via `paper_style.py`; eccezione full-width: fig_separatrix_gallery)
- PaperFigures/: fig_indicatrici, fig_brachistocrone_confronto
- FLRWfigures/: fig_flrw_cinematica, _worldlines, _variazionale
- Vaidyafigures/: fig_vaidya_orizzonti, _kerr_a0, _traiettorie,
  _variazionale, _bounce, _timing, _penetration_map,
  fig_kodama_conservazione (R12), fig_verifica_minimo_brachi (R14),
  fig_vaidya_plunge_t_tau (R15),
  fig_vaidya_no_inversione_evaporazione (R16, risultato negativo:
  l'evaporazione fisica non inverte; superata la vecchia figura ora
  rinominata DEPRECATED_fig_vaidya_inversione_evaporazione, ERRATA)
- Thakurtafigures/: fig_thakurta_cattura, _kerr_superfici,
  _kerr_residui, _kerr_rami, _kerr_plunge_map,
  fig_thakurta_kerr_inversione_AJ (R10b, sostituisce _kerr_inversione_t_tau
  ora DEPRECATED: J-range troppo stretto nascondeva l'inversione),
  _kerr_plunge_t_tau, _kerr_quasicostante,
  _kerr_drift_map, _kerr_drift_map_A,
  fig_thakurta_tricotomia_Jneg (R12f),
  fig_thakurta_cuspide_ergosfera (R12g)
- (opzionale, da fare) SdS: figura "scatola a doppia barriera"
