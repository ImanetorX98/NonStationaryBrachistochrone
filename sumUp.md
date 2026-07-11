# Riassunto di sessione — brachistocrone in spaziotempi non stazionari

Sintesi dei risultati prodotti in questa sessione, con puntatori ai file
di dettaglio, agli script e alle figure, per la stesura dell'articolo.

File di riferimento:
- `FLRWmetric/FLRWresults.md` — caso base FLRW
- `VaidyaMetric/VaidyaResults.md` — Vaidya (buco nero dinamico)
- `ThakurtaMetric/ThakurtaResults.md` — Thakurta / Thakurta-Kerr conforme
- `SdSMetric/SdSresults.md` — Schwarzschild–de Sitter
- `paperOutline.md` — struttura dell'articolo + inventario figure
- `paper_style.py` — stile figure single-column (radice del progetto)

---

## 0. Infrastruttura figure (fatto in questa sessione)

Tutte le ~22 figure del paper rigenerate in **inglese** e **layout a
colonna singola** (larghezza 3.4″, font 7–8pt, PDF vettoriale + PNG) via
il modulo condiviso `paper_style.py` (`COL`, `DCOL`, `set_style`,
`savefig`). Multi-pannello riorganizzati in verticale per stare in
colonna; unica eccezione full-width: `fig_separatrix_gallery`.
Inventario completo in `paperOutline.md` §"Inventario figure".

---

## 1. FLRW — degenerazione dei rami (già in FLRWresults.md)

**Risultato concettuale ripreso**: in FLRW i tre tempi (η conforme, t
cosmico, τ proprio) danno la **STESSA curva** (retta comovente). Causa:
omogeneità spaziale — `a(η)` non dipende dalla posizione, quindi i tre
funzionali di tempo sono pesi di η soltanto e hanno lo stesso minimo.
La scissione t/τ è un effetto di **gradiente gravitazionale**; FLRW è il
caso degenere di base. → Ottima diagnostica per §3 dell'articolo.

---

## 2. Vaidya — fondazione, leggi di plunge, inversione
(dettagli in `VaidyaMetric/VaidyaResults.md`, §3f e R12–R16)

### R12 — Energia di Kodama conservata lungo la brachistocrona
Script `kodama_conservation.py`, fig `fig_kodama_conservazione`.
- Miracolo di Kodama (simbolico): `K=∂_v`, `∇·K=0`, carica = massa di
  Misner–Sharp `=m(v)`.
- Identità della rotaia: `−u·K = E` cade ESATTA dal vincolo (non
  imposta). Verifica numerica `|−u·K−E|<6.7e-16`.
- Non-banalità: per geodetica `−u·K` driftrebbe `∝m′` (controllo, non
  Noether). **Test discriminante**: solo il Kodama NON normalizzato ha
  costo nullo nel limite statico ⟹ `W=∂_v` fissato da "statico ⟺
  costo-zero".

### R13 — Segno fisso di r_min^t − r_min^τ
Ai punti di svolta: τ: `(E²−f)J²=f r²`, t: `(E²−f)J²=E²r²/f`.
Rapporto locale `E²/f² > 1` ⟹ **r_min^t < r_min^τ** (t più fondo).
Robusto anche in dinamica (R16): l'evaporazione fisica non inverte.

### R14 — Validazione diretta del principio di minimo (tre test)
Script `verifica_minimo_brachi.py`, fig `fig_verifica_minimo_brachi`.
Perturbazione a estremi fissi: minimo esatto sulla soluzione. Tre test
indipendenti: FLRW (`fig_flrw_variazionale`), Vaidya non autonoma
(`fig_vaidya_variazionale`), t-vs-τ statico. Ogni ramo minimizza il SUO
tempo e non l'altro ⟹ curve genuinamente distinte + minimi veri.

### R15 — Legge del raggio minimo di plunge t vs τ (dinamico)
Script `plunge_vaidya_t_tau.py`, fig `fig_vaidya_plunge_t_tau`.
Ramo t = minimizza tempo avanzato v (≡ min t a estremi radiali fissi).
Legge di svolta = forma statica R13 con `f_*=1−2m(v_peri)/r_min`.
- Coerenza m′=0: riproduce le radici statiche R13 (residuo ~1e-9).
- Segno robusto per μ∈[−0.02,0.02]; accrescimento allarga il gap.
- **Correzione dinamica ∝m′** dal momento teleologico `p_v` (memoria):
  la legge NON è la statica a m(v_peri), il residuo cresce ~lineare in
  m′. La svolta ricorda l'intera storia.

### R16 — NON c'è inversione evaporativa fisica [RETTIFICATO]
Script `inversione_fisica.py`, fig `fig_vaidya_no_inversione_evaporazione`.
Rimpiazza la versione errata (`inversione_evaporazione.py`,
`fig_vaidya_inversione_evaporazione` — DA NON USARE).
- Una prima analisi (lineare `m=1+μv`, v1=40) sembrava mostrare
  inversione a `μ_inv=−0.058`: **artefatto della massa negativa**. A
  quel μ la svolta era a m≈0.67 ma l'ancoraggio a v1=40 era a m=−1.31
  (singolarità nuda). La serie perturbativa ha `R~|μ_inv|` perché la
  singolarità È l'ingresso in m≤0.
- **Controlli**: scan 2D lineare (ogni gap<0 ha m<0; ogni m>0 ha gap>0)
  + esponenziale `m=e^{−λv}` (m>0 sempre) → gap plateau positivo, **mai
  inverte** (fino a λ=1, m(40)=4e-18).
- **Conclusione**: l'evaporazione fisica restringe il gap ma non lo
  inverte; `r_min^t < r_min^τ` robusto. **Un solo meccanismo di
  inversione, ROTAZIONALE** (Kerr, R12c). "Serve rotazione" (R13) è
  CORRETTO.

---

## 3. Thakurta-Kerr — separatrice chiusa e tricotomia
(dettagli in `ThakurtaMetric/ThakurtaResults.md`, R10–R12f)

### R10b — Teorema dell'inversione conformale [NUOVO]
Script `inversione_conformale_AJ.py`, fig `fig_thakurta_kerr_inversione_AJ`
(sostituisce la deprecata `fig_thakurta_kerr_inversione_t_tau`, J-range
troppo stretto). Il fattore conforme A INVERTE il plunge (τ più fondo).
Condizione chiusa `cond(r,A)=A²bQ−P̄(Ê−A²f)=0`, Q=b·v̄²+√(v̄²Δ),
J_inv=P̄/Q. Teorema: `cond(1)<0`, `cond(A_freeze)=ÊrΔ(Ê−1)/(r−2M)>0`
⟹ per IVT `A_inv∈(1,A_freeze)`: inversione **sempre esistente e
raggiungibile** (far-field A_inv→√Ê). Vive a **J grande (~16–24)**:
il vecchio scan `J∈[3.9,5.5]` la mancava (cadeva TRA l'inversione
rotazionale a J~1.6 e quella conformale a J~24). Verificato: integrazione
diretta r_τ=r_t a ~1e-4; mappa (J,A) 1029/3271 celle, cond≡numerico.

### R12e — Galleria: la sovrapposizione a 3 metodi è generale
Script `../KerrMetric/kerr_separatrix_gallery.py`, fig
`fig_separatrix_gallery` (full-width). La forma chiusa Weierstrass della
separatrice (R12a) verificata su 4 parametri (a,E), reticoli diversi:
|W−quad| ~1e-10 in tutti; reticolo sempre rombico (E>1); grazing a
r_e=2M strutturale. Non è un caso speciale.

### R12f — Tricotomia retrograda: J negativo di cattura [NUOVO]
Script `tricotomia_equatoriale.py`, `tricotomia_figura.py`, fig
`fig_thakurta_tricotomia_Jneg`.
Banda di cattura (raggiungere l'ergosfera r_e=2M) al variare di J,
incluso J<0:
- **Ramo τ**: banda SIMMETRICA `[−sA²/Ê, +sA²/Ê]` (A=1: ±0.75). La
  simmetria viene dalla struttura `J²` di doranTau.
- **Ramo t (η)**: banda larga e ASIMMETRICA, `J_neg^t≈−8.05`
  (retrogrado forte); il ramo t è molto più capture-prone (coerente con
  R13/R15).

### PRECISAZIONE importante sulla "penetrazione" (da chiarire nel paper)
Test (barriera a r_e): `min_{p_r}H = 0` ESATTO a r_e per ogni J≤J_c, con
barriera `>0` appena sotto. Integrando senza fermare a r_e, TUTTA la
banda svolta esattamente a r_e (r_min=2.0). Conclusioni:
- **Nessun J penetra davvero in senso hamiltoniano/PMP** — nemmeno J_c:
  tutte svoltano a r_e per la barriera.
- La spirale su r_+ è la **continuazione analitica** (doranTau, genus-1),
  definita a `J=±J_c` (fattorizzazione), come worldline vincolato non
  come brachistocrona PMP.
- Sotto r_e esistono regioni di nuovo permesse ⟹ **rami interni
  disconnessi**, irraggiungibili da fuori.
Quindi la "tricotomia" va enunciata con cura: scatter (|J|>J_c) /
svolta-a-r_e (|J|≤J_c) / separatrice-analitica (J=±J_c).

---

## 4. Struttura per l'articolo (riferimento paperOutline.md)

Filo conduttore suggerito dai risultati di sessione:
1. Rotaia lecita come invariante controllato; W dalla gerarchia
   Killing → CKV → **Kodama** (Vaidya, R12) → convenzione.
2. FLRW degenere (omogeneità) → il buco nero rompe la degenerazione t/τ.
3. Leggi di plunge: statica chiusa (R13), dinamica con memoria (R15).
4. **Un solo meccanismo di inversione**: quello ROTAZIONALE (Kerr,
   chiuso, R12c). L'evaporazione NON inverte in regime fisico (R16):
   punto forte come risultato negativo (falso allarme dell'artefatto
   a massa negativa, poi smentito da lineare-2D + esponenziale).
5. Separatrice in forma chiusa (Weierstrass, R12a–e) e tricotomia
   equatoriale con la precisazione sulla barriera a r_e (R12f).
6. Validazione del principio di minimo su tre spaziotempi (R14).

**Attenzione a non confondere due fenomeni distinti** (bordi diversi,
spaziotempi diversi, carte diverse):
- **Penetrazione dell'ORIZZONTE — Vaidya (non rotante), DICOTOMIA**:
  fig_vaidya_penetration_map. Bordo = orizzonte apparente `r=2m(v)`;
  niente ergosfera (a=0); soglia `J_c(v₀)` cattura/riflette. Coordinate
  EF avanzate `(v,r)` GIÀ regolari all'orizzonte → **nessun Doran**.
- **Tricotomia dell'ERGOSFERA — Thakurta-Kerr (rotante)**:
  fig_thakurta_tricotomia_Jneg (R12f). Bordo = ergosfera `r_e=2M`;
  scatter / grazing-a-r_e / separatrice-analitica; struttura da
  trascinamento (2Mas/r) + barriera a r_e. Le **coordinate di Doran**
  (R12c) servono QUI, per la forma chiusa `φ(r)` attraverso l'ergosfera.

Aperti / da fare: figura SdS "scatola a doppia barriera"; enunciato
finale pulito della tricotomia con la distinzione grazing/continuazione.

---

## 5. Tabella figura → risultato → script → file .md

Cartelle: FLRW = `FLRWmetric/FLRWfigures/`, VA = `VaidyaMetric/Vaidyafigures/`,
TH = `ThakurtaMetric/Thakurtafigures/`, PF = `PaperFigures/`,
KM = `KerrMetric/`. Tutte in `.pdf` (paper) + `.png` (anteprima).

| figura | risultato | script | .md di dettaglio |
|---|---|---|---|
| FLRW/fig_flrw_cinematica | cinematica rotaia, congelamento a=Ê | genera_figure_flrw.py | FLRWresults.md |
| FLRW/fig_flrw_worldlines | worldline comoventi, Δx_max (R5) | genera_figure_flrw.py | FLRWresults.md |
| FLRW/fig_flrw_variazionale | min. variazionale t,τ (R2/R14) | genera_figure_flrw.py | FLRWresults.md |
| VA/fig_vaidya_orizzonti | EH vs AH; raggi autosimilari | genera_figure_vaidya.py | VaidyaResults.md §1,R6 |
| VA/fig_vaidya_traiettorie | orbite τ Schw vs Vaidya + p_v | genera_figure_vaidya.py | VaidyaResults.md R4 |
| VA/fig_vaidya_penetration_map | soglia J_c(v₀), asimmetria (R8) | vaidya_penetration_map.py | VaidyaResults.md §3c |
| VA/fig_vaidya_timing | r_min vs istante d'arrivo (R7) | vaidya_brachistochrone_vparam.py | VaidyaResults.md §3b |
| VA/fig_vaidya_bounce | rimbalzo al periasse (parametro v) | vaidya_brachistochrone_vparam.py | VaidyaResults.md §3b |
| VA/fig_vaidya_variazionale | min. EL non autonomo (R14) | genera_figure_vaidya.py | VaidyaResults.md R14 |
| VA/fig_vaidya_kerr_a0 | residui limite Kerr a=0 (R3) | genera_figure_vaidya.py | VaidyaResults.md R3 |
| VA/fig_kodama_conservazione | −u·K=E lungo brachi (R12) | kodama_conservation.py | VaidyaResults.md R12 |
| VA/fig_verifica_minimo_brachi | t/τ minimi veri per perturbaz. (R14) | verifica_minimo_brachi.py | VaidyaResults.md R14 |
| VA/fig_vaidya_plunge_t_tau | legge plunge t vs τ, gap vs μ (R15) | plunge_vaidya_t_tau.py | VaidyaResults.md R15 |
| VA/fig_vaidya_no_inversione_evaporazione | NO inversione fisica (R16) | inversione_fisica.py | VaidyaResults.md R16 |
| TH/fig_thakurta_cattura | cattura da espansione, freezing | genera_figure_thakurta.py | ThakurtaResults.md |
| TH/fig_thakurta_kerr_superfici | orizzonte rigido vs freezing | genera_figure_thakurta.py | ThakurtaResults.md |
| TH/fig_thakurta_kerr_residui | flusso H vs forme chiuse Kerr | genera_figure_thakurta.py | ThakurtaResults.md |
| TH/fig_thakurta_kerr_rami | tre rami dall'indicatrice ellittica | genera_figure_thakurta_rami.py | ThakurtaResults.md R8a |
| TH/fig_thakurta_kerr_plunge_map | mappa cattura/fuga (J,A) | thakurta_kerr_plunge_inversion.py | ThakurtaResults.md |
| TH/fig_thakurta_kerr_inversione_AJ | inversione conformale (J,A), teorema cond=0 (R10b) | inversione_conformale_AJ.py | ThakurtaResults.md R10b |
| TH/fig_thakurta_kerr_plunge_t_tau | traiettorie t,τ profondità | thakurta_kerr_plunge_t_tau.py | ThakurtaResults.md |
| TH/fig_thakurta_kerr_quasicostante | K_τ NLO e drift | thakurta_kerr_quasicostanti.py | ThakurtaResults.md R11 |
| TH/fig_thakurta_kerr_drift_map | drift (θ,r) log, ergosfera | genera_drift_map.py | ThakurtaResults.md R11e |
| TH/fig_thakurta_kerr_drift_map_A | drift (A,r), muro r_w(A) | genera_drift_map_A.py | ThakurtaResults.md R11e |
| TH/fig_thakurta_tricotomia_Jneg | J<0 di cattura, banda t/τ (R12f) | tricotomia_figura.py | ThakurtaResults.md R12f |
| TH/fig_thakurta_cuspide_ergosfera | cuspide a r_e, potenza 2/3, tricotomia (R12g) | cuspide_ergosfera.py | ThakurtaResults.md R12g |
| KM/fig_separatrix_3traiettorie | 3 metodi separatrice (R12d) | kerr_separatrix_trajectories.py | ThakurtaResults.md R12d |
| KM/fig_separatrix_gallery | galleria 4 parametri (R12e, full-width) | kerr_separatrix_gallery.py | ThakurtaResults.md R12e |
| KM/fig_bvp_estremi_fissi | BVP estremi fissi sim/asim, τ più fondo (R10c) | bvp_estremi_fissi.py | ThakurtaResults.md R10c |
| KM/fig_bvp_kerr_inversione | Δr(a) estremi fissi, no inversione spin (R10c) | bvp_kerr_simmetrico.py | ThakurtaResults.md R10c |
| KM/fig_bvp_conforme_inversione | Δr(E_eff) estremi fissi, no inversione A (R10c) | bvp_conforme_inversione.py | ThakurtaResults.md R10c |
| KM/fig_colormap_spin_estremi_fissi | colormap Δr(J,a) estremi fissi, Δr>0 (R10c) | colormap_spin_estremi_fissi.py | ThakurtaResults.md R10c |
| KM/fig_colormap_conforme_estremi_fissi | colormap Δr(J,A) estremi fissi, Δr>0 (R10c) | colormap_conforme_estremi_fissi.py | ThakurtaResults.md R10c |
| KM/fig_colormap_spin_asimmetrico | colormap Δr(J,a) estremi ASIMMETRICI, Δr>0 (R10c) | colormap_spin_asimmetrico.py | ThakurtaResults.md R10c |
| KM/fig_colormap_conforme_asimmetrico | colormap Δr(J,A) estremi ASIMMETRICI, Δr>0 (R10c) | colormap_conforme_asimmetrico.py | ThakurtaResults.md R10c |
| TH/fig_colormaps_inversione | colormap stesso-lancio (J,a) e (J,A) invertono | colormaps_inversione_aJ_AJ.py | ThakurtaResults.md R10c |
| PF/fig_indicatrici | indicatrici ellittiche (concettuale) | genera_figure_thakurta.py | ThakurtaResults.md |

Script di sola analisi (senza figura dedicata):
`tricotomia_equatoriale.py` (R12f, tabella soglie J),
`kerr_separatrix_weierstrass.py` / `kerr_separatrix_validation.py` /
`kerr_doran_shift_weierstrass.py` (forma chiusa separatrice, R12a–c).
