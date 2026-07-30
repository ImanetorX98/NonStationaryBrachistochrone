# Split map — master `paper/main.tex` → Paper I / Paper II (I.0 change-log)

Sorgenti: `paper1/paper1.tex` (rail/Kodama/Vaidya sferico), `paper2/paper2.tex` (Thakurta-Kerr conforme).
Risorse condivise (via percorso relativo, nessuna duplicazione): `../paper/Immagini/` (figure),
`../paper/refs.bib` (bibliografia). Classe: `iopart` (CQG) per entrambi.

## Mappa sezioni del master (righe indicative)
| Master (main.tex) | Paper I | Paper II | Note |
|---|---|---|---|
| §1 Introduction (84) | ✓ (adattata, sferico) | ✓ (adattata, cita Paper I) | due intro distinte |
| §2 Controlled-rail principle (171) | **✓ (fondamenta)** | cita Paper I | §2.1 indicatrix+problema controllo; §2.2 esist./norm./HJB; §2.3 W-hierarchy |
| §3 FLRW (520) | ✓ | — | caso degenere |
| §4 Vaidya (576) | ✓ (tutte le sottosezioni) | — | Kodama energy, penetrazione/timing/bounce, plunge law, adiabatica+dilog, sorgente `[r p_r]−λ` |
| §5 Thakurta-Kerr (1040) | — | ✓ | indicatrix/H_eta/H_t/H_tau, η-brachistocrona t≡η, Randers, adiabatica |
| §5.1 Quasi-constants & drift (1190) | — | **✗ RIMOSSA** | off-equatoriale/Carter/O(a²) → materiale paper futuro 3D |
| §5.2 Semi-analytic first-order (1241) | — | ✓ | |
| §6 Equatorial closed forms & ergosphere (1689) | — | ✓ | separatrix Weierstrass, tricotomia/cuspide, inversione, fixed-endpoint tiered, breathing |
| §7 Conclusions (2141) | ✓ (Vaidya scope) | ✓ (TK scope) | 3 blocchi ciascuna |
| App A Validation (2228) | ✓ (Vaidya/shared) | ✓ (TK) | split per geometria |
| App B genus-degeneration (2812) | ✓ (Vaidya J_deg) | ✓ (TK sep + J_deg) | |
| App C first-order (3004) | ✓ (Vaidya source) | ✓ (TK canonical + conformal source) | |

## Contenuti condivisi (in Paper I, citati da Paper II)
Formalismo rail `−u·W=Ê`; enunciato problema controllo ottimo; esistenza/normalità/HJB verification;
gerarchia W Killing→CKV→Kodama; equivalenza Perlick (limite stazionario). Paper II li richiama.

## Regole di rietichettatura durante la migrazione
- Ogni figura/eq/sez ri-numerata: rigenerare i cross-reference DOPO la migrazione.
- Terminologia già armonizzata nel master (main15): mantenere J_deg/J_sep/J_pen, negative-rate,
  axial control costate, weight-1/weight-2, HJB non-autonoma, App C canonica.
- Data availability + DOI Zenodo identici in entrambi.

## Stato I.0
- [x] Creati `paper1/paper1.tex` e `paper2/paper2.tex` (scheletri con preambolo, graphicspath+bib
      condivisi, placeholder di sezione con \todo di migrazione).
- [x] Questa mappa (SPLIT_MAP.md).
- [ ] Prossimo: I.1/I.2 migrazione contenuti sezione per sezione.
