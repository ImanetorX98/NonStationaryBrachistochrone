# Inventario della bibliografia locale e audit delle fonti

**Data dell'audit:** 25 agosto 2026  
**Bibliografia controllata:** `paper/refs.bib`  
**Cartella delle fonti:** `Fonti/`

## Esito sintetico

La raccolta locale è stata ampliata con **15 PDF** provenienti esclusivamente da
Zenodo, arXiv, repository d'autore/universitari, dall'editore open access PeerJ e
dall'archivio ufficiale JETP. Non sono stati usati Sci-Hub, mirror pirata o copie
di provenienza giuridicamente incerta.

| Voce di controllo | Esito |
|---|---:|
| Voci in `paper/refs.bib` | 98 |
| Riferimenti software/web, per i quali non è atteso un PDF | 5 |
| Pubblicazioni per le quali è sensato cercare un testo completo | 93 |
| PDF presenti in `Fonti/` dopo l'intervento | 83 |
| PDF validi e leggibili da `pdfinfo` | 83/83 |
| Voci coperte dall'opera citata o da una sua versione d'autore/preprint | 72/93 |
| Voci con solo materiale sostitutivo, non identico alla fonte citata | 2/93 |
| Voci senza testo completo locale, neppure sostitutivo | 19/93 |
| Pubblicazioni per cui manca ancora la versione esatta citata | 21/93 |

La copertura **stretta** è quindi 72/93 (77,4%). La copertura utile per lo studio,
includendo i due sostituti esplicitamente marcati, è 74/93 (79,6%). I nove PDF
rimanenti sono letture di contesto o duplicati/preprint aggiuntivi non associati
univocamente a una voce BibTeX.

Tutti i PDF sono ignorati da Git e restano locali. Non sono stati modificati né
`paper/refs.bib` né i due paper.

## Nuovi documenti acquisiti e verificati

Ogni file seguente inizia con la firma `%PDF-`, viene aperto correttamente da
`pdfinfo` e ha un numero di pagine plausibile. I PDF arXiv sono versioni d'autore
della stessa opera citata e non vanno confusi con l'impaginato editoriale finale.

| Chiave BibTeX | File locale | Pagine | Provenienza autorizzata |
|---|---|---:|---|
| `RosignoliKerr` | [Rosignoli-2026-Kerr-coordinate-proper-time-brachistochrones.pdf](Fonti/Rosignoli-2026-Kerr-coordinate-proper-time-brachistochrones.pdf) | 15 | [Zenodo, record 21139490](https://zenodo.org/records/21139490), CC BY 4.0 |
| `RosignoliIntSchw` | [Rosignoli-2026-interior-Schwarzschild-NLO-NNLO.pdf](Fonti/Rosignoli-2026-interior-Schwarzschild-NLO-NNLO.pdf) | 8 | [Zenodo, record 20678106](https://zenodo.org/records/20678106), CC BY 4.0 |
| `RosignoliTOV` | [Rosignoli-2026-TOV-coordinate-proper-time-brachistochrones.pdf](Fonti/Rosignoli-2026-TOV-coordinate-proper-time-brachistochrones.pdf) | 11 | [Zenodo, record 20708060](https://zenodo.org/records/20708060), CC BY 4.0 |
| `CaponioJavaloyes2026` | [CaponioJavaloyesSanchez-2026-Fermat-principle-formulations.pdf](Fonti/CaponioJavaloyesSanchez-2026-Fermat-principle-formulations.pdf) | 22 | [arXiv:2605.01532](https://arxiv.org/abs/2605.01532) |
| `SchmidtSmith` | [FrittelliNewman-1999-exact-universal-gravitational-lensing-equation.pdf](Fonti/FrittelliNewman-1999-exact-universal-gravitational-lensing-equation.pdf) | 13 | [arXiv:gr-qc/9810017](https://arxiv.org/abs/gr-qc/9810017) |
| `Liberzon` | [Liberzon-2012-calculus-variations-optimal-control.pdf](Fonti/Liberzon-2012-calculus-variations-optimal-control.pdf) | 200 | [pagina d'autore, University of Illinois](https://liberzon.csl.illinois.edu/teaching/cvoc.pdf); **precursore, non edizione finale** |
| `NielsenYeom2009` | [NielsenYoon-2008-dynamical-surface-gravity.pdf](Fonti/NielsenYoon-2008-dynamical-surface-gravity.pdf) | 17 | [arXiv:0711.1445](https://arxiv.org/abs/0711.1445) |
| `Kottler1918` | [Kottler-1918-physikalischen-Grundlagen-Gravitationstheorie.pdf](Fonti/Kottler-1918-physikalischen-Grundlagen-Gravitationstheorie.pdf) | 62 | [Zenodo, record 1424336](https://zenodo.org/records/1424336), articolo storico open |
| `ClaudelViracVirbhadra2001` | [ClaudelVirbhadraEllis-2001-geometry-photon-surfaces.pdf](Fonti/ClaudelVirbhadraEllis-2001-geometry-photon-surfaces.pdf) | 24 | [arXiv:gr-qc/0005050](https://arxiv.org/abs/gr-qc/0005050) |
| `NumPy` | [HarrisEtAl-2020-array-programming-NumPy.pdf](Fonti/HarrisEtAl-2020-array-programming-NumPy.pdf) | 19 | [arXiv:2006.10256](https://arxiv.org/abs/2006.10256), preprint dell'articolo OA |
| `SciPy` | [VirtanenEtAl-2020-SciPy-1.0.pdf](Fonti/VirtanenEtAl-2020-SciPy-1.0.pdf) | 22 | [arXiv:1907.10121](https://arxiv.org/abs/1907.10121), preprint dell'articolo OA |
| `SymPy` | [MeurerEtAl-2017-SymPy-symbolic-computing.pdf](Fonti/MeurerEtAl-2017-SymPy-symbolic-computing.pdf) | 27 | PDF editoriale PeerJ, CC BY 4.0 |
| `BuchstaberEnolskii1997` | [BuchstaberEnolskiiLeykin-1997-hyperelliptic-Kleinian-functions.pdf](Fonti/BuchstaberEnolskiiLeykin-1997-hyperelliptic-Kleinian-functions.pdf) | 24 | [arXiv:solv-int/9603005](https://arxiv.org/abs/solv-int/9603005) |
| `Zagier1990` | [Zagier-1990-Bloch-Wigner-Ramakrishnan-polylogarithm.pdf](Fonti/Zagier-1990-Bloch-Wigner-Ramakrishnan-polylogarithm.pdf) | 12 | [pagina d'autore, MPIM](https://people.mpim-bonn.mpg.de/zagier/files/doi/10.1007/BF01453591/fulltext.pdf) |
| `Timofeev1978` | [Timofeev-1978-constancy-adiabatic-invariant.pdf](Fonti/Timofeev-1978-constancy-adiabatic-invariant.pdf) | 4 | [archivio ufficiale JETP](https://www.jetp.ras.ru/cgi-bin/dn/e_048_04_0656.pdf) |

### Due sostituti da non scambiare per l'opera citata

1. `Liberzon`: il PDF locale è una versione autorizzata delle lecture notes. La
   prima pagina avverte espressamente che differisce dall'edizione Princeton per
   correzioni, editing e numerazione. È utile per controllare gli argomenti, ma
   **non certifica una citazione a pagina dell'edizione pubblicata**.
2. `Bishop1972`: [SinghGupta-undated-Clairaut-submersions-Bishop-theorem.pdf](Fonti/SinghGupta-undated-Clairaut-submersions-Bishop-theorem.pdf)
   è un capitolo moderno che espone e cita il teorema di Bishop. Non è
   “Clairaut submersions”, pp. 21–31, nel volume del 1972.

## Pubblicazioni per cui manca ancora la versione esatta citata

### Articoli o capitoli

| Chiave | Stato al 25 agosto 2026 | Azione lecita consigliata |
|---|---|---|
| `Bishop1972` | L'originale nel volume in onore di K. Yano non è stato trovato in un repository autorevole; il PDF locale è solo secondario. | Prestito interbibliotecario o scansione personale ottenuta dalla biblioteca. |
| `Thakurta1981` | SciSpace lo indica come open access, ma non è stato individuato un file primario/editoriale con provenienza verificabile. | Cercare nell'archivio dell'Indian Journal of Physics o chiedere all'autore/editore; importare solo una copia con provenienza chiara. |
| `DyerHonig1979` | DOI `10.1063/1.524078`; nessuna copia OA di repository verificata. | Accesso AIP tramite biblioteca o richiesta agli autori. |
| `SultanaDyer2004` | DOI `10.1063/1.1814417`; nessuna copia OA di repository verificata. | Accesso AIP tramite biblioteca o richiesta agli autori. |
| `Matplotlib` | DOI `10.1109/MCSE.2007.55`; nessuna copia OA verificata. | Accesso IEEE tramite biblioteca o richiesta di author manuscript. |
| `BeilinsonLevin1994` | Capitolo in *Motives*, Proc. Symp. Pure Math. 55(2), pp. 126–196; nessuna copia autorizzata verificata. | AMS/biblioteca o richiesta agli autori. |
| `Ree1958` | DOI `10.2307/1970243`; nessuna copia OA verificata. | JSTOR/biblioteca. |
| `Neishtadt1986` | L'articolo esatto del *Soviet Journal of Plasma Physics* non è stato trovato in un archivio ufficiale; sono locali lavori correlati e ora anche Timofeev 1978. | Richiesta all'autore o ricerca nell'archivio della traduzione sovietica. |

### Libri

Non sono state importate copie non autorizzate di questi volumi:

- `Pontryagin` — *The Mathematical Theory of Optimal Processes*.
- `Cesari1983` — *Optimization—Theory and Applications*.
- `Chandrasekhar` — *The Mathematical Theory of Black Holes*.
- `FaraoniBook` — *Cosmological and Black Hole Apparent Horizons*.
- `WhittakerWatson` — *A Course of Modern Analysis*, 4a edizione.
- `Fay1973` — *Theta Functions on Riemann Surfaces*, DOI `10.1007/BFb0060090`.
- `BakerAbelian1897` — scansioni statunitensi dichiarate public domain esistono,
  ma nell'agosto 2026 la situazione territoriale UE richiede cautela; una copia
  da biblioteca è la via più pulita.
- `MumfordTata` — *Tata Lectures on Theta I, II*.
- `Bloch2000` — *Higher Regulators, Algebraic K-Theory, and Zeta Functions of
  Elliptic Curves*.
- `BenderOrszag1999` — *Advanced Mathematical Methods for Scientists and
  Engineers*.
- `GriffithsHarris1978` — *Principles of Algebraic Geometry*.
- `BairdWood2003` — *Harmonic Morphisms Between Riemannian Manifolds*, DOI
  `10.1093/acprof:oso/9780198503620.001.0001`.

Per questi libri la soluzione consigliata è l'accesso istituzionale, il prestito
interbibliotecario o una scansione personale nei limiti consentiti dalla biblioteca.

## Riferimenti che non richiedono un PDF

Queste cinque voci sono correttamente trattate come risorse digitali, non come
“articoli mancanti”:

- `DLMF`: risorsa web versionata del NIST.
- `Rosignoli2026Code`: release software Zenodo/GitHub.
- `Rosignoli2026CodeJMP`: release software Zenodo/GitHub.
- `SageMath`: software e relativa documentazione.
- `abelfunctions`: repository software; il tree esatto è archiviato su Zenodo.

Per la riproducibilità può essere utile conservare in futuro anche gli ZIP delle
due release Zenodo, ma non aumenterebbero la copertura degli articoli scientifici.

## Correzioni bibliografiche consigliate, non ancora applicate

1. La chiave `SchmidtSmith` è fuorviante: i campi autore sono correttamente
   Frittelli e Newman. Aggiungere DOI `10.1103/PhysRevD.59.124001` e, se lo stile
   lo consente, `eprint = {gr-qc/9810017}`. Una rinomina della chiave va fatta solo
   con sostituzione globale delle citazioni nei `.tex`.
2. La chiave `NielsenYeom2009` è fuorviante: i campi sono correttamente Nielsen e
   **Yoon**, anno 2008. La rinomina naturale sarebbe `NielsenYoon2008`, sempre con
   sostituzione globale.
3. La chiave `ClaudelViracVirbhadra2001` contiene probabilmente il refuso `Virac`.
   I campi autore sono corretti. Aggiungere DOI `10.1063/1.1308507` e
   `eprint = {gr-qc/0005050}`.
4. `Liberzon` riporta `year = {2011}`, mentre cataloghi Princeton/JSTOR indicano
   l'edizione pubblicata nel **2012**. Verificare la copia effettivamente citata e
   uniformare anno, ISBN/DOI e citazioni di pagina.
5. La presenza di un PDF locale non giustifica automaticamente una citazione:
   per ogni affermazione centrale va controllato il passaggio preciso e, quando
   disponibile, va citata la numerazione della versione editoriale.

## Riproducibilità dell'audit

Gli strumenti temporanei, deliberatamente ignorati da Git, sono:

- `scratch_codex/inventory_bibliography_codex.py` — associazione indicativa tra
  BibTeX e nomi file;
- `scratch_codex/download_open_sources_codex.py` — manifest dei soli download
  autorizzati e controllo PDF;
- `scratch_codex/download_open_sources_log_codex.json` — dimensioni, pagine e
  SHA-256 dei quindici nuovi documenti.

L'associazione automatica è stata riesaminata manualmente proprio per eliminare
il falso positivo `Bishop1972`. Il conteggio finale di questo documento, non il
solo output dello script, è quello da usare come stato ufficiale.
