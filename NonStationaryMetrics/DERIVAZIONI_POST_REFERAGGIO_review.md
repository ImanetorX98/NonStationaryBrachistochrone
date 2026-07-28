# Review di validità — `derivazioni_post_referaggio.pdf` (28 lug 2026)

Verifica sezione per sezione della nota tecnica "Derivazioni e dimostrazioni successive al
referaggio di main11" (10 pp). **Verdetto: derivazioni e dimostrazioni VALIDE sotto le ipotesi
dichiarate; nessun errore matematico; onesta sui limiti.** La §4 (adiabatica) coincide con la
verifica numerica indipendente di questa serie di sessioni.

## Stato per sezione

| Sez | Contenuto | Validità |
|-----|-----------|----------|
| §2 | Indicatrice / gauge di Minkowski `F`, tempo d'arrivo ODE | **valida**. Eq (5) `(1/R)|v|≤F≤(1/ρ)|v|` corretta (verificate entrambe direzioni da `ρB₁⊂B⊂RB₁`). Eq (8) risolve Issue 6 (bordo reale, non convessificazione) |
| §3.1 | Normalità (Prop 3.1) | **valida, pulita**. `p₀=0`⟹`max p·v=0` via `H*=h_B+p₀=0`; ma `h_B(p)≥ρ|p|>0`⟹`p=0`, contraddizione. Chiude Issue 6 (normalità/non-annullamento costato) |
| §3.2 | Limite stazionario (Prop 3.2) | **valida**. `∂_sH=0`⟹`H`cost⟹`H≡0` via `H(s_f)=0` |
| §3.3 | Recupero Perlick | **INCOMPLETA (dichiarato)**. Eq (15) `F_rail=F_Perlick` da mostrare, non provata → Issue 1/2 referee |
| §4 | Correzione adiabatica: `−(Ȧ/A)P_r`, `D=P_r∂_Pr+E∂_E+J∂_J`, Eq (28) | **valida, coerente col numerico**. Verificato algebricamente che usare `D` uniformemente = mia macchina con `Θ` (l'extra `λP_r` in `δP_r` cancella l'extra `P_r∂_Pr` in `DG`). Oss 4.1 = il bug di validazione circolare che avevo trovato |
| §5 | Esistenza (Teorema 5.1, metodo diretto) | **valida**. Tonelli: coercività(5)+Arzelà–Ascoli+semicontinuità da convessità in `v`; `s`-coupling via stabilità ODE+Grönwall. Oss 5.2 onesta sui limiti |
| §6.1 | FLRW piatta | **valida**. Triangolare, uguaglianza⟺direzione costante⟹retta comovente unica |
| §6.2 | Schwarzschild radiale | **valida**. `√(ṙ²+fr²φ̇²)≥|ṙ|`, `=`⟺`φ̇=0`, classe radiale monotona |
| §7 | HJB + Teorema 7.1 (verifica) | **valida**. DPP→HJB corretto; sub-soluzione Lipschitz⟹bound inferiore⟹ottimalità globale su saturazione. Fornisce il certificato globale (PMP necessario ≠ globale) |
| §8 | Coniugati/Maxwell/cut locus | **valida**. Jacobi (44-46), Maxwell come shock HJB (47-48); nota che convessità stretta non esclude Maxwell |
| §9-11 | Proposizioni + frontiera + testo proposto | **oneste**. Prop 9.1-9.3 seguono; lista aperta corretta |

## Punti INCOMPLETI del documento (da chiudere)
1. **§3.3 — Equivalenza esplicita con Perlick** e **relazione `p_φ` (costato) ↔ momento meccanico
   `L_mech`**: il documento asserisce `E_eff=Ê/A, J_eff=J/A` per peso conforme e lascia
   `F_rail=F_Perlick` (oppure `X_{H_rail}∥X_{H_Perlick}` su `H=0`) da mostrare. Il referee (Issue 1/2)
   vuole la mappa di Legendre DERIVATA, non asserita. **Aperto.**
2. **§9 lista aperta** — parzialmente superata: l'"assemblaggio speciale chiuso del termine off-shell"
   è stato FATTO per il t-branch generico (verificato 1e-14, ancorato al flusso, commit e903544/f15c2c4).
   Restano: Vaidya generico + theta-naming.
3. Esclusioni globali (dal documento stesso §9): funzione valore esplicita Kerr/Vaidya/TK; esclusione
   punti di Maxwell tra winding diversi; minimalità dentro l'ergosfera (selettore spacelike); matching
   uniforme attraverso turning points/separatrici; prova globale equivalenza Perlick a livello di
   funzionale.

## Azione consigliata sul paper
Il "Testo proposto" (§10) è prudente e corretto — integrabile come blocco control-theory (esistenza +
normalità + HJB + coniugati/Maxwell), che chiude gran parte degli Issues 1/2/6 del referee a livello di
FORMULAZIONE (la globalità resta condizionale, come giustamente dichiarato).
