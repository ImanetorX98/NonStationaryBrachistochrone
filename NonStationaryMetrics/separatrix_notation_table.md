# Issue 3 — Tabella di notazione delle soglie (autoritativa)

Il referee (Issue 3) chiede una notazione UNICA e coerente: `\Jc` era usato per oggetti inequivalenti.
Classificazione autoritativa: `separatrix_classification.py` (regression test per la terminologia).

## Notazione proposta
| Simbolo | Significato | Come identificarlo |
|---|---|---|
| `J_deg` | radice del discriminante con `r_d` NON fisico (negativo o interno a r_+) | doppio root sextica ma `r_d<r_+` o `r_d<0` |
| `J_sep^phys` | separatrice fisica ESTERNA | doppio root con `r_d > r_+` |
| `J_ergo` | soglia ergosferica/grazing (`≈ a/E`) | collasso ovale / ingresso ergosfera |
| `J_pen(v_0)` | soglia dinamica di penetrazione Vaidya | finestra `(J_pen^-,J_pen^+)` |

## Classificazione verificata (separatrix_classification.py)
- **Schwarzschild τ, E=1.4**: `J=7.0266, r_d=−3.3637` → **`J_deg`** (raggio NEGATIVO, degenerazione algebrica).
  NON è una separatrice fisica esterna. **Va etichettato `J_deg/m`, non `\Jc`** (righe main.tex ~737, ~1829).
- **TK t-branch, a=0.9, E=1.2**:
  - `J=−8.0535, r_d=+3.5139` → **`J_sep^phys`** (retrograda esterna)
  - `J=+2.9364, r_d=+1.5123` → **`J_sep^phys`** (prograda, in ergosfera; la separatrice off-shell/tracking)
  - `J=+1.2666, r_d=+0.5509` (interno), `J=+19.089, r_d=−6.62`, `J=−18.67, r_d=−6.59` → **`J_deg`**

## Fix richiesti nel testo (per propagazione)
1. `7.0266` (Schwarzschild/Vaidya τ, E=1.4) → etichettare `J_deg/m` ovunque; NON "separatrice".
2. Rimuovere l'interpretazione capture/escape da `J_deg` (raggio negativo): l'annullamento algebrico
   del triplo polo resta valido, ma NON implica una transizione esterna cattura/fuga.
3. Ogni doppio root nel paper/figura: stampare `r_d`, `r_+`/`r_ergo`, e dire se `r_d` è nel dominio fisico
   (usare `separatrix_classification.py` come riferimento).
4. `\Jc` → uno solo dei quattro simboli per ogni occorrenza.

Nota: TK = oggetto compatto (non BH); `Δ=0` = superficie null seme Kerr, non orizzonte (memoria
thakurta-kerr-not-black-hole). Linguaggio BH/orizzonte solo per Vaidya.
