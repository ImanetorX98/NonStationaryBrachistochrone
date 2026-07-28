# Chiusura off-shell generico — Vaidya (ramo τ), verificata math + fisica

## Esito
Il wrap off-shell generico (`J≠Jc`) del ramo **τ di Vaidya** è chiuso in forma esplicita, con
**coefficienti simbolici in tutti i (m,E,J)**, e verificato **contro la matematica E contro la
fisica**. Il ramo **v** ha struttura derivata ma non ancora assemblata (blocco elementare extra).

## Forma pulita (niente sqrt-hell)
Hamiltoniana τ NON-surrogata (`−f/E`, non `−1`): la shell è pulita,
```
S = r(r-2m)·DE·Q3,   DE=(E²-1)r+2m,   Q3=r²(r-2m)-J²·DE,   p_r0 = -√S/(r(r-2m)DE)  (ingoing)
```
Semplificazioni analitiche sulla shell (verificate vs numerico full-Hamiltoniano a 1e-16):
```
G = J/(Δ p_r),   H_pr = (r-2m)DE·p_r/(E r²),   H_m = C0(r) + C2(r) p_r²   (il √(Δv) si cancella)
kernel  d_pr G/H_pr = A_V/√S ,  A_V = J E r³ DE / Q3
inner   m H_m/H_pr  = P_inner/√S ,  P_inner = -m·N4 / ((r-2m)DE)
```
`N4 = E⁴J²r² + 4E²J²mr − 2E²J²r² + 4J²m² − 4J²mr + J²r² + 4m²r² − 4mr³ + r⁴`.

## Assemblaggio A+B+C
Sorgente `Θ = m∂_m` (modulo metrico, **niente dilatazione** — a differenza di TK).
Wrap `Φ = -∫ (A_V/√S)·Σ dr`, `Σ(r)=∫ P_inner/√S dz`:
- **A (2ª specie)**: kernel-poly × `Σ`-poly → `Σ Q_jk W_jk` (Kleiniane ζ,σ + dilog genus-2).
- **B (terza specie)**: `Σ` ha poli a **r=2m (orizzonte Schwarzschild)** e DE=0 (fuori arco);
  residui `ρ(2m)=-2E²J²m²`, `ρ(DE)=8E²m⁴/…` → dilog genus-2 all'orizzonte.
- **C (Hermite)**: `A_V=A_poly + N3/Q3` (poli a Q3=0, branch points); riduzione di Hermite
  `N3/(Q3√S)=d[P√S/Q3] + Σremₖ rᵏ/√S`, remainder esatto = 0 SIMBOLICO; termine algebrico → elementare.
Kernel effettivo `g_k = p_j + remₖ`; coefficienti tutti razionali in (m,E,J).

## Verifiche (math + fisica)
| Check | Risultato |
|---|---|
| kernel `A_V`, inner `P_inner` vs numerico on-shell | 1e-16 |
| Hermite exact-division remainder | 0 (simbolico) |
| **A+B+C == wrap diretto** (config 1: m=1,E=1.4,J=2.5) | **1e-16** |
| **A+B+C == wrap diretto** (config 2: m=6/5,E=3/2,J=3) | **1e-15** |
| **PHYSICS**: closed form == sub-pezzo off-shell del VERO flusso τ | **6.9e-10** (floor flusso) |
| **PHYSICS**: totale (partA+partB) vs vero flusso non-autonomo τ | **slope 2.03** (~2, O(ε²)) |

Nota metodologica: il physics anchor ha inizialmente rivelato **slope 1** — dovuto alla convenzione
di segno (Vaidya `m` CRESCE, opposto a TK dove E,J calano: `δp_r=(S_D−λΘH)/H_pr`). Corretto il segno,
slope → 2.03. Conferma il valore del confronto con la fisica (non solo math↔math).

## Scripts
- `VaidyaMetric/vaidya_offshell_shellpoly.py` — forme pulite A_V, P_inner (verifica vs numerico)
- `VaidyaMetric/vaidya_offshell_FULL_assembly.py` — A+B+C, verifica vs diretto (2 config)
- `VaidyaMetric/vaidya_tau_physics_anchor.py` — ancora fisica (slope 2.03, closed==partB)

## Aperto: ramo v
Il ramo v (tempo avanzato, costo "−1" con termine lineare in p_r) ha `p_r0` con parte additiva
razionale (curva spettrale `D_v`=discriminante). Kernel `= A_part/D_v^{3/2} (2ª specie) + B_part/D_v
(ELEMENTARE)`. Chiudibile con A+B+C + blocco elementare extra; struttura derivata, assemblaggio TODO.
