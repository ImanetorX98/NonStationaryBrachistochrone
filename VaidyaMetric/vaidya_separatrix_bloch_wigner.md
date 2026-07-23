# Chiusura Bloch–Wigner del dilog d'orizzonte D_k sulla separatrice

Su `|J|=Jc` la sestica `S` acquista radice doppia, il genere scende a 1 e il dilog
iperellittico `D_k = ∫U_k/(r-2m) dr` chiude nel **dilogaritmo ellittico
(Bloch–Wigner)**, forma tabulata e riutilizzabile.

## Dati verificati (M=1, E=1.4)
- **Jc = 7.02662374** (radice doppia del bracket `B=r²(r-2m)-J²DE`, via discriminante).
- **r_d = -3.36371118** (radice doppia; `√S = (r-r_d)√Q4`, `Q4` quartica).
- Curva ellittica **E: w² = Q4(r)**, `a4 = E²-1 = 0.96`, punti di dirama­zione
  (radici Q4) = **{-2.08333, 0, 2, 8.72742}**.
  - `-2.08333 = -2m/(E²-1)` (radice di DE); `2 = 2m` (**orizzonte**); `8.72742`
    turning fisico; `0`.
- **Modulo ellittico k² = 0.60671566**;  periodo reale = 0.66912986,
  periodo immag = 0.60621385.
- **τ = 0.90597338 i** (reticolo RETTANGOLARE). Verificato DUE modi indipendenti:
  Legendre `τ=iK'(k)/K(k)` **=** Sage `EllipticCurve.period_lattice()` a 15 cifre.

## Riduzioni verificate
1. **U_k ellittici** (Stage A): `U_k = ∫r^k/[(r-r_d)√Q4]dr` vs `∫r^k/√S` diretto
   genus-2 @Jc → coincidono a **1e-14…1e-17** (precisione macchina). Split
   `r^k/(r-r_d)=poly_{k-1}+r_d^k/(r-r_d)` ⇒ ellittico 1a/2a specie + 3a specie
   (polo r_d).
2. **Mappa di Abel** `z(r)=∫dr/√Q4`: `dz/dr = 1/√Q4` verificato (1e-7).
3. **Orizzonte r=2m = punto 2-TORSIONE**: essendo `r=2m` branch point di Q4, la sua
   immagine di Abel è a mezzo-periodo immaginario `z_h = ω_im/2` (verificato,
   diff 7e-13). Lettera dlog al 2-torsione.

## Forma chiusa — MATCH NUMERICO (vaidya_ell_dilog_match.py)
`D_k = U_k ln(r-2m) - G_k`,  `G_k = ∫ln(r-2m) r^k/√S dr`.

### PESO-1: U_0 in Weierstrass σ,ζ (3a specie) — verificato 1e-14
Toolkit σ,ζ da θ₁ (mpmath), reticolo auto-consistente (om1=0.66913, om2=0.60621 i,
τ=0.90597 i). Base al PUNTO DI WEIERSTRASS e4 (z=0) ⇒ involuzione z→−z ⇒ le due
preimmagini di r_d sono ±z_d.

    U_0(r) = ρ[ ln σ(z−z_d) − ln σ(z+z_d) ] + C z + cost
    ρ = 1/√Q4(r_d) = 0.061069   (RESIDUO derivato, NON fittato)
    C = 1/(e4−r_d) + 2ρ ζ(z_d) = 0.333272   (parte olomorfa, derivata)

Verifica derivata dU_0/dz = 1/(r−r_d) e primitiva U_0: diff **1e-14** (macchina).
z_d = 0.46104 (immagine di Abel di r_d, via ovale reale per infinito).

### PESO-2: dilog d'orizzonte = dilogaritmo ellittico
`G_0 = ∫ ln(r−2m) dU_0 = ∫ ln(r(z)−2m)[ρ(ζ(z−z_d)−ζ(z+z_d))+C] dz`.
`ln(r−2m)` ha zero doppio al 2-torsione z_h e poli doppi ai due infiniti ⇒ G_0 =
combinazione di `∫ ln σ(z−a) ζ(z−b) dz` = **dilogaritmo ellittico**.
Verifica cerchio `D_0 = U_0 ln(r−2m) − G_0` (con U_0 in σ,ζ): diff **1e-14**.

### Funzione tabulata (Zagier)
`D^E(z) = Σ_{n≥0} D(qⁿζ) − Σ_{n≥1} D(qⁿ/ζ)`, `q=e^{2πiτ}`, `ζ=e^{2πi z/2ω1}`,
D=Bloch–Wigner. Implementato, converge; D^E(0.2+0.3i)=0.59621, D^E(z_d)=0.67494.
NB (correzione): z reale ⇒ ζ sul cerchio unitario (fase) ⇒ qⁿζ complessi ⇒ D^E≠0:
il dilog ellittico singola-valore È il contenuto (NON si annulla su config. reale).

### Cosa resta (unico punto aperto)
La normalizzazione ESATTA del regulator che lega l'integrale fisico (reale,
indefinito) D_0 a una combinazione specifica di D^E — è teoria dei polilog ellittici
(regulator reale+immag di Beilinson–Levin). Struttura + riduzione ellittica: fatte a
macchina. Funzione tabulata: in mano.

Script: `vaidya_separatrix_ell.py` (Stage A), `vaidya_bloch_wigner.py` (τ, Abel,
2-torsione), `vaidya_ell_dilog_match.py` (σ,ζ closed form + D^E q-serie),
`/tmp/sep_abel2.sage` (Sage cross-check τ).
