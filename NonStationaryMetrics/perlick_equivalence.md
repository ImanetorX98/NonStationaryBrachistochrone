# Equivalenza controlled-rail ↔ Perlick (limite stazionario) — DIMOSTRATA

Chiude l'Issue 1 del referee ("recovery of Perlick formalism") a livello funzionale: la funzione di
Finsler tempo-d'arrivo del vincolo controlled-rail **è** la metrica ottica di Randers di Perlick per
il problema a energia fissata in uno spazio-tempo stazionario. §3.3 della nota post-referaggio la
lasciava "da mostrare" — ora è mostrata e verificata a precisione macchina.

## Enunciato
Metrica stazionaria `g_tt=−f, g_ta, g_ab` (t-indipendenti), rail `−u·ξ=Ê` (ξ=∂_t Killing), `g(u,u)=−1`.
Minimizzando il tempo d'arrivo `∫dt`, la Finsler `F_rail=dt/dσ` risolta dalla shell quadratica è
```
F_rail = β_a v^a + √(a_ab v^a v^b)          (RANDERS)
  β_a  = g_ta/f
  a_ab = Ê²(g_ta g_tb + f g_ab)/(f²(Ê²−f))
```

## Forma threading (Perlick esplicito)
Con `g = −f(dt − ω_a dx^a)² + h_ab dx^a dx^b`: `g_tt=−f`, `g_ta=f ω_a`, `g_ab=h_ab−f ω_a ω_b`, quindi
`g_ta g_tb + f g_ab = f h_ab`. Segue
```
β_a  = ω_a                    (1-forma frame-dragging / shift)
a_ab = Ê² h_ab / (f(Ê²−f))    (metrica Jacobi–Maupertuis a energia fissata)
```
- **Limite null `Ê→∞`**: `a_ab → h_ab/f` = optical metric di Fermat/Perlick standard per raggi di luce.
- **Massivo (Ê finito)**: `a_ab = h_ab/(f(1−f/Ê²))` = metrica ottica a energia fissata di Perlick.

## Verifica
`KerrSessionScripts/perlick_equivalence.py`: Kerr equatoriale, 16 casi (r, v^r, v^φ, Ê diversi):
`F_rail` (risolta da rail+mass-shell) == forma Randers a **0–4.4e-16** (precisione macchina).

## Conseguenza per il paper
Chiude l'Issue 1/§3.3: il gauge ottico del controlled-rail = intersezione mass-shell ∩ livello-energia
È la metrica di Randers/Perlick a energia fissata. Da aggiungere come proposizione autonoma:
"*Nel limite stazionario (W=ξ Killing), F_rail è la metrica di Randers β_a v^a+√(a_ab v^a v^b) con
β_a=ω_a, a_ab=Ê² h_ab/(f(Ê²−f)); ossia il formalismo controlled-rail recupera l'ottica di Perlick a
energia fissata, con limite null Ê→∞ = Fermat.*" Distinto dall'equivalenza al livello di flusso
Hamiltoniano (X_{H_rail}∥X_{H_Perlick} su H=0), qui più forte (uguaglianza delle Finsler).
