# Brachistocrone in Schwarzschild–de Sitter (Kottler) — risultati

**Script**: `sds_brachistochrone_sympy.py` (tutte le verifiche passano).
Contesto: `../VaidyaMetric/VaidyaResults.md` (formalismo EF/Hamiltoniane),
`../FLRWmetric/FLRWresults.md`, `../KerrMetric/doranT.md` (dicotomia t).

```
f(r) = 1 − 2M/r − Λr²/3      (EF avanzate: ds² = −f dv² + 2dv dr + r²dΩ²)
```

## Stazionarietà — il posizionamento rispetto a Perlick 1991

`∂_t` è Killing **genuino**, ma timelike solo nella **vasca**
`r_b < r < r_c` (f > 0). Oltre r_c (come sotto r_b) è spacelike: lo
spaziotempo è dinamico lì. Perlick 1991 vale solo nella vasca; la sua
riduzione ottica degenera (`1/F`) a **entrambi** gli orizzonti. Il
worldline vincolato/Hamiltoniane EF li attraversa: SdS = "Perlick
recuperato dentro, oltrepassato ai bordi".

## R1. Struttura (simbolico, esatto)
- orizzonti r_b < r_c ⇔ `0 < 9ΛM² < 1`;
- `f_max = 1 − (9ΛM²)^{1/3}` a `r* = (3M/Λ)^{1/3}`;
- sfera fotonica `r_ph = 3M` **esatta** (Λ si cancella da `rf′ = 2f`);
- congelamento statico: `v = √w/E = 0` dove `f = E²` (`w = E²−f`):
  per `E < √f_max` **due barriere spaziali** — l'analogo statico del
  congelamento FLRW (`a = Ê`). Per viaggiare ovunque nella vasca serve
  `E > √f_max`.

## R2. Soglie ai due orizzonti (simbolico)
- **ramo τ**: `g_τ = r²f − J²w = −J²E² < 0` a f=0 ⇒ **riflette** per ogni
  `J ≠ 0` a *entrambi* gli orizzonti (J_c = 0, misura nulla, doppia);
- **ramo t/v**: `𝒦 = fJ/E → 0` agli orizzonti (auto-sintonizzazione,
  doranT con a=0); `∂_r g_t|_{f=0} = r_h²f′(r_h)`: `> 0` a r_b, `< 0` a
  r_c ⇒ regione permessa confina con entrambi i bordi ⇒ il ramo t
  **attraversa marginalmente r_b E r_c, per ogni J**.

**Dicotomia massimale nella vasca**: τ rimbalza in una scatola
`(r_in, r_out)` strettamente interna; t tocca entrambi gli orizzonti.

## R3. Numerico (M=1, Λ=0.03, E=1.2, J=1.3; Hamilton, H≡0)
- r_b = 2.091488, r_c = 8.788851, f_max = 0.354;
- svolte τ = radici di g_τ a **1e-12** (interna 2.947238, esterna
  8.568457 — riflessione dall'orizzonte cosmologico);
- `dφ/dr` lungo il flusso = forma chiusa `𝒦r√(wf)/(Δ√(Δ−𝒦²w))` con
  `Δ = r²f_SdS` a 10 cifre — **le forme chiuse Kerr a=0 valgono per f
  generica**;
- ramo v: raggiunge r_b e r_c a 1e-6 (f ~ 1e-7): attraversamento
  marginale confermato.

## Uso nel programma
- limite `a → cost` di Thakurta e limite tardo di McVittie (Λ>0);
- laboratorio a doppia barriera per la fenomenologia "scatola";
- prossimo: Vaidya–de Sitter (`f = 1−2m(v)/r−Λr²/3`) = corsa
  orbita-orizzonte con secondo orizzonte.
