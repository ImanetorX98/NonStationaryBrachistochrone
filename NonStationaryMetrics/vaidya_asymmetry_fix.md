# Issue 4 (Vaidya accrescimento/evaporazione) — FIX

## (a) Segno dell'asimmetria (referee Issue 4)
Coefficienti (δφ/ṁ con segno): `accr = A+B`, `evap = A−B`, con
`A = ∫dM F·E U_3` (infinito, polilog), `B = ∫dM F·r_*` (orizzonte, dilog tortoise).
FISICA (`δφ = ṁ·[·]`; ṁ>0 accr, ṁ<0 evap): `δφ_accr=+|ṁ|(A+B)`, `δφ_evap=−|ṁ|(A−B)`. Quindi
- **Differenza** `accr−evap = 2A = +17.455` (portata da A, NON B)
- **Somma** `accr+evap = 2B = +12.225` = **fallimento dell'antisimmetria** (portata da B, orizzonte)

La vecchia dicitura "asimmetria netta = 2B" confondeva la differenza dei COEFFICIENTI `(A+B)−(A−B)=2B`
con la differenza FISICA `accr−evap=2A`. Corretto in `vaidya_asymmetry.py` (verificato).

## (b) Clock uscente DERIVATO (non flip ad-hoc)
Metrica outgoing `ds²=−f du²−2 du dr`. `u=t−r_*`, `v=t+r_*` (t = tempo di Killing condiviso) ⟹ lungo
qualsiasi traiettoria `v−u=2r_*`, quindi `du/dr = dv/dr − 2 dr_*/dr`. Con `dv/dr=E r³/√S+r/(r−2m)`
(v-branch, metrica ingoing) e `dr_*/dr=r/(r−2m)`:
```
du/dr = E r³/√S − r/(r−2m)   =>   u = E U_3 − r_*
```
DERIVATO dalla metrica outgoing via l'identità del tortoise. Verificato `du/dr==dv/dr−2dr_*/dr` a 2.2e-16.

## Nomenclatura (referee Issue 4)
Accrescimento = INGOING (v, ṁ>0). Evaporazione = OUTGOING (u, ṁ<0). NON presentare massa-decrescente
in chart ingoing come evaporazione. La riduzione `dM F` (S congelata) è comune: ingoing/outgoing
Schwarzschild = stessa geometria; cambiano solo clock (v↔u) e segno di ṁ.

Script: `VaidyaMetric/vaidya_asymmetry.py`.
