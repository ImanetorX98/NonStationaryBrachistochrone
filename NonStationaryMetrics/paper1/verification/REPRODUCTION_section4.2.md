# Reproduction log — section 4.2 phenomenology and quoted numerics

Section 4.2 makes claims about integrated trajectories rather than closed-form
identities, so it cannot be settled by the symbolic suite in this directory. It
was instead reproduced by re-running the archived scripts. Every number below
came out of those scripts as shipped; nothing was reimplemented.

Parameters throughout: `E = 1.2`, `J = 1.3`, `r_1 = 10`, `m(v) = 1 + mu v`.

## Penetration threshold (`VaidyaMetric/vaidya_penetration_map.py`)

Claim: accretion opens a finite window `J_pen > 0`; the negative-rate
continuation closes it to zero measure.

| `mu` | `v_0 = 0` | `v_0 = 30` | `v_0 = 60` |
|---|---|---|---|
| `+0.01` | `J_c = 0.91851` | `1.10642` | `1.28634` |
| `-0.01` | `0.00000` | `0.00000` | `0.00000` |

Confirmed. The script also reports the small-rate scaling `J_c/sqrt(mu)` drifting
from 8.33 to 10.49 over `mu` in `[0.0025, 0.02]`, and notes that the static case
has `J_c = 0` up to a numerical floor — consistent with the text's statement that
the static threshold is of zero measure.

## Bounce (`VaidyaMetric/vaidya_brachistochrone_vparam.py`)

Claim: the `v`-parametrisation integrates through the turning point, reaching
periapsis exactly where the `r`-parametrisation diverges; `r_min = 3.0105` at
`v = 7.29` for `mu = 0.01`.

Reproduced: `r_min = 3.010492` at `v = 7.2857`, with the script reporting that
the `r`-parametrised flow "died at r = 3.0105". Supporting checks from the same
run: the static periapsis integrates to `2.72713516` against the exact root
`2.72713516` (difference `9.0e-12`), and the `v`- and `r`-parametrisations agree
on `v(3.5)`, `Dphi` and `T_tau` to between `5e-8` and `5e-7`.

## Timing (`VaidyaMetric/vaidya_brachistochrone_vparam.py`, section V4)

Claim: absolute and relative penetration have opposite winners.

The script prints `r_min` but not the ratio; the ratio was formed by calling its
own `orbita()` with `2m(v_peri)`.

| `v_1` | `mu = -0.01`: `r_min` / ratio | `mu = +0.01`: `r_min` / ratio |
|---|---|---|
| 25 | `2.6816` / `1.2656` | `2.7628` / `1.4799` |
| 40 | `2.4153` / `1.3431` | `3.0105` / `1.4030` |
| 55 | `2.1522` / `1.4549` | `3.2598` / `1.3450` |
| 70 | `1.8907` / `1.6267` | `3.5112` / `1.2999` |

Confirmed: under the negative rate the absolute periapsis descends while the
ratio rises (penetration worsens); under accretion the periapsis rises while the
ratio falls (late arrivals graze the relative horizon).

Note: for the earliest arrivals the periapsis occurs at negative advanced time
(`v_peri = -5.94` and `-6.66` at `v_1 = 25`), so the linear mass model is
evaluated slightly before `v = 0`. It remains positive there (`m = 0.933` at
worst), so no unphysical `m <= 0` is reached, but the trajectories do extend
outside `[0, v_1]`.

## Quoted numerics elsewhere

`VaidyaMetric/universal_SD_source_vaidya.py`:

- `m H_m + (r H_r + J H_J) = 0` and `(J d_J + p_r d_pr) H - (H+1) = 0`, both
  exactly zero — the two identities proved symbolically as `eq:selfsimilar` and
  `eq:finsler-euler`
- `max |S_D - ([r p_r] - lambda)| = 3.24e-9` on the frozen orbit, the residual
  being quadrature error; the identity itself is proved in appendix C

`VaidyaMetric/vaidya_first_order_offshell.py`:

- `SLOPE full = 2.00`, `SLOPE on-shell only = 1.00`, matching the text's
  `O(eps^2)` closure for the complete term and `O(eps)` for the on-shell piece
  alone

## Verdict

All three phenomenological claims of section 4.2 reproduce, as do the quoted
numerics. No discrepancy was found. These are reproductions of the author's own
integrations, not independent derivations: they confirm that the text reports
what the scripts produce, not that the scripts model what they intend to.
