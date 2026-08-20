# Symbolic verification of the formal core of Paper I

Two independent computer-algebra checks of every algebraic step that
`paper1/paper1_JMP.tex` states as proved in sections 2 and in appendices A–C.
The two scripts are deliberately redundant: they verify the same claims in
different systems, so an artefact of one simplifier cannot pass unnoticed.

```
wolframscript -file verify_paper1_core.wls     # 55 checks
python3       verify_paper1_core.py            # 22 checks, exit 1 on failure
wolframscript -file verify_appB_residues.wls   # 9 checks (appendix B residues)
wolframscript -file verify_appB_blocks.wls     # 4 checks (appendix B block form)
wolframscript -file verify_perlick_recovery.wls # 6 checks (stationary limit)
python3       verify_section44_closedform.py   # end-to-end quadrature of eq:vaidya-full
```

Both print `[OK]` per claim, and every `[OK]` is an exact symbolic zero — not a
numerical tolerance. `verify_paper1_core.py` exits non-zero if any check fails,
so it can be dropped into a regression run.

## What is checked, and where it appears in the manuscript

| Claim | Manuscript | Check |
|---|---|---|
| Rail locus is empty / a point / a sphere as `Ehat < = > \|W\|` | Lemma I.A(i)–(iii) | rest-frame decomposition; `\|wbar\|^2 = Ehat^2/\|W\|^2 - 1` |
| Vaidya indicatrix is the ellipse `(V^r+(E^2-f))^2 + E^2 (V^phi)^2 = E^2 w` | Lemma I.A, Remark | eliminate `u^v` from `g(u,u)=-1` and `-u_v=Ehat` |
| That ellipse is exactly `eq:vaidya-ind` | §4 | paper's parametrization substituted into the relation |
| Support function of the ellipse is the branch Hamiltonian `eq:Hv` | Theorem I.1(ii) | `max_theta (p_r r' + J phi') - 1 - H_v` |
| `det M = -32 (Ehat^2 - 1) disc_r S` | Lemma I.4, `eq:det-disc` | symbolic 11×11 determinant (Mathematica); sample points (SymPy) |
| Kernel at the degeneration locus is one-dimensional, spanned by `A = (r-r_d) Q_4`, `C = A' - Q_4` | Lemma I.4, proof | `rank M = 10`; `2A'S - AS'` divisible by `S`, `deg C = 4` |
| `J_deg` matches the closed form of Appendix B | Appendix B | `J_deg - 5 sqrt(3011/3072 + 581 sqrt(249)/9216) = 0` |
| `S` is jointly homogeneous of degree six in `(r,m,J)` | Lemma I.6, proof | Euler identity and explicit rescaling |
| `N_tot(r_d) = 0` while `N(r_d)`, `N_J(r_d)` are separately nonzero | Lemma I.6, `eq:track-cancel` | evaluation at the exact double root |
| Randers–Euler `(J d_J + p_r d_pr) H = H + 1` | `eq:finsler-euler` | direct |
| Self-similarity `r H_r + m H_m + J H_J = 0` — **identically, not only on shell** | `eq:selfsimilar` | direct and by explicit rescaling |
| `Theta H_v` equals `eq:ThetaHv` | Appendix C | direct |
| `d(r p_r)/dlambda - (1 + Theta H) = H`, hence `S_D = [r p_r] - lambda` on shell | `eq:SD-vaidya` | direct |

### Appendix B (`verify_appB_residues.wls`)

| Claim | Manuscript | Check |
|---|---|---|
| `r - r_d = s_d eps (1 + a_1 eps + a_2 eps^2 + ...)`, `a_1 = Q_4'(r_d)/(4 s_d)`, `a_2 = Q_4''(r_d)/12` | Lemma I.E(i) | order matching in `(dr/dz)^2 = Q_4(r)` |
| Laurent coefficients `b_3, b_2, b_1` at the triple pole `z_d` | Lemma I.E(ii) | series product, generic `Q_4` and `F` |
| `r(z)` is even about a two-torsion point: `k_3 = k_5 = 0`, `k_4 != 0` | Lemma I.E(iii) | order matching at a simple root of `Q_4` |
| `b_2^{e_i} = 4 N_m(e)/((e-r_d)^3 Q_4'(e)^2)` and no residue there | Lemma I.E(iii) | double-pole coefficient |

| `r(z) = -1/(sqrt(a_4)(z-z_inf)) + B + ...` with `B = -a_3/(4 a_4)` | Lemma I.F | order matching at the pole |
| `e_2^{z_inf} = 1/a_4` and `e_1^{z_inf} = -(2B + r_d - 2m)/sqrt(a_4)` | Lemma I.F | double-pole coefficient and residue of the clock integrand |
| symmetrised by parts `Int A'B = AB/2 + Int(A'B - AB')/2` | Lemma I.F | differentiate both sides |

Nothing Vaidya-specific enters these: they hold for any quartic `Q_4` and any
numerator `F`, which is why they are stated as a lemma rather than as a table.

`B` was used but never defined in the manuscript, and two coefficients were
written with `m` set to 1 (`-2` where `-2m` belongs). Both are fixed in the
current text; the scripts carry `m` symbolically so the omission cannot recur.

### Perlick recovery (`verify_perlick_recovery.wls`)

Section 2.3, `eq:perlick-randers`. Entirely symbolic, generic stationary metric:

| Claim | Check |
|---|---|
| `F_rail = beta.nu + sqrt(a.nu.nu)` solves the rail | substitute into the quadratic; rational part and the coefficient of the radical vanish separately |
| it is the future-directed branch | `f F_rail - g_ta nu^a = f sqrt(...) > 0`, so the unsquared rail selects it |
| `beta_a = omega_a` in threading form | direct |
| `g_ta g_tb + f g_ab = f h_ab` | direct |
| `a_ab = Ehat^2 h_ab/(f(Ehat^2-f))` — Perlick's fixed-energy optical metric | direct |
| `Ehat -> infinity` gives `h_ab/f` — Fermat-Randers | limit |

Substituting into the quadratic avoids comparing nested radicals, which no
simplifier resolves without positivity assumptions; the identity is then exact
rather than numerical.

## Notes

- The genus-degeneration double root is **negative**, `r_d = -3.3637111847...`
  at `m = 1`, `Ehat = 7/5`. This is not an anomaly: the manuscript states it
  ("the double root sits at `r_d < 0`", off the physical arc). The scripts use a
  plain symbol for `r` so the root is not silently filtered out — declaring `r`
  positive makes the locus appear empty.
- `det M` vanishes also at `Ehat = 1`, where `D_E = (Ehat^2-1) r + 2m` stops
  contributing a root and `deg S` drops from six: the reduction is genuinely
  singular there, which is why Lemma I.4 carries the hypothesis `Ehat != 1`.
- In SymPy the support-function identity needs `powsimp(..., force=True)` to
  combine `sqrt(A) sqrt(B)` into `sqrt(AB)`. On the regular domain
  (`r > 2m > 0`, `Ehat > 1`) both radicands are positive, so the rewriting is
  valid; Mathematica performs it directly from the assumptions.

## Not covered

These scripts verify the formal core only. They do **not** check the closed
forms of §4.4 (`eq:vaidya-full`, the `W_3k` blocks, the horizon dilogarithm),
the genus-two machinery of Appendix B beyond Lemma I.6, or §3 (FLRW). The
numerical validation of those lives in the reproducibility package
(`VaidyaMetric/`, `reproduce_reductions.py`, `run_regression.py`).
