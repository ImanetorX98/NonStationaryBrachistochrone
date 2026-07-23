# TODO — closed-form of the exact first-order adiabatic correction

Status of the non-autonomous (slow-`A`) first-order correction to the optical
brachistochrone shape `phi(r)`. Context: referee 4.6 (off-shell term) / 4.7 / 4.8.

## DONE (consolidated)
- **First order CLOSED to O(eps^2)** by exact canonical perturbation theory of the
  non-autonomous optical Hamiltonian flow `H2` (E_eff=Ehat/A, J_eff=J/A, A=exp(eps*lambda)):
  ```
  delta p_r/eps = [ lambda*Euler_H - S ] / H_pr ,   Euler_H = E H_E + J H_J ,  S = int_0^lambda Euler_H dlambda
  delta phi/eps = int [ G_pr*(delta p_r/eps) - lambda*(E G_E + J G_J) ] dr ,    G = H_J/H_pr = dphi/dr
  ```
  Verified: coeff -0.3298 vs true -0.328; robust slope ~1.90 (=> O(eps^2)).
  Script: `ThakurtaMetric/adiabatic_first_order_exact.py`. Steps 1-3 numerically checked.
- **Decomposition** (both LARGE, ~half each; they nearly cancel):
  - `partA = int lambda*[ (G_pr/H_pr)Euler_H - Euler_G ] dr`  ~ -0.671  (single integral, ~ -full Euler)
  - `partB = -int (G_pr/H_pr) * S dr`  ~ +0.341  (OFF-SHELL; S nested => length-2 iterated integral)
  - The paper's `1/2`-Euler is the APPROXIMATION `partB ~= -1/2 partA` (accurate to ~2%);
    the true off-shell term is `partB`, not a 2% afterthought.
- **Structure identified** (via IBP, no heavy simplify): `partA = [Lambda*lambda] - int Lambda*lambda' dr`
  with `lambda' = 1/H_pr` = the clock differential (paper Table: `rho_t/sqrt(R6)`), so
  `partA` lives in the SAME genus-2 2nd/3rd-kind class as the paper's on-shell `C + psi`.
- **Honest validation figure**: `ThakurtaMetric/fig_phi_validation_corrected.py` adds the true
  dynamic-test panel (~2% physical accuracy) vs the algebraic IBP floor (~1e-6). Referee 4.6/4.14.

## TODO (task a) — explicit closed form with symbolic coefficients
1. **Frame conversion** `partA`, `partB` from the Hamiltonian frame (`G_pr, H_pr, Euler_H,
   Euler_G`, with the nested-sqrt `H2`) to the ALGEBRAIC frame (`F = K/sqrt(R)`, `R` sextic,
   E/J derivatives). Relations: `G = F` on-shell; `H_pr = 1/lambda'` (clock density);
   the tricky piece is `G_pr` (off-shell response). NB: a direct `sympy.simplify` of `Pi`
   HANGS (nested sqrt) -- substitute the on-shell `p_r` (via `sqrt(R)`) EARLY, then reduce
   rational-in-`r` times powers of `sqrt(R)`. Verify numerically against `adiabatic_first_order_exact.py`.
2. **Genus-2 reduction** of `partA` = int (algebraic) x (clock differential): run the existing
   pipeline (`KerrSessionScripts/`) -- principal parts -> `c_k, Q_kj`, second/third-kind
   differentials -> Weierstrass `zeta/wp` (genus-1) + genus-2 `theta` (the `W_kj` basis).
   Expect the same special functions as the paper's `C + psi`.
3. **partB (off-shell) closed form**: `S = int Euler_H dlambda` is a length-2 iterated
   integral; `partB = -int (G_pr/H_pr) S dr` is therefore a depth-2 iterated integral on the
   genus-2 curve -> a genus-2 DILOGARITHM (the `W_kj` / Kronecker-Eisenstein structure). This
   is exactly the object the paper's genus-2 machinery handles -- assemble in that basis with
   ALGEBRAIC (symbolic, not fitted) coefficients.
4. **Cross-check** the assembled closed form numerically (residual ~machine precision) AND
   independently in Mathematica, per the repo convention; then it can replace the `1/2`-Euler
   in the paper as the exact first-order term (closing 4.6 fully).
5. NB: this overlaps the genus-2 rigor flagged by **referee 4.12** (base points, single-valued
   completion, irreducibility) -- keep the algebraic reduction (rigorous) separate from any
   special-function-basis / irreducibility claims (conjectural).
