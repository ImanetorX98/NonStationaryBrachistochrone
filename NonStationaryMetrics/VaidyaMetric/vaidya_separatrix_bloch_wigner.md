# Bloch–Wigner closure of the horizon dilog D_k on the separatrix

On `|J|=Jc` the sextic `S` acquires a double root, the genus drops to 1 and the
hyperelliptic dilog `D_k = ∫U_k/(r-2m) dr` closes into the **elliptic dilogarithm
(Bloch–Wigner)**, a tabulated and reusable form.

## Verified data (M=1, E=1.4)
- **Jc = 7.02662374** (double root of the bracket `B=r²(r-2m)-J²DE`, via the discriminant).
- **r_d = -3.36371118** (double root; `√S = (r-r_d)√Q4`, `Q4` quartic).
- Elliptic curve **E: w² = Q4(r)**, `a4 = E²-1 = 0.96`, branch points
  (roots of Q4) = **{-2.08333, 0, 2, 8.72742}**.
  - `-2.08333 = -2m/(E²-1)` (root of DE); `2 = 2m` (**horizon**); `8.72742`
    physical turning point; `0`.
- **Elliptic modulus k² = 0.60671566**;  real period = 0.66912986,
  imaginary period = 0.60621385.
- **τ = 0.90597338 i** (RECTANGULAR lattice). Verified in TWO independent ways:
  Legendre `τ=iK'(k)/K(k)` **=** Sage `EllipticCurve.period_lattice()` to 15 digits.

## Verified reductions
1. **Elliptic U_k** (Stage A): `U_k = ∫r^k/[(r-r_d)√Q4]dr` against the direct genus-2
   `∫r^k/√S` at Jc → they agree to **1e-14…1e-17** (machine precision). The split
   `r^k/(r-r_d)=poly_{k-1}+r_d^k/(r-r_d)` ⇒ elliptic of the first/second kind + third kind
   (pole r_d).
2. **Abel map** `z(r)=∫dr/√Q4`: `dz/dr = 1/√Q4` verified (1e-7).
3. **Horizon r=2m = 2-TORSION point**: since `r=2m` is a branch point of Q4, its
   Abel image is at the imaginary half-period `z_h = ω_im/2` (verified,
   diff 7e-13). dlog letter at the 2-torsion point.

## Closed form — NUMERICAL MATCH (vaidya_ell_dilog_match.py)
`D_k = U_k ln(r-2m) - G_k`,  `G_k = ∫ln(r-2m) r^k/√S dr`.

### WEIGHT 1: U_0 in Weierstrass σ,ζ (third kind) — verified to 1e-14
σ,ζ toolkit from θ₁ (mpmath), self-consistent lattice (om1=0.66913, om2=0.60621 i,
τ=0.90597 i). Base point at the WEIERSTRASS POINT e4 (z=0) ⇒ involution z→−z ⇒ the two
preimages of r_d are ±z_d.

    U_0(r) = ρ[ ln σ(z−z_d) − ln σ(z+z_d) ] + C z + const
    ρ = 1/√Q4(r_d) = 0.061069   (derived RESIDUE, NOT fitted)
    C = 1/(e4−r_d) + 2ρ ζ(z_d) = 0.333272   (holomorphic part, derived)

Check of the derivative dU_0/dz = 1/(r−r_d) and of the antiderivative U_0: diff **1e-14** (machine precision).
z_d = 0.46104 (Abel image of r_d, through the real oval at infinity).

### WEIGHT 2: horizon dilog = elliptic dilogarithm
`G_0 = ∫ ln(r−2m) dU_0 = ∫ ln(r(z)−2m)[ρ(ζ(z−z_d)−ζ(z+z_d))+C] dz`.
`ln(r−2m)` has a double zero at the 2-torsion point z_h and double poles at the two infinities ⇒ G_0 =
a combination of `∫ ln σ(z−a) ζ(z−b) dz` = **elliptic dilogarithm**.
Closed-loop check `D_0 = U_0 ln(r−2m) − G_0` (with U_0 in σ,ζ): diff **1e-14**.

### Tabulated function (Zagier)
`D^E(z) = Σ_{n≥0} D(qⁿζ) − Σ_{n≥1} D(qⁿ/ζ)`, `q=e^{2πiτ}`, `ζ=e^{2πi z/2ω1}`,
D = Bloch–Wigner. Implemented, it converges; D^E(0.2+0.3i)=0.59621, D^E(z_d)=0.67494.
NB (correction): real z ⇒ ζ on the unit circle (a phase) ⇒ the qⁿζ are complex ⇒ D^E≠0:
the single-valued elliptic dilog IS the content (it does NOT vanish on a real configuration).

### What remains (the only open point)
The EXACT normalisation of the regulator that relates the physical (real,
indefinite) integral D_0 to a specific combination of D^E — this is the theory of elliptic polylogarithms
(real+imaginary Beilinson–Levin regulator). Structure + elliptic reduction: done to
machine precision. Tabulated function: in hand.

Scripts: `vaidya_separatrix_ell.py` (Stage A), `vaidya_bloch_wigner.py` (τ, Abel map,
2-torsion), `vaidya_ell_dilog_match.py` (σ,ζ closed form + D^E q-series); the Sage cross-check of τ
was a working script that is not distributed.
