import ViaB.DirectionGeodesicGerm

open Filter
open scoped Topology

namespace ViaB

/-- Connection contraction in the fixed three-coordinate chart. -/
noncomputable def connectionQuadratic
    (Γ : CoordinateIndex → CoordinateIndex → CoordinateIndex → ℝ)
    (u : CoordinateIndex → ℝ) (i : CoordinateIndex) : ℝ :=
  ∑ j, ∑ k, Γ i j k * u j * u k

/-- First variation before using torsion freedom. DG is the derivative of
connection coefficients along the parameter family, not an assumed curvature. -/
noncomputable def connectionVariation
    (Γ DG : CoordinateIndex → CoordinateIndex → CoordinateIndex → ℝ)
    (u U : CoordinateIndex → ℝ) (i : CoordinateIndex) : ℝ :=
  ∑ j, ∑ k, (DG i j k * u j * u k +
    Γ i j k * U j * u k + Γ i j k * u j * U k)

/-- Differentiate the connection contraction with respect to a family parameter. -/
theorem connectionQuadratic_hasDerivAt {a : ℝ}
    {Γ : ℝ → CoordinateIndex → CoordinateIndex → CoordinateIndex → ℝ}
    {u : ℝ → CoordinateIndex → ℝ}
    {DG : CoordinateIndex → CoordinateIndex → CoordinateIndex → ℝ}
    {U : CoordinateIndex → ℝ}
    (hΓ : ∀ i j k, HasDerivAt (fun e => Γ e i j k) (DG i j k) a)
    (hu : ∀ j, HasDerivAt (fun e => u e j) (U j) a) (i : CoordinateIndex) :
    HasDerivAt (fun e => connectionQuadratic (Γ e) (u e) i)
      (connectionVariation (Γ a) DG (u a) U i) a := by
  unfold connectionQuadratic connectionVariation
  apply HasDerivAt.fun_sum
  intro j _
  apply HasDerivAt.fun_sum
  intro k _
  convert ((hΓ i j k).mul (hu j)).mul (hu k) using 1
  dsimp
  ring

/-- Torsion freedom combines the two velocity-variation terms. -/
theorem connectionVariation_torsion_free
    (Γ DG : CoordinateIndex → CoordinateIndex → CoordinateIndex → ℝ)
    (u U : CoordinateIndex → ℝ)
    (hs : ∀ i j k, Γ i j k = Γ i k j) (i : CoordinateIndex) :
    connectionVariation Γ DG u U i =
      (∑ j, ∑ k, DG i j k * u j * u k) +
        2 * (∑ j, ∑ k, Γ i j k * u j * U k) := by
  have swap : (∑ j, ∑ k, Γ i j k * U j * u k) =
      ∑ j, ∑ k, Γ i j k * u j * U k := by
    rw [Finset.sum_comm]
    apply Finset.sum_congr rfl
    intro j _
    apply Finset.sum_congr rfl
    intro k _
    rw [hs i k j]
    ring
  unfold connectionVariation
  simp only [Finset.sum_add_distrib]
  rw [swap]
  ring

/-- A family satisfying the geodesic equation on a neighborhood has the
linearized acceleration given by the differentiated connection contraction.
No mixed-derivative commutation is claimed here. -/
theorem geodesic_family_acceleration_variation {a : ℝ}
    {Γ : ℝ → CoordinateIndex → CoordinateIndex → CoordinateIndex → ℝ}
    {u b : ℝ → CoordinateIndex → ℝ}
    {DG : CoordinateIndex → CoordinateIndex → CoordinateIndex → ℝ}
    {U : CoordinateIndex → ℝ}
    (hΓ : ∀ i j k, HasDerivAt (fun e => Γ e i j k) (DG i j k) a)
    (hu : ∀ j, HasDerivAt (fun e => u e j) (U j) a)
    (hgeo : ∀ᶠ e in 𝓝 a, ∀ i, b e i = -connectionQuadratic (Γ e) (u e) i)
    (i : CoordinateIndex) :
    HasDerivAt (fun e => b e i)
      (-connectionVariation (Γ a) DG (u a) U i) a := by
  apply ((connectionQuadratic_hasDerivAt hΓ hu i).neg).congr_of_eventuallyEq
  filter_upwards [hgeo] with e he
  exact he i

/-- Conditional coordinate Jacobi equation: the derivatives of velocity and
acceleration with respect to the family parameter must have been identified
with J' and J''. Those identifications remain separate geometric obligations. -/
theorem coordinate_jacobi_of_commuted_variation {a : ℝ}
    {Γ : ℝ → CoordinateIndex → CoordinateIndex → CoordinateIndex → ℝ}
    {u b : ℝ → CoordinateIndex → ℝ}
    {DG : CoordinateIndex → CoordinateIndex → CoordinateIndex → ℝ}
    {Jprime Jsecond : CoordinateIndex → ℝ}
    (hΓ : ∀ i j k, HasDerivAt (fun e => Γ e i j k) (DG i j k) a)
    (hu : ∀ j, HasDerivAt (fun e => u e j) (Jprime j) a)
    (hb : ∀ i, HasDerivAt (fun e => b e i) (Jsecond i) a)
    (hgeo : ∀ᶠ e in 𝓝 a, ∀ i, b e i = -connectionQuadratic (Γ e) (u e) i)
    (hs : ∀ i j k, Γ a i j k = Γ a i k j) (i : CoordinateIndex) :
    Jsecond i + 2 * (∑ j, ∑ k, Γ a i j k * u a j * Jprime k) +
      (∑ j, ∑ k, DG i j k * u a j * u a k) = 0 := by
  have he := (hb i).unique (geodesic_family_acceleration_variation hΓ hu hgeo i)
  rw [connectionVariation_torsion_free (Γ a) DG (u a) Jprime hs i] at he
  linarith

/-- The already verified Fermat contraction is the quadratic contraction used
by the variation theorem; its coefficients satisfy the required symmetry. -/
theorem fermat_connectionQuadratic (q m μ r : ℝ)
    (u : CoordinateIndex → ℝ) (i : CoordinateIndex) :
    connectionQuadratic (fermatChristoffel q m μ r) u i =
      fermatChristoffelContraction q m μ r (u 0) (u 1) (u 2) i := by
  have hv : (![u 0, u 1, u 2] : CoordinateIndex → ℝ) = u := by
    funext j
    fin_cases j <;> rfl
  simp only [connectionQuadratic, fermatChristoffelContraction, hv]

/-- Angular translations preserve the actual affine null geodesic equation.
This provides a concrete variation without an existence theorem for families. -/
theorem coordinate_null_geodesic_angle_shift {q s c : ℝ}
    {m μ v R Φ : ℝ → ℝ} (h : CoordinateNullGeodesicAt q m μ v R Φ s) :
    CoordinateNullGeodesicAt q m μ v R (fun t => Φ t + c) s := by
  constructor
  · simpa only [deriv_add_const] using h.equation
  · simpa only [deriv_add_const] using h.null
  · exact h.future

/-- The variation vector of the angular-translation family is the angular
coordinate vector. It does not vanish at the launch and is not the screen
Jacobi field used in the scalar Sturm comparison. -/
theorem angle_shift_family_position_derivative (v R Φ : ℝ → ℝ) (s a : ℝ) :
    HasDerivAt (fun e => (v s, R s, Φ s + e)) (0, 0, 1) a := by
  simpa using (hasDerivAt_const a (v s)).prodMk
    ((hasDerivAt_const a (R s)).prodMk ((hasDerivAt_const a (Φ s)).add (hasDerivAt_id a)))

end ViaB
