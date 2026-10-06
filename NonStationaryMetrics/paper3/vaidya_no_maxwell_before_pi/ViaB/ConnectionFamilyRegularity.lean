import ViaB.GeodesicVariation

open Filter
open scoped Topology

namespace ViaB

/-- Differentiability of the inverse metric along a regular parameter family. -/
theorem fermatInverseMetric_family_differentiable {q a : ℝ} {m r : ℝ → ℝ}
    (hm : DifferentiableAt ℝ m a) (hrd : DifferentiableAt ℝ r a)
    (hr : r a ≠ 0) (hw : vaidyaW q (m a) (r a) ≠ 0) (i j : CoordinateIndex) :
    DifferentiableAt ℝ (fun e => fermatInverseMetric q (m e) (r e) i j) a := by
  have hwd : DifferentiableAt ℝ (fun e => vaidyaW q (m e) (r e)) a := by
    unfold vaidyaW
    fun_prop (disch := assumption)
  fin_cases i <;> fin_cases j
  · exact (differentiableAt_const (-1)).div hwd hw
  · exact differentiableAt_const 1
  · exact differentiableAt_const 0
  · exact differentiableAt_const 1
  · exact (differentiableAt_const q).sub hwd
  · exact differentiableAt_const 0
  · exact differentiableAt_const 0
  · exact differentiableAt_const 0
  · exact (differentiableAt_const 1).div (hrd.pow 2) (pow_ne_zero 2 hr)

/-- Differentiating metric partials requires differentiability of the mass
rate as well as the mass. No C1-to-C2 regularity upgrade is assumed. -/
theorem fermatMetricPartial_family_differentiable {q a : ℝ} {m μ r : ℝ → ℝ}
    (hq : q ≠ 0) (hm : DifferentiableAt ℝ m a)
    (hμ : DifferentiableAt ℝ μ a) (hrd : DifferentiableAt ℝ r a)
    (hr : r a ≠ 0) (k i j : CoordinateIndex) :
    DifferentiableAt ℝ (fun e => fermatMetricPartial q (m e) (μ e) (r e) k i j) a := by
  have hwd : DifferentiableAt ℝ (fun e => vaidyaW q (m e) (r e)) a := by
    unfold vaidyaW
    fun_prop (disch := assumption)
  fin_cases k <;> fin_cases i <;> fin_cases j <;> dsimp [fermatMetricPartial]
  all_goals
    have hqr : q * r a ≠ 0 := mul_ne_zero hq hr
    have hr2 : r a ^ 2 ≠ 0 := pow_ne_zero 2 hr
    have hqr2 : q * r a ^ 2 ≠ 0 := mul_ne_zero hq hr2
    fun_prop (disch := assumption)

/-- The actual Fermat coefficients are differentiable along any family with
regular radius and differentiable mass and mass rate. -/
theorem fermatChristoffel_family_differentiable {q a : ℝ} {m μ r : ℝ → ℝ}
    (hq : q ≠ 0) (hm : DifferentiableAt ℝ m a)
    (hμ : DifferentiableAt ℝ μ a) (hrd : DifferentiableAt ℝ r a)
    (hr : r a ≠ 0) (hw : vaidyaW q (m a) (r a) ≠ 0) (i j k : CoordinateIndex) :
    DifferentiableAt ℝ (fun e => fermatChristoffel q (m e) (μ e) (r e) i j k) a := by
  unfold fermatChristoffel
  apply DifferentiableAt.fun_sum
  intro h _
  apply (fermatInverseMetric_family_differentiable hm hrd hr hw i h).mul
  unfold coordinateKoszul
  exact (((fermatMetricPartial_family_differentiable hq hm hμ hrd hr j h k).add
    (fermatMetricPartial_family_differentiable hq hm hμ hrd hr k h j)).sub
    (fermatMetricPartial_family_differentiable hq hm hμ hrd hr h j k)).div_const 2

/-- C2 physical mass supplies differentiability of both sampled mass and
sampled mass rate along a differentiable time-coordinate family. -/
theorem fermatChristoffel_sampled_C2_mass_differentiable {q a : ℝ}
    {m v r : ℝ → ℝ} (hq : q ≠ 0) (hm : ContDiffAt ℝ 2 m (v a))
    (hv : DifferentiableAt ℝ v a) (hrd : DifferentiableAt ℝ r a)
    (hr : r a ≠ 0) (hw : vaidyaW q (m (v a)) (r a) ≠ 0)
    (i j k : CoordinateIndex) :
    DifferentiableAt ℝ (fun e => fermatChristoffel q (m (v e))
      (deriv m (v e)) (r e) i j k) a := by
  have hdm : ContDiffAt ℝ 1 (deriv m) (v a) := by
    have hf := (hm.fderiv_right (m := 1) (by norm_num)).clm_apply
      (contDiffAt_const (c := (1 : ℝ)))
    simpa [← deriv_fderiv] using hf
  exact fermatChristoffel_family_differentiable hq
    ((hm.differentiableAt (by norm_num)).comp a hv)
    ((hdm.differentiableAt le_rfl).comp a hv) hrd hr hw i j k

/-- Conditional Fermat Jacobi equation with the coefficient derivative
computed as the actual derivative of the connection along the family.
Velocity/acceleration parameter derivatives are still explicit premises. -/
theorem fermat_coordinate_jacobi_of_family {q a : ℝ} {m μ r : ℝ → ℝ}
    {u b : ℝ → CoordinateIndex → ℝ} {Jprime Jsecond : CoordinateIndex → ℝ}
    (hq : q ≠ 0) (hm : DifferentiableAt ℝ m a)
    (hμ : DifferentiableAt ℝ μ a) (hrd : DifferentiableAt ℝ r a)
    (hr : r a ≠ 0) (hw : vaidyaW q (m a) (r a) ≠ 0)
    (hu : ∀ j, HasDerivAt (fun e => u e j) (Jprime j) a)
    (hb : ∀ i, HasDerivAt (fun e => b e i) (Jsecond i) a)
    (hgeo : ∀ᶠ e in 𝓝 a, ∀ i, b e i =
      -fermatChristoffelContraction q (m e) (μ e) (r e) (u e 0) (u e 1) (u e 2) i)
    (i : CoordinateIndex) :
    Jsecond i + 2 * (∑ j, ∑ k,
      fermatChristoffel q (m a) (μ a) (r a) i j k * u a j * Jprime k) +
      (∑ j, ∑ k, deriv (fun e => fermatChristoffel q (m e) (μ e) (r e) i j k) a *
        u a j * u a k) = 0 := by
  apply coordinate_jacobi_of_commuted_variation
    (Γ := fun e => fermatChristoffel q (m e) (μ e) (r e))
    (fun i j k => (fermatChristoffel_family_differentiable hq hm hμ hrd hr hw i j k).hasDerivAt)
    hu hb
  · filter_upwards [hgeo] with e he
    intro i
    simpa only [fermat_connectionQuadratic] using he i
  · exact fermatChristoffel_torsion_free q (m a) (μ a) (r a)

end ViaB
