import ViaB.MixedGeodesicDerivatives

open Filter
open scoped Topology

namespace ViaB

/-- Actual coordinate Jacobi equation for a supplied jointly C3 Fermat
geodesic family with C2 physical mass. Neither mixed derivatives nor a
curvature coefficient are assumed. Existence of this family is not claimed. -/
theorem fermat_C3_family_actual_coordinate_jacobi {q a s : ℝ}
    {m : ℝ → ℝ} {X : CoordinateIndex → ℝ × ℝ → ℝ}
    (hq : q ≠ 0) (hX : ∀ i, ContDiffAt ℝ 3 (X i) (a, s))
    (hm : ContDiffAt ℝ 2 m (X 0 (a, s)))
    (hr : X 1 (a, s) ≠ 0)
    (hw : vaidyaW q (m (X 0 (a, s))) (X 1 (a, s)) ≠ 0)
    (hgeo : ∀ᶠ e in 𝓝 a, ∀ i,
      affinePartial (affinePartial (X i)) (e, s) =
        -fermatChristoffelContraction q (m (X 0 (e, s)))
          (deriv m (X 0 (e, s))) (X 1 (e, s))
          (affinePartial (X 0) (e, s)) (affinePartial (X 1) (e, s))
          (affinePartial (X 2) (e, s)) i)
    (i : CoordinateIndex) :
    deriv (deriv (fun t => launchPartial (X i) (a, t))) s +
      2 * (∑ j, ∑ k,
        fermatChristoffel q (m (X 0 (a, s))) (deriv m (X 0 (a, s))) (X 1 (a, s)) i j k *
          affinePartial (X j) (a, s) * deriv (fun t => launchPartial (X k) (a, t)) s) +
      (∑ j, ∑ k,
        deriv (fun e => fermatChristoffel q (m (X 0 (e, s)))
          (deriv m (X 0 (e, s))) (X 1 (e, s)) i j k) a *
          affinePartial (X j) (a, s) * affinePartial (X k) (a, s)) = 0 := by
  have hd : ∀ j, DifferentiableAt ℝ (fun e => X j (e, s)) a :=
    fun j => (launch_slice_hasDerivAt ((hX j).differentiableAt (by norm_num))).differentiableAt
  apply C3_geodesic_family_coordinate_jacobi hX
    (fun i j k => (fermatChristoffel_sampled_C2_mass_differentiable hq hm
      (hd 0) (hd 1) hr hw i j k).hasDerivAt)
  · exact fermatChristoffel_torsion_free q (m (X 0 (a, s)))
      (deriv m (X 0 (a, s))) (X 1 (a, s))
  · filter_upwards [hgeo] with e he
    intro j
    simpa only [fermat_connectionQuadratic] using he j

/-- Apply the C3 variation theorem to genuine coordinate null-geodesic
families, eliminating a separate assumption about the symbolic acceleration.
The jointly C3 family still has to be supplied or constructed. -/
theorem fermat_C3_null_geodesic_family_coordinate_jacobi {q a s : ℝ}
    {m : ℝ → ℝ} {X : CoordinateIndex → ℝ × ℝ → ℝ}
    (hq : q ≠ 0) (hX : ∀ i, ContDiffAt ℝ 3 (X i) (a, s))
    (hm : ContDiffAt ℝ 2 m (X 0 (a, s)))
    (hr : X 1 (a, s) ≠ 0)
    (hw : vaidyaW q (m (X 0 (a, s))) (X 1 (a, s)) ≠ 0)
    (hgeo : ∀ᶠ e in 𝓝 a, CoordinateNullGeodesicAt q m (deriv m)
      (fun t => X 0 (e, t)) (fun t => X 1 (e, t)) (fun t => X 2 (e, t)) s)
    (i : CoordinateIndex) :
    deriv (deriv (fun t => launchPartial (X i) (a, t))) s +
      2 * (∑ j, ∑ k,
        fermatChristoffel q (m (X 0 (a, s))) (deriv m (X 0 (a, s))) (X 1 (a, s)) i j k *
          affinePartial (X j) (a, s) * deriv (fun t => launchPartial (X k) (a, t)) s) +
      (∑ j, ∑ k,
        deriv (fun e => fermatChristoffel q (m (X 0 (e, s)))
          (deriv m (X 0 (e, s))) (X 1 (e, s)) i j k) a *
          affinePartial (X j) (a, s) * affinePartial (X k) (a, s)) = 0 := by
  apply fermat_C3_family_actual_coordinate_jacobi hq hX hm hr hw
  have hc : ContinuousAt (fun e : ℝ => (e, s)) a :=
    continuousAt_id.prodMk continuousAt_const
  have hnear : ∀ᶠ e in 𝓝 a, ∀ j, ContDiffAt ℝ 3 (X j) (e, s) := by
    exact Filter.eventually_all.mpr (fun j => hc.eventually ((hX j).eventually (by simp)))
  filter_upwards [hgeo, hnear] with e he hregular
  intro j
  have hvel : ∀ k, deriv (fun t => X k (e, t)) s = affinePartial (X k) (e, s) :=
    fun k => (affine_slice_hasDerivAt ((hregular k).differentiableAt (by norm_num))).deriv
  have hacc : HasDerivAt (deriv (fun t => X j (e, t)))
      (-fermatChristoffelContraction q (m (X 0 (e, s)))
        (deriv m (X 0 (e, s))) (X 1 (e, s))
        (deriv (fun t => X 0 (e, t)) s) (deriv (fun t => X 1 (e, t)) s)
        (deriv (fun t => X 2 (e, t)) s) j) s := by
    fin_cases j
    · simpa using he.equation 0
    · simpa using he.equation 1
    · simpa using he.equation 2
  have heq := (affine_slice_second_hasDerivAt ((hregular j).of_le (by norm_num))).unique hacc
  simpa only [hvel] using heq

end ViaB
