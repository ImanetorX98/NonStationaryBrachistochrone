import ViaB.ConnectionFamilyRegularity
import Mathlib.Analysis.Calculus.FDeriv.Symmetric

open Filter
open scoped Topology

namespace ViaB

/-- A directional partial of a family on (launch parameter, affine time). -/
noncomputable def familyPartial (d : ℝ × ℝ) (F : ℝ × ℝ → ℝ)
    (p : ℝ × ℝ) : ℝ := fderiv ℝ F p d

noncomputable def launchPartial := familyPartial (1, 0)
noncomputable def affinePartial := familyPartial (0, 1)

/-- Regularity loss by one derivative, with no smooth-flow premise hidden. -/
theorem familyPartial_contDiffAt {n : WithTop ℕ∞} {F : ℝ × ℝ → ℝ}
    {p : ℝ × ℝ} (hF : ContDiffAt ℝ (n + 1) F p) (d : ℝ × ℝ) :
    ContDiffAt ℝ n (familyPartial d F) p := by
  exact (hF.fderiv_right le_rfl).clm_apply contDiffAt_const

/-- The affine coordinate slice has the actual derivative specified by its partial. -/
theorem affine_slice_hasDerivAt {F : ℝ × ℝ → ℝ} {a s : ℝ}
    (hF : DifferentiableAt ℝ F (a, s)) :
    HasDerivAt (fun t => F (a, t)) (affinePartial F (a, s)) s := by
  exact hF.hasFDerivAt.comp_hasDerivAt s
    ((hasDerivAt_const s a).prodMk (hasDerivAt_id s))

/-- The launch coordinate slice has the actual derivative specified by its partial. -/
theorem launch_slice_hasDerivAt {F : ℝ × ℝ → ℝ} {a s : ℝ}
    (hF : DifferentiableAt ℝ F (a, s)) :
    HasDerivAt (fun e => F (e, s)) (launchPartial F (a, s)) a := by
  exact hF.hasFDerivAt.comp_hasDerivAt a
    ((hasDerivAt_id a).prodMk (hasDerivAt_const a s))

/-- Derivative of a directional partial comes from the actual second Frechet derivative. -/
theorem familyPartial_hasFDerivAt {F : ℝ × ℝ → ℝ} {p : ℝ × ℝ}
    (hF : ContDiffAt ℝ 2 F p) (d : ℝ × ℝ) :
    HasFDerivAt (familyPartial d F) ((fderiv ℝ (fderiv ℝ F) p).flip d) p := by
  have hd := ((hF.fderiv_right (m := 1) (by norm_num)).differentiableAt le_rfl).hasFDerivAt
  simpa [familyPartial] using hd.clm_apply (hasFDerivAt_const d p)

/-- The second directional partials commute for C2 families. -/
theorem familyPartial_commute {F : ℝ × ℝ → ℝ} {p : ℝ × ℝ}
    (hF : ContDiffAt ℝ 2 F p) (d e : ℝ × ℝ) :
    familyPartial d (familyPartial e F) p = familyPartial e (familyPartial d F) p := by
  change fderiv ℝ (familyPartial e F) p d = fderiv ℝ (familyPartial d F) p e
  rw [(familyPartial_hasFDerivAt hF e).fderiv, (familyPartial_hasFDerivAt hF d).fderiv]
  exact (hF.isSymmSndFDerivAt (by simp)).eq d e

/-- Actual derivative of velocity with respect to launch equals the affine
 derivative of the actual launch-variation field. -/
theorem velocity_launch_variation_hasDerivAt {F : ℝ × ℝ → ℝ} {a s : ℝ}
    (hF : ContDiffAt ℝ 2 F (a, s)) :
    HasDerivAt (fun e => affinePartial F (e, s))
      (affinePartial (launchPartial F) (a, s)) a := by
  have h := launch_slice_hasDerivAt
    ((familyPartial_contDiffAt (n := 1) hF (0, 1)).differentiableAt le_rfl)
  have hc := familyPartial_commute hF (1, 0) (0, 1)
  simpa only [launchPartial, affinePartial, hc] using h

/-- Actual derivative of acceleration with respect to launch equals the
second affine partial of the launch-variation field for C3 families. -/
theorem acceleration_launch_variation_hasDerivAt {F : ℝ × ℝ → ℝ} {a s : ℝ}
    (hF : ContDiffAt ℝ 3 F (a, s)) :
    HasDerivAt (fun e => affinePartial (affinePartial F) (e, s))
      (affinePartial (affinePartial (launchPartial F)) (a, s)) a := by
  have hG : ContDiffAt ℝ 2 (affinePartial F) (a, s) :=
    familyPartial_contDiffAt (n := 2) hF (0, 1)
  have h := velocity_launch_variation_hasDerivAt hG
  have heq : launchPartial (affinePartial F) =ᶠ[𝓝 (a, s)]
      affinePartial (launchPartial F) := by
    filter_upwards [hF.eventually (by simp)] with p hp
    exact familyPartial_commute (hp.of_le (by norm_num)) (1, 0) (0, 1)
  have hdEq : fderiv ℝ (launchPartial (affinePartial F)) (a, s) =
      fderiv ℝ (affinePartial (launchPartial F)) (a, s) := heq.fderiv_eq
  have hvalue := congrArg (fun f : (ℝ × ℝ) →L[ℝ] ℝ => f (0, 1)) hdEq
  change affinePartial (launchPartial (affinePartial F)) (a, s) =
    affinePartial (affinePartial (launchPartial F)) (a, s) at hvalue
  rw [hvalue] at h
  exact h

/-- The launch field has an actual affine derivative, without identifying it
by name alone. -/
theorem launch_field_affine_hasDerivAt {F : ℝ × ℝ → ℝ} {a s : ℝ}
    (hF : ContDiffAt ℝ 2 F (a, s)) :
    HasDerivAt (fun t => launchPartial F (a, t))
      (affinePartial (launchPartial F) (a, s)) s := by
  exact affine_slice_hasDerivAt
    ((familyPartial_contDiffAt (n := 1) hF (1, 0)).differentiableAt le_rfl)

/-- The second derivative is the derivative of the actual first derivative
of the launch field. Equality of germs is used before differentiating. -/
theorem launch_field_affine_second_hasDerivAt {F : ℝ × ℝ → ℝ} {a s : ℝ}
    (hF : ContDiffAt ℝ 3 F (a, s)) :
    HasDerivAt (deriv (fun t => launchPartial F (a, t)))
      (affinePartial (affinePartial (launchPartial F)) (a, s)) s := by
  have hJ : ContDiffAt ℝ 2 (launchPartial F) (a, s) :=
    familyPartial_contDiffAt (n := 2) hF (1, 0)
  have hJprime : ContDiffAt ℝ 1 (affinePartial (launchPartial F)) (a, s) :=
    familyPartial_contDiffAt (n := 1) hJ (0, 1)
  have hd := affine_slice_hasDerivAt (hJprime.differentiableAt le_rfl)
  have hc : ContinuousAt (fun t : ℝ => (a, t)) s :=
    continuousAt_const.prodMk continuousAt_id
  have heq : deriv (fun t => launchPartial F (a, t)) =ᶠ[𝓝 s]
      (fun t => affinePartial (launchPartial F) (a, t)) := by
    filter_upwards [hc.eventually (hF.eventually (by simp))] with t ht
    exact (launch_field_affine_hasDerivAt (ht.of_le (by norm_num))).deriv
  exact hd.congr_of_eventuallyEq heq

/-- Coordinate Jacobi identity for a C3 geodesic family. Mixed derivatives
are conclusions here, not independent premises. The family itself is supplied. -/
theorem C3_geodesic_family_coordinate_jacobi {a s : ℝ}
    {X : CoordinateIndex → ℝ × ℝ → ℝ}
    {Γ : ℝ → CoordinateIndex → CoordinateIndex → CoordinateIndex → ℝ}
    {DG : CoordinateIndex → CoordinateIndex → CoordinateIndex → ℝ}
    (hX : ∀ i, ContDiffAt ℝ 3 (X i) (a, s))
    (hΓ : ∀ i j k, HasDerivAt (fun e => Γ e i j k) (DG i j k) a)
    (hs : ∀ i j k, Γ a i j k = Γ a i k j)
    (hgeo : ∀ᶠ e in 𝓝 a, ∀ i,
      affinePartial (affinePartial (X i)) (e, s) =
        -connectionQuadratic (Γ e) (fun j => affinePartial (X j) (e, s)) i)
    (i : CoordinateIndex) :
    deriv (deriv (fun t => launchPartial (X i) (a, t))) s +
      2 * (∑ j, ∑ k, Γ a i j k * affinePartial (X j) (a, s) *
        deriv (fun t => launchPartial (X k) (a, t)) s) +
      (∑ j, ∑ k, DG i j k * affinePartial (X j) (a, s) *
        affinePartial (X k) (a, s)) = 0 := by
  have hc := coordinate_jacobi_of_commuted_variation hΓ
    (fun j => velocity_launch_variation_hasDerivAt ((hX j).of_le (by norm_num)))
    (fun j => acceleration_launch_variation_hasDerivAt (hX j)) hgeo hs i
  have hJ : ∀ j, deriv (fun t => launchPartial (X j) (a, t)) s =
      affinePartial (launchPartial (X j)) (a, s) :=
    fun j => (launch_field_affine_hasDerivAt ((hX j).of_le (by norm_num))).deriv
  rw [(launch_field_affine_second_hasDerivAt (hX i)).deriv]
  simp_rw [hJ]
  exact hc

/-- Launching all neighboring curves from one event forces the actual
variation field to vanish there. This is deduced from fixed-position data. -/
theorem fixed_launch_event_variation_zero {F : ℝ × ℝ → ℝ} {a s x0 : ℝ}
    (hF : DifferentiableAt ℝ F (a, s))
    (hfixed : ∀ᶠ e in 𝓝 a, F (e, s) = x0) :
    launchPartial F (a, s) = 0 := by
  have hd : HasDerivAt (fun e => F (e, s)) 0 a :=
    (hasDerivAt_const a x0).congr_of_eventuallyEq hfixed
  exact (launch_slice_hasDerivAt hF).unique hd

/-- Actual second derivative of an affine coordinate slice under C2 regularity. -/
theorem affine_slice_second_hasDerivAt {F : ℝ × ℝ → ℝ} {a s : ℝ}
    (hF : ContDiffAt ℝ 2 F (a, s)) :
    HasDerivAt (deriv (fun t => F (a, t)))
      (affinePartial (affinePartial F) (a, s)) s := by
  have hU : ContDiffAt ℝ 1 (affinePartial F) (a, s) :=
    familyPartial_contDiffAt (n := 1) hF (0, 1)
  have hd := affine_slice_hasDerivAt (hU.differentiableAt le_rfl)
  have hc : ContinuousAt (fun t : ℝ => (a, t)) s :=
    continuousAt_const.prodMk continuousAt_id
  have heq : deriv (fun t => F (a, t)) =ᶠ[𝓝 s]
      (fun t => affinePartial F (a, t)) := by
    filter_upwards [hc.eventually (hF.eventually (by simp))] with t ht
    exact (affine_slice_hasDerivAt (ht.differentiableAt (by norm_num))).deriv
  exact hd.congr_of_eventuallyEq heq

end ViaB
