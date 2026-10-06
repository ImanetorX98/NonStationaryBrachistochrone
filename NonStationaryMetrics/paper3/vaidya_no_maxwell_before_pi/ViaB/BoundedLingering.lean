import ViaB.RadialConvergence
import ViaB.UnboundedEscape

open Filter Topology

namespace ViaB

theorem derivative_limit_zero_on_bounded_positive_radius {r dr : ℝ → ℝ} {a ρ R L : ℝ}
    (hρ : 0 < ρ) (hr : ∀ t ≥ a, ρ ≤ r t ∧ r t ≤ R)
    (hder : ∀ t ≥ a, HasDerivAt r (dr t) t)
    (hlim : Tendsto dr atTop (𝓝 L)) : L = 0 := by
  by_contra hne
  rcases lt_or_gt_of_ne hne with hneg | hpos
  · have hev : ∀ᶠ t in atTop, dr t < L / 2 :=
      hlim.eventually (eventually_lt_nhds (by linarith))
    obtain ⟨T, hT⟩ := eventually_atTop.1 hev
    let A := max a T
    have haA : a ≤ A := le_max_left _ _
    have hTA : T ≤ A := le_max_right _ _
    exact positive_future_incompatible_with_negative_derivative (C := -L / 2) (by linarith)
      (fun t ht => hρ.trans_le (hr t (haA.trans ht)).1)
      (fun t ht => hder t (haA.trans ht)) (fun t ht => by have hh := hT t (hTA.trans ht); linarith)
  · have hev : ∀ᶠ t in atTop, L / 2 < dr t :=
      hlim.eventually (eventually_gt_nhds (by linarith))
    obtain ⟨T, hT⟩ := eventually_atTop.1 hev
    let A := max a T
    have haA : a ≤ A := le_max_left _ _
    have hTA : T ≤ A := le_max_right _ _
    exact positive_future_incompatible_with_negative_derivative (r := fun t => R + 1 - r t)
      (dr := fun t => -dr t) (C := L / 2) (by linarith)
      (fun t ht => by have hh := (hr t (haA.trans ht)).2; linarith)
      (fun t ht => (hder t (haA.trans ht)).const_sub (R + 1))
      (fun t ht => by have hh := hT t (hTA.trans ht); linarith)

theorem zero_outgoing_speed_forces_charge {q W : ℝ} (hq : 1 < q)
    (hW : q - 1 ≤ W) (hzero : radialSpeed q W 1 = 0) : W = q := by
  have hqpos : 0 < q := by linarith
  have hWpos : 0 < W := by linarith
  have hroot : Real.sqrt (q * W) = W := by unfold radialSpeed at hzero; linarith
  have hs := Real.sq_sqrt (mul_pos hqpos hWpos).le
  rw [hroot] at hs
  have heq : W * (q - W) = 0 := by nlinarith
  rcases mul_eq_zero.mp heq with hwzero | hdiff
  · exact False.elim ((ne_of_gt hWpos) hwzero)
  · linarith

/-- Full bounded-lingering exclusion for the declared future ODE solution.
An angular primitive cannot remain bounded when radius stays in a compact
positive interval, the solution is exterior, and the bounded mass is monotone. -/
theorem no_bounded_exterior_finite_angle {q M a ρ R U : ℝ} (flow : DirectionFlow q M a)
    (φ : ℝ → ℝ) (hq : 1 < q) (hρ : 0 < ρ)
    (hr : ∀ t ≥ a, ρ ≤ flow.radius t ∧ flow.radius t ≤ R)
    (hexterior : ∀ t ≥ a, 2 * flow.mass t ≤ flow.radius t)
    (hmmono : MonotoneOn flow.mass (Set.Ici a))
    (hφode : ∀ t ≥ a, HasDerivAt φ
      (Real.sqrt (q - 1 + 2 * flow.mass t / flow.radius t) / flow.radius t * Real.sin (flow.angle t)) t)
    (hφupper : ∀ t ≥ a, φ t ≤ U) : False := by
  let w := fun t => q - 1 + 2 * flow.mass t / flow.radius t
  have hw : ∀ t ≥ a, q - 1 ≤ w t ∧ w t ≤ q := by
    intro t ht
    have hm := (flow.mass_range t ht).1
    have hrt := flow.radius_pos t ht
    have hdivn : 0 ≤ 2 * flow.mass t / flow.radius t := by positivity
    have hdivu : 2 * flow.mass t / flow.radius t ≤ 1 := (div_le_iff₀ hrt).2 (by simpa using hexterior t ht)
    dsimp [w]
    constructor <;> linarith
  have hαlim := bounded_exterior_finite_angle_outgoing_limit hq hρ hr hw flow.angle_range
    flow.angular_ode hφode flow.radial_ode hφupper
  obtain ⟨Lr, hrlim⟩ := bounded_radius_converges_from_outgoing_limit hq hρ hr hw
    flow.angle_range hφode flow.radial_ode hφupper hαlim
  have hLr : ρ ≤ Lr := ge_of_tendsto hrlim (by
    filter_upwards [eventually_ge_atTop a] with t ht
    exact (hr t ht).1)
  have hLrpos : 0 < Lr := hρ.trans_le hLr
  obtain ⟨Lm, hmlim⟩ := bounded_monotone_future_converges hmmono
    (fun t ht => (flow.mass_range t ht).2)
  let W := q - 1 + 2 * Lm / Lr
  have hwlim : Tendsto w atTop (𝓝 W) :=
    tendsto_const_nhds.add ((hmlim.const_mul 2).div hrlim (ne_of_gt hLrpos))
  have hWlow : q - 1 ≤ W := ge_of_tendsto hwlim (by
    filter_upwards [eventually_ge_atTop a] with t ht
    exact (hw t ht).1)
  have hcos : Tendsto (fun t => Real.cos (flow.angle t)) atTop (𝓝 1) := by
    simpa only [Real.cos_zero] using Real.continuous_cos.continuousAt.tendsto.comp hαlim
  have hroot : Tendsto (fun t => Real.sqrt (q * w t)) atTop (𝓝 (Real.sqrt (q * W))) :=
    Real.continuous_sqrt.continuousAt.tendsto.comp (hwlim.const_mul q)
  have hdrLim : Tendsto (fun t => radialSpeed q (w t) (Real.cos (flow.angle t))) atTop
      (𝓝 (radialSpeed q W 1)) := hwlim.neg.add (hroot.mul hcos)
  have hzero := derivative_limit_zero_on_bounded_positive_radius hρ hr flow.radial_ode hdrLim
  have hWq := zero_outgoing_speed_forces_charge hq hWlow hzero
  have hwq : Tendsto w atTop (𝓝 q) := by rwa [hWq] at hwlim
  exact no_radial_horizon_limit (by linarith : 0 < q) flow.radius_pos flow.angle_range
    flow.angular_ode hwq hαlim

/-- The bounded-lingering lemma in its positive-mass exterior form: the swept
angle diverges, rather than merely failing to have a finite limit. -/
theorem bounded_exterior_flow_angle_diverges {q M a R : ℝ} (flow : DirectionFlow q M a)
    (φ : ℝ → ℝ) (hq : 1 < q) (hmpositive : 0 < flow.mass a)
    (hrupper : ∀ t ≥ a, flow.radius t ≤ R)
    (hexterior : ∀ t ≥ a, 2 * flow.mass t ≤ flow.radius t)
    (hmmono : MonotoneOn flow.mass (Set.Ici a))
    (hφode : ∀ t ≥ a, HasDerivAt φ
      (Real.sqrt (q - 1 + 2 * flow.mass t / flow.radius t) / flow.radius t * Real.sin (flow.angle t)) t) :
    Tendsto φ atTop atTop := by
  have hρ : 0 < 2 * flow.mass a := by positivity
  have hr : ∀ t ≥ a, 2 * flow.mass a ≤ flow.radius t ∧ flow.radius t ≤ R := by
    intro t ht
    have hm := hmmono (by simp : a ∈ Set.Ici a)
      (Set.mem_Ici.mpr ht) ht
    exact ⟨by linarith [hexterior t ht], hrupper t ht⟩
  have hunbounded : ∀ K : ℝ, ∃ t ≥ a, K < φ t := by
    intro K
    by_contra hn
    have hupper : ∀ t ≥ a, φ t ≤ K := by
      intro t ht
      by_contra hnot
      exact hn ⟨t, ht, lt_of_not_ge hnot⟩
    exact no_bounded_exterior_finite_angle flow φ hq hρ hr hexterior hmmono hφode hupper
  apply tendsto_atTop.2
  intro K
  obtain ⟨T, haT, hKT⟩ := hunbounded K
  filter_upwards [eventually_ge_atTop T] with t ht
  have hg := increment_ge_of_derivative_lower (C := 0) ht
    (fun x hx => hφode x (haT.trans hx.1)) (fun x hx => ?_)
  · linarith
  · exact mul_nonneg
      (div_nonneg (Real.sqrt_nonneg _) (flow.radius_pos x (haT.trans hx.1.le)).le)
      (Real.sin_pos_of_pos_of_lt_pi (flow.angle_range x (haT.trans hx.1.le)).1
        (flow.angle_range x (haT.trans hx.1.le)).2).le

end ViaB
