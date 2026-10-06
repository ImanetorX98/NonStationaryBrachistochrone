import ViaB.BoundedAngle
import ViaB.LingeringLimits
import Mathlib.Analysis.SpecialFunctions.Trigonometric.Inverse

open Filter Topology

namespace ViaB

theorem sine_decay_forces_radial_limit {α : ℝ → ℝ} {a : ℝ}
    (hcont : ∀ t ≥ a, ContinuousAt α t)
    (hrange : ∀ t ≥ a, 0 < α t ∧ α t < Real.pi)
    (hsin : Tendsto (fun t => Real.sin (α t)) atTop (𝓝 0)) :
    Tendsto α atTop (𝓝 0) ∨ Tendsto α atTop (𝓝 Real.pi) := by
  have hev : ∀ᶠ t in atTop, Real.sin (α t) < 1 :=
    hsin.eventually (eventually_lt_nhds (by norm_num : (0 : ℝ) < 1))
  obtain ⟨T, hT⟩ := eventually_atTop.1 hev
  let S := max a T
  have haS : a ≤ S := le_max_left _ _
  have hTS : T ≤ S := le_max_right _ _
  have hne : α S ≠ Real.pi / 2 := by
    intro hh
    have hval := hT S hTS
    rw [hh, Real.sin_pi_div_two] at hval
    exact (lt_irrefl _) hval
  have hasin : Tendsto (fun t => Real.arcsin (Real.sin (α t))) atTop (𝓝 0) := by
    simpa only [Real.arcsin_zero] using Real.continuous_arcsin.continuousAt.tendsto.comp hsin
  by_cases hlow : α S < Real.pi / 2
  · have hstay : ∀ t ≥ S, α t < Real.pi / 2 := by
      intro t hSt
      by_contra hn
      have hc : ContinuousOn α (Set.Icc S t) :=
        fun x hx => (hcont x (haS.trans hx.1)).continuousWithinAt
      obtain ⟨x, hx, hval⟩ := intermediate_value_Icc hSt hc ⟨hlow.le, le_of_not_gt hn⟩
      have hb := hT x (hTS.trans hx.1)
      rw [hval, Real.sin_pi_div_two] at hb
      exact (lt_irrefl _) hb
    left
    have heq : (fun t => Real.arcsin (Real.sin (α t))) =ᶠ[atTop] α := by
      filter_upwards [eventually_ge_atTop S] with t ht
      exact Real.arcsin_sin (by linarith [(hrange t (haS.trans ht)).1, Real.pi_pos]) (hstay t ht).le
    exact (tendsto_congr' heq).mp hasin
  · have hhigh : Real.pi / 2 < α S := lt_of_le_of_ne (le_of_not_gt hlow) (Ne.symm hne)
    have hstay : ∀ t ≥ S, Real.pi / 2 < α t := by
      intro t hSt
      by_contra hn
      have hc : ContinuousOn α (Set.Icc S t) :=
        fun x hx => (hcont x (haS.trans hx.1)).continuousWithinAt
      obtain ⟨x, hx, hval⟩ := intermediate_value_Icc' hSt hc ⟨le_of_not_gt hn, hhigh.le⟩
      have hb := hT x (hTS.trans hx.1)
      rw [hval, Real.sin_pi_div_two] at hb
      exact (lt_irrefl _) hb
    right
    have hlim : Tendsto (fun t => Real.pi - Real.arcsin (Real.sin (α t))) atTop (𝓝 Real.pi) := by
      simpa using hasin.const_sub Real.pi
    have heq : (fun t => Real.pi - Real.arcsin (Real.sin (α t))) =ᶠ[atTop] α := by
      filter_upwards [eventually_ge_atTop S] with t ht
      have hval := Real.arcsin_sin (x := Real.pi - α t)
        (by linarith [(hrange t (haS.trans ht)).2, Real.pi_pos]) (by linarith [hstay t ht])
      rw [Real.sin_pi_sub] at hval
      linarith
    exact (tendsto_congr' heq).mp hlim

theorem bounded_exterior_finite_angle_outgoing_limit {r w α φ : ℝ → ℝ} {q a ρ R U : ℝ}
    (hq : 1 < q) (hρ : 0 < ρ)
    (hr : ∀ t ≥ a, ρ ≤ r t ∧ r t ≤ R)
    (hw : ∀ t ≥ a, q - 1 ≤ w t ∧ w t ≤ q)
    (hα : ∀ t ≥ a, 0 < α t ∧ α t < Real.pi)
    (hαode : ∀ t ≥ a, HasDerivAt α
      (Real.sin (α t) / r t * directionB q (w t) (Real.cos (α t))) t)
    (hφode : ∀ t ≥ a, HasDerivAt φ (Real.sqrt (w t) / r t * Real.sin (α t)) t)
    (hrode : ∀ t ≥ a, HasDerivAt r (radialSpeed q (w t) (Real.cos (α t))) t)
    (hφupper : ∀ t ≥ a, φ t ≤ U) : Tendsto α atTop (𝓝 0) := by
  have hs := bounded_exterior_finite_angle_sine_decay hq hρ hr hw hα hαode hφode hφupper
  rcases sine_decay_forces_radial_limit (fun t ht => (hαode t ht).continuousAt) hα hs with hzero | hpi
  · exact hzero
  · exact False.elim (no_incoming_direction_limit hq
      (fun t ht => hρ.trans_le (hr t ht).1) (fun t ht => (hw t ht).1) hrode hpi)

end ViaB
