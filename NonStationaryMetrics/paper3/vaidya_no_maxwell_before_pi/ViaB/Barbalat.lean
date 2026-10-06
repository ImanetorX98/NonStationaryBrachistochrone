import ViaB.Dynamics
import Mathlib.Topology.Order.MonotoneConvergence

open Filter Topology

namespace ViaB

theorem bounded_nondecreasing_primitive_converges {P dP : ℝ → ℝ} {a U : ℝ}
    (hder : ∀ t ≥ a, HasDerivAt P (dP t) t)
    (hpos : ∀ t ≥ a, 0 ≤ dP t) (hupper : ∀ t ≥ a, P t ≤ U) :
    ∃ L : ℝ, Tendsto P atTop (𝓝 L) := by
  let g := fun t => P (max a t)
  have hmono : Monotone g := by
    intro x y hxy
    have hh := increment_ge_of_derivative_lower (C := 0) (max_le_max_left a hxy)
      (fun t ht => hder t ((le_max_left _ _).trans ht.1))
      (fun t ht => hpos t ((le_max_left _ _).trans ht.1.le))
    dsimp [g]
    linarith
  have hbdd : BddAbove (Set.range g) := by
    refine ⟨U, ?_⟩
    rintro _ ⟨t, rfl⟩
    exact hupper _ (le_max_left _ _)
  have hlim := tendsto_atTop_ciSup hmono hbdd
  have heq : g =ᶠ[atTop] P := by
    filter_upwards [eventually_ge_atTop a] with t ht
    simp [g, max_eq_right ht]
  exact ⟨⨆ t, g t, (tendsto_congr' heq).mp hlim⟩

/-- Barbalat in a form suited to an angle primitive: nonnegative f, a bounded
derivative, and a convergent primitive with P' ≥ c f imply f → 0. No assumption
that the derivative of f is integrable is made. -/
theorem decay_from_convergent_primitive {f df P dP : ℝ → ℝ} {a C c L : ℝ}
    (hC : 0 < C) (hc : 0 < c)
    (hf : ∀ t ≥ a, HasDerivAt f (df t) t) (hP : ∀ t ≥ a, HasDerivAt P (dP t) t)
    (hfn : ∀ t ≥ a, 0 ≤ f t) (hdf : ∀ t ≥ a, |df t| ≤ C)
    (hdom : ∀ t ≥ a, c * f t ≤ dP t) (hlim : Tendsto P atTop (𝓝 L)) :
    Tendsto f atTop (𝓝 0) := by
  apply tendsto_order.2
  constructor
  · intro l hl
    filter_upwards [eventually_ge_atTop a] with t ht
    exact hl.trans_le (hfn t ht)
  · intro ε hε
    let h := ε / (2 * C)
    have hh : 0 < h := div_pos hε (by positivity)
    have hcancel : C * h = ε / 2 := by dsimp [h]; field_simp
    let D := (c * ε / 2) * h
    have hD : 0 < D := by dsimp [D]; positivity
    have hshift : Tendsto (fun t : ℝ => t + h) atTop atTop :=
      tendsto_atTop_add_const_right atTop h tendsto_id
    have hinc : Tendsto (fun t => P (t + h) - P t) atTop (𝓝 0) := by
      simpa using (hlim.comp hshift).sub hlim
    have hevent : ∀ᶠ t in atTop, P (t + h) - P t < D :=
      hinc.eventually (eventually_lt_nhds hD)
    filter_upwards [hevent, eventually_ge_atTop a] with t ht hat
    by_contra hn
    have hpeak : ε ≤ f t := le_of_not_gt hn
    have hlow : ∀ s ∈ Set.Ioo t (t + h), ε / 2 ≤ f s := by
      intro s hs
      have hg := increment_ge_of_derivative_lower (C := -C) hs.1.le
        (fun x hx => hf x (hat.trans hx.1))
        (fun x hx => (abs_le.mp (hdf x (hat.trans hx.1.le))).1)
      have htime : s - t ≤ h := by linarith [hs.2]
      have hmult := mul_le_mul_of_nonneg_left htime hC.le
      linarith
    have hPgrowth := increment_ge_of_derivative_lower (C := c * ε / 2)
      (show t ≤ t + h by linarith)
      (fun s hs => hP s (hat.trans hs.1)) (fun s hs => ?_)
    · dsimp [D] at ht
      nlinarith
    · have hmult := mul_le_mul_of_nonneg_left (hlow s hs) hc.le
      have hdoms := hdom s (hat.trans hs.1.le)
      nlinarith

theorem decay_from_bounded_primitive {f df P dP : ℝ → ℝ} {a C c U : ℝ}
    (hC : 0 < C) (hc : 0 < c)
    (hf : ∀ t ≥ a, HasDerivAt f (df t) t) (hP : ∀ t ≥ a, HasDerivAt P (dP t) t)
    (hfn : ∀ t ≥ a, 0 ≤ f t) (hdf : ∀ t ≥ a, |df t| ≤ C)
    (hdom : ∀ t ≥ a, c * f t ≤ dP t) (hupper : ∀ t ≥ a, P t ≤ U) :
    Tendsto f atTop (𝓝 0) := by
  have hnonneg : ∀ t ≥ a, 0 ≤ dP t :=
    fun t ht => (mul_nonneg hc.le (hfn t ht)).trans (hdom t ht)
  obtain ⟨L, hlim⟩ := bounded_nondecreasing_primitive_converges hP hnonneg hupper
  exact decay_from_convergent_primitive hC hc hf hP hfn hdf hdom hlim

end ViaB
