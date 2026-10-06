import ViaB.Dynamics

open Filter Topology

namespace ViaB

theorem no_zero_limit_of_positive_nondecreasing {f df : ℝ → ℝ} {a : ℝ}
    (hpos : 0 < f a) (hder : ∀ t ≥ a, HasDerivAt f (df t) t)
    (hbound : ∀ t ≥ a, 0 ≤ df t) : ¬ Tendsto f atTop (𝓝 0) := by
  intro hlim
  have hupper : ∀ᶠ t in atTop, f t < f a := hlim.eventually (eventually_lt_nhds hpos)
  have hlower : ∀ᶠ t in atTop, f a ≤ f t := by
    filter_upwards [eventually_ge_atTop a] with t hat
    have hh := increment_ge_of_derivative_lower (C := 0) hat
      (fun x hx => hder x hx.1) (fun x hx => hbound x hx.1.le)
    linarith
  have hfalse : ∀ᶠ (t : ℝ) in atTop, False := by
    filter_upwards [hupper, hlower] with t hu hl
    exact (not_lt_of_ge hl) hu
  obtain ⟨t, ht⟩ := hfalse.exists
  exact ht

/-- Critical lingering limit excluded directly by the direction ODE.
The preceding analytic steps producing w → q and alpha → 0 are not assumed proved. -/
theorem no_radial_horizon_limit {r w α : ℝ → ℝ} {q a : ℝ} (hq : 0 < q)
    (hr : ∀ t ≥ a, 0 < r t)
    (hα : ∀ t ≥ a, 0 < α t ∧ α t < Real.pi)
    (hode : ∀ t ≥ a, HasDerivAt α
      (Real.sin (α t) / r t * directionB q (w t) (Real.cos (α t))) t)
    (hw : Tendsto w atTop (𝓝 q)) : ¬ Tendsto α atTop (𝓝 0) := by
  intro hlim
  have hc : Tendsto (fun t => Real.cos (α t)) atTop (𝓝 1) := by
    simpa only [Real.cos_zero] using Real.continuous_cos.continuousAt.tendsto.comp hlim
  have hB : Tendsto (fun t => directionB q (w t) (Real.cos (α t))) atTop (𝓝 (1 / 2)) := by
    have hh := (directionB_continuousAt (q := q) (w := q) (c := 1) (ne_of_gt hq)).tendsto.comp
      (hw.prodMk_nhds hc)
    simpa only [Function.comp_def, directionB_at_horizon hq] using hh
  have hBevent : ∀ᶠ t in atTop, 0 < directionB q (w t) (Real.cos (α t)) :=
    hB.eventually (eventually_gt_nhds (by norm_num : (0 : ℝ) < 1 / 2))
  obtain ⟨T, hT⟩ := eventually_atTop.1 hBevent
  let A := max a T
  have haA : a ≤ A := le_max_left _ _
  have hTA : T ≤ A := le_max_right _ _
  have hnonneg : ∀ t ≥ A, 0 ≤ Real.sin (α t) / r t * directionB q (w t) (Real.cos (α t)) := by
    intro t ht
    have hat := haA.trans ht
    exact (mul_pos (div_pos (Real.sin_pos_of_pos_of_lt_pi (hα t hat).1 (hα t hat).2)
      (hr t hat)) (hT t (hTA.trans ht))).le
  exact no_zero_limit_of_positive_nondecreasing (hα A haA).1
    (fun t ht => hode t (haA.trans ht)) hnonneg hlim

theorem positive_future_incompatible_with_negative_derivative {r dr : ℝ → ℝ} {a C : ℝ}
    (hC : 0 < C) (hr : ∀ t ≥ a, 0 < r t)
    (hder : ∀ t ≥ a, HasDerivAt r (dr t) t)
    (hbound : ∀ t ≥ a, dr t ≤ -C) : False := by
  let b := a + r a / C + 1
  have ha : 0 < r a := hr a le_rfl
  have hab : a ≤ b := by dsimp [b]; linarith [div_pos ha hC]
  have hg := increment_le_of_derivative_upper hab
    (fun t ht => hder t ht.1) (fun t ht => hbound t ht.1.le)
  have hb := hr b hab
  have heq : C * (b - a) = r a + C := by dsimp [b]; field_simp; ring
  nlinarith

/-- A positive-radius future solution cannot converge to the incoming radial direction. -/
theorem no_incoming_direction_limit {r w α : ℝ → ℝ} {q a : ℝ} (hq : 1 < q)
    (hr : ∀ t ≥ a, 0 < r t) (hw : ∀ t ≥ a, q - 1 ≤ w t)
    (hode : ∀ t ≥ a, HasDerivAt r (radialSpeed q (w t) (Real.cos (α t))) t) :
    ¬ Tendsto α atTop (𝓝 Real.pi) := by
  intro hlim
  have hc : Tendsto (fun t => Real.cos (α t)) atTop (𝓝 (-1)) := by
    simpa only [Real.cos_pi] using Real.continuous_cos.continuousAt.tendsto.comp hlim
  have hevent : ∀ᶠ t in atTop, Real.cos (α t) < 0 :=
    hc.eventually (eventually_lt_nhds (by norm_num : (-1 : ℝ) < 0))
  obtain ⟨T, hT⟩ := eventually_atTop.1 hevent
  let A := max a T
  have haA : a ≤ A := le_max_left _ _
  have hTA : T ≤ A := le_max_right _ _
  apply positive_future_incompatible_with_negative_derivative (C := q - 1) (by linarith)
    (fun t ht => hr t (haA.trans ht)) (fun t ht => hode t (haA.trans ht))
  intro t ht
  have hmul := mul_nonpos_of_nonneg_of_nonpos (Real.sqrt_nonneg (q * w t))
    (hT t (hTA.trans ht)).le
  have hh := hw t (haA.trans ht)
  unfold radialSpeed
  linarith

end ViaB
