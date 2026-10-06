import ViaB.AngleLimits

open Filter Topology

namespace ViaB

theorem angular_speed_lower {q w r R : ℝ} (hr : 0 < r)
    (hrR : r ≤ R) (hw : q - 1 ≤ w) : Real.sqrt (q - 1) / R ≤ Real.sqrt w / r := by
  have hR : 0 < R := hr.trans_le hrR
  have hsqrt := Real.sqrt_le_sqrt hw
  have hs1 := mul_le_mul_of_nonneg_left hrR (Real.sqrt_nonneg (q - 1))
  have hs2 := mul_le_mul_of_nonneg_right hsqrt hR.le
  apply (div_le_div_iff₀ hR hr).2
  nlinarith

theorem radialSpeed_lower_near_outgoing {q w α : ℝ} (hq : 0 < q)
    (hw : 0 ≤ w) (hwq : w ≤ q) (hα : 0 ≤ α ∧ α ≤ Real.pi / 2) :
    -q * Real.sin α ≤ radialSpeed q w (Real.cos α) := by
  have hsnonneg := Real.sqrt_nonneg (q * w)
  have hssq := Real.sq_sqrt (mul_nonneg hq.le hw)
  have hrootup : Real.sqrt (q * w) ≤ q := by nlinarith
  have hrootlow : w ≤ Real.sqrt (q * w) :=
    (sq_le_sq₀ hw hsnonneg).mp (by nlinarith [mul_nonneg hw (sub_nonneg.mpr hwq)])
  have hsin : 0 ≤ Real.sin α := Real.sin_nonneg_of_nonneg_of_le_pi hα.1 (by linarith [Real.pi_pos])
  have hcos : 0 ≤ Real.cos α := by
    have hh := Real.cos_le_cos_of_nonneg_of_le_pi hα.1
      (show Real.pi / 2 ≤ Real.pi by linarith [Real.pi_pos]) hα.2
    simpa only [Real.cos_pi_div_two] using hh
  have hsquares := Real.sin_sq_add_cos_sq α
  have hsum : 1 ≤ Real.sin α + Real.cos α := by nlinarith [mul_nonneg hsin hcos]
  have hprod1 := mul_nonneg hsnonneg (show 0 ≤ Real.sin α + Real.cos α - 1 by linarith)
  have hprod2 := mul_nonneg (sub_nonneg.mpr hrootup) hsin
  unfold radialSpeed
  nlinarith

theorem bounded_radius_converges_from_outgoing_limit {r w α φ : ℝ → ℝ} {q a ρ R U : ℝ}
    (hq : 1 < q) (hρ : 0 < ρ)
    (hr : ∀ t ≥ a, ρ ≤ r t ∧ r t ≤ R)
    (hw : ∀ t ≥ a, q - 1 ≤ w t ∧ w t ≤ q)
    (hα : ∀ t ≥ a, 0 < α t ∧ α t < Real.pi)
    (hφode : ∀ t ≥ a, HasDerivAt φ (Real.sqrt (w t) / r t * Real.sin (α t)) t)
    (hrode : ∀ t ≥ a, HasDerivAt r (radialSpeed q (w t) (Real.cos (α t))) t)
    (hφupper : ∀ t ≥ a, φ t ≤ U) (hαlim : Tendsto α atTop (𝓝 0)) :
    ∃ L : ℝ, Tendsto r atTop (𝓝 L) := by
  have hqpos : 0 < q := by linarith
  have hqm : 0 < q - 1 := by linarith
  have hR : 0 < R := hρ.trans_le ((hr a le_rfl).1.trans (hr a le_rfl).2)
  let c := Real.sqrt (q - 1) / R
  have hc : 0 < c := div_pos (Real.sqrt_pos.2 hqm) hR
  have hφnonneg : ∀ t ≥ a, 0 ≤ Real.sqrt (w t) / r t * Real.sin (α t) := by
    intro t ht
    exact mul_nonneg (div_nonneg (Real.sqrt_nonneg _) (hρ.trans_le (hr t ht).1).le)
      (Real.sin_pos_of_pos_of_lt_pi (hα t ht).1 (hα t ht).2).le
  obtain ⟨Lφ, hφlim⟩ := bounded_nondecreasing_primitive_converges hφode hφnonneg hφupper
  have hev : ∀ᶠ t in atTop, α t < Real.pi / 2 :=
    hαlim.eventually (eventually_lt_nhds (half_pos Real.pi_pos))
  obtain ⟨T, hT⟩ := eventually_atTop.1 hev
  let A := max a T
  have haA : a ≤ A := le_max_left _ _
  have hTA : T ≤ A := le_max_right _ _
  let Y := fun t => r t + (q / c) * φ t
  let dY := fun t => radialSpeed q (w t) (Real.cos (α t)) +
    (q / c) * (Real.sqrt (w t) / r t * Real.sin (α t))
  have hderY : ∀ t ≥ A, HasDerivAt Y (dY t) t :=
    fun t ht => (hrode t (haA.trans ht)).add ((hφode t (haA.trans ht)).const_mul (q / c))
  have hnonnegY : ∀ t ≥ A, 0 ≤ dY t := by
    intro t ht
    have hat := haA.trans ht
    have hrt : 0 < r t := hρ.trans_le (hr t hat).1
    have hwt : 0 ≤ w t := by linarith [(hw t hat).1]
    have hrbound := radialSpeed_lower_near_outgoing hqpos hwt (hw t hat).2
      ⟨(hα t hat).1.le, (hT t (hTA.trans ht)).le⟩
    have hcpoint : c ≤ Real.sqrt (w t) / r t := angular_speed_lower hrt (hr t hat).2 (hw t hat).1
    have hsin : 0 ≤ Real.sin (α t) := (Real.sin_pos_of_pos_of_lt_pi (hα t hat).1 (hα t hat).2).le
    have hm := mul_le_mul_of_nonneg_left (mul_le_mul_of_nonneg_right hcpoint hsin)
      (div_pos hqpos hc).le
    have heq : (q / c) * (c * Real.sin (α t)) = q * Real.sin (α t) := by field_simp
    rw [heq] at hm
    dsimp [dY]
    linarith
  have hupperY : ∀ t ≥ A, Y t ≤ R + (q / c) * U := by
    intro t ht
    have hrt := (hr t (haA.trans ht)).2
    have hφt := mul_le_mul_of_nonneg_left (hφupper t (haA.trans ht)) (div_pos hqpos hc).le
    dsimp [Y]
    linarith
  obtain ⟨LY, hYlim⟩ := bounded_nondecreasing_primitive_converges hderY hnonnegY hupperY
  have hh := hYlim.sub (hφlim.const_mul (q / c))
  have heq : (fun t => Y t - (q / c) * φ t) = r := by funext t; dsimp [Y]; ring
  rw [heq] at hh
  exact ⟨LY - (q / c) * Lφ, hh⟩

theorem bounded_monotone_future_converges {f : ℝ → ℝ} {a U : ℝ}
    (hmono : MonotoneOn f (Set.Ici a)) (hupper : ∀ t ≥ a, f t ≤ U) :
    ∃ L : ℝ, Tendsto f atTop (𝓝 L) := by
  let g := fun t => f (max a t)
  have hmg : Monotone g := fun x y hxy =>
    hmono (show max a x ∈ Set.Ici a from le_max_left a x)
      (show max a y ∈ Set.Ici a from le_max_left a y) (max_le_max_left a hxy)
  have hbdd : BddAbove (Set.range g) := by
    refine ⟨U, ?_⟩
    rintro _ ⟨t, rfl⟩
    exact hupper _ (le_max_left _ _)
  have heq : g =ᶠ[atTop] f := by
    filter_upwards [eventually_ge_atTop a] with t ht
    simp [g, max_eq_right ht]
  exact ⟨⨆ t, g t, (tendsto_congr' heq).mp (tendsto_atTop_ciSup hmg hbdd)⟩

end ViaB
