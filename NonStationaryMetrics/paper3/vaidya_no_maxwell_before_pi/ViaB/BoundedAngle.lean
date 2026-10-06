import ViaB.Barbalat
import Mathlib.Analysis.SpecialFunctions.Trigonometric.Deriv

open Filter Topology

namespace ViaB

noncomputable def directionBound (q : ℝ) : ℝ :=
  q + 1 + (1 / 2) * Real.sqrt (q / (q - 1)) * (2 * q + 1)

theorem directionBound_pos {q : ℝ} (hq : 1 < q) : 0 < directionBound q := by
  have hqp : 0 < q + 1 := by linarith
  have hcoef : 0 < 2 * q + 1 := by linarith
  unfold directionBound
  positivity

theorem directionB_abs_bound {q w c : ℝ} (hq : 1 < q)
    (hwl : q - 1 ≤ w) (hwu : w ≤ q) (hc : -1 ≤ c ∧ c ≤ 1) :
    |directionB q w c| ≤ directionBound q := by
  have hqm : 0 < q - 1 := by linarith
  have hqpos : 0 < q := by linarith
  have hwpos : 0 < w := hqm.trans_le hwl
  have hsqrt := Real.sqrt_le_sqrt (div_le_div_of_nonneg_left hqpos.le hqm hwl)
  have hroot := Real.sqrt_nonneg (q / w)
  have hrootMax := Real.sqrt_nonneg (q / (q - 1))
  have hA : 0 ≤ 2 * w - q + 1 := by linarith
  have hAup : 2 * w - q + 1 ≤ q + 1 := by linarith
  have hc1 := mul_nonneg (show 0 ≤ c + 1 by linarith [hc.1]) hA
  have hc2 := mul_nonneg (sub_nonneg.mpr hc.2) hA
  have hT : 0 ≤ 3 * w - q + 1 := by linarith
  have hTup : 3 * w - q + 1 ≤ 2 * q + 1 := by linarith
  have hprod := mul_le_mul hsqrt hTup hT hrootMax
  have hprodpos := mul_nonneg hroot hT
  apply abs_le.2
  unfold directionB directionBound
  constructor <;> nlinarith

theorem bounded_exterior_finite_angle_sine_decay {r w α φ : ℝ → ℝ} {q a ρ R U : ℝ}
    (hq : 1 < q) (hρ : 0 < ρ)
    (hr : ∀ t ≥ a, ρ ≤ r t ∧ r t ≤ R)
    (hw : ∀ t ≥ a, q - 1 ≤ w t ∧ w t ≤ q)
    (hα : ∀ t ≥ a, 0 < α t ∧ α t < Real.pi)
    (hαode : ∀ t ≥ a, HasDerivAt α
      (Real.sin (α t) / r t * directionB q (w t) (Real.cos (α t))) t)
    (hφode : ∀ t ≥ a, HasDerivAt φ (Real.sqrt (w t) / r t * Real.sin (α t)) t)
    (hφupper : ∀ t ≥ a, φ t ≤ U) :
    Tendsto (fun t => Real.sin (α t)) atTop (𝓝 0) := by
  let A := directionBound q
  let C := A / ρ + 1
  let c := Real.sqrt (q - 1) / R
  have hA : 0 < A := directionBound_pos hq
  have hR : 0 < R := hρ.trans_le ((hr a le_rfl).1.trans (hr a le_rfl).2)
  have hqm : 0 < q - 1 := by linarith
  have hc : 0 < c := div_pos (Real.sqrt_pos.2 hqm) hR
  have hC : 0 < C := by dsimp [C]; positivity
  apply decay_from_bounded_primitive (a := a) (C := C) (c := c) hC hc
    (fun t ht => (hαode t ht).sin) hφode
    (fun t ht => (Real.sin_pos_of_pos_of_lt_pi (hα t ht).1 (hα t ht).2).le) ?_ ?_ hφupper
  · intro t ht
    have hrpos : 0 < r t := hρ.trans_le (hr t ht).1
    have hsin : 0 ≤ Real.sin (α t) := (Real.sin_pos_of_pos_of_lt_pi (hα t ht).1 (hα t ht).2).le
    have hfac : 0 ≤ Real.sin (α t) / r t := div_nonneg hsin hrpos.le
    have hfacup : Real.sin (α t) / r t ≤ 1 / ρ := (div_le_div_iff₀ hrpos hρ).2 (by
      have hm := mul_le_mul_of_nonneg_right (Real.sin_le_one (α t)) hρ.le
      nlinarith [(hr t ht).1])
    have hB := directionB_abs_bound (c := Real.cos (α t)) hq (hw t ht).1 (hw t ht).2
      ⟨Real.neg_one_le_cos (α t), Real.cos_le_one (α t)⟩
    have hmul := mul_le_mul hfacup hB (abs_nonneg _) (show 0 ≤ 1 / ρ by positivity)
    have hαbound : |Real.sin (α t) / r t * directionB q (w t) (Real.cos (α t))| ≤ A / ρ := by
      rw [abs_mul, abs_of_nonneg hfac]
      simpa only [A, one_div, div_eq_mul_inv, mul_comm, one_mul, mul_one] using hmul
    have hcos : |Real.cos (α t)| ≤ 1 := abs_le.2 ⟨Real.neg_one_le_cos _, Real.cos_le_one _⟩
    rw [abs_mul]
    have hm := mul_le_mul_of_nonneg_right hcos
      (abs_nonneg (Real.sin (α t) / r t * directionB q (w t) (Real.cos (α t))))
    dsimp [C]
    nlinarith
  · intro t ht
    have hrpos : 0 < r t := hρ.trans_le (hr t ht).1
    have hsqrt : Real.sqrt (q - 1) ≤ Real.sqrt (w t) := Real.sqrt_le_sqrt (hw t ht).1
    have hs1 := mul_le_mul_of_nonneg_left (hr t ht).2 (Real.sqrt_nonneg (q - 1))
    have hs2 := mul_le_mul_of_nonneg_right hsqrt hR.le
    have hdiv : c ≤ Real.sqrt (w t) / r t := (div_le_div_iff₀ hR hrpos).2 (by nlinarith)
    exact mul_le_mul_of_nonneg_right hdiv
      (Real.sin_pos_of_pos_of_lt_pi (hα t ht).1 (hα t ht).2).le

end ViaB
