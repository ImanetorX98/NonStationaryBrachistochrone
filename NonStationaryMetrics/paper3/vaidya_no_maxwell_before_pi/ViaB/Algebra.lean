import Mathlib.Analysis.SpecialFunctions.Sqrt
import Mathlib.Tactic.Linarith
import Mathlib.Tactic.Positivity
import Mathlib.Tactic.Ring

/-! Scalar lemmas of Via B. Here q denotes the constant rail charge Ê².
These lemmas do not assert that a metric has the stated curvature formula. -/
namespace ViaB

def curvaturePolynomial (q m r : ℝ) : ℝ :=
  3 * (2 * q - 3) * (q - 1) * r ^ 2 + (33 * q - 38) * m * r + 40 * m ^ 2

theorem curvaturePolynomial_pos {q m r : ℝ}
    (hq : 3 / 2 ≤ q) (hm : 0 < m) (hr : 0 < r) :
    0 < curvaturePolynomial q m r := by
  have hq1 : 0 < q - 1 := by linarith
  have hq2 : 0 ≤ 2 * q - 3 := by linarith
  have hq3 : 0 < 33 * q - 38 := by linarith
  have hfirst : 0 ≤ 3 * (2 * q - 3) * (q - 1) * r ^ 2 := by positivity
  have hsecond : 0 < (33 * q - 38) * m * r := by positivity
  have hthird : 0 < 40 * m ^ 2 := by positivity
  unfold curvaturePolynomial
  linarith

noncomputable def screenGap (q m r massRate angularMomentum : ℝ) : ℝ :=
  angularMomentum ^ 2 / (r ^ 4 * (q * r + 2 * m - r) ^ 2) *
    (m * curvaturePolynomial q m r / r + massRate * r * (3 * (q - 1) * r + 4 * m))

theorem screenGap_pos {q m r massRate angularMomentum : ℝ}
    (hq : 3 / 2 ≤ q) (hm : 0 < m) (hr : 0 < r)
    (hmassRate : 0 ≤ massRate) (hL : angularMomentum ≠ 0) :
    0 < screenGap q m r massRate angularMomentum := by
  have hq1 : 0 < q - 1 := by linarith
  have hden : 0 < q * r + 2 * m - r := by
    nlinarith [mul_pos hq1 hr]
  have hP := curvaturePolynomial_pos hq hm hr
  have hLsq : 0 < angularMomentum ^ 2 := sq_pos_of_ne_zero hL
  unfold screenGap
  positivity

noncomputable def directionB (q w c : ℝ) : ℝ :=
  c * (2 * w - q + 1) - (1 / 2) * Real.sqrt (q / w) * (3 * w - q + 1)

theorem directionB_at_horizon {q : ℝ} (hq : 0 < q) :
    directionB q q 1 = 1 / 2 := by
  simp only [directionB, div_self (ne_of_gt hq), Real.sqrt_one]
  ring

/-- Uniformly negative limiting direction coefficient, even without c ≥ -1. -/
theorem directionB_at_infinity_negative {q c : ℝ} (hq : 1 < q) (hc : c ≤ 1) :
    directionB q (q - 1) c < 0 := by
  have hqm : 0 < q - 1 := by linarith
  have hratio : 1 < q / (q - 1) := (lt_div_iff₀ hqm).2 (by linarith)
  have hsnonneg : 0 ≤ Real.sqrt (q / (q - 1)) := Real.sqrt_nonneg _
  have hssq := Real.sq_sqrt (show 0 ≤ q / (q - 1) by linarith)
  have hs : 1 < Real.sqrt (q / (q - 1)) := by nlinarith
  have hform : directionB q (q - 1) c =
      (q - 1) * (c - Real.sqrt (q / (q - 1))) := by
    unfold directionB
    ring
  rw [hform]
  exact mul_neg_of_pos_of_neg hqm (by linarith)

end ViaB
