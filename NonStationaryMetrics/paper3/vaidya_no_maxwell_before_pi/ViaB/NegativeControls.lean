import ViaB.Algebra
import Mathlib.Tactic.NormNum

namespace ViaB

/-- Exact counterexample to extending the polynomial sign to all q > 1.
This is a sign counterexample, not a Vaidya Maxwell counterexample. -/
theorem curvaturePolynomial_negative_below_threshold :
    curvaturePolynomial (5 / 4) 1 20 = -45 := by
  norm_num [curvaturePolynomial]

theorem screenGap_negative_below_threshold :
    screenGap (5 / 4) 1 20 0 1 < 0 := by
  norm_num [screenGap, curvaturePolynomial]

/-- Strict positivity cannot be claimed on a radial ray using this formula. -/
theorem screenGap_radial (q m r massRate : ℝ) :
    screenGap q m r massRate 0 = 0 := by
  simp [screenGap]

end ViaB
