import ViaB.FlowDependence
import Mathlib.Tactic.FieldSimp

namespace ViaB

/-- Coordinate quadratic form of the equatorial Vaidya Fermat metric. -/
noncomputable def metricQuad (q w r V R P : ℝ) : ℝ :=
  -w * V ^ 2 + (R + w * V) ^ 2 / q + r ^ 2 * P ^ 2

/-- Coordinate quadratic form of the inverse metric. -/
noncomputable def cometricQuad (q w r Pv Pr L : ℝ) : ℝ :=
  -Pv ^ 2 / w + 2 * Pv * Pr + (q - w) * Pr ^ 2 + L ^ 2 / r ^ 2

noncomputable def lowerMomentum (q w r : ℝ) (x : ℝ × ℝ × ℝ) : ℝ × ℝ × ℝ :=
  (-w * x.1 + w * (x.2.1 + w * x.1) / q,
    (x.2.1 + w * x.1) / q, r ^ 2 * x.2.2)

noncomputable def raiseMomentum (q w r : ℝ) (p : ℝ × ℝ × ℝ) : ℝ × ℝ × ℝ :=
  (-p.1 / w + p.2.1, p.1 + (q - w) * p.2.1, p.2.2 / r ^ 2)

theorem raise_lower_identity {q w r : ℝ} (hq : q ≠ 0) (hw : w ≠ 0) (hr : r ≠ 0)
    (x : ℝ × ℝ × ℝ) : raiseMomentum q w r (lowerMomentum q w r x) = x := by
  rcases x with ⟨V, R, P⟩
  ext <;> simp only [raiseMomentum, lowerMomentum]
  all_goals field_simp
  all_goals ring

theorem lower_raise_identity {q w r : ℝ} (hq : q ≠ 0) (hw : w ≠ 0) (hr : r ≠ 0)
    (p : ℝ × ℝ × ℝ) : lowerMomentum q w r (raiseMomentum q w r p) = p := by
  rcases p with ⟨Pv, Pr, L⟩
  ext <;> simp only [raiseMomentum, lowerMomentum]
  all_goals field_simp
  all_goals ring

/-- The inverse quadratic form is obtained by raising the covector. -/
theorem metric_of_raised_momentum {q w r : ℝ}
    (hq : q ≠ 0) (hw : w ≠ 0) (hr : r ≠ 0) (Pv Pr L : ℝ) :
    metricQuad q w r (-Pv / w + Pr) (Pv + (q - w) * Pr) (L / r ^ 2) =
      cometricQuad q w r Pv Pr L := by
  unfold metricQuad cometricQuad
  field_simp
  ring

/-- Nullness gives the normalized spatial velocity ellipse without any
geodesic equation or assumed direction dynamics. -/
theorem null_velocity_ellipse {q w r V R P : ℝ}
    (hq : q ≠ 0) (hw : w ≠ 0) (hV : V ≠ 0)
    (hnull : metricQuad q w r V R P = 0) :
    (R / V + w) ^ 2 / (q * w) + (r * (P / V)) ^ 2 / w = 1 := by
  unfold metricQuad at hnull
  have hn : (R + w * V) ^ 2 + q * r ^ 2 * P ^ 2 = q * w * V ^ 2 := by
    field_simp at hnull
    nlinarith
  field_simp
  nlinarith [hn]

/-- Twice the Hamiltonian, completed as a square in the time momentum. -/
theorem cometric_complete_square {q w r Pv Pr L : ℝ} (hw : w ≠ 0) :
    cometricQuad q w r Pv Pr L =
      -(Pv - w * Pr) ^ 2 / w + q * Pr ^ 2 + L ^ 2 / r ^ 2 := by
  unfold cometricQuad
  field_simp
  ring

noncomputable def nullAmplitude (q r Pr L : ℝ) : ℝ :=
  Real.sqrt (q * Pr ^ 2 + L ^ 2 / r ^ 2)

theorem nullAmplitude_pos {q r Pr L : ℝ} (hq : 0 < q) (hr : r ≠ 0) (hL : L ≠ 0) :
    0 < nullAmplitude q r Pr L := by
  unfold nullAmplitude
  apply Real.sqrt_pos.2
  have hterm : 0 < L ^ 2 / r ^ 2 := div_pos (sq_pos_of_ne_zero hL) (sq_pos_of_ne_zero hr)
  positivity

theorem nullAmplitude_sq {q r Pr L : ℝ} (hq : 0 ≤ q) :
    nullAmplitude q r Pr L ^ 2 = q * Pr ^ 2 + L ^ 2 / r ^ 2 := by
  exact Real.sq_sqrt (by positivity)

/-- The future null sheet selected by V>0. -/
theorem future_null_sheet {q w r Pr L : ℝ} (hq : 0 < q) (hw : 0 < w)
    (hr : r ≠ 0) (hL : L ≠ 0) :
    cometricQuad q w r (w * Pr - Real.sqrt w * nullAmplitude q r Pr L) Pr L = 0 ∧
    0 < -(w * Pr - Real.sqrt w * nullAmplitude q r Pr L) / w + Pr := by
  have hS := nullAmplitude_pos (Pr := Pr) hq hr hL
  have hs := nullAmplitude_sq (r := r) (Pr := Pr) (L := L) hq.le
  have hroot := Real.sq_sqrt hw.le
  have heq : -(w * Pr - Real.sqrt w * nullAmplitude q r Pr L) / w + Pr =
      nullAmplitude q r Pr L / Real.sqrt w := by
    field_simp
    nlinarith [hroot]
  constructor
  · rw [cometric_complete_square (ne_of_gt hw)]
    have hh : w * Pr - Real.sqrt w * nullAmplitude q r Pr L - w * Pr =
        -Real.sqrt w * nullAmplitude q r Pr L := by ring
    rw [hh, mul_pow, neg_sq, hroot]
    field_simp at hs
    field_simp
    nlinarith [hs]
  · rw [heq]
    exact div_pos hS (Real.sqrt_pos.2 hw)

/-- Nullness and positive time velocity uniquely select the future sheet. -/
theorem future_null_sheet_unique {q w r Pv Pr L : ℝ}
    (hq : 0 < q) (hw : 0 < w) (hr : r ≠ 0) (hL : L ≠ 0)
    (hnull : cometricQuad q w r Pv Pr L = 0) (hfuture : 0 < -Pv / w + Pr) :
    Pv = w * Pr - Real.sqrt w * nullAmplitude q r Pr L := by
  have hS := nullAmplitude_pos (Pr := Pr) hq hr hL
  have hs := nullAmplitude_sq (r := r) (Pr := Pr) (L := L) hq.le
  have hroot := Real.sq_sqrt hw.le
  rw [cometric_complete_square (ne_of_gt hw)] at hnull
  simp only [neg_div] at hnull
  have hdiv : (Pv - w * Pr) ^ 2 / w = nullAmplitude q r Pr L ^ 2 := by linarith
  have hsq := (div_eq_iff (ne_of_gt hw)).mp hdiv
  have hv : -Pv / w + Pr = -(Pv - w * Pr) / w := by field_simp; ring
  rw [hv] at hfuture
  have hneg := (div_pos_iff_of_pos_right hw).mp hfuture
  have hpos := mul_pos (Real.sqrt_pos.2 hw) hS
  have hsq' : (Pv - w * Pr) ^ 2 = (Real.sqrt w * nullAmplitude q r Pr L) ^ 2 := by
    rw [mul_pow, hroot]
    nlinarith [hsq]
  nlinarith [hsq']

end ViaB
