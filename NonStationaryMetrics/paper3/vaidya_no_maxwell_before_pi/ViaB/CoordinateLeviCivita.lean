import ViaB.DirectionHamiltonianLift
import Mathlib.LinearAlgebra.Matrix.Notation

namespace ViaB

abbrev CoordinateIndex := Fin 3

noncomputable def fermatMetric (q m r : ℝ) : Matrix CoordinateIndex CoordinateIndex ℝ :=
  let w := vaidyaW q m r
  !![-w + w ^ 2 / q, w / q, 0; w / q, 1 / q, 0; 0, 0, r ^ 2]

noncomputable def fermatInverseMetric (q m r : ℝ) : Matrix CoordinateIndex CoordinateIndex ℝ :=
  !![-1 / vaidyaW q m r, 1, 0; 1, q - vaidyaW q m r, 0; 0, 0, 1 / r ^ 2]

/-- Coordinate order is (v,r,phi); mu is the instantaneous mass derivative. -/
noncomputable def fermatMetricPartial (q m μ r : ℝ)
    (k : CoordinateIndex) : Matrix CoordinateIndex CoordinateIndex ℝ :=
  fun i j =>
    let w := vaidyaW q m r
    match k.val, i.val, j.val with
    | 0, 0, 0 => (2 * w / q - 1) * (2 * μ / r)
    | 0, 0, 1 => 2 * μ / (q * r)
    | 0, 1, 0 => 2 * μ / (q * r)
    | 1, 0, 0 => (2 * w / q - 1) * (-2 * m / r ^ 2)
    | 1, 0, 1 => -2 * m / (q * r ^ 2)
    | 1, 1, 0 => -2 * m / (q * r ^ 2)
    | 1, 2, 2 => 2 * r
    | _, _, _ => 0

noncomputable def coordinateKoszul (D : CoordinateIndex → Matrix CoordinateIndex CoordinateIndex ℝ)
    (i j k : CoordinateIndex) : ℝ := (D j i k + D k i j - D i j k) / 2

noncomputable def fermatChristoffel (q m μ r : ℝ) (i j k : CoordinateIndex) : ℝ :=
  ∑ h : CoordinateIndex, fermatInverseMetric q m r i h *
    coordinateKoszul (fermatMetricPartial q m μ r) h j k

theorem fermatMetric_symmetric (q m r : ℝ) : ∀ i j, fermatMetric q m r i j = fermatMetric q m r j i := by
  intro i j
  fin_cases i <;> fin_cases j <;> simp [fermatMetric]

theorem fermatMetricPartial_symmetric (q m μ r : ℝ) (k : CoordinateIndex) :
    ∀ i j, fermatMetricPartial q m μ r k i j = fermatMetricPartial q m μ r k j i := by
  intro i j
  fin_cases k <;> fin_cases i <;> fin_cases j <;> simp [fermatMetricPartial]

theorem fermatMetric_mul_inverse {q m r : ℝ}
    (hq : q ≠ 0) (hr : r ≠ 0) (hw : vaidyaW q m r ≠ 0) :
    fermatMetric q m r * fermatInverseMetric q m r = 1 := by
  ext i j
  fin_cases i <;> fin_cases j <;>
    simp [Matrix.mul_apply, Fin.sum_univ_three, fermatMetric, fermatInverseMetric]
  all_goals field_simp
  all_goals ring

theorem fermatInverseMetric_mul_metric {q m r : ℝ}
    (hq : q ≠ 0) (hr : r ≠ 0) (hw : vaidyaW q m r ≠ 0) :
    fermatInverseMetric q m r * fermatMetric q m r = 1 := by
  ext i j
  fin_cases i <;> fin_cases j <;>
    simp [Matrix.mul_apply, Fin.sum_univ_three, fermatMetric, fermatInverseMetric]
  all_goals field_simp
  all_goals ring

/-- These are actual time derivatives of the metric, with r held fixed. -/
theorem fermatMetric_time_derivative {q r t μ : ℝ} {m : ℝ → ℝ}
    (hm : HasDerivAt m μ t) (i j : CoordinateIndex) :
    HasDerivAt (fun v => fermatMetric q (m v) r i j)
      (fermatMetricPartial q (m t) μ r 0 i j) t := by
  have hdw : HasDerivAt (fun v => vaidyaW q (m v) r) (2 * μ / r) t :=
    ((hm.const_mul 2).div_const r).const_add (q - 1)
  have hA := hdw.neg.add ((hdw.pow 2).div_const q)
  have hB := hdw.div_const q
  fin_cases i <;> fin_cases j
  · convert hA using 1
    dsimp [fermatMetricPartial]
    ring
  · convert hB using 1
    dsimp [fermatMetricPartial]
    ring
  · exact hasDerivAt_const _ _
  · convert hB using 1
    dsimp [fermatMetricPartial]
    ring
  all_goals exact hasDerivAt_const _ _

/-- These are actual radial derivatives of the metric, with m held fixed. -/
theorem fermatMetric_radius_derivative {q m μ r : ℝ} (hr : r ≠ 0)
    (i j : CoordinateIndex) :
    HasDerivAt (fun x => fermatMetric q m x i j)
      (fermatMetricPartial q m μ r 1 i j) r := by
  have hdw : HasDerivAt (fun x => vaidyaW q m x) (-2 * m / r ^ 2) r := by
    convert ((hasDerivAt_const r (2 * m)).div (hasDerivAt_id r) hr).const_add (q - 1) using 1
    dsimp [vaidyaW]
    ring
  have hA := hdw.neg.add ((hdw.pow 2).div_const q)
  have hB := hdw.div_const q
  fin_cases i <;> fin_cases j
  · convert hA using 1
    dsimp [fermatMetricPartial]
    ring
  · convert hB using 1
    dsimp [fermatMetricPartial]
    ring
  · exact hasDerivAt_const _ _
  · convert hB using 1
    dsimp [fermatMetricPartial]
    ring
  · exact hasDerivAt_const _ _
  · exact hasDerivAt_const _ _
  · exact hasDerivAt_const _ _
  · exact hasDerivAt_const _ _
  · change HasDerivAt (fun x : ℝ => x ^ 2) (2 * r) r
    convert (hasDerivAt_id r).pow 2 using 1
    norm_num

theorem fermatMetric_angle_derivative (q m μ r φ : ℝ) (i j : CoordinateIndex) :
    HasDerivAt (fun _ : ℝ => fermatMetric q m r i j)
      (fermatMetricPartial q m μ r 2 i j) φ := by
  simpa [fermatMetricPartial] using hasDerivAt_const φ (fermatMetric q m r i j)

theorem coordinateKoszul_torsion_free {D : CoordinateIndex → Matrix CoordinateIndex CoordinateIndex ℝ}
    (hD : ∀ k i j, D k i j = D k j i) (i j k : CoordinateIndex) :
    coordinateKoszul D i j k = coordinateKoszul D i k j := by
  have hs : D i j k = D i k j := hD i j k
  unfold coordinateKoszul
  rw [hs]
  ring

theorem coordinateKoszul_metric_compatible {D : CoordinateIndex → Matrix CoordinateIndex CoordinateIndex ℝ}
    (hD : ∀ k i j, D k i j = D k j i) (i j k : CoordinateIndex) :
    D k i j = coordinateKoszul D i k j + coordinateKoszul D j k i := by
  have h₁ : D k j i = D k i j := hD k j i
  have h₂ : D j k i = D j i k := hD j k i
  have h₃ : D i k j = D i j k := hD i k j
  unfold coordinateKoszul
  rw [h₁, h₂, h₃]
  ring

/-- Torsion freedom and metric compatibility determine the lowered
connection uniquely; the Koszul formula is derived from those properties. -/
theorem coordinateKoszul_unique {D : CoordinateIndex → Matrix CoordinateIndex CoordinateIndex ℝ}
    {C : CoordinateIndex → CoordinateIndex → CoordinateIndex → ℝ}
    (hT : ∀ i j k, C i j k = C i k j)
    (hM : ∀ i j k, D k i j = C i k j + C j k i) (i j k : CoordinateIndex) :
    C i j k = coordinateKoszul D i j k := by
  unfold coordinateKoszul
  rw [hM i k j, hM i j k, hM j k i,
    hT i k j, hT k j i, hT j k i]
  ring

theorem fermatChristoffel_lowered {q m μ r : ℝ}
    (hq : q ≠ 0) (hr : r ≠ 0) (hw : vaidyaW q m r ≠ 0)
    (i j k : CoordinateIndex) :
    (∑ h : CoordinateIndex, fermatMetric q m r i h * fermatChristoffel q m μ r h j k) =
      coordinateKoszul (fermatMetricPartial q m μ r) i j k := by
  calc
    _ = ∑ o : CoordinateIndex, (fermatMetric q m r * fermatInverseMetric q m r) i o *
        coordinateKoszul (fermatMetricPartial q m μ r) o j k := by
      simp only [fermatChristoffel, Matrix.mul_apply, Finset.mul_sum, Finset.sum_mul]
      rw [Finset.sum_comm]
      apply Finset.sum_congr rfl
      intro o _
      apply Finset.sum_congr rfl
      intro h _
      ring
    _ = _ := by
      rw [fermatMetric_mul_inverse hq hr hw]
      simp only [Matrix.one_apply, ite_mul, one_mul, zero_mul, Fintype.sum_ite_eq]

theorem fermatChristoffel_torsion_free (q m μ r : ℝ) (i j k : CoordinateIndex) :
    fermatChristoffel q m μ r i j k = fermatChristoffel q m μ r i k j := by
  unfold fermatChristoffel
  apply Finset.sum_congr rfl
  intro h _
  rw [coordinateKoszul_torsion_free (fermatMetricPartial_symmetric q m μ r)]

theorem fermatChristoffel_metric_compatible {q m μ r : ℝ}
    (hq : q ≠ 0) (hr : r ≠ 0) (hw : vaidyaW q m r ≠ 0)
    (i j k : CoordinateIndex) :
    fermatMetricPartial q m μ r k i j =
      (∑ h : CoordinateIndex, fermatMetric q m r i h * fermatChristoffel q m μ r h k j) +
      (∑ h : CoordinateIndex, fermatMetric q m r j h * fermatChristoffel q m μ r h k i) := by
  rw [fermatChristoffel_lowered hq hr hw, fermatChristoffel_lowered hq hr hw]
  exact coordinateKoszul_metric_compatible (fermatMetricPartial_symmetric q m μ r) i j k

/-- In this chart the computed connection is the unique torsion-free,
metric-compatible connection, with uniqueness verified also for upper indices. -/
theorem fermatChristoffel_unique {q m μ r : ℝ}
    {C : CoordinateIndex → CoordinateIndex → CoordinateIndex → ℝ}
    (hq : q ≠ 0) (hr : r ≠ 0) (hw : vaidyaW q m r ≠ 0)
    (hT : ∀ i j k, C i j k = C i k j)
    (hM : ∀ i j k, fermatMetricPartial q m μ r k i j =
      (∑ h : CoordinateIndex, fermatMetric q m r i h * C h k j) +
      (∑ h : CoordinateIndex, fermatMetric q m r j h * C h k i))
    (i j k : CoordinateIndex) : C i j k = fermatChristoffel q m μ r i j k := by
  let lower : CoordinateIndex → CoordinateIndex → CoordinateIndex → ℝ :=
    fun i j k => ∑ h : CoordinateIndex, fermatMetric q m r i h * C h j k
  have hLT : ∀ i j k, lower i j k = lower i k j := by
    intro i j k
    apply Finset.sum_congr rfl
    intro h _
    rw [hT h j k]
  have hLM : ∀ i j k, fermatMetricPartial q m μ r k i j = lower i k j + lower j k i := hM
  have hK : ∀ i j k, lower i j k = coordinateKoszul (fermatMetricPartial q m μ r) i j k :=
    coordinateKoszul_unique hLT hLM
  calc
    C i j k = ∑ o : CoordinateIndex,
        (fermatInverseMetric q m r * fermatMetric q m r) i o * C o j k := by
      rw [fermatInverseMetric_mul_metric hq hr hw]
      simp only [Matrix.one_apply, ite_mul, one_mul, zero_mul, Fintype.sum_ite_eq]
    _ = ∑ h : CoordinateIndex, fermatInverseMetric q m r i h * lower h j k := by
      simp only [Matrix.mul_apply, lower, Finset.sum_mul, Finset.mul_sum]
      rw [Finset.sum_comm]
      apply Finset.sum_congr rfl
      intro h _
      apply Finset.sum_congr rfl
      intro o _
      ring
    _ = _ := by
      unfold fermatChristoffel
      apply Finset.sum_congr rfl
      intro h _
      rw [hK h j k]

noncomputable def fermatChristoffelContraction (q m μ r V R P : ℝ)
    (i : CoordinateIndex) : ℝ :=
  ∑ j : CoordinateIndex, ∑ k : CoordinateIndex,
    fermatChristoffel q m μ r i j k * (![V, R, P] j) * (![V, R, P] k)

noncomputable def loweredConnectionContraction (q m μ r V R P : ℝ) : ℝ × ℝ × ℝ :=
  let w := vaidyaW q m r
  ((2 * w / q - 1) * μ / r * V ^ 2 -
    2 * m * (2 * w / q - 1) / r ^ 2 * V * R - 2 * m / (q * r ^ 2) * R ^ 2,
    (2 * μ / (q * r) + m * (2 * w / q - 1) / r ^ 2) * V ^ 2 - r * P ^ 2,
    2 * r * R * P)

theorem lowered_connection_contraction_formula {q m μ r : ℝ}
    (hq : q ≠ 0) (hr : r ≠ 0) (V R P : ℝ) (i : CoordinateIndex) :
    (∑ j : CoordinateIndex, ∑ k : CoordinateIndex,
      coordinateKoszul (fermatMetricPartial q m μ r) i j k *
        (![V, R, P] j) * (![V, R, P] k)) =
      let C := loweredConnectionContraction q m μ r V R P
      ![C.1, C.2.1, C.2.2] i := by
  fin_cases i <;>
    simp [Fin.sum_univ_three, coordinateKoszul,
      fermatMetricPartial, loweredConnectionContraction]
  all_goals field_simp
  all_goals ring

/-- The connection contraction is computed from metric derivatives and
the Koszul formula, rather than declared as the desired acceleration. -/
theorem fermatChristoffel_contraction_formula {q m μ r : ℝ}
    (hq : q ≠ 0) (hr : r ≠ 0) (V R P : ℝ)
    (i : CoordinateIndex) :
    fermatChristoffelContraction q m μ r V R P i =
      let C := raiseMomentum q (vaidyaW q m r) r
        (loweredConnectionContraction q m μ r V R P)
      ![C.1, C.2.1, C.2.2] i := by
  have hsum : fermatChristoffelContraction q m μ r V R P i =
      ∑ h : CoordinateIndex, fermatInverseMetric q m r i h *
        (∑ j : CoordinateIndex, ∑ k : CoordinateIndex,
          coordinateKoszul (fermatMetricPartial q m μ r) h j k *
            (![V, R, P] j) * (![V, R, P] k)) := by
    simp only [fermatChristoffelContraction, fermatChristoffel, Finset.sum_mul]
    calc
      _ = ∑ j : CoordinateIndex, ∑ h : CoordinateIndex, ∑ k : CoordinateIndex,
          fermatInverseMetric q m r i h *
            coordinateKoszul (fermatMetricPartial q m μ r) h j k *
            (![V, R, P] j) * (![V, R, P] k) := by
        apply Finset.sum_congr rfl
        intro j _
        rw [Finset.sum_comm]
      _ = _ := by
        rw [Finset.sum_comm]
        simp only [Finset.mul_sum]
        apply Finset.sum_congr rfl
        intro h _
        apply Finset.sum_congr rfl
        intro j _
        apply Finset.sum_congr rfl
        intro k _
        ring
  rw [hsum]
  simp_rw [lowered_connection_contraction_formula hq hr]
  fin_cases i <;> simp [Fin.sum_univ_three, fermatInverseMetric, raiseMomentum] <;> ring

end ViaB
