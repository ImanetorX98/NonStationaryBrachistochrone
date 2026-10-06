import ViaB.MetricHamiltonian

namespace ViaB

noncomputable def coordinateHamiltonian (q w r Pv Pr L : ℝ) : ℝ :=
  (1 / 2) * cometricQuad q w r Pv Pr L

theorem hamiltonian_time_momentum_derivative (q w r Pv Pr L : ℝ) :
    HasDerivAt (fun p => coordinateHamiltonian q w r p Pr L) (-Pv / w + Pr) Pv := by
  have hh := ((((hasDerivAt_id Pv).pow 2).neg.div_const w).add
    ((hasDerivAt_id Pv).const_mul (2 * Pr))).add
      (hasDerivAt_const Pv ((q - w) * Pr ^ 2 + L ^ 2 / r ^ 2))
  convert hh.const_mul (1 / 2) using 1
  · funext p; dsimp [coordinateHamiltonian, cometricQuad]; ring
  · dsimp [coordinateHamiltonian, cometricQuad]; ring

theorem hamiltonian_radial_momentum_derivative (q w r Pv Pr L : ℝ) :
    HasDerivAt (fun p => coordinateHamiltonian q w r Pv p L) (Pv + (q - w) * Pr) Pr := by
  have hh := ((hasDerivAt_const Pr (-Pv ^ 2 / w)).add
    ((hasDerivAt_id Pr).const_mul (2 * Pv))).add
      (((hasDerivAt_id Pr).pow 2).const_mul (q - w))
  convert (hh.add (hasDerivAt_const Pr (L ^ 2 / r ^ 2))).const_mul (1 / 2) using 1
  dsimp [coordinateHamiltonian, cometricQuad]
  ring

theorem hamiltonian_angular_momentum_derivative (q w r Pv Pr L : ℝ) :
    HasDerivAt (fun p => coordinateHamiltonian q w r Pv Pr p) (L / r ^ 2) L := by
  have hh := (hasDerivAt_const L (-Pv ^ 2 / w + 2 * Pv * Pr + (q - w) * Pr ^ 2)).add
    (((hasDerivAt_id L).pow 2).div_const (r ^ 2))
  convert hh.const_mul (1 / 2) using 1
  dsimp [coordinateHamiltonian, cometricQuad]
  ring

noncomputable def directionCos (q r Pr L : ℝ) : ℝ :=
  Real.sqrt q * Pr / nullAmplitude q r Pr L

noncomputable def directionSin (q r Pr L : ℝ) : ℝ :=
  L / (r * nullAmplitude q r Pr L)

theorem momentum_direction_circle {q r Pr L : ℝ}
    (hq : 0 < q) (hr : r ≠ 0) (hL : L ≠ 0) :
    directionCos q r Pr L ^ 2 + directionSin q r Pr L ^ 2 = 1 := by
  have hS := nullAmplitude_pos (Pr := Pr) hq hr hL
  have hs := nullAmplitude_sq (r := r) (Pr := Pr) (L := L) hq.le
  have hQ := Real.sq_sqrt hq.le
  unfold directionCos directionSin
  rw [div_pow, mul_pow, hQ, div_pow, mul_pow]
  field_simp at hs ⊢
  nlinarith [hs]

theorem momentum_direction_sin_pos {q r Pr L : ℝ}
    (hq : 0 < q) (hr : 0 < r) (hL : 0 < L) : 0 < directionSin q r Pr L := by
  exact div_pos hL (mul_pos hr (nullAmplitude_pos (Pr := Pr) hq (ne_of_gt hr) (ne_of_gt hL)))

/-- Reparametrizing the future Hamiltonian velocity gives the radial ODE
already used in DirectionFlow. -/
theorem null_radial_velocity {q w r Pr L : ℝ}
    (hq : 0 < q) (hw : 0 < w) (hr : r ≠ 0) (hL : L ≠ 0) :
    ((w * Pr - Real.sqrt w * nullAmplitude q r Pr L) + (q - w) * Pr) /
      (nullAmplitude q r Pr L / Real.sqrt w) =
      radialSpeed q w (directionCos q r Pr L) := by
  have hS := nullAmplitude_pos (Pr := Pr) hq hr hL
  have hrootw := Real.sq_sqrt hw.le
  have hrootq := Real.sq_sqrt hq.le
  unfold radialSpeed directionCos
  rw [Real.sqrt_mul hq.le]
  field_simp
  ring_nf
  rw [hrootw, hrootq]

theorem null_angular_velocity {q w r Pr L : ℝ}
    (hq : 0 < q) (hw : 0 < w) (hr : r ≠ 0) (hL : L ≠ 0) :
    (L / r ^ 2) / (nullAmplitude q r Pr L / Real.sqrt w) =
      Real.sqrt w / r * directionSin q r Pr L := by
  have hS := nullAmplitude_pos (Pr := Pr) hq hr hL
  unfold directionSin
  field_simp

noncomputable def vaidyaW (q m r : ℝ) : ℝ := q - 1 + 2 * m / r

theorem vaidyaW_pos {q m r : ℝ} (hq : 1 < q) (hm : 0 ≤ m) (hr : 0 < r) :
    0 < vaidyaW q m r := by
  unfold vaidyaW
  have hqm : 0 < q - 1 := by linarith
  positivity

theorem exterior_iff_vaidyaW_lt_charge {q m r : ℝ} (hr : 0 < r) :
    vaidyaW q m r < q ↔ 2 * m < r := by
  unfold vaidyaW
  have hh : 2 * m / r < 1 ↔ 2 * m < r := by rw [div_lt_iff₀ hr]; simp
  constructor
  · intro h; apply hh.mp; linarith
  · intro h; have hdiv := hh.mpr h; linarith

theorem receiver_timelike {q w r : ℝ} (hq : 0 < q) (hw : 0 < w) (hqw : w < q) :
    metricQuad q w r 1 0 0 < 0 := by
  unfold metricQuad
  simp only [one_pow, mul_one, zero_add, zero_pow (by norm_num : (2 : ℕ) ≠ 0), mul_zero, add_zero]
  have hdiv : w ^ 2 / q < w := (div_lt_iff₀ hq).mpr
    (by nlinarith [mul_pos hw (sub_pos.mpr hqw)])
  nlinarith

theorem time_covector_timelike {q w r : ℝ} (hw : 0 < w) :
    cometricQuad q w r 1 0 0 < 0 := by
  simp [cometricQuad]
  exact div_neg_of_neg_of_pos (by norm_num) hw

noncomputable def radialMomentumAffine (q m r Pv Pr L : ℝ) : ℝ :=
  m / r ^ 2 * (Pv ^ 2 / vaidyaW q m r ^ 2 - Pr ^ 2) + L ^ 2 / r ^ 3

/-- The radial force is derived from the Hamiltonian at fixed mass, not
assumed as a direction equation. -/
theorem hamiltonian_radius_derivative {q m r Pv Pr L : ℝ}
    (hr : r ≠ 0) (hw : vaidyaW q m r ≠ 0) :
    HasDerivAt (fun x => coordinateHamiltonian q (vaidyaW q m x) x Pv Pr L)
      (-radialMomentumAffine q m r Pv Pr L) r := by
  have hdw : HasDerivAt (fun x => vaidyaW q m x) (-2 * m / r ^ 2) r := by
    have hh := ((hasDerivAt_const r (2 * m)).div (hasDerivAt_id r) hr).const_add (q - 1)
    convert hh using 1
    dsimp [vaidyaW]
    ring
  have ha := (hasDerivAt_const r (-Pv ^ 2)).div hdw hw
  have hb := (hdw.neg.const_add q).const_mul (Pr ^ 2)
  have hc := (hasDerivAt_const r (L ^ 2)).div ((hasDerivAt_id r).pow 2) (pow_ne_zero 2 hr)
  have hh := (((ha.add (hasDerivAt_const r (2 * Pv * Pr))).add hb).add hc).const_mul (1 / 2)
  convert hh using 1
  · funext x
    dsimp [coordinateHamiltonian, cometricQuad]
    ring
  · dsimp [radialMomentumAffine]
    field_simp
    ring

noncomputable def radialMomentumV (q m r Pr L : ℝ) : ℝ :=
  m * nullAmplitude q r Pr L / (r ^ 2 * Real.sqrt (vaidyaW q m r)) -
    2 * m * Pr / r ^ 2 + Real.sqrt (vaidyaW q m r) * L ^ 2 / (r ^ 3 * nullAmplitude q r Pr L)

/-- Reduction of the radial force on the future null sheet. -/
theorem radial_momentum_future_reduction {q m r Pr L : ℝ}
    (hq : 0 < q) (hr : r ≠ 0) (hL : L ≠ 0) (hw : 0 < vaidyaW q m r) :
    radialMomentumAffine q m r
      (vaidyaW q m r * Pr - Real.sqrt (vaidyaW q m r) * nullAmplitude q r Pr L) Pr L /
      (nullAmplitude q r Pr L / Real.sqrt (vaidyaW q m r)) = radialMomentumV q m r Pr L := by
  have hS := nullAmplitude_pos (Pr := Pr) hq hr hL
  have hroot := Real.sq_sqrt hw.le
  unfold radialMomentumAffine radialMomentumV
  field_simp
  ring_nf
  have h3 : Real.sqrt (vaidyaW q m r) ^ 3 = vaidyaW q m r * Real.sqrt (vaidyaW q m r) := by
    calc _ = Real.sqrt (vaidyaW q m r) ^ 2 * Real.sqrt (vaidyaW q m r) := by ring
         _ = _ := by rw [hroot]
  have h4 : Real.sqrt (vaidyaW q m r) ^ 4 = vaidyaW q m r ^ 2 := by
    calc _ = (Real.sqrt (vaidyaW q m r) ^ 2) ^ 2 := by ring
         _ = _ := by rw [hroot]
  rw [h3, h4]
  ring

/-- The combination entering the derivative of the direction cosine has
no mass-rate term. It remains regular at Pr=0 and at radial turning points (r_v=0). -/
theorem radial_momentum_direction_combination {q m r Pr L : ℝ}
    (hq : 0 < q) (hr : r ≠ 0) (hL : L ≠ 0) (hw : 0 < vaidyaW q m r) :
    radialMomentumV q m r Pr L + Pr / r * radialSpeed q (vaidyaW q m r) (directionCos q r Pr L) =
      nullAmplitude q r Pr L / r *
        ((3 * vaidyaW q m r - q + 1) / (2 * Real.sqrt (vaidyaW q m r)) -
          Pr / nullAmplitude q r Pr L * (2 * vaidyaW q m r - q + 1)) := by
  have hS := nullAmplitude_pos (Pr := Pr) hq hr hL
  have hs := nullAmplitude_sq (r := r) (Pr := Pr) (L := L) hq.le
  have hQ := Real.sq_sqrt hq.le
  have hW := Real.sq_sqrt hw.le
  unfold radialMomentumV radialSpeed directionCos
  rw [Real.sqrt_mul hq.le]
  field_simp at hs ⊢
  ring_nf
  simp only [hQ, hW]
  have hws : vaidyaW q m r * r = (q - 1) * r + 2 * m := by
    unfold vaidyaW
    field_simp
  linear_combination
    (2 * nullAmplitude q r Pr L * r * Pr * Real.sqrt (vaidyaW q m r) -
      nullAmplitude q r Pr L ^ 2 * r) * hws - 2 * vaidyaW q m r * hs

end ViaB
