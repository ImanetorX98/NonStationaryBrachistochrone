import ViaB.DirectionFromHamiltonian
import ViaB.TimeDependentHamiltonian

namespace ViaB

/-- Recover the radial covector from a nonradial direction and a positive
angular-momentum normalization. No radial velocity is inverted. -/
noncomputable def reconstructedPr (q L r α : ℝ) : ℝ :=
  L * Real.cos α / (Real.sqrt q * r * Real.sin α)

noncomputable def reconstructedPv (q L m r α : ℝ) : ℝ :=
  vaidyaW q m r * reconstructedPr q L r α -
    Real.sqrt (vaidyaW q m r) * L / (r * Real.sin α)

noncomputable def reconstructedTimeSpeed (L w r α : ℝ) : ℝ :=
  L / (r * Real.sqrt w * Real.sin α)

theorem reconstructed_nullAmplitude {q L r α : ℝ}
    (hq : 0 < q) (hL : 0 < L) (hr : 0 < r)
    (hα : 0 < α ∧ α < Real.pi) :
    nullAmplitude q r (reconstructedPr q L r α) L = L / (r * Real.sin α) := by
  have hs : 0 < Real.sin α := Real.sin_pos_of_pos_of_lt_pi hα.1 hα.2
  have hQ := Real.sq_sqrt hq.le
  have hamp := nullAmplitude_pos (Pr := reconstructedPr q L r α)
    hq (ne_of_gt hr) (ne_of_gt hL)
  have htarget : 0 < L / (r * Real.sin α) := div_pos hL (mul_pos hr hs)
  have heq : nullAmplitude q r (reconstructedPr q L r α) L ^ 2 =
      (L / (r * Real.sin α)) ^ 2 := by
    rw [nullAmplitude_sq hq.le]
    unfold reconstructedPr
    field_simp
    rw [hQ]
    nlinarith [Real.sin_sq_add_cos_sq α]
  nlinarith

theorem reconstructed_direction {q L r α : ℝ}
    (hq : 0 < q) (hL : 0 < L) (hr : 0 < r)
    (hα : 0 < α ∧ α < Real.pi) :
    directionCos q r (reconstructedPr q L r α) L = Real.cos α ∧
    directionSin q r (reconstructedPr q L r α) L = Real.sin α := by
  have hs := Real.sin_pos_of_pos_of_lt_pi hα.1 hα.2
  have hQ := Real.sqrt_pos.2 hq
  unfold directionCos directionSin
  rw [reconstructed_nullAmplitude hq hL hr hα]
  unfold reconstructedPr
  constructor <;> field_simp

/-- The reconstructed covector is future null; the time speed is positive.
This is an algebraic reconstruction, not yet the full affine ODE. -/
theorem reconstructed_future_null {q L m r α : ℝ}
    (hq : 0 < q) (hL : 0 < L) (hr : 0 < r)
    (hα : 0 < α ∧ α < Real.pi) (hw : 0 < vaidyaW q m r) :
    cometricQuad q (vaidyaW q m r) r
      (reconstructedPv q L m r α) (reconstructedPr q L r α) L = 0 ∧
    0 < -reconstructedPv q L m r α / vaidyaW q m r + reconstructedPr q L r α := by
  have hh := future_null_sheet (Pr := reconstructedPr q L r α)
    hq hw (ne_of_gt hr) (ne_of_gt hL)
  rw [reconstructed_nullAmplitude hq hL hr hα] at hh
  simpa only [reconstructedPv, mul_div_assoc] using hh

theorem reconstructed_time_speed_identity {q L m r α : ℝ}
    (hr : r ≠ 0) (hα : Real.sin α ≠ 0) (hw : 0 < vaidyaW q m r) :
    -reconstructedPv q L m r α / vaidyaW q m r + reconstructedPr q L r α =
      reconstructedTimeSpeed L (vaidyaW q m r) r α := by
  have hW := Real.sq_sqrt hw.le
  unfold reconstructedPv reconstructedTimeSpeed
  field_simp
  ring_nf
  rw [hW]

theorem reconstructed_radial_momentum_hasDerivAt {q L t dr dα : ℝ}
    {r α : ℝ → ℝ} (hq : 0 < q) (hr : r t ≠ 0) (hs : Real.sin (α t) ≠ 0)
    (hrder : HasDerivAt r dr t) (hαder : HasDerivAt α dα t) :
    HasDerivAt (fun v => reconstructedPr q L (r v) (α v))
      (-reconstructedPr q L (r t) (α t) / r t * dr -
        L / (Real.sqrt q * r t * Real.sin (α t) ^ 2) * dα) t := by
  have hQ := Real.sqrt_pos.2 hq
  have hden : Real.sqrt q * r t * Real.sin (α t) ≠ 0 := by positivity
  have hd := (hαder.cos.const_mul L).div
    ((hrder.const_mul (Real.sqrt q)).mul hαder.sin) hden
  convert hd using 1
  dsimp [reconstructedPr]
  field_simp
  linear_combination L * dα * r t * (Real.sin_sq_add_cos_sq (α t))

/-- The angular equation recovers the reduced radial Hamilton equation.
The mass rate is not set to zero: this equation contains its instantaneous
mass value, while the separate time-momentum equation carries the mass rate. -/
theorem reconstructed_radial_momentum_dynamics {q L t : ℝ} {m r α : ℝ → ℝ}
    (hq : 0 < q) (hL : 0 < L) (hr : 0 < r t)
    (hα : 0 < α t ∧ α t < Real.pi) (hw : 0 < vaidyaW q (m t) (r t))
    (hrder : HasDerivAt r (radialSpeed q (vaidyaW q (m t) (r t)) (Real.cos (α t))) t)
    (hαder : HasDerivAt α (Real.sin (α t) / r t *
      directionB q (vaidyaW q (m t) (r t)) (Real.cos (α t))) t) :
    HasDerivAt (fun v => reconstructedPr q L (r v) (α v))
      (radialMomentumV q (m t) (r t) (reconstructedPr q L (r t) (α t)) L) t := by
  have hs := Real.sin_pos_of_pos_of_lt_pi hα.1 hα.2
  have hQ := Real.sqrt_pos.2 hq
  have hW := Real.sqrt_pos.2 hw
  have hd := reconstructed_radial_momentum_hasDerivAt (L := L) hq (ne_of_gt hr) (ne_of_gt hs)
    hrder hαder
  have hcombo := radial_momentum_direction_combination
    (Pr := reconstructedPr q L (r t) (α t)) hq (ne_of_gt hr) (ne_of_gt hL) hw
  rw [(reconstructed_direction hq hL hr hα).1,
    reconstructed_nullAmplitude hq hL hr hα] at hcombo
  have hB : -(L / (Real.sqrt q * r t * Real.sin (α t) ^ 2) *
      (Real.sin (α t) / r t * directionB q (vaidyaW q (m t) (r t)) (Real.cos (α t)))) =
      (L / (r t * Real.sin (α t))) / r t *
        ((3 * vaidyaW q (m t) (r t) - q + 1) /
          (2 * Real.sqrt (vaidyaW q (m t) (r t))) -
          reconstructedPr q L (r t) (α t) / (L / (r t * Real.sin (α t))) *
          (2 * vaidyaW q (m t) (r t) - q + 1)) := by
    unfold directionB reconstructedPr
    rw [Real.sqrt_div hq.le]
    field_simp
    ring_nf
  apply hd.congr_deriv
  rw [show -reconstructedPr q L (r t) (α t) / r t *
        radialSpeed q (vaidyaW q (m t) (r t)) (Real.cos (α t)) -
        L / (Real.sqrt q * r t * Real.sin (α t) ^ 2) *
          (Real.sin (α t) / r t * directionB q (vaidyaW q (m t) (r t)) (Real.cos (α t))) =
      -(reconstructedPr q L (r t) (α t) / r t *
        radialSpeed q (vaidyaW q (m t) (r t)) (Real.cos (α t))) +
      -(L / (Real.sqrt q * r t * Real.sin (α t) ^ 2) *
        (Real.sin (α t) / r t * directionB q (vaidyaW q (m t) (r t)) (Real.cos (α t)))) by ring,
    hB, ← hcombo]
  ring

/-- All three normalized Hamiltonian velocities/forces match the direction
variables after reconstructing the future null covector. -/
theorem reconstructed_normalized_hamiltonian_identities {q L m r α : ℝ}
    (hq : 0 < q) (hL : 0 < L) (hr : 0 < r)
    (hα : 0 < α ∧ α < Real.pi) (hw : 0 < vaidyaW q m r) :
    ((reconstructedPv q L m r α + (q - vaidyaW q m r) * reconstructedPr q L r α) /
      reconstructedTimeSpeed L (vaidyaW q m r) r α =
      radialSpeed q (vaidyaW q m r) (Real.cos α)) ∧
    ((L / r ^ 2) / reconstructedTimeSpeed L (vaidyaW q m r) r α =
      Real.sqrt (vaidyaW q m r) / r * Real.sin α) ∧
    (radialMomentumAffine q m r (reconstructedPv q L m r α) (reconstructedPr q L r α) L /
      reconstructedTimeSpeed L (vaidyaW q m r) r α =
      radialMomentumV q m r (reconstructedPr q L r α) L) := by
  have hs := Real.sin_pos_of_pos_of_lt_pi hα.1 hα.2
  have htime : reconstructedTimeSpeed L (vaidyaW q m r) r α =
      nullAmplitude q r (reconstructedPr q L r α) L / Real.sqrt (vaidyaW q m r) := by
    rw [reconstructed_nullAmplitude hq hL hr hα]
    unfold reconstructedTimeSpeed
    field_simp
  have hpv : reconstructedPv q L m r α =
      vaidyaW q m r * reconstructedPr q L r α -
      Real.sqrt (vaidyaW q m r) * nullAmplitude q r (reconstructedPr q L r α) L := by
    rw [reconstructed_nullAmplitude hq hL hr hα]
    simp only [reconstructedPv, mul_div_assoc]
  rw [htime, hpv]
  exact ⟨by rw [null_radial_velocity hq hw (ne_of_gt hr) (ne_of_gt hL),
      (reconstructed_direction hq hL hr hα).1],
    by rw [null_angular_velocity hq hw (ne_of_gt hr) (ne_of_gt hL),
      (reconstructed_direction hq hL hr hα).2],
    radial_momentum_future_reduction hq (ne_of_gt hr) (ne_of_gt hL) hw⟩

end ViaB
