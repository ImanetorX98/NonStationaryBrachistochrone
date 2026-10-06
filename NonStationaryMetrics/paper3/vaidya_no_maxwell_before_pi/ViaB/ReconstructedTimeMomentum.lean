import ViaB.DirectionMomentumReconstruction

open Filter
open scoped Topology

namespace ViaB

/-- On the null level, the radial Hamilton equations and differentiability
force the missing time-momentum equation. The null level is an algebraic
premise here, not a consequence obtained from that missing equation. -/
theorem time_momentum_force_from_null_radial_equations {q L t μ : ℝ}
    {m r pv p : ℝ → ℝ} (hr0 : r t ≠ 0) (hw0 : vaidyaW q (m t) (r t) ≠ 0)
    (hV : -pv t / vaidyaW q (m t) (r t) + p t ≠ 0)
    (hm : HasDerivAt m μ t) (hpv : DifferentiableAt ℝ pv t)
    (hr : HasDerivAt r ((pv t + (q - vaidyaW q (m t) (r t)) * p t) /
      (-pv t / vaidyaW q (m t) (r t) + p t)) t)
    (hp : HasDerivAt p (radialMomentumAffine q (m t) (r t) (pv t) (p t) L /
      (-pv t / vaidyaW q (m t) (r t) + p t)) t)
    (hnull : ∀ᶠ v in 𝓝 t,
      coordinateHamiltonian q (vaidyaW q (m v) (r v)) (r v) (pv v) (p v) L = 0) :
    HasDerivAt pv (-μ * massHamiltonianGradient q (m t) (r t) (pv t) (p t) /
      (-pv t / vaidyaW q (m t) (r t) + p t)) t := by
  have hd := hamiltonian_along_path_derivative hm hr hpv.hasDerivAt hp
    (hasDerivAt_const t L) hr0 hw0
  have hz : HasDerivAt (fun v => coordinateHamiltonian q (vaidyaW q (m v) (r v))
      (r v) (pv v) (p v) L) 0 t :=
    (hasDerivAt_const t (0 : ℝ)).congr_of_eventuallyEq hnull
  have heq := hd.unique hz
  apply hpv.hasDerivAt.congr_deriv
  apply (eq_div_iff hV).mpr
  linear_combination heq

theorem reconstructedPv_differentiableAt {q L t : ℝ} {m r α : ℝ → ℝ}
    (hq : 0 < q) (hr : r t ≠ 0) (hs : Real.sin (α t) ≠ 0)
    (hw : 0 < vaidyaW q (m t) (r t))
    (hm : DifferentiableAt ℝ m t) (hrder : DifferentiableAt ℝ r t)
    (hαder : DifferentiableAt ℝ α t) :
    DifferentiableAt ℝ (fun v => reconstructedPv q L (m v) (r v) (α v)) t := by
  have hQ : Real.sqrt q ≠ 0 := ne_of_gt (Real.sqrt_pos.2 hq)
  unfold reconstructedPv reconstructedPr vaidyaW at *
  fun_prop (disch := positivity)

/-- The reconstructed time momentum obeys the nonstationary force law.
Its derivative is derived from algebraic nullness and the already proved
radial equation; no conservation of Pv or frozen mass is assumed. -/
theorem reconstructed_time_momentum_dynamics {q L t μ : ℝ} {m r α : ℝ → ℝ}
    (hq : 0 < q) (hL : 0 < L) (hr : 0 < r t)
    (hα : 0 < α t ∧ α t < Real.pi) (hw : 0 < vaidyaW q (m t) (r t))
    (hm : HasDerivAt m μ t)
    (hrder : HasDerivAt r (radialSpeed q (vaidyaW q (m t) (r t)) (Real.cos (α t))) t)
    (hαder : HasDerivAt α (Real.sin (α t) / r t *
      directionB q (vaidyaW q (m t) (r t)) (Real.cos (α t))) t) :
    HasDerivAt (fun v => reconstructedPv q L (m v) (r v) (α v))
      (-μ * massHamiltonianGradient q (m t) (r t)
        (reconstructedPv q L (m t) (r t) (α t)) (reconstructedPr q L (r t) (α t)) /
        reconstructedTimeSpeed L (vaidyaW q (m t) (r t)) (r t) (α t)) t := by
  let p : ℝ → ℝ := fun v => reconstructedPr q L (r v) (α v)
  let pv : ℝ → ℝ := fun v => reconstructedPv q L (m v) (r v) (α v)
  have hs := Real.sin_pos_of_pos_of_lt_pi hα.1 hα.2
  have htime := reconstructed_time_speed_identity (q := q) (L := L)
    (ne_of_gt hr) (ne_of_gt hs) hw
  have hfuture := (reconstructed_future_null hq hL hr hα hw).2
  have hsheet : pv t = vaidyaW q (m t) (r t) * p t -
      Real.sqrt (vaidyaW q (m t) (r t)) * nullAmplitude q (r t) (p t) L := by
    dsimp [pv, p, reconstructedPv]
    rw [reconstructed_nullAmplitude hq hL hr hα, mul_div_assoc]
  have hamp := reconstructed_nullAmplitude hq hL hr hα
  have hVamp : -pv t / vaidyaW q (m t) (r t) + p t =
      nullAmplitude q (r t) (p t) L / Real.sqrt (vaidyaW q (m t) (r t)) := by
    rw [htime]
    dsimp [p, reconstructedTimeSpeed]
    rw [hamp]
    field_simp
  have hr' : HasDerivAt r ((pv t + (q - vaidyaW q (m t) (r t)) * p t) /
      (-pv t / vaidyaW q (m t) (r t) + p t)) t := by
    rw [hVamp, hsheet, null_radial_velocity hq hw (ne_of_gt hr) (ne_of_gt hL),
      (reconstructed_direction hq hL hr hα).1]
    exact hrder
  have hp' : HasDerivAt p (radialMomentumAffine q (m t) (r t) (pv t) (p t) L /
      (-pv t / vaidyaW q (m t) (r t) + p t)) t := by
    rw [hVamp, hsheet, radial_momentum_future_reduction hq (ne_of_gt hr) (ne_of_gt hL) hw]
    exact reconstructed_radial_momentum_dynamics hq hL hr hα hw hrder hαder
  have hcw : ContinuousAt (fun v => vaidyaW q (m v) (r v)) t :=
    continuousAt_const.add ((hm.continuousAt.const_mul 2).div hrder.continuousAt (ne_of_gt hr))
  have hnull : ∀ᶠ v in 𝓝 t,
      coordinateHamiltonian q (vaidyaW q (m v) (r v)) (r v) (pv v) (p v) L = 0 := by
    filter_upwards [continuousAt_const.eventually_lt hrder.continuousAt hr,
      continuousAt_const.eventually_lt hαder.continuousAt hα.1,
      hαder.continuousAt.eventually_lt continuousAt_const hα.2,
      continuousAt_const.eventually_lt hcw hw] with v hrv hαv hαπv hwv
    have hn := (reconstructed_future_null hq hL hrv ⟨hαv, hαπv⟩ hwv).1
    simp only [coordinateHamiltonian, p, pv, hn, mul_zero]
  have hpvd := reconstructedPv_differentiableAt (L := L) hq (ne_of_gt hr) (ne_of_gt hs) hw
    hm.differentiableAt hrder.differentiableAt hαder.differentiableAt
  have hd := time_momentum_force_from_null_radial_equations (ne_of_gt hr) (ne_of_gt hw)
    (ne_of_gt hfuture) hm hpvd hr' hp' hnull
  simpa only [htime, pv, p] using hd

end ViaB
