import ViaB.DirectionFromHamiltonian

namespace ViaB

/-- Coordinate Hamiltonian equations transform into the direction system
under a supplied differentiable time reparametrization. This theorem does
not assume the angle ODE and does not assert existence of the inverse map. -/
theorem affine_hamiltonian_to_direction {q m L Pv t : ℝ}
    {r p φ ℓ : ℝ → ℝ} (hq : 0 < q) (hr : 0 < r (ℓ t)) (hL : 0 < L)
    (hw : 0 < vaidyaW q m (r (ℓ t)))
    (hnull : cometricQuad q (vaidyaW q m (r (ℓ t))) (r (ℓ t)) Pv (p (ℓ t)) L = 0)
    (hfuture : 0 < -Pv / vaidyaW q m (r (ℓ t)) + p (ℓ t))
    (hℓ : HasDerivAt ℓ (1 / (-Pv / vaidyaW q m (r (ℓ t)) + p (ℓ t))) t)
    (hr_affine : HasDerivAt r (Pv + (q - vaidyaW q m (r (ℓ t))) * p (ℓ t)) (ℓ t))
    (hp_affine : HasDerivAt p (radialMomentumAffine q m (r (ℓ t)) Pv (p (ℓ t)) L) (ℓ t))
    (hphi_affine : HasDerivAt φ (L / (r (ℓ t)) ^ 2) (ℓ t)) :
    HasDerivAt (fun v => r (ℓ v))
      (radialSpeed q (vaidyaW q m (r (ℓ t)))
        (Real.cos (momentumAngle q (r (ℓ t)) (p (ℓ t)) L))) t ∧
    HasDerivAt (fun v => momentumAngle q (r (ℓ v)) (p (ℓ v)) L)
      (Real.sin (momentumAngle q (r (ℓ t)) (p (ℓ t)) L) / r (ℓ t) *
        directionB q (vaidyaW q m (r (ℓ t)))
          (Real.cos (momentumAngle q (r (ℓ t)) (p (ℓ t)) L))) t ∧
    HasDerivAt (fun v => φ (ℓ v))
      (Real.sqrt (vaidyaW q m (r (ℓ t))) / r (ℓ t) *
        Real.sin (momentumAngle q (r (ℓ t)) (p (ℓ t)) L)) t := by
  have hpv := future_null_sheet_unique hq hw (ne_of_gt hr) (ne_of_gt hL) hnull hfuture
  have htime : -Pv / vaidyaW q m (r (ℓ t)) + p (ℓ t) =
      nullAmplitude q (r (ℓ t)) (p (ℓ t)) L / Real.sqrt (vaidyaW q m (r (ℓ t))) := by
    rw [hpv]
    have hroot := Real.sq_sqrt hw.le
    field_simp
    ring_nf
    rw [hroot]
  have hrd : HasDerivAt (fun v => r (ℓ v))
      (radialSpeed q (vaidyaW q m (r (ℓ t))) (directionCos q (r (ℓ t)) (p (ℓ t)) L)) t := by
    apply (hr_affine.comp t hℓ).congr_deriv
    rw [htime, hpv]
    simpa only [one_div, div_eq_mul_inv, one_mul] using null_radial_velocity
      (Pr := p (ℓ t)) hq hw (ne_of_gt hr) (ne_of_gt hL)
  have hpd : HasDerivAt (fun v => p (ℓ v)) (radialMomentumV q m (r (ℓ t)) (p (ℓ t)) L) t := by
    apply (hp_affine.comp t hℓ).congr_deriv
    rw [htime, hpv]
    simpa only [one_div, div_eq_mul_inv, one_mul] using radial_momentum_future_reduction
      (Pr := p (ℓ t)) hq (ne_of_gt hr) (ne_of_gt hL) hw
  have hφd : HasDerivAt (fun v => φ (ℓ v))
      (Real.sqrt (vaidyaW q m (r (ℓ t))) / r (ℓ t) * directionSin q (r (ℓ t)) (p (ℓ t)) L) t := by
    apply (hphi_affine.comp t hℓ).congr_deriv
    rw [htime]
    simpa only [one_div, div_eq_mul_inv, one_mul] using null_angular_velocity
      (Pr := p (ℓ t)) hq hw (ne_of_gt hr) (ne_of_gt hL)
  have hangle := momentumAngle_dynamics
    (r := fun v => r (ℓ v)) (p := fun v => p (ℓ v)) (t := t) hq hr hL hw hrd hpd
  have hprops := momentumAngle_properties (Pr := p (ℓ t)) hq hr hL
  rw [← hprops.2.1] at hrd
  rw [← hprops.2.2] at hφd
  exact ⟨hrd, hangle, hφd⟩

end ViaB
