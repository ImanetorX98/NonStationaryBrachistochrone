import ViaB.CoordinateLeviCivita

open Filter
open scoped Topology

namespace ViaB

noncomputable def hamiltonianCoordinateAcceleration (q m μ r Pv Pr L : ℝ) : ℝ × ℝ × ℝ :=
  let w := vaidyaW q m r
  let V := -Pv / w + Pr
  let R := Pv + (q - w) * Pr
  let F := radialMomentumAffine q m r Pv Pr L
  let G := massHamiltonianGradient q m r Pv Pr
  let Wdot := 2 * μ * V / r - 2 * m * R / r ^ 2
  (F + μ * G / w + Pv * Wdot / w ^ 2,
    -μ * G + (q - w) * F - Wdot * Pr,
    -2 * L * R / r ^ 3)

/-- Differentiating raised momenta gives the coordinate acceleration.
The mass rate along the affine curve is derived by the chain rule. -/
theorem hamiltonian_raised_velocity_derivative {q L s : ℝ}
    {m μ v R Φ pv p : ℝ → ℝ}
    (h : HamiltonianLiftAt q L m μ v R Φ pv p s)
    (hm : HasDerivAt m (μ (v s)) (v s)) (hr : R s ≠ 0)
    (hw : vaidyaW q (m (v s)) (R s) ≠ 0) :
    HasDerivAt (fun t => raiseMomentum q (vaidyaW q (m (v t)) (R t)) (R t)
      (pv t, p t, L))
      (hamiltonianCoordinateAcceleration q (m (v s)) (μ (v s)) (R s) (pv s) (p s) L) s := by
  have hdw : HasDerivAt (fun t => vaidyaW q (m (v t)) (R t))
      ((2 * (μ (v s) * (-pv s / vaidyaW q (m (v s)) (R s) + p s)) * R s -
        2 * m (v s) * (pv s + (q - vaidyaW q (m (v s)) (R s)) * p s)) / (R s) ^ 2) s :=
    (((hm.comp s h.time).const_mul 2).div h.radial hr).const_add (q - 1)
  have hV := (h.time_momentum.neg.div hdw hw).add h.radial_momentum
  have hR := h.time_momentum.add ((hdw.neg.const_add q).mul h.radial_momentum)
  have hP := (hasDerivAt_const s L).div (h.radial.pow 2) (pow_ne_zero 2 hr)
  have hd := hV.prodMk (hR.prodMk hP)
  convert hd using 1
  ext <;> dsimp [hamiltonianCoordinateAcceleration]
  all_goals field_simp
  all_goals ring

/-- The Hamiltonian acceleration is minus the contraction of the connection
computed independently from the metric derivatives and the Koszul formula.
The identity holds without a null constraint. -/
theorem hamiltonian_acceleration_equals_connection {q m μ r Pv Pr L : ℝ}
    (hq : q ≠ 0) (hr : r ≠ 0) (hw : vaidyaW q m r ≠ 0) :
    hamiltonianCoordinateAcceleration q m μ r Pv Pr L =
      let u := raiseMomentum q (vaidyaW q m r) r (Pv, Pr, L)
      let C := raiseMomentum q (vaidyaW q m r) r
        (loweredConnectionContraction q m μ r u.1 u.2.1 u.2.2)
      (-C.1, -C.2.1, -C.2.2) := by
  ext <;> dsimp [hamiltonianCoordinateAcceleration, loweredConnectionContraction,
    raiseMomentum, radialMomentumAffine, massHamiltonianGradient]
  all_goals field_simp
  all_goals ring

/-- An actual Hamiltonian germ has second coordinate derivatives satisfying
the affine geodesic equation for the verified coordinate connection. -/
theorem hamiltonian_germ_coordinate_geodesic {q L s : ℝ}
    {m μ v R Φ pv p : ℝ → ℝ} (hq : q ≠ 0) (hr : R s ≠ 0)
    (hw : vaidyaW q (m (v s)) (R s) ≠ 0)
    (hm : HasDerivAt m (μ (v s)) (v s))
    (h : ∀ᶠ t in 𝓝 s, HamiltonianLiftAt q L m μ v R Φ pv p t) :
    ∀ i : CoordinateIndex,
      HasDerivAt (fun t => ![deriv v t, deriv R t, deriv Φ t] i)
        (-fermatChristoffelContraction q (m (v s)) (μ (v s)) (R s)
          (deriv v s) (deriv R s) (deriv Φ s) i) s := by
  have hs := h.self_of_nhds
  have hd := hamiltonian_raised_velocity_derivative hs hm hr hw
  have hvelocity : (fun t => (deriv v t, deriv R t, deriv Φ t)) =ᶠ[𝓝 s]
      (fun t => raiseMomentum q (vaidyaW q (m (v t)) (R t)) (R t) (pv t, p t, L)) := by
    filter_upwards [h] with t ht
    simp only [ht.time.deriv, ht.radial.deriv, ht.angular.deriv, raiseMomentum]
  have hd' := hd.congr_of_eventuallyEq hvelocity
  have ha := hamiltonian_acceleration_equals_connection
    (μ := μ (v s)) (Pv := pv s) (Pr := p s) (L := L) hq hr hw
  rw [ha] at hd'
  intro i
  have hC := fermatChristoffel_contraction_formula (m := m (v s)) (μ := μ (v s)) hq hr
    (-pv s / vaidyaW q (m (v s)) (R s) + p s)
    (pv s + (q - vaidyaW q (m (v s)) (R s)) * p s) (L / (R s) ^ 2) i
  rw [hs.time.deriv, hs.radial.deriv, hs.angular.deriv]
  rw [hC]
  fin_cases i
  · simpa only [raiseMomentum] using
      (ContinuousLinearMap.fst ℝ ℝ (ℝ × ℝ)).hasFDerivAt.comp_hasDerivAt s hd'
  · simpa only [raiseMomentum] using (ContinuousLinearMap.fst ℝ ℝ ℝ).hasFDerivAt.comp_hasDerivAt s
      ((ContinuousLinearMap.snd ℝ ℝ (ℝ × ℝ)).hasFDerivAt.comp_hasDerivAt s hd')
  · simpa only [raiseMomentum] using (ContinuousLinearMap.snd ℝ ℝ ℝ).hasFDerivAt.comp_hasDerivAt s
      ((ContinuousLinearMap.snd ℝ ℝ (ℝ × ℝ)).hasFDerivAt.comp_hasDerivAt s hd')

theorem hamiltonian_lift_null_coordinate_velocity {q L s : ℝ}
    {m μ v R Φ pv p : ℝ → ℝ} (hq : q ≠ 0) (hr : R s ≠ 0)
    (hw : vaidyaW q (m (v s)) (R s) ≠ 0)
    (h : HamiltonianLiftAt q L m μ v R Φ pv p s) :
    metricQuad q (vaidyaW q (m (v s)) (R s)) (R s)
      (deriv v s) (deriv R s) (deriv Φ s) = 0 := by
  rw [h.time.deriv, h.radial.deriv, h.angular.deriv,
    metric_of_raised_momentum hq hw hr]
  exact h.null

end ViaB
