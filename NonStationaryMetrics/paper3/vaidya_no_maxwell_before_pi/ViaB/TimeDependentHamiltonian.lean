import ViaB.NullHamiltonianDynamics
import Mathlib.Analysis.Calculus.MeanValue

namespace ViaB

noncomputable def massHamiltonianGradient (q m r Pv Pr : ℝ) : ℝ :=
  (Pv ^ 2 / vaidyaW q m r ^ 2 - Pr ^ 2) / r

/-- The affine mass-rate premise follows from M=m∘v by the chain rule. -/
theorem mass_along_affine_derivative {m v : ℝ → ℝ} {a μ V : ℝ}
    (hm : HasDerivAt m μ (v a)) (hv : HasDerivAt v V a) :
    HasDerivAt (fun x => m (v x)) (μ * V) a :=
  hm.comp a hv

/-- The mass derivative of H, with coordinates and momenta held fixed. -/
theorem hamiltonian_mass_derivative {q m r Pv Pr L : ℝ}
    (hr : r ≠ 0) (hw : vaidyaW q m r ≠ 0) :
    HasDerivAt (fun x => coordinateHamiltonian q (vaidyaW q x r) r Pv Pr L)
      (massHamiltonianGradient q m r Pv Pr) m := by
  have hdw : HasDerivAt (fun x => vaidyaW q x r) (2 / r) m := by
    convert ((hasDerivAt_id m).const_mul 2).div_const r |>.const_add (q - 1) using 1
    dsimp [vaidyaW]; ring
  have ha := (hasDerivAt_const m (-Pv ^ 2)).div hdw hw
  have hb := (hdw.neg.const_add q).const_mul (Pr ^ 2)
  have hh := (((ha.add (hasDerivAt_const m (2 * Pv * Pr))).add hb).add
    (hasDerivAt_const m (L ^ 2 / r ^ 2))).const_mul (1 / 2)
  convert hh using 1
  · funext x; dsimp [coordinateHamiltonian, cometricQuad]; ring
  · dsimp [massHamiltonianGradient]; field_simp; ring

/-- Explicit time dependence produces a time-momentum force; Pv is not
assumed conserved. -/
theorem hamiltonian_time_coordinate_derivative {q r Pv Pr L t μ : ℝ} {m : ℝ → ℝ}
    (hm : HasDerivAt m μ t) (hr : r ≠ 0) (hw : vaidyaW q (m t) r ≠ 0) :
    HasDerivAt (fun v => coordinateHamiltonian q (vaidyaW q (m v) r) r Pv Pr L)
      (μ * massHamiltonianGradient q (m t) r Pv Pr) t := by
  simpa only [mul_comm] using (hamiltonian_mass_derivative hr hw).comp t hm

/-- Chain rule for H along an arbitrary differentiable coordinate path.
Here M is the mass evaluated along the affine path, not a constant. -/
theorem hamiltonian_along_path_derivative {q a dm dr dpv dp dL : ℝ}
    {M r pv p L : ℝ → ℝ}
    (hm : HasDerivAt M dm a) (hr : HasDerivAt r dr a)
    (hpv : HasDerivAt pv dpv a) (hp : HasDerivAt p dp a)
    (hL : HasDerivAt L dL a) (hr0 : r a ≠ 0)
    (hw0 : vaidyaW q (M a) (r a) ≠ 0) :
    HasDerivAt (fun x => coordinateHamiltonian q (vaidyaW q (M x) (r x))
        (r x) (pv x) (p x) (L x))
      (massHamiltonianGradient q (M a) (r a) (pv a) (p a) * dm -
        radialMomentumAffine q (M a) (r a) (pv a) (p a) (L a) * dr +
        (-pv a / vaidyaW q (M a) (r a) + p a) * dpv +
        (pv a + (q - vaidyaW q (M a) (r a)) * p a) * dp +
        L a / (r a) ^ 2 * dL) a := by
  have hdw : HasDerivAt (fun x => vaidyaW q (M x) (r x))
      ((2 * dm * r a - 2 * M a * dr) / (r a) ^ 2) a :=
    ((hm.const_mul 2).div hr hr0).const_add (q - 1)
  have ha := (hpv.pow 2).neg.div hdw hw0
  have hb := (hpv.mul hp).const_mul 2
  have hc := (hdw.neg.const_add q).mul (hp.pow 2)
  have hd := (hL.pow 2).div (hr.pow 2) (pow_ne_zero 2 hr0)
  have hh := (((ha.add hb).add hc).add hd).const_mul (1 / 2)
  convert hh using 1
  · funext x
    dsimp [coordinateHamiltonian, cometricQuad, vaidyaW]
    ring
  · dsimp [massHamiltonianGradient, radialMomentumAffine]
    field_simp
    ring

/-- Full Hamiltonian equations cancel all terms in dH/dλ, including the
mass-rate and time-momentum terms. No null constraint is an input. -/
theorem hamiltonian_flow_energy_derivative_zero {q a μ : ℝ}
    {M r pv p L : ℝ → ℝ}
    (hr0 : r a ≠ 0) (hw0 : vaidyaW q (M a) (r a) ≠ 0)
    (hm : HasDerivAt M (μ * (-pv a / vaidyaW q (M a) (r a) + p a)) a)
    (hr : HasDerivAt r (pv a + (q - vaidyaW q (M a) (r a)) * p a) a)
    (hpv : HasDerivAt pv (-μ * massHamiltonianGradient q (M a) (r a) (pv a) (p a)) a)
    (hp : HasDerivAt p (radialMomentumAffine q (M a) (r a) (pv a) (p a) (L a)) a)
    (hL : HasDerivAt L 0 a) :
    HasDerivAt (fun x => coordinateHamiltonian q (vaidyaW q (M x) (r x))
      (r x) (pv x) (p x) (L x)) 0 a := by
  convert hamiltonian_along_path_derivative hm hr hpv hp hL hr0 hw0 using 1
  ring

/-- A zero derivative on a convex domain propagates the initial null level.
The derivative premise can be discharged by the full Hamiltonian equations. -/
theorem hamiltonian_null_level_propagates {H : ℝ → ℝ} {s : Set ℝ} {a b : ℝ}
    (hs : Convex ℝ s) (ha : a ∈ s) (hb : b ∈ s)
    (hder : ∀ x ∈ s, HasDerivAt H 0 x) (hnull : H a = 0) : H b = 0 := by
  have hh := hs.norm_image_sub_le_of_norm_hasFDerivWithin_le (C := (0 : ℝ))
    (fun x hx => (hder x hx).hasFDerivAt.hasFDerivWithinAt)
    (fun x hx => by simp) ha hb
  have heq : H b = H a := by
    simpa only [zero_mul, norm_le_zero_iff, sub_eq_zero] using hh
  exact heq.trans hnull

/-- H is independent of the angular coordinate, so its angular-coordinate
partial derivative is zero. -/
theorem hamiltonian_angular_coordinate_derivative (q m r Pv Pr L φ : ℝ) :
    HasDerivAt (fun _ : ℝ => coordinateHamiltonian q (vaidyaW q m r) r Pv Pr L) 0 φ :=
  hasDerivAt_const _ _

/-- The angular Hamilton equation L_λ=0 preserves L on a convex interval. -/
theorem angular_momentum_constant {L : ℝ → ℝ} {s : Set ℝ} {a b : ℝ}
    (hs : Convex ℝ s) (ha : a ∈ s) (hb : b ∈ s)
    (hL : ∀ x ∈ s, HasDerivAt L 0 x) : L b = L a := by
  have hd : ∀ x ∈ s, HasDerivAt (fun y => L y - L a) 0 x := by
    intro x hx
    simpa only [sub_zero] using (hL x hx).sub_const (L a)
  have hh := hamiltonian_null_level_propagates hs ha hb hd (sub_self (L a))
  exact sub_eq_zero.mp hh

/-- Nullness is propagated from the initial event by the full nonstationary
Hamiltonian equations. It is not assumed at every event. -/
theorem hamiltonian_flow_null_propagates {q a b : ℝ} {s : Set ℝ}
    {M r pv p L μ : ℝ → ℝ}
    (hs : Convex ℝ s) (ha : a ∈ s) (hb : b ∈ s)
    (hr0 : ∀ x ∈ s, r x ≠ 0)
    (hw0 : ∀ x ∈ s, vaidyaW q (M x) (r x) ≠ 0)
    (hm : ∀ x ∈ s, HasDerivAt M (μ x * (-pv x / vaidyaW q (M x) (r x) + p x)) x)
    (hr : ∀ x ∈ s, HasDerivAt r (pv x + (q - vaidyaW q (M x) (r x)) * p x) x)
    (hpv : ∀ x ∈ s, HasDerivAt pv
      (-μ x * massHamiltonianGradient q (M x) (r x) (pv x) (p x)) x)
    (hp : ∀ x ∈ s, HasDerivAt p
      (radialMomentumAffine q (M x) (r x) (pv x) (p x) (L x)) x)
    (hL : ∀ x ∈ s, HasDerivAt L 0 x)
    (hnull : cometricQuad q (vaidyaW q (M a) (r a)) (r a) (pv a) (p a) (L a) = 0) :
    cometricQuad q (vaidyaW q (M b) (r b)) (r b) (pv b) (p b) (L b) = 0 := by
  have hH : coordinateHamiltonian q (vaidyaW q (M b) (r b)) (r b) (pv b) (p b) (L b) = 0 :=
    hamiltonian_null_level_propagates
      (H := fun x => coordinateHamiltonian q (vaidyaW q (M x) (r x))
        (r x) (pv x) (p x) (L x)) (a := a) (b := b) hs ha hb
      (fun x hx => hamiltonian_flow_energy_derivative_zero (a := x) (hr0 x hx) (hw0 x hx)
        (hm x hx) (hr x hx) (hpv x hx) (hp x hx) (hL x hx))
      (by simp [coordinateHamiltonian, hnull])
  dsimp [coordinateHamiltonian] at hH
  linarith

end ViaB
