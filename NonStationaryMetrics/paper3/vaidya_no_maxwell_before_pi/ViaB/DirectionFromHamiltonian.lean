import ViaB.NullHamiltonianDynamics
import Mathlib.Analysis.SpecialFunctions.Trigonometric.InverseDeriv

namespace ViaB

/-- Chain rule for the momentum direction cosine, with no division by Pr. -/
theorem directionCos_hasDerivAt {q L t dr dp : ℝ} {r p : ℝ → ℝ}
    (hq : 0 < q) (hr : 0 < r t) (hL : L ≠ 0)
    (hrder : HasDerivAt r dr t) (hpder : HasDerivAt p dp t) :
    HasDerivAt (fun v => directionCos q (r v) (p v) L)
      (Real.sqrt q * L ^ 2 / ((r t) ^ 2 * nullAmplitude q (r t) (p t) L ^ 3) *
        (dp + p t / r t * dr)) t := by
  have hS := nullAmplitude_pos (Pr := p t) hq (ne_of_gt hr) hL
  have hs := nullAmplitude_sq (r := r t) (Pr := p t) (L := L) hq.le
  have hbase : 0 < q * (p t) ^ 2 + L ^ 2 / (r t) ^ 2 := by
    have hterm := div_pos (sq_pos_of_ne_zero hL) (sq_pos_of_ne_zero (ne_of_gt hr))
    positivity
  have hi := ((hpder.pow 2).const_mul q).add
    ((hasDerivAt_const t (L ^ 2)).div (hrder.pow 2) (pow_ne_zero 2 (ne_of_gt hr)))
  have hd := (hpder.const_mul (Real.sqrt q)).div (hi.sqrt (ne_of_gt hbase)) (ne_of_gt hS)
  convert hd using 1
  dsimp
  simp only [pow_one]
  change _ = (Real.sqrt q * dp * nullAmplitude q (r t) (p t) L -
    Real.sqrt q * p t *
      ((q * (2 * p t * dp) + (0 * (r t) ^ 2 - L ^ 2 * (2 * r t * dr)) / ((r t) ^ 2) ^ 2) /
        (2 * nullAmplitude q (r t) (p t) L))) / nullAmplitude q (r t) (p t) L ^ 2
  field_simp at hs ⊢
  linear_combination -2 * r t * dp * hs

/-- The cosine dynamics is deduced from the two reduced Hamilton equations.
The direction ODE is a conclusion, not an input. -/
theorem directionCos_dynamics {q m L t : ℝ} {r p : ℝ → ℝ}
    (hq : 0 < q) (hr : 0 < r t) (hL : L ≠ 0) (hw : 0 < vaidyaW q m (r t))
    (hrder : HasDerivAt r (radialSpeed q (vaidyaW q m (r t)) (directionCos q (r t) (p t) L)) t)
    (hpder : HasDerivAt p (radialMomentumV q m (r t) (p t) L) t) :
    HasDerivAt (fun v => directionCos q (r v) (p v) L)
      (-(directionSin q (r t) (p t) L) ^ 2 / r t *
        directionB q (vaidyaW q m (r t)) (directionCos q (r t) (p t) L)) t := by
  have hc := directionCos_hasDerivAt hq hr hL hrder hpder
  have hS := nullAmplitude_pos (Pr := p t) hq (ne_of_gt hr) hL
  have hcombo := radial_momentum_direction_combination (Pr := p t) hq (ne_of_gt hr) hL hw
  rw [hcombo] at hc
  convert hc using 1
  unfold directionSin directionB directionCos
  rw [Real.sqrt_div hq.le]
  field_simp
  ring

/-- Recover the angle equation by the chain rule for cosine on (0,π).
No assumption about radial velocity or cos α being nonzero is used. -/
theorem angle_derivative_from_cosine {α c : ℝ → ℝ} {t dα w q r : ℝ}
    (hα : 0 < α t ∧ α t < Real.pi)
    (hder : HasDerivAt α dα t)
    (heq : (fun v => Real.cos (α v)) = c)
    (hcder : HasDerivAt c
      (-(Real.sin (α t)) ^ 2 / r * directionB q w (Real.cos (α t))) t) :
    dα = Real.sin (α t) / r * directionB q w (Real.cos (α t)) := by
  have hcos := hder.cos
  rw [heq] at hcos
  have hs : Real.sin (α t) ≠ 0 := ne_of_gt (Real.sin_pos_of_pos_of_lt_pi hα.1 hα.2)
  have hh := hcos.unique hcder
  have he : -Real.sin (α t) * dα =
      -Real.sin (α t) * (Real.sin (α t) / r * directionB q w (Real.cos (α t))) := by
    calc _ = _ := hh
         _ = _ := by ring
  exact (mul_left_cancel₀ (neg_ne_zero.mpr hs)) he

noncomputable def momentumAngle (q r Pr L : ℝ) : ℝ :=
  Real.arccos (directionCos q r Pr L)

theorem momentumAngle_properties {q r Pr L : ℝ}
    (hq : 0 < q) (hr : 0 < r) (hL : 0 < L) :
    (0 < momentumAngle q r Pr L ∧ momentumAngle q r Pr L < Real.pi) ∧
      Real.cos (momentumAngle q r Pr L) = directionCos q r Pr L ∧
      Real.sin (momentumAngle q r Pr L) = directionSin q r Pr L := by
  have hc := momentum_direction_circle (Pr := Pr) hq (ne_of_gt hr) (ne_of_gt hL)
  have hs := momentum_direction_sin_pos (Pr := Pr) hq hr hL
  have hs2 : 0 < directionSin q r Pr L ^ 2 := sq_pos_of_pos hs
  have hlow : -1 < directionCos q r Pr L := by nlinarith [sq_nonneg (directionCos q r Pr L + 1)]
  have hupp : directionCos q r Pr L < 1 := by nlinarith [sq_nonneg (directionCos q r Pr L - 1)]
  have hne : Real.arccos (directionCos q r Pr L) ≠ Real.pi := by
    intro heq
    have hh := Real.arccos_eq_pi.mp heq
    linarith
  have heq : 1 - directionCos q r Pr L ^ 2 = directionSin q r Pr L ^ 2 := by linarith
  unfold momentumAngle
  refine ⟨⟨Real.arccos_pos.mpr hupp, lt_of_le_of_ne (Real.arccos_le_pi _) hne⟩,
    Real.cos_arccos hlow.le hupp.le, ?_⟩
  rw [Real.sin_arccos, heq, Real.sqrt_sq hs.le]

/-- Constructing the angle by arccos gives its derivative and range;
the angular equation is not assumed in the Hamiltonian input. -/
theorem momentumAngle_dynamics {q m L t : ℝ} {r p : ℝ → ℝ}
    (hq : 0 < q) (hr : 0 < r t) (hL : 0 < L) (hw : 0 < vaidyaW q m (r t))
    (hrder : HasDerivAt r (radialSpeed q (vaidyaW q m (r t)) (directionCos q (r t) (p t) L)) t)
    (hpder : HasDerivAt p (radialMomentumV q m (r t) (p t) L) t) :
    HasDerivAt (fun v => momentumAngle q (r v) (p v) L)
      (Real.sin (momentumAngle q (r t) (p t) L) / r t *
        directionB q (vaidyaW q m (r t)) (Real.cos (momentumAngle q (r t) (p t) L))) t := by
  have hc := directionCos_dynamics hq hr (ne_of_gt hL) hw hrder hpder
  have hp := momentumAngle_properties (Pr := p t) hq hr hL
  have hcircle := momentum_direction_circle (Pr := p t) hq (ne_of_gt hr) (ne_of_gt hL)
  have hs := momentum_direction_sin_pos (Pr := p t) hq hr hL
  have hlow : -1 < directionCos q (r t) (p t) L := by nlinarith [sq_pos_of_pos hs]
  have hupp : directionCos q (r t) (p t) L < 1 := by nlinarith [sq_pos_of_pos hs]
  have heq : Real.sqrt (1 - directionCos q (r t) (p t) L ^ 2) = directionSin q (r t) (p t) L := by
    rw [← Real.sin_arccos]
    exact hp.2.2
  have hd := (Real.hasDerivAt_arccos (ne_of_gt hlow) (ne_of_lt hupp)).comp t hc
  convert hd using 1
  rw [hp.2.1, hp.2.2]
  simp only [heq]
  field_simp

end ViaB
