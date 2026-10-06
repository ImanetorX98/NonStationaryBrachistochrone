import ViaB.AffineReparametrization
import Mathlib.Analysis.Calculus.InverseFunctionTheorem.Deriv
import Mathlib.Analysis.Calculus.MeanValue

open Filter
open scoped Topology

namespace ViaB

/-- A continuous positive time velocity constructs a genuine local inverse.
Neither existence of an inverse nor its derivative is an input. -/
theorem local_time_inverse_exists {v V : ℝ → ℝ} {a : ℝ}
    (hder : ∀ᶠ x in 𝓝 a, HasDerivAt v (V x) x)
    (hcont : ContinuousAt V a) (hpos : 0 < V a) :
    ∃ ℓ : ℝ → ℝ, ℓ (v a) = a ∧
      (∀ᶠ x in 𝓝 a, ℓ (v x) = x) ∧
      (∀ᶠ t in 𝓝 (v a), v (ℓ t) = t) ∧
      HasDerivAt ℓ (1 / V a) (v a) := by
  have hs := hasStrictDerivAt_of_hasDerivAt_of_continuousAt hder hcont
  have he := hs.hasStrictFDerivAt_equiv (ne_of_gt hpos)
  refine ⟨hs.localInverse v (V a) a (ne_of_gt hpos),
    he.localInverse_apply_image, he.eventually_left_inverse,
    he.eventually_right_inverse, ?_⟩
  simpa only [one_div] using (hs.to_localInverse (ne_of_gt hpos)).hasDerivAt

noncomputable def affineTimeSpeed (q : ℝ) (m r pv p : ℝ → ℝ) (x : ℝ) : ℝ :=
  -pv x / vaidyaW q (m x) (r x) + p x

/-- Continuity of the Hamiltonian time velocity follows from continuity of
the state and mass, away from r=0 and w=0. -/
theorem affineTimeSpeed_continuousAt {q a : ℝ} {m r pv p : ℝ → ℝ}
    (hm : ContinuousAt m a) (hr : ContinuousAt r a)
    (hpv : ContinuousAt pv a) (hp : ContinuousAt p a)
    (hr0 : r a ≠ 0) (hw : vaidyaW q (m a) (r a) ≠ 0) :
    ContinuousAt (affineTimeSpeed q m r pv p) a := by
  have hcw : ContinuousAt (fun x => vaidyaW q (m x) (r x)) a :=
    continuousAt_const.add ((continuousAt_const.mul hm).div hr hr0)
  exact (hpv.neg.div hcw hw).add hp

/-- The local inverse is constructed from a regular future Hamiltonian time
equation. The mass here is evaluated along the affine trajectory. -/
theorem hamiltonian_time_inverse_exists {q a : ℝ} {v m r pv p : ℝ → ℝ}
    (hm : ContinuousAt m a) (hr : ContinuousAt r a)
    (hpv : ContinuousAt pv a) (hp : ContinuousAt p a)
    (hr0 : r a ≠ 0) (hw : vaidyaW q (m a) (r a) ≠ 0)
    (hder : ∀ᶠ x in 𝓝 a, HasDerivAt v (affineTimeSpeed q m r pv p x) x)
    (hfuture : 0 < affineTimeSpeed q m r pv p a) :
    ∃ ℓ : ℝ → ℝ, ℓ (v a) = a ∧
      (∀ᶠ x in 𝓝 a, ℓ (v x) = x) ∧
      (∀ᶠ t in 𝓝 (v a), v (ℓ t) = t) ∧
      HasDerivAt ℓ (1 / affineTimeSpeed q m r pv p a) (v a) :=
  local_time_inverse_exists hder
    (affineTimeSpeed_continuousAt hm hr hpv hp hr0 hw) hfuture

/-- Composition with the coordinate derivation removes the supplied-inverse
premise. All Hamiltonian and null hypotheses remain explicit. -/
theorem hamiltonian_local_direction_exists {q L a : ℝ}
    {v m r pv p φ : ℝ → ℝ}
    (hq : 0 < q) (hrpos : 0 < r a) (hL : 0 < L)
    (hwpos : 0 < vaidyaW q (m a) (r a))
    (hm : ContinuousAt m a) (hr : ContinuousAt r a)
    (hpv : ContinuousAt pv a) (hp : ContinuousAt p a)
    (htime : ∀ᶠ x in 𝓝 a, HasDerivAt v (affineTimeSpeed q m r pv p x) x)
    (hfuture : 0 < affineTimeSpeed q m r pv p a)
    (hnull : cometricQuad q (vaidyaW q (m a) (r a)) (r a) (pv a) (p a) L = 0)
    (hr_affine : HasDerivAt r (pv a + (q - vaidyaW q (m a) (r a)) * p a) a)
    (hp_affine : HasDerivAt p (radialMomentumAffine q (m a) (r a) (pv a) (p a) L) a)
    (hphi_affine : HasDerivAt φ (L / (r a) ^ 2) a) :
    ∃ ℓ : ℝ → ℝ, ℓ (v a) = a ∧
      (∀ᶠ x in 𝓝 a, ℓ (v x) = x) ∧
      (∀ᶠ t in 𝓝 (v a), v (ℓ t) = t) ∧
      HasDerivAt (fun t => r (ℓ t))
        (radialSpeed q (vaidyaW q (m a) (r a))
          (Real.cos (momentumAngle q (r a) (p a) L))) (v a) ∧
      HasDerivAt (fun t => momentumAngle q (r (ℓ t)) (p (ℓ t)) L)
        (Real.sin (momentumAngle q (r a) (p a) L) / r a *
          directionB q (vaidyaW q (m a) (r a))
            (Real.cos (momentumAngle q (r a) (p a) L))) (v a) ∧
      HasDerivAt (fun t => φ (ℓ t))
        (Real.sqrt (vaidyaW q (m a) (r a)) / r a *
          Real.sin (momentumAngle q (r a) (p a) L)) (v a) := by
  obtain ⟨ℓ, hbase, hleft, hright, hℓ⟩ := hamiltonian_time_inverse_exists
    hm hr hpv hp (ne_of_gt hrpos) (ne_of_gt hwpos) htime hfuture
  have hd := affine_hamiltonian_to_direction (q := q) (m := m a) (L := L)
    (Pv := pv a) (t := v a) (r := r) (p := p) (φ := φ) (ℓ := ℓ)
  simp only [hbase] at hd
  have hℓ' : HasDerivAt ℓ (1 / (-pv a / vaidyaW q (m a) (r a) + p a)) (v a) := hℓ
  exact ⟨ℓ, hbase, hleft, hright,
    hd hq hrpos hL hwpos hnull hfuture hℓ' hr_affine hp_affine hphi_affine⟩

end ViaB
