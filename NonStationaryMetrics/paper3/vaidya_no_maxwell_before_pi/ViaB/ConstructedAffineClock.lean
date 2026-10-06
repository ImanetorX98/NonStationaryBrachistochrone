import ViaB.LocalTimeInverse
import ViaB.ReconstructedTimeMomentum
import Mathlib.MeasureTheory.Integral.IntervalIntegral.FundThmCalculus

open Filter Set
open scoped Topology

namespace ViaB

/-- Construct a primitive on an interval by extending its continuous slope
using a clamp. The clamp only supplies values outside the interval. -/
theorem continuous_interval_primitive_exists {F : ℝ → ℝ} {l b c : ℝ}
    (hlb : l ≤ b) (hF : ContinuousOn F (Icc l b)) :
    ∃ P : ℝ → ℝ, P c = 0 ∧
      ∀ t ∈ Ioo l b, HasDerivAt P (F t) t := by
  let G : ℝ → ℝ := fun t => F (max l (min b t))
  have hG : Continuous G := hF.comp_continuous
    (continuous_const.max (continuous_const.min continuous_id))
    (fun t => ⟨le_max_left _ _, max_le hlb (min_le_left _ _)⟩)
  refine ⟨fun t => ∫ v in c..t, G v, by simp, ?_⟩
  intro t ht
  have hd := intervalIntegral.integral_hasDerivAt_right
    (hG.intervalIntegrable c t)
    hG.aestronglyMeasurable.stronglyMeasurableAtFilter hG.continuousAt
  simpa only [G, min_eq_right ht.2.le, max_eq_right ht.1.le] using hd

/-- The derivative of one chosen inverse holds throughout a neighborhood,
not only at the launch. Local inverse identities and strict derivatives
are transported to nearby events. -/
theorem local_inverse_derivative_on_neighborhood {F ℓ v : ℝ → ℝ} {a : ℝ}
    (hder : ∀ᶠ t in 𝓝 a, HasDerivAt ℓ (F t) t)
    (hF : ∀ᶠ t in 𝓝 a, ContinuousAt F t) (hpos : 0 < F a)
    (hva : v 0 = a) (hvc : ContinuousAt v 0)
    (hleft : ∀ᶠ t in 𝓝 a, v (ℓ t) = t)
    (hright : ∀ᶠ s in 𝓝 (0 : ℝ), ℓ (v s) = s) :
    ∀ᶠ s in 𝓝 (0 : ℝ), HasDerivAt v (1 / F (v s)) s := by
  have hnear : ∀ᶠ t in 𝓝 a,
      (∀ᶠ u in 𝓝 t, HasDerivAt ℓ (F u) u) ∧ ContinuousAt F t ∧
      0 < F t ∧ (∀ᶠ u in 𝓝 t, v (ℓ u) = u) := by
    filter_upwards [eventually_eventually_nhds.mpr hder, hF,
      continuousAt_const.eventually_lt hF.self_of_nhds hpos,
      eventually_eventually_nhds.mpr hleft] with t hd hc hp hl
    exact ⟨hd, hc, hp, hl⟩
  have hv : Tendsto v (𝓝 (0 : ℝ)) (𝓝 a) := by simpa only [hva] using hvc.tendsto
  filter_upwards [hv.eventually hnear, hright] with s hs hr
  have hstrict := hasStrictDerivAt_of_hasDerivAt_of_continuousAt hs.1 hs.2.1
  have hd := (hstrict.to_local_left_inverse (ne_of_gt hs.2.2.1) hs.2.2.2).hasDerivAt
  rw [hr] at hd
  simpa only [one_div] using hd

/-- A locally continuous positive clock rate supplies both its primitive
and a genuine local inverse. Neither clock nor inverse is a premise. -/
theorem positive_clock_with_inverse_exists {F : ℝ → ℝ} {a : ℝ}
    (hF : ∀ᶠ t in 𝓝 a, ContinuousAt F t) (hpos : 0 < F a) :
    ∃ ℓ v : ℝ → ℝ, ℓ a = 0 ∧ v 0 = a ∧
      (∀ᶠ t in 𝓝 a, HasDerivAt ℓ (F t) t) ∧
      (∀ᶠ t in 𝓝 a, v (ℓ t) = t) ∧
      (∀ᶠ s in 𝓝 (0 : ℝ), ℓ (v s) = s) ∧
      HasDerivAt v (1 / F a) 0 ∧
      (∀ᶠ s in 𝓝 (0 : ℝ), HasDerivAt v (1 / F (v s)) s) := by
  obtain ⟨ε, hε, hball⟩ := Metric.mem_nhds_iff.mp hF
  have hcont : ContinuousOn F (Icc (a - ε / 2) (a + ε / 2)) := by
    intro t ht
    apply (hball ?_).continuousWithinAt
    rw [Real.ball_eq_Ioo]
    constructor <;> linarith [ht.1, ht.2]
  obtain ⟨ℓ, hℓa, hder⟩ := continuous_interval_primitive_exists
    (c := a) (by linarith : a - ε / 2 ≤ a + ε / 2) hcont
  have hder' : ∀ᶠ t in 𝓝 a, HasDerivAt ℓ (F t) t := by
    filter_upwards [Ioo_mem_nhds (by linarith : a - ε / 2 < a)
      (by linarith : a < a + ε / 2)] with t ht
    exact hder t ht
  obtain ⟨v, hv, hleft, hright, hdv⟩ := local_time_inverse_exists hder' hF.self_of_nhds hpos
  rw [hℓa] at hv hright hdv
  exact ⟨ℓ, v, hℓa, hv, hder', hleft, hright, hdv,
    local_inverse_derivative_on_neighborhood hder' hF hpos hv hdv.continuousAt hleft hright⟩

noncomputable def directionAffineClockRate (q L : ℝ) (m r α : ℝ → ℝ) (t : ℝ) : ℝ :=
  r t * Real.sqrt (vaidyaW q (m t) (r t)) * Real.sin (α t) / L

theorem directionAffineClockRate_pos {q L t : ℝ} {m r α : ℝ → ℝ}
    (hL : 0 < L) (hr : 0 < r t) (hα : 0 < α t ∧ α t < Real.pi)
    (hw : 0 < vaidyaW q (m t) (r t)) : 0 < directionAffineClockRate q L m r α t := by
  have hs := Real.sin_pos_of_pos_of_lt_pi hα.1 hα.2
  unfold directionAffineClockRate
  positivity

theorem directionAffineClockRate_continuousAt {q L t : ℝ} {m r α : ℝ → ℝ}
    (hm : ContinuousAt m t) (hr : ContinuousAt r t) (hα : ContinuousAt α t)
    (hr0 : r t ≠ 0) : ContinuousAt (directionAffineClockRate q L m r α) t := by
  have hcw : ContinuousAt (fun v => vaidyaW q (m v) (r v)) t :=
    continuousAt_const.add ((hm.const_mul 2).div hr hr0)
  exact ((hr.mul hcw.sqrt).mul (Real.continuous_sin.continuousAt.comp hα)).div_const L

/-- The inverse clock derivative is exactly the reconstructed Hamiltonian
time speed, with the desired future orientation. -/
theorem direction_affine_clock_exists {q L a : ℝ} {m r α : ℝ → ℝ}
    (hL : 0 < L) (hr : 0 < r a) (hα : 0 < α a ∧ α a < Real.pi)
    (hw : 0 < vaidyaW q (m a) (r a))
    (hm : ∀ᶠ t in 𝓝 a, ContinuousAt m t)
    (hrc : ∀ᶠ t in 𝓝 a, ContinuousAt r t)
    (hαc : ∀ᶠ t in 𝓝 a, ContinuousAt α t) :
    ∃ ℓ v : ℝ → ℝ, ℓ a = 0 ∧ v 0 = a ∧
      (∀ᶠ t in 𝓝 a, HasDerivAt ℓ (directionAffineClockRate q L m r α t) t) ∧
      (∀ᶠ t in 𝓝 a, v (ℓ t) = t) ∧
      (∀ᶠ s in 𝓝 (0 : ℝ), ℓ (v s) = s) ∧
      HasDerivAt v (reconstructedTimeSpeed L (vaidyaW q (m a) (r a)) (r a) (α a)) 0 ∧
      (∀ᶠ s in 𝓝 (0 : ℝ), HasDerivAt v
        (reconstructedTimeSpeed L (vaidyaW q (m (v s)) (r (v s))) (r (v s)) (α (v s))) s) := by
  have hrate : ∀ᶠ t in 𝓝 a, ContinuousAt (directionAffineClockRate q L m r α) t := by
    filter_upwards [hm, hrc, hαc,
      continuousAt_const.eventually_lt hrc.self_of_nhds hr] with t hmt hrt hαt hrtpos
    exact directionAffineClockRate_continuousAt hmt hrt hαt (ne_of_gt hrtpos)
  obtain ⟨ℓ, v, hℓa, hva, hder, hleft, hright, hdv, hdvnear⟩ :=
    positive_clock_with_inverse_exists hrate (directionAffineClockRate_pos hL hr hα hw)
  have hrecip : ∀ t, 1 / directionAffineClockRate q L m r α t =
      reconstructedTimeSpeed L (vaidyaW q (m t) (r t)) (r t) (α t) := by
    intro t
    unfold directionAffineClockRate reconstructedTimeSpeed
    simp only [one_div, inv_div]
  exact ⟨ℓ, v, hℓa, hva, hder, hleft, hright,
    by simpa only [hrecip] using hdv, by simpa only [hrecip] using hdvnear⟩

end ViaB
