import ViaB.EscapeStability
import Mathlib.Analysis.ODE.Gronwall
import Mathlib.Analysis.Calculus.ContDiff.Bounds

open Filter Topology
open scoped NNReal

namespace ViaB

noncomputable def directionField (q m : ℝ) (x : ℝ × ℝ) : ℝ × ℝ :=
  (radialSpeed q (q - 1 + 2 * m / x.1) (Real.cos x.2),
    Real.sin x.2 / x.1 * directionB q (q - 1 + 2 * m / x.1) (Real.cos x.2))

theorem directionField_contDiffAt {q : ℝ} {z : (ℝ × ℝ) × ℝ}
    (hq : 1 < q) (hr : 0 < z.1.1) (hm : 0 ≤ z.2) :
    ContDiffAt ℝ 1 (fun y : (ℝ × ℝ) × ℝ => directionField q y.2 y.1) z := by
  have hqpos : 0 < q := by linarith
  have hqm : 0 < q - 1 := by linarith
  have hw : 0 < q - 1 + 2 * z.2 / z.1.1 := by positivity
  have hqw : 0 < q * (q - 1 + 2 * z.2 / z.1.1) := mul_pos hqpos hw
  have hdiv : 0 < q / (q - 1 + 2 * z.2 / z.1.1) := div_pos hqpos hw
  unfold directionField radialSpeed directionB
  fun_prop (disch := positivity)

/-- The state and mass field is locally Lipschitz at every regular point.
This does not itself construct a family of globally existing solutions. -/
theorem directionField_locally_lipschitz {q : ℝ} {z : (ℝ × ℝ) × ℝ}
    (hq : 1 < q) (hr : 0 < z.1.1) (hm : 0 ≤ z.2) :
    ∃ K : ℝ≥0, ∃ U ∈ 𝓝 z,
      LipschitzOnWith K (fun y : (ℝ × ℝ) × ℝ => directionField q y.2 y.1) U := by
  exact (directionField_contDiffAt hq hr hm).exists_lipschitzOnWith

/-- Along a finite interval, a derivative gap bounded by K times the state
gap yields an explicit dependence estimate for existing trajectories. -/
theorem trajectory_gap_bound {E : Type*} [NormedAddCommGroup E] [NormedSpace ℝ E]
    {f g df dg : ℝ → E} {a b : ℝ} {K : ℝ}
    (hf : ∀ t ∈ Set.Icc a b, HasDerivAt f (df t) t)
    (hg : ∀ t ∈ Set.Icc a b, HasDerivAt g (dg t) t)
    (hgap : ∀ t ∈ Set.Ico a b, ‖df t - dg t‖ ≤ K * ‖f t - g t‖) :
    ∀ t ∈ Set.Icc a b,
      dist (f t) (g t) ≤ dist (f a) (g a) * Real.exp (K * (t - a)) := by
  have hc : ContinuousOn (fun t => f t - g t) (Set.Icc a b) :=
    fun t ht => ((hf t ht).sub (hg t ht)).continuousAt.continuousWithinAt
  have hd : ∀ t ∈ Set.Ico a b,
      HasDerivWithinAt (fun t => f t - g t) (df t - dg t) (Set.Ici t) t :=
    fun t ht => ((hf t ⟨ht.1, ht.2.le⟩).sub (hg t ⟨ht.1, ht.2.le⟩)).hasDerivWithinAt
  have hh := norm_le_gronwallBound_of_norm_deriv_right_le
    (δ := ‖f a - g a‖) (K := K) (ε := 0) hc hd (le_refl _)
    (fun t ht => by simpa only [add_zero] using hgap t ht)
  intro t ht
  simpa only [gronwallBound_ε0, dist_eq_norm] using hh t ht

/-- Continuity of the initial state propagates to any fixed finite time when
nearby existing trajectories obey a common derivative-gap estimate. -/
theorem trajectory_eval_continuousAt {P E : Type*} [TopologicalSpace P]
    [NormedAddCommGroup E] [NormedSpace ℝ E]
    (f df : P → ℝ → E) (p : P) {a b : ℝ} (hab : a ≤ b) (K : ℝ)
    (hder : ∀ x t, t ∈ Set.Icc a b → HasDerivAt (f x) (df x t) t)
    (hgap : ∀ᶠ x in 𝓝 p, ∀ t ∈ Set.Ico a b,
      ‖df x t - df p t‖ ≤ K * ‖f x t - f p t‖)
    (hinitial : ContinuousAt (fun x => f x a) p) :
    ContinuousAt (fun x => f x b) p := by
  apply Metric.tendsto_nhds.2
  intro ε hε
  have hC : 0 < Real.exp (K * (b - a)) := Real.exp_pos _
  have hev := Metric.tendsto_nhds.1 hinitial (ε / Real.exp (K * (b - a))) (div_pos hε hC)
  filter_upwards [hgap, hev] with x hx hxa
  have hh := trajectory_gap_bound (fun t ht => hder x t ht)
    (fun t ht => hder p t ht) hx b ⟨hab, le_rfl⟩
  have hsmall := (lt_div_iff₀ hC).mp hxa
  exact hh.trans_lt hsmall

theorem directionFlow_state_derivative {q M a : ℝ} (flow : DirectionFlow q M a)
    {t : ℝ} (ht : a ≤ t) :
    HasDerivAt (fun v => (flow.radius v, flow.angle v))
      (directionField q (flow.mass t) (flow.radius t, flow.angle t)) t := by
  exact (flow.radial_ode t ht).prodMk (flow.angular_ode t ht)

/-- For one prescribed mass history, a common Lipschitz tube and initial
continuity imply finite-time continuity of the radius/direction pair. -/
theorem directionFlow_eval_continuousAt {P : Type*} [TopologicalSpace P]
    {q M a b : ℝ} (flows : P → DirectionFlow q M a) (p : P) (hab : a ≤ b)
    (m : ℝ → ℝ) (tube : ℝ → Set (ℝ × ℝ)) (K : ℝ≥0)
    (hm : ∀ x t, t ∈ Set.Icc a b → (flows x).mass t = m t)
    (hLip : ∀ t ∈ Set.Ico a b, LipschitzOnWith K (directionField q (m t)) (tube t))
    (hp : ∀ t ∈ Set.Ico a b, ((flows p).radius t, (flows p).angle t) ∈ tube t)
    (hnear : ∀ᶠ x in 𝓝 p, ∀ t ∈ Set.Ico a b,
      ((flows x).radius t, (flows x).angle t) ∈ tube t)
    (hinitial : ContinuousAt (fun x => ((flows x).radius a, (flows x).angle a)) p) :
    ContinuousAt (fun x => ((flows x).radius b, (flows x).angle b)) p := by
  apply trajectory_eval_continuousAt
    (f := fun x t => ((flows x).radius t, (flows x).angle t))
    (df := fun x t => directionField q (m t) ((flows x).radius t, (flows x).angle t))
    p hab (K : ℝ)
  · intro x t ht
    have hh := directionFlow_state_derivative (flows x) ht.1
    rwa [hm x t ht] at hh
  · filter_upwards [hnear] with x hx
    intro t ht
    simpa only [dist_eq_norm] using (hLip t ht).dist_le_mul _ (hx t ht) _ (hp t ht)
  · exact hinitial

/-- The finite-time dependence estimate feeds directly into permanent
escape of nearby initial data, subject to the common finite-time tube. -/
theorem cone_escape_stable_from_initial_data {P : Type*} [TopologicalSpace P]
    {q M a s β R κ δ : ℝ} (flows : P → DirectionFlow q M a) (p : P)
    (has : a ≤ s) (hβ : 0 < β ∧ β < Real.pi / 2) (hκ : 0 < κ) (hδ : 0 < δ)
    (hmargins : ∀ r m c : ℝ, R ≤ r → 0 ≤ m → m ≤ M →
      (c ≤ 1 → directionB q (q - 1 + 2 * m / r) c ≤ -κ) ∧
      (Real.cos β ≤ c → δ ≤ radialSpeed q (q - 1 + 2 * m / r) c))
    (hrp : R < (flows p).radius s) (hαp : (flows p).angle s < β)
    (m : ℝ → ℝ) (tube : ℝ → Set (ℝ × ℝ)) (K : ℝ≥0)
    (hm : ∀ x t, t ∈ Set.Icc a s → (flows x).mass t = m t)
    (hLip : ∀ t ∈ Set.Ico a s, LipschitzOnWith K (directionField q (m t)) (tube t))
    (hp : ∀ t ∈ Set.Ico a s, ((flows p).radius t, (flows p).angle t) ∈ tube t)
    (hnear : ∀ᶠ x in 𝓝 p, ∀ t ∈ Set.Ico a s,
      ((flows x).radius t, (flows x).angle t) ∈ tube t)
    (hinitial : ContinuousAt (fun x => ((flows x).radius a, (flows x).angle a)) p) :
    ∀ᶠ x in 𝓝 p,
      (∀ t ≥ s, δ * (t - s) ≤ (flows x).radius t - (flows x).radius s) ∧
      Tendsto (flows x).radius atTop atTop := by
  have hc := directionFlow_eval_continuousAt flows p has m tube K hm hLip hp hnear hinitial
  exact cone_entry_stable_in_family flows p has hβ hκ hδ hmargins hc.fst hc.snd hrp hαp

end ViaB
