import ViaB.FlowDependence

open Filter Topology

namespace ViaB

/-- A local derivative-gap estimate suffices: a small initial gap cannot
reach the tube boundary during the prescribed finite interval. -/
theorem trajectory_stays_in_tube {E : Type*} [NormedAddCommGroup E] [NormedSpace ℝ E]
    {f g df dg : ℝ → E} {a b K ε : ℝ} (hab : a ≤ b) (hK : 0 ≤ K)
    (hf : ∀ t ∈ Set.Icc a b, HasDerivAt f (df t) t)
    (hg : ∀ t ∈ Set.Icc a b, HasDerivAt g (dg t) t)
    (hgap : ∀ t ∈ Set.Ico a b, dist (f t) (g t) ≤ ε →
      ‖df t - dg t‖ ≤ K * ‖f t - g t‖)
    (hsmall : dist (f a) (g a) * Real.exp (K * (b - a)) < ε) :
    ∀ t ∈ Set.Icc a b, dist (f t) (g t) < ε := by
  let D := fun t => dist (f t) (g t)
  have hcont : ContinuousOn D (Set.Icc a b) := by
    simpa only [D, dist_eq_norm] using continuous_norm.comp_continuousOn
      ((HasDerivAt.continuousOn hf).sub (HasDerivAt.continuousOn hg))
  have hexp : 1 ≤ Real.exp (K * (b - a)) :=
    Real.one_le_exp_iff.mpr (mul_nonneg hK (sub_nonneg.mpr hab))
  have hDa : D a < ε := by
    have hh := mul_le_mul_of_nonneg_left hexp (dist_nonneg (x := f a) (y := g a))
    dsimp [D]
    linarith
  intro t ht
  by_contra hn
  have hDt : ε ≤ D t := le_of_not_gt hn
  have hct := hcont.mono (Set.Icc_subset_Icc le_rfl ht.2)
  let Z := Set.Icc a t ∩ D ⁻¹' {ε}
  have hclosed : IsClosed Z := hct.preimage_isClosed_of_isClosed isClosed_Icc isClosed_singleton
  have hcompact : IsCompact Z := isCompact_Icc.of_isClosed_subset hclosed Set.inter_subset_left
  have hnonempty : Z.Nonempty := by
    obtain ⟨d, hd, hDd⟩ := intermediate_value_Icc ht.1 hct ⟨hDa.le, hDt⟩
    exact ⟨d, hd, hDd⟩
  obtain ⟨d, hdZ, hdleast⟩ := hcompact.exists_isLeast hnonempty
  have hDd : D d = ε := hdZ.2
  have hbefore : ∀ x ∈ Set.Ico a d, D x < ε := by
    intro x hx
    by_contra hn'
    have hDx : ε ≤ D x := le_of_not_gt hn'
    have hcx := hcont.mono (Set.Icc_subset_Icc le_rfl (hx.2.le.trans (hdZ.1.2.trans ht.2)))
    obtain ⟨y, hy, hDy⟩ := intermediate_value_Icc hx.1 hcx ⟨hDa.le, hDx⟩
    have hyZ : y ∈ Z := ⟨⟨hy.1, hy.2.trans (hx.2.le.trans hdZ.1.2)⟩, hDy⟩
    exact (not_lt_of_ge (hdleast hyZ)) (hy.2.trans_lt hx.2)
  have hdb : d ≤ b := hdZ.1.2.trans ht.2
  have hbnd := trajectory_gap_bound
    (fun x hx => hf x ⟨hx.1, hx.2.trans hdb⟩)
    (fun x hx => hg x ⟨hx.1, hx.2.trans hdb⟩)
    (fun x hx => hgap x ⟨hx.1, hx.2.trans_le hdb⟩ (hbefore x hx).le)
    d ⟨hdZ.1.1, le_rfl⟩
  have he := Real.exp_le_exp.mpr (mul_le_mul_of_nonneg_left (sub_le_sub_right hdb a) hK)
  have hm := mul_le_mul_of_nonneg_left he (dist_nonneg (x := f a) (y := g a))
  change D d ≤ _ at hbnd
  rw [hDd] at hbnd
  linarith

/-- Initial continuity plus a uniform local gap bound gives finite-time
continuity without assuming nearby trajectories remain in the tube. -/
theorem trajectory_eval_continuousAt_of_local_gap {P E : Type*} [TopologicalSpace P]
    [NormedAddCommGroup E] [NormedSpace ℝ E]
    (f df : P → ℝ → E) (p : P) {a b K ε : ℝ}
    (hab : a ≤ b) (hK : 0 ≤ K) (hε : 0 < ε)
    (hder : ∀ x t, t ∈ Set.Icc a b → HasDerivAt (f x) (df x t) t)
    (hgap : ∀ᶠ x in 𝓝 p, ∀ t ∈ Set.Ico a b,
      dist (f x t) (f p t) ≤ ε → ‖df x t - df p t‖ ≤ K * ‖f x t - f p t‖)
    (hinitial : ContinuousAt (fun x => f x a) p) :
    ContinuousAt (fun x => f x b) p := by
  have hC : 0 < Real.exp (K * (b - a)) := Real.exp_pos _
  have hev := Metric.tendsto_nhds.1 hinitial (ε / Real.exp (K * (b - a))) (div_pos hε hC)
  apply trajectory_eval_continuousAt f df p hab K hder _ hinitial
  filter_upwards [hgap, hev] with x hx hxa
  have hsmall := (lt_div_iff₀ hC).mp hxa
  have htube := trajectory_stays_in_tube hab hK (hder x) (hder p) hx hsmall
  intro t ht
  exact hx t ht (htube t ⟨ht.1, ht.2.le⟩).le

/-- Stability from initial data with a local field estimate around the
reference trajectory; no nearby-trajectory tube membership is assumed. -/
theorem cone_escape_stable_of_local_field_bound {P : Type*} [TopologicalSpace P]
    {q M a s β R κ δ ε : ℝ} (flows : P → DirectionFlow q M a) (p : P)
    (has : a ≤ s) (hβ : 0 < β ∧ β < Real.pi / 2) (hκ : 0 < κ) (hδ : 0 < δ)
    (hmargins : ∀ r m c : ℝ, R ≤ r → 0 ≤ m → m ≤ M →
      (c ≤ 1 → directionB q (q - 1 + 2 * m / r) c ≤ -κ) ∧
      (Real.cos β ≤ c → δ ≤ radialSpeed q (q - 1 + 2 * m / r) c))
    (hrp : R < (flows p).radius s) (hαp : (flows p).angle s < β)
    (m : ℝ → ℝ) (K : ℝ) (hK : 0 ≤ K) (hε : 0 < ε)
    (hm : ∀ x t, t ∈ Set.Icc a s → (flows x).mass t = m t)
    (hfield : ∀ t ∈ Set.Ico a s, ∀ y : ℝ × ℝ,
      dist y ((flows p).radius t, (flows p).angle t) ≤ ε →
      dist (directionField q (m t) y)
        (directionField q (m t) ((flows p).radius t, (flows p).angle t)) ≤
        K * dist y ((flows p).radius t, (flows p).angle t))
    (hinitial : ContinuousAt (fun x => ((flows x).radius a, (flows x).angle a)) p) :
    ∀ᶠ x in 𝓝 p,
      (∀ t ≥ s, δ * (t - s) ≤ (flows x).radius t - (flows x).radius s) ∧
      Tendsto (flows x).radius atTop atTop := by
  have hc : ContinuousAt (fun x => ((flows x).radius s, (flows x).angle s)) p := by
    apply trajectory_eval_continuousAt_of_local_gap
      (f := fun x t => ((flows x).radius t, (flows x).angle t))
      (df := fun x t => directionField q (m t) ((flows x).radius t, (flows x).angle t))
      p has hK hε
    · intro x t ht
      have hh := directionFlow_state_derivative (flows x) ht.1
      rwa [hm x t ht] at hh
    · exact Filter.Eventually.of_forall (fun x t ht hclose => by
        simpa only [dist_eq_norm] using
          hfield t ht ((flows x).radius t, (flows x).angle t) hclose)
    · exact hinitial
  exact cone_entry_stable_in_family flows p has hβ hκ hδ hmargins hc.fst hc.snd hrp hαp

end ViaB
