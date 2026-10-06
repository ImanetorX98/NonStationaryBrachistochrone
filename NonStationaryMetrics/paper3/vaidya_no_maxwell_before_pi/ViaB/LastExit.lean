import ViaB.Excursion

open Filter Topology

namespace ViaB

theorem exists_last_exit {r : ℝ → ℝ} {a b R : ℝ} (hab : a ≤ b)
    (hcont : ContinuousOn r (Set.Icc a b)) (hra : r a < R) (hrb : R < r b) :
    ∃ d ∈ Set.Ico a b, r d = R ∧ ∀ t ∈ Set.Ioc d b, R < r t := by
  let K := Set.Icc a b ∩ r ⁻¹' {R}
  have hclosed : IsClosed K :=
    hcont.preimage_isClosed_of_isClosed isClosed_Icc isClosed_singleton
  have hcompact : IsCompact K := isCompact_Icc.of_isClosed_subset hclosed Set.inter_subset_left
  have hnonempty : K.Nonempty := by
    obtain ⟨d, hd, hrd⟩ := intermediate_value_Icc hab hcont ⟨hra.le, hrb.le⟩
    exact ⟨d, hd, hrd⟩
  obtain ⟨d, hdK, hdgreatest⟩ := hcompact.exists_isGreatest hnonempty
  have hrd : r d = R := hdK.2
  have hdb : d < b := by
    have hne : d ≠ b := by intro heq; rw [heq] at hrd; linarith
    exact lt_of_le_of_ne hdK.1.2 hne
  refine ⟨d, ⟨hdK.1.1, hdb⟩, hrd, ?_⟩
  intro t ht
  by_contra hn
  have hrt : r t ≤ R := le_of_not_gt hn
  have hcont' : ContinuousOn r (Set.Icc t b) :=
    hcont.mono (Set.Icc_subset_Icc (hdK.1.1.trans ht.1.le) le_rfl)
  obtain ⟨s, hs, hrs⟩ := intermediate_value_Icc ht.2 hcont' ⟨hrt, hrb.le⟩
  have hsK : s ∈ K := ⟨⟨hdK.1.1.trans (ht.1.le.trans hs.1), hs.2⟩, hrs⟩
  have hsd : s ≤ d := hdgreatest hsK
  exact (not_lt_of_ge (hs.1.trans hsd)) ht.1

theorem derivative_nonneg_at_right_exit {r : ℝ → ℝ} {d b R dr : ℝ}
    (hdb : d < b) (hder : HasDerivAt r dr d) (hrd : r d = R)
    (hafter : ∀ t ∈ Set.Ioo d b, R < r t) : 0 ≤ dr := by
  have ht : Tendsto (slope r d) (𝓝[>] d) (𝓝 dr) :=
    (hasDerivWithinAt_iff_tendsto_slope' (show d ∉ Set.Ioi d by simp)).mp hder.hasDerivWithinAt
  have hnonneg : ∀ᶠ t in 𝓝[>] d, 0 ≤ slope r d t := by
    filter_upwards [self_mem_nhdsWithin,
      (eventually_lt_nhds hdb).filter_mono nhdsWithin_le_nhds] with t hdt htb
    simp only [slope_def_field, hrd]
    exact (div_pos (sub_pos.mpr (hafter t ⟨hdt, htb⟩)) (sub_pos.mpr hdt)).le
  exact ge_of_tendsto ht hnonneg

def DirectionFlow.restrict {q M a : ℝ} (flow : DirectionFlow q M a) (s : ℝ) (has : a ≤ s) :
    DirectionFlow q M s where
  radius := flow.radius
  angle := flow.angle
  mass := flow.mass
  radius_pos := fun t ht => flow.radius_pos t (has.trans ht)
  angle_range := fun t ht => flow.angle_range t (has.trans ht)
  mass_range := fun t ht => flow.mass_range t (has.trans ht)
  radial_ode := fun t ht => flow.radial_ode t (has.trans ht)
  angular_ode := fun t ht => flow.angular_ode t (has.trans ht)

theorem strict_cone_entry {q M a β R κ δ : ℝ} (flow : DirectionFlow q M a)
    (hβ : 0 < β ∧ β < Real.pi / 2) (hδ : 0 < δ)
    (hmargins : ∀ r m c : ℝ, R ≤ r → 0 ≤ m → m ≤ M →
      (c ≤ 1 → directionB q (q - 1 + 2 * m / r) c ≤ -κ) ∧
      (Real.cos β ≤ c → δ ≤ radialSpeed q (q - 1 + 2 * m / r) c))
    (hra : R ≤ flow.radius a) (hαa : flow.angle a < β) :
    ∃ s > a, R < flow.radius s ∧ flow.angle s < β := by
  have hc := Real.cos_le_cos_of_nonneg_of_le_pi (flow.angle_range a le_rfl).1.le
    (show β ≤ Real.pi by linarith [Real.pi_pos]) hαa.le
  have hm := flow.mass_range a le_rfl
  have hspeed := (hmargins _ _ _ hra hm.1 hm.2).2 hc
  have hspeedpos := hδ.trans_le hspeed
  obtain ⟨cr, hcr, hr⟩ := positive_right_of_positive_derivative
    ((flow.radial_ode a le_rfl).sub_const (flow.radius a)) (by simp) hspeedpos
  have hcα : ContinuousAt (fun t => β - flow.angle t) a :=
    continuousAt_const.sub (flow.angular_ode a le_rfl).continuousAt
  obtain ⟨cα, hcαa, hα⟩ := positive_right_of_positive_value hcα (sub_pos.mpr hαa)
  obtain ⟨s, has, hs⟩ := exists_between (lt_min hcr hcαa)
  have hrs := hr s ⟨has, hs.trans_le (min_le_left _ _)⟩
  have hαs := hα s ⟨has, hs.trans_le (min_le_right _ _)⟩
  exact ⟨s, has, by linarith, by linarith⟩

theorem cone_escape_from_nonstrict_radius {q M a β R κ δ : ℝ} (flow : DirectionFlow q M a)
    (hβ : 0 < β ∧ β < Real.pi / 2) (hκ : 0 < κ) (hδ : 0 < δ)
    (hmargins : ∀ r m c : ℝ, R ≤ r → 0 ≤ m → m ≤ M →
      (c ≤ 1 → directionB q (q - 1 + 2 * m / r) c ≤ -κ) ∧
      (Real.cos β ≤ c → δ ≤ radialSpeed q (q - 1 + 2 * m / r) c))
    (hra : R ≤ flow.radius a) (hαa : flow.angle a < β) : Tendsto flow.radius atTop atTop := by
  obtain ⟨s, has, hrs, hαs⟩ := strict_cone_entry flow hβ hδ hmargins hra hαa
  exact radial_growth_tendsto_atTop hδ
    (outgoing_cone_linear_growth (flow.restrict s has.le) hβ hκ hδ hmargins hrs hαs)

end ViaB
