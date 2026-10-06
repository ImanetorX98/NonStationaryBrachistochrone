import ViaB.FirstZero
import Mathlib.Topology.Order.Compact
import Mathlib.Topology.Order.IntermediateValue

open Filter Topology

namespace ViaB

theorem positive_right_of_positive_derivative {z : ℝ → ℝ} {a dz : ℝ}
    (hz : HasDerivAt z dz a) (hza : z a = 0) (hd : 0 < dz) :
    ∃ c > a, ∀ x ∈ Set.Ioo a c, 0 < z x := by
  have ht : Tendsto (slope z a) (𝓝[>] a) (𝓝 dz) :=
    (hasDerivWithinAt_iff_tendsto_slope' (show a ∉ Set.Ioi a by simp)).mp
      hz.hasDerivWithinAt
  have hp : ∀ᶠ x in 𝓝[>] a, 0 < z x := by
    filter_upwards [ht.eventually (eventually_gt_nhds hd), self_mem_nhdsWithin] with x hs hx
    have hdiv : 0 < z x / (x - a) := by simpa only [slope_def_field, hza, sub_zero] using hs
    have hm := mul_pos hdiv (sub_pos.mpr hx)
    rwa [div_mul_cancel₀ _ (ne_of_gt (sub_pos.mpr hx))] at hm
  obtain ⟨c, hac, hc⟩ := mem_nhdsGT_iff_exists_Ioo_subset.mp hp
  exact ⟨c, hac, hc⟩

/-- Construct the first zero, rather than assuming that such a first zero exists. -/
theorem exists_first_zero_of_nonpos {z : ℝ → ℝ} {a b c₀ : ℝ}
    (hab : a < b) (hac₀ : a < c₀) (hcont : ContinuousOn z (Set.Icc a b))
    (hstart : ∀ x ∈ Set.Ioo a c₀, 0 < z x) (hzb : z b ≤ 0) :
    ∃ d ∈ Set.Ioc a b, z d = 0 ∧ ∀ x ∈ Set.Ioo a d, 0 < z x := by
  obtain ⟨c, hac, hcmin⟩ := exists_between (lt_min hac₀ hab)
  have hcc₀ : c < c₀ := hcmin.trans_le (min_le_left _ _)
  have hcb : c < b := hcmin.trans_le (min_le_right _ _)
  have hzc : 0 < z c := hstart c ⟨hac, hcc₀⟩
  have hcontcb : ContinuousOn z (Set.Icc c b) :=
    hcont.mono (Set.Icc_subset_Icc hac.le le_rfl)
  let K := Set.Icc c b ∩ z ⁻¹' {0}
  have hKclosed : IsClosed K :=
    hcontcb.preimage_isClosed_of_isClosed isClosed_Icc (isClosed_singleton)
  have hKcompact : IsCompact K :=
    isCompact_Icc.of_isClosed_subset hKclosed Set.inter_subset_left
  have hKnonempty : K.Nonempty := by
    obtain ⟨s, hs, hsz⟩ := intermediate_value_Icc' hcb.le hcontcb ⟨hzb, hzc.le⟩
    exact ⟨s, hs, hsz⟩
  obtain ⟨d, hdK, hdleast⟩ := hKcompact.exists_isLeast hKnonempty
  have hzd : z d = 0 := hdK.2
  have hcd : c < d := by
    have hne : c ≠ d := by intro heq; rw [← heq] at hzd; linarith
    exact lt_of_le_of_ne hdK.1.1 hne
  refine ⟨d, ⟨hac.trans hcd, hdK.1.2⟩, hzd, ?_⟩
  intro x hx
  by_cases hxc : x < c
  · exact hstart x ⟨hx.1, hxc.trans hcc₀⟩
  · by_contra hnot
    have hzx : z x ≤ 0 := le_of_not_gt hnot
    have hcx : c ≤ x := le_of_not_gt hxc
    have hcontcx : ContinuousOn z (Set.Icc c x) :=
      hcont.mono (Set.Icc_subset_Icc hac.le (hx.2.le.trans hdK.1.2))
    obtain ⟨s, hs, hsz⟩ := intermediate_value_Icc' hcx hcontcx ⟨hzx, hzc.le⟩
    have hsK : s ∈ K := ⟨⟨hs.1, hs.2.trans (hx.2.le.trans hdK.1.2)⟩, hsz⟩
    have hds : d ≤ s := hdleast hsK
    exact (not_lt_of_ge (hds.trans hs.2)) hx.2

/-- Full scalar Sturm comparison before the first zero of the reference mode.
Unlike `no_first_zero`, this theorem does not assume positivity of z anywhere:
it derives it from the launch derivative and the two ODEs. -/
theorem sturm_positive {u du z dz kperp kpar : ℝ → ℝ} {a b : ℝ}
    (hab : a < b)
    (hu : ∀ x ∈ Set.Icc a b, HasDerivAt u (du x) x)
    (hz : ∀ x ∈ Set.Icc a b, HasDerivAt z (dz x) x)
    (hdu : ∀ x ∈ Set.Icc a b, HasDerivAt du (-kperp x * u x) x)
    (hdz : ∀ x ∈ Set.Icc a b, HasDerivAt dz (-kpar x * z x) x)
    (hgap : ∀ x ∈ Set.Ioo a b, 0 < kperp x - kpar x)
    (hupos : ∀ x ∈ Set.Ioc a b, 0 < u x)
    (hua : u a = 0) (hza : z a = 0) (hdza : 0 < dz a) :
    ∀ x ∈ Set.Ioc a b, 0 < z x := by
  have ha : a ∈ Set.Icc a b := ⟨le_rfl, hab.le⟩
  obtain ⟨c₀, hac₀, hstart⟩ := positive_right_of_positive_derivative (hz a ha) hza hdza
  intro x hx
  by_contra hnpos
  have hsub : Set.Icc a x ⊆ Set.Icc a b := Set.Icc_subset_Icc le_rfl hx.2
  have hcont : ContinuousOn z (Set.Icc a x) :=
    fun t ht => (hz t (hsub ht)).continuousAt.continuousWithinAt
  obtain ⟨d, hd, hzd, hbefore⟩ :=
    exists_first_zero_of_nonpos hx.1 hac₀ hcont hstart (le_of_not_gt hnpos)
  have hsubd : Set.Icc a d ⊆ Set.Icc a b := Set.Icc_subset_Icc le_rfl (hd.2.trans hx.2)
  have hno := no_first_zero hd.1
    (fun t ht => hu t (hsubd ht)) (fun t ht => hz t (hsubd ht))
    (fun t ht => hdu t (hsubd ht)) (fun t ht => hdz t (hsubd ht))
    (fun t ht => hgap t ⟨ht.1, ht.2.trans_le (hd.2.trans hx.2)⟩)
    (fun t ht => hupos t ⟨ht.1, ht.2.le.trans (hd.2.trans hx.2)⟩)
    hbefore hua hza (hupos d ⟨hd.1, hd.2.trans hx.2⟩)
  exact hno hzd

/-- At the first zero of the reference mode the comparison mode stays positive,
provided the reference crosses with negative derivative. This is the scalar
ingredient used at the antipode; no geometric identification is assumed here. -/
theorem sturm_positive_at_reference_zero {u du z dz kperp kpar : ℝ → ℝ} {a b : ℝ}
    (hab : a < b)
    (hu : ∀ x ∈ Set.Icc a b, HasDerivAt u (du x) x)
    (hz : ∀ x ∈ Set.Icc a b, HasDerivAt z (dz x) x)
    (hdu : ∀ x ∈ Set.Icc a b, HasDerivAt du (-kperp x * u x) x)
    (hdz : ∀ x ∈ Set.Icc a b, HasDerivAt dz (-kpar x * z x) x)
    (hgap : ∀ x ∈ Set.Ioo a b, 0 < kperp x - kpar x)
    (hupos : ∀ x ∈ Set.Ioo a b, 0 < u x)
    (hua : u a = 0) (hza : z a = 0) (hdza : 0 < dz a)
    (hub : u b = 0) (hdub : du b < 0) : 0 < z b := by
  have hzpos : ∀ x ∈ Set.Ioo a b, 0 < z x := by
    intro x hx
    have hsub : Set.Icc a x ⊆ Set.Icc a b := Set.Icc_subset_Icc le_rfl hx.2.le
    have hpositive := sturm_positive hx.1
      (fun t ht => hu t (hsub ht)) (fun t ht => hz t (hsub ht))
      (fun t ht => hdu t (hsub ht)) (fun t ht => hdz t (hsub ht))
      (fun t ht => hgap t ⟨ht.1, ht.2.trans hx.2⟩)
      (fun t ht => hupos t ⟨ht.1, ht.2.trans_lt hx.2⟩)
      hua hza hdza
    exact hpositive x ⟨hx.1, le_rfl⟩
  have hmono := wronskian_strictMonoOn hu hz hdu hdz hgap hupos hzpos
  have hwa : wronskian u du z dz a = 0 := by simp [wronskian, hua, hza]
  have hwb : 0 < wronskian u du z dz b := by
    have hlt := hmono (show a ∈ Set.Icc a b from ⟨le_rfl, hab.le⟩)
      (show b ∈ Set.Icc a b from ⟨hab.le, le_rfl⟩) hab
    rwa [hwa] at hlt
  simp only [wronskian, hub, mul_zero, zero_sub] at hwb
  by_contra hnpos
  have hprod := mul_nonpos_of_nonpos_of_nonneg (le_of_not_gt hnpos) (le_of_lt (neg_pos.mpr hdub))
  nlinarith

end ViaB
