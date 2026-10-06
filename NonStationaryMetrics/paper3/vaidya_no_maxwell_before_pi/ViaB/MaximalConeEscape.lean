import ViaB.MaximalExteriorCapture

open Filter Set
open scoped Topology

namespace ViaB

/-- The cone barrier only needs a solution on the finite interval being
tested. In particular, it does not require a future-global DirectionFlow. -/
theorem outgoing_cone_invariant_on_finite_interval {q M c t β R κ δ : ℝ}
    {m : ℝ → ℝ} {f : ℝ → ℝ × ℝ}
    (hβ : 0 < β ∧ β < Real.pi / 2) (hκ : 0 < κ) (hδ : 0 < δ)
    (hode : ∀ v ∈ Icc c t, HasDerivAt f (directionField q (m v) (f v)) v)
    (hr : ∀ v ∈ Icc c t, 0 < (f v).1)
    (hα : ∀ v ∈ Icc c t, 0 < (f v).2 ∧ (f v).2 < Real.pi)
    (hm : ∀ v ∈ Icc c t, 0 ≤ m v ∧ m v ≤ M)
    (hmargins : ∀ r mass z : ℝ, R ≤ r → 0 ≤ mass → mass ≤ M →
      (z ≤ 1 → directionB q (q - 1 + 2 * mass / r) z ≤ -κ) ∧
      (Real.cos β ≤ z → δ ≤ radialSpeed q (q - 1 + 2 * mass / r) z))
    (hrc : R < (f c).1) (hαc : (f c).2 < β) :
    ∀ v ∈ Icc c t, R < (f v).1 ∧ (f v).2 < β := by
  have hb := positive_pair_invariant
    (a := c) (b := t)
    (f := fun v => (f v).1 - R) (g := fun v => β - (f v).2)
    (fun v hv => (direction_radial_derivative (hode v hv)).sub_const R)
    (fun v hv => (direction_angular_derivative (hode v hv)).const_sub β)
    (sub_pos.mpr hrc) (sub_pos.mpr hαc)
    (fun v hv hzero hother => ?_) (fun v hv hzero hother => ?_)
  · intro v hv
    obtain ⟨h₁, h₂⟩ := hb v hv
    exact ⟨by linarith, by linarith⟩
  · have hvr : (f v).1 = R := by linarith
    have hangle : (f v).2 ≤ β := by linarith
    have hz := Real.cos_le_cos_of_nonneg_of_le_pi (hα v hv).1.le
      (show β ≤ Real.pi by linarith [Real.pi_pos]) hangle
    have hmargin := (hmargins (f v).1 (m v) (Real.cos (f v).2)
      (by linarith) (hm v hv).1 (hm v hv).2).2 hz
    exact hδ.trans_le hmargin
  · have hangle : (f v).2 = β := by linarith
    have hrR : R ≤ (f v).1 := by linarith
    have hB := (hmargins (f v).1 (m v) (Real.cos (f v).2)
      hrR (hm v hv).1 (hm v hv).2).1 (Real.cos_le_one (f v).2)
    have hBneg : directionB q (q - 1 + 2 * m v / (f v).1)
        (Real.cos (f v).2) < 0 := by linarith
    have hsin : 0 < Real.sin (f v).2 :=
      Real.sin_pos_of_pos_of_lt_pi (hα v hv).1 (hα v hv).2
    exact neg_pos.mpr (mul_neg_of_pos_of_neg (div_pos hsin (hr v hv)) hBneg)

/-- Strict cone entry rules out a bounded maximal future domain. This closes
the finite-endpoint gap before applying any future-global flow theorem. -/
theorem maximal_exterior_endpoints_unbounded_after_cone_entry
    {q M a l c β R κ δ : ℝ} {m : ℝ → ℝ} {x : ℝ × ℝ}
    (hq : 1 < q) (hla : l < a)
    (hmcont : ∀ v ∈ maximalExteriorDomain q m a l x, ContinuousAt m v)
    (hmC1 : ∀ v ∈ Ici a, ContDiffAt ℝ 1 m v)
    (hmono : MonotoneOn m (Ici a)) (hmM : ∀ v ∈ Ici a, m v ≤ M)
    (s : ExteriorSegment q m a l x)
    (hc : c ∈ maximalExteriorDomain q m a l x) (hac : a ≤ c)
    (hβ : 0 < β ∧ β < Real.pi / 2) (hκ : 0 < κ) (hδ : 0 < δ)
    (hR : 2 * M < R)
    (hmargins : ∀ r mass z : ℝ, R ≤ r → 0 ≤ mass → mass ≤ M →
      (z ≤ 1 → directionB q (q - 1 + 2 * mass / r) z ≤ -κ) ∧
      (Real.cos β ≤ z → δ ≤ radialSpeed q (q - 1 + 2 * mass / r) z))
    (hrc : R < (maximalExteriorCurve q m a l x c).1)
    (hαc : (maximalExteriorCurve q m a l x c).2 < β) :
    ¬ BddAbove (range (ExteriorSegment.endpoint (q := q) (m := m)
      (a := a) (l := l) (x := x))) := by
  intro hbdd
  let b := sSup (range (ExteriorSegment.endpoint (q := q) (m := m)
    (a := a) (l := l) (x := x)))
  let f := maximalExteriorCurve q m a l x
  have hdomain : maximalExteriorDomain q m a l x = Ioo l b :=
    maximalExteriorDomain_eq_Ioo s hbdd
  have hcb : c < b := (hdomain ▸ hc).2
  have hab : a < b := hac.trans_lt hcb
  obtain ⟨h, ε, hε, heq, hodeh, hcapture⟩ :=
    maximal_exterior_finite_endpoint_is_capture hq hla hmcont hmono s hbdd
      (hmC1 b hab.le)
  have hbarrier : ∀ t ∈ Ico c b, R < (f t).1 := by
    intro t ht
    have hreg : ∀ v ∈ Icc c t,
        HasDerivAt f (directionField q (m v) (f v)) v ∧
        0 < m v ∧ 2 * m v < (f v).1 ∧ 0 < (f v).2 ∧ (f v).2 < Real.pi := by
      intro v hv
      apply maximalExteriorCurve_solves hq hla hmcont
      rw [hdomain]
      exact ⟨hla.trans_le (hac.trans hv.1), hv.2.trans_lt ht.2⟩
    have hinv := outgoing_cone_invariant_on_finite_interval hβ hκ hδ
      (fun v hv => (hreg v hv).1)
      (fun v hv => by have hh := (hreg v hv).2; linarith [hh.1, hh.2.1])
      (fun v hv => (hreg v hv).2.2.2)
      (fun v hv => ⟨(hreg v hv).2.1.le, hmM v (hac.trans hv.1)⟩)
      hmargins hrc hαc
    exact (hinv t ⟨ht.1, le_rfl⟩).1
  have hlim : Tendsto (fun v => (h v).1) (𝓝[<] b) (𝓝 (h b).1) :=
    (hodeh b ⟨hab, by linarith⟩).continuousAt.fst.tendsto.mono_left nhdsWithin_le_nhds
  have hRfinal : R ≤ (h b).1 := by
    apply ge_of_tendsto hlim
    filter_upwards [Ioo_mem_nhdsLT hcb] with v hv
    have hav : a < v := hac.trans_lt hv.1
    rw [heq ⟨hav, hv.2⟩]
    exact (hbarrier v ⟨hv.1.le, hv.2⟩).le
  have hmass := hmM b hab.le
  rw [hcapture] at hRfinal
  linarith

/-- Build the old DirectionFlow interface from the constructed maximal
curve when its admissible endpoints are unbounded. -/
theorem maximal_directionFlow_of_unbounded_endpoints {q M a l : ℝ}
    {m : ℝ → ℝ} {x : ℝ × ℝ} (hq : 1 < q) (hla : l < a)
    (hmcont : ∀ v ∈ maximalExteriorDomain q m a l x, ContinuousAt m v)
    (hmM : ∀ v ∈ Ici a, m v ≤ M) (s : ExteriorSegment q m a l x)
    (hbdd : ¬ BddAbove (range (ExteriorSegment.endpoint (q := q) (m := m)
      (a := a) (l := l) (x := x)))) :
    ∃ flow : DirectionFlow q M a,
      (flow.radius a, flow.angle a) = x ∧ flow.mass = m ∧
      (∀ t ≥ a, (flow.radius t, flow.angle t) = maximalExteriorCurve q m a l x t) ∧
      ∀ t ≥ a, 2 * flow.mass t < flow.radius t := by
  let f := maximalExteriorCurve q m a l x
  have hreg : ∀ t ≥ a,
      HasDerivAt f (directionField q (m t) (f t)) t ∧
      0 < m t ∧ 2 * m t < (f t).1 ∧ 0 < (f t).2 ∧ (f t).2 < Real.pi := by
    intro t ht
    apply maximalExteriorCurve_solves hq hla hmcont
    rw [maximalExteriorDomain_eq_Ioi hbdd]
    exact hla.trans_le ht
  let flow : DirectionFlow q M a := {
    radius := fun t => (f t).1
    angle := fun t => (f t).2
    mass := m
    radius_pos := fun t ht => by
      have hh := (hreg t ht).2
      linarith [hh.1, hh.2.1]
    angle_range := fun t ht => (hreg t ht).2.2.2
    mass_range := fun t ht => ⟨(hreg t ht).2.1.le, hmM t ht⟩
    radial_ode := fun t ht => direction_radial_derivative (hreg t ht).1
    angular_ode := fun t ht => direction_angular_derivative (hreg t ht).1 }
  refine ⟨flow, ?_, rfl, ?_, ?_⟩
  · exact maximalExteriorCurve_preserves_launch hq hla hmcont s
  · intro t _
    exact Prod.eta (f t)
  · exact fun t ht => (hreg t ht).2.2.1

/-- Cone entry produces a future-global exterior flow and radial escape.
Global existence is a conclusion, rather than a DirectionFlow premise. -/
theorem directionFlow_exists_and_escapes_after_cone_entry
    {q M a l c β R κ δ : ℝ} {m : ℝ → ℝ} {x : ℝ × ℝ}
    (hq : 1 < q) (hla : l < a)
    (hmcont : ∀ v ∈ maximalExteriorDomain q m a l x, ContinuousAt m v)
    (hmC1 : ∀ v ∈ Ici a, ContDiffAt ℝ 1 m v)
    (hmono : MonotoneOn m (Ici a)) (hmM : ∀ v ∈ Ici a, m v ≤ M)
    (s : ExteriorSegment q m a l x)
    (hc : c ∈ maximalExteriorDomain q m a l x) (hac : a ≤ c)
    (hβ : 0 < β ∧ β < Real.pi / 2) (hκ : 0 < κ) (hδ : 0 < δ)
    (hR : 2 * M < R)
    (hmargins : ∀ r mass z : ℝ, R ≤ r → 0 ≤ mass → mass ≤ M →
      (z ≤ 1 → directionB q (q - 1 + 2 * mass / r) z ≤ -κ) ∧
      (Real.cos β ≤ z → δ ≤ radialSpeed q (q - 1 + 2 * mass / r) z))
    (hrc : R < (maximalExteriorCurve q m a l x c).1)
    (hαc : (maximalExteriorCurve q m a l x c).2 < β) :
    ∃ flow : DirectionFlow q M a,
      (flow.radius a, flow.angle a) = x ∧ flow.mass = m ∧
      (∀ t ≥ a, (flow.radius t, flow.angle t) = maximalExteriorCurve q m a l x t) ∧
      (∀ t ≥ a, 2 * flow.mass t < flow.radius t) ∧
      Tendsto flow.radius atTop atTop ∧
      ∀ t ≥ c, δ * (t - c) ≤ flow.radius t - flow.radius c := by
  have hbdd := maximal_exterior_endpoints_unbounded_after_cone_entry hq hla hmcont
    hmC1 hmono hmM s hc hac hβ hκ hδ hR hmargins hrc hαc
  obtain ⟨flow, hinit, hmass, hstate, hext⟩ :=
    maximal_directionFlow_of_unbounded_endpoints hq hla hmcont hmM s hbdd
  have hrflow : R < flow.radius c := by
    have heq := congrArg Prod.fst (hstate c hac)
    simpa only [← heq] using hrc
  have hαflow : flow.angle c < β := by
    have heq := congrArg Prod.snd (hstate c hac)
    simpa only [← heq] using hαc
  have hgrowth := outgoing_cone_linear_growth (flow.restrict c hac)
    hβ hκ hδ hmargins hrflow hαflow
  exact ⟨flow, hinit, hmass, hstate, hext,
    radial_growth_tendsto_atTop hδ hgrowth, hgrowth⟩

/-- End-to-end existence and escape from a strict cone launch. The seed
segment and uniform cone margins are constructed internally. Global C1 mass
is a sufficient explicit regularity hypothesis for this packaged theorem. -/
theorem directionFlow_exists_from_strict_cone_launch {q M a β : ℝ} {m : ℝ → ℝ}
    (hq : 1 < q) (hmpos : 0 < m a) (hmC1 : ContDiff ℝ 1 m)
    (hmono : MonotoneOn m (Ici a)) (hmM : ∀ v ∈ Ici a, m v ≤ M)
    (hβ : 0 < β ∧ β < Real.pi / 2) (hcone : q - 1 < q * Real.cos β ^ 2) :
    ∃ R > 2 * M, ∀ x : ℝ × ℝ,
      R < x.1 → 0 < x.2 → x.2 < β →
      ∃ flow : DirectionFlow q M a,
        (flow.radius a, flow.angle a) = x ∧ flow.mass = m ∧
        (∀ t ≥ a, 2 * flow.mass t < flow.radius t) ∧
        Tendsto flow.radius atTop atTop := by
  have hM : 0 < M := hmpos.trans_le (hmM a (by simp))
  obtain ⟨R, hR, κ, hκ, δ, hδ, hmargins⟩ := outgoing_cone_margins hq hM hβ hcone
  refine ⟨R, hR, ?_⟩
  intro x hrx hαx hαβ
  have hext : 2 * m a < x.1 := by linarith [hmM a (by simp)]
  have hαπ : x.2 < Real.pi := by linarith [Real.pi_pos]
  obtain ⟨l, hla, ⟨s⟩⟩ :=
    exterior_segment_exists hq hmpos hext ⟨hαx, hαπ⟩ hmC1.contDiffAt
  have hmcont : ∀ v ∈ maximalExteriorDomain q m a l x, ContinuousAt m v :=
    fun _ _ => hmC1.continuous.continuousAt
  have hlaunch := maximalExteriorCurve_preserves_launch hq hla hmcont s
  have ha : a ∈ maximalExteriorDomain q m a l x := mem_iUnion.mpr ⟨s, hla, s.launch_lt⟩
  obtain ⟨flow, hinit, hmass, _, hextflow, hescape, _⟩ :=
    directionFlow_exists_and_escapes_after_cone_entry hq hla hmcont
      (fun _ _ => hmC1.contDiffAt) hmono hmM s ha le_rfl hβ hκ hδ hR hmargins
      (by simpa only [hlaunch] using hrx) (by simpa only [hlaunch] using hαβ)
  exact ⟨flow, hinit, hmass, hextflow, hescape⟩

end ViaB
