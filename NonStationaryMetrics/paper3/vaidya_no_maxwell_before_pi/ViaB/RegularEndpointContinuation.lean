import ViaB.FiniteEndpointLimit
import Mathlib.Analysis.Calculus.FDeriv.Extend

open Filter Set
open scoped Topology NNReal

namespace ViaB

/-- The ODE extends to a regular finite endpoint as a left derivative. -/
theorem direction_endpoint_left_derivative {q a b : ℝ} {m : ℝ → ℝ}
    {f : ℝ → ℝ × ℝ} (hq : 1 < q) (hab : a < b)
    (hr : 0 < (f b).1) (hm0 : 0 ≤ m b)
    (hm : ContinuousAt m b) (hc : ContinuousAt f b)
    (hode : ∀ t ∈ Ioo a b, HasDerivAt f (directionField q (m t) (f t)) t) :
    HasDerivWithinAt f (directionField q (m b) (f b)) (Iic b) b := by
  have houter : ContinuousAt (fun z : (ℝ × ℝ) × ℝ => directionField q z.2 z.1) (f b, m b) :=
    (directionField_contDiffAt (z := (f b, m b)) hq hr hm0).continuousAt
  have hfield : ContinuousAt (fun t => directionField q (m t) (f t)) b :=
    houter.comp (f := fun t => (f t, m t)) (hc.prodMk hm)
  have heq : (fun t => deriv f t) =ᶠ[𝓝[<] b] (fun t => directionField q (m t) (f t)) := by
    filter_upwards [Ioo_mem_nhdsLT hab] with t ht
    exact (hode t ht).deriv
  exact hasDerivWithinAt_Iic_of_tendsto_deriv
    (fun t ht => (hode t ht).differentiableAt.differentiableWithinAt)
    hc.continuousWithinAt (Ioo_mem_nhdsLT hab)
    ((tendsto_congr' heq).mpr hfield.continuousWithinAt)

/-- Backwards local uniqueness identifies a restarted solution with the
incoming path, including a merely one-sided ODE at the endpoint. -/
theorem direction_endpoint_backward_agreement {q a b : ℝ} {m : ℝ → ℝ}
    {f g : ℝ → ℝ × ℝ} (hq : 1 < q) (hab : a < b)
    (hr : 0 < (f b).1) (hm0 : 0 ≤ m b) (hm : ContinuousAt m b)
    (hc : Continuous f)
    (hf : ∀ t ∈ Ioo a b, HasDerivAt f (directionField q (m t) (f t)) t)
    (hg : ∀ᶠ t in 𝓝 b, HasDerivAt g (directionField q (m t) (g t)) t)
    (hbase : f b = g b) :
    ∃ c ∈ Ioo a b, EqOn f g (Icc c b) := by
  have hleft := direction_endpoint_left_derivative hq hab hr hm0 hm hc.continuousAt hf
  obtain ⟨K, U, hU, hLip⟩ := directionField_locally_lipschitz (z := (f b, m b)) hq hr hm0
  have hgc := hg.self_of_nhds.continuousAt
  have hfe : ∀ᶠ t in 𝓝 b, (f t, m t) ∈ U := (hc.continuousAt.prodMk hm).eventually hU
  have hge : ∀ᶠ t in 𝓝 b, (g t, m t) ∈ U := by
    have hUg : U ∈ 𝓝 (g b, m b) := by rwa [← hbase]
    exact (hgc.prodMk hm).eventually hUg
  obtain ⟨δ, hδ, hball⟩ := Metric.mem_nhds_iff.mp (hg.and (hfe.and hge))
  let d := min (δ / 2) ((b - a) / 2)
  have hd : 0 < d := lt_min (by positivity) (by linarith)
  have hdδ : d < δ := (min_le_left _ _).trans_lt (by linarith)
  have hda : d < b - a := (min_le_right _ _).trans_lt (by linarith)
  have hab' : a < b - d := by linarith
  have hball' : ∀ t ∈ Icc (b - d) b, t ∈ Metric.ball b δ := by
    intro t ht
    rw [Real.ball_eq_Ioo]
    constructor <;> linarith [ht.1, ht.2]
  have hl : ∀ t : ℝ, LipschitzOnWith K (directionField q (m t))
      {y : ℝ × ℝ | (y, m t) ∈ U} := by
    intro t
    apply LipschitzOnWith.of_dist_le_mul
    intro y hy z hz
    have hdist : dist (y, m t) (z, m t) = dist y z := by
      rw [Prod.dist_eq, dist_self, max_eq_left dist_nonneg]
    simpa only [hdist] using hLip.dist_le_mul (y, m t) hy (z, m t) hz
  refine ⟨b - d, ⟨hab', by linarith⟩, ?_⟩
  apply ODE_solution_unique_of_mem_Icc_left (fun t _ => hl t) hc.continuousOn
  · intro t ht
    by_cases htb : t = b
    · subst t; exact hleft
    · exact (hf t ⟨hab'.trans ht.1, lt_of_le_of_ne ht.2 htb⟩).hasDerivWithinAt
  · intro t ht
    exact (hball (hball' t ⟨ht.1.le, ht.2⟩)).2.1
  · intro t ht
    exact (hball (hball' t ht)).1.continuousAt.continuousWithinAt
  · intro t ht
    exact (hball (hball' t ⟨ht.1.le, ht.2⟩)).1.hasDerivWithinAt
  · intro t ht
    exact (hball (hball' t ⟨ht.1.le, ht.2⟩)).2.2
  · exact hbase

/-- A solution confined to a regular compact state-mass tube continues past
a finite endpoint. The outgoing solution is constructed and glued by backwards
uniqueness, not assumed to be an extension of the incoming trajectory. -/
theorem compact_direction_continues {q a b : ℝ} {m : ℝ → ℝ}
    {f : ℝ → ℝ × ℝ} {C : Set ((ℝ × ℝ) × ℝ)}
    (hq : 1 < q) (hab : a < b) (hC : IsCompact C)
    (hr : ∀ z ∈ C, 0 < z.1.1) (hm0 : ∀ z ∈ C, 0 ≤ z.2)
    (hm : ContDiffAt ℝ 1 m b)
    (hode : ∀ t ∈ Ioo a b, HasDerivAt f (directionField q (m t) (f t)) t)
    (htube : ∀ t ∈ Ioo a b, (f t, m t) ∈ C) :
    ∃ h : ℝ → ℝ × ℝ, ∃ ε > 0, EqOn h f (Ioo a b) ∧
      ∀ t ∈ Ioo a (b + ε), HasDerivAt h (directionField q (m t) (h t)) t := by
  obtain ⟨K, hK⟩ := compact_direction_speed_bound hq hC hr hm0
  have hLip := (convex_Ioo a b).lipschitzOnWith_of_nnnorm_hasDerivWithin_le
    (fun t ht => (hode t ht).hasDerivWithinAt)
    (fun t ht => hK (f t, m t) (htube t ht))
  obtain ⟨F, hFLip, heq⟩ := planar_lipschitz_extension hLip
  have hFode : ∀ t ∈ Ioo a b, HasDerivAt F (directionField q (m t) (F t)) t := by
    intro t ht
    have hloc : F =ᶠ[𝓝 t] f := by
      filter_upwards [isOpen_Ioo.mem_nhds ht] with u hu
      exact (heq hu).symm
    have hd := (hode t ht).congr_of_eventuallyEq hloc
    simpa only [heq ht] using hd
  have hb : b ∈ closure (Ioo a b) := by
    rw [closure_Ioo hab.ne]; exact ⟨hab.le, le_rfl⟩
  haveI : NeBot (𝓝[Ioo a b] b) := mem_closure_iff_nhdsWithin_neBot.mp hb
  have hzm : Tendsto (fun t => (F t, m t)) (𝓝[Ioo a b] b) (𝓝 (F b, m b)) :=
    (hFLip.continuous.continuousAt.prodMk hm.continuousAt).continuousWithinAt
  have hzc : (F b, m b) ∈ C := hC.isClosed.mem_of_tendsto hzm (by
    filter_upwards [self_mem_nhdsWithin] with t ht
    rw [← heq ht]; exact htube t ht)
  obtain ⟨g, hgbase, ε, hε, hgode⟩ :=
    direction_local_solution_exists hq (hr _ hzc) (hm0 _ hzc) hm
  have hge : ∀ᶠ t in 𝓝 b, HasDerivAt g (directionField q (m t) (g t)) t := by
    filter_upwards [Ioo_mem_nhds (show b - ε < b by linarith)
      (show b < b + ε by linarith)] with t ht
    exact hgode t ht
  obtain ⟨c, hc, hagree⟩ := direction_endpoint_backward_agreement hq hab
    (hr _ hzc) (hm0 _ hzc) hm.continuousAt hFLip.continuous hFode hge hgbase.symm
  let h : ℝ → ℝ × ℝ := fun t => if t < b then f t else g t
  refine ⟨h, ε, hε, ?_, ?_⟩
  · intro t ht; exact if_pos ht.2
  · intro t ht
    by_cases htb : t < b
    · have hloc : h =ᶠ[𝓝 t] f := by
        filter_upwards [Iio_mem_nhds htb] with u hu
        exact if_pos hu
      have hd := (hode t ⟨ht.1, htb⟩).congr_of_eventuallyEq hloc
      simpa only [h, if_pos htb] using hd
    · have hbt : b ≤ t := le_of_not_gt htb
      have hloc : h =ᶠ[𝓝 t] g := by
        filter_upwards [Ioo_mem_nhds (hc.2.trans_le hbt) ht.2] with u hu
        change (if u < b then f u else g u) = g u
        split_ifs with hub
        · exact (heq ⟨hc.1.trans hu.1, hub⟩).trans (hagree ⟨hu.1.le, hub.le⟩)
        · rfl
      have hd := (hgode t ⟨by linarith, ht.2⟩).congr_of_eventuallyEq hloc
      simpa only [h, if_neg htb] using hd

end ViaB
