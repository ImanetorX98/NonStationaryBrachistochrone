import ViaB.IntervalUniqueness

open Filter Set
open scoped Topology

namespace ViaB

/-- An exterior nonradial solution on an open interval containing the launch.
The left endpoint is fixed; only the future endpoint varies. -/
structure ExteriorSegment (q : ℝ) (m : ℝ → ℝ) (a l : ℝ) (x : ℝ × ℝ) where
  endpoint : ℝ
  curve : ℝ → ℝ × ℝ
  launch_lt : a < endpoint
  initial : curve a = x
  regular : ∀ t ∈ Ioo l endpoint,
    HasDerivAt curve (directionField q (m t) (curve t)) t ∧
    0 < m t ∧ 2 * m t < (curve t).1 ∧ 0 < (curve t).2 ∧ (curve t).2 < Real.pi

/-- The union includes every admissible future extension, without presuming
that the set of admissible endpoints is bounded or unbounded. -/
def maximalExteriorDomain (q : ℝ) (m : ℝ → ℝ) (a l : ℝ) (x : ℝ × ℝ) : Set ℝ :=
  ⋃ s : ExteriorSegment q m a l x, Ioo l s.endpoint

noncomputable def maximalExteriorCurve (q : ℝ) (m : ℝ → ℝ)
    (a l : ℝ) (x : ℝ × ℝ) (t : ℝ) : ℝ × ℝ := by
  classical
  exact if h : ∃ s : ExteriorSegment q m a l x, t ∈ Ioo l s.endpoint
  then (Classical.choose h).curve t else x

theorem exterior_segments_agree {q a l : ℝ} {m : ℝ → ℝ} {x : ℝ × ℝ}
    (hq : 1 < q) (hla : l < a)
    (hm : ∀ t ∈ maximalExteriorDomain q m a l x, ContinuousAt m t)
    (s u : ExteriorSegment q m a l x) :
    EqOn s.curve u.curve (Ioo l (min s.endpoint u.endpoint)) := by
  have hsub : Ioo l (min s.endpoint u.endpoint) ⊆ Ioo l s.endpoint :=
    fun _ ht => ⟨ht.1, lt_of_lt_of_le ht.2 (min_le_left _ _)⟩
  have husub : Ioo l (min s.endpoint u.endpoint) ⊆ Ioo l u.endpoint :=
    fun _ ht => ⟨ht.1, lt_of_lt_of_le ht.2 (min_le_right _ _)⟩
  apply direction_solution_unique_on_interval hq ⟨hla, lt_min s.launch_lt u.launch_lt⟩
  · intro t ht
    exact (hm t (mem_iUnion.mpr ⟨s, hsub ht⟩)).continuousWithinAt
  · intro t ht
    exact (s.regular t (hsub ht)).2.1.le
  · intro t ht
    have hs := s.regular t (hsub ht)
    linarith [hs.2.1, hs.2.2.1]
  · exact fun t ht => (s.regular t (hsub ht)).1
  · exact fun t ht => (u.regular t (husub ht)).1
  · exact s.initial.trans u.initial.symm

theorem maximalExteriorDomain_isOpen (q : ℝ) (m : ℝ → ℝ) (a l : ℝ) (x : ℝ × ℝ) :
    IsOpen (maximalExteriorDomain q m a l x) :=
  isOpen_iUnion (fun _ => isOpen_Ioo)

theorem maximalExteriorDomain_isPreconnected {q a l : ℝ} {m : ℝ → ℝ} {x : ℝ × ℝ}
    (hla : l < a) : IsPreconnected (maximalExteriorDomain q m a l x) :=
  isPreconnected_iUnion ⟨a, mem_iInter.mpr (fun s => ⟨hla, s.launch_lt⟩)⟩
    (fun s => (convex_Ioo l s.endpoint).isPreconnected)

/-- The chosen curve agrees with every admissible segment, so its values do
not depend on the local choice used in its definition. -/
theorem maximalExteriorCurve_agrees {q a l : ℝ} {m : ℝ → ℝ} {x : ℝ × ℝ}
    (hq : 1 < q) (hla : l < a)
    (hm : ∀ t ∈ maximalExteriorDomain q m a l x, ContinuousAt m t)
    (s : ExteriorSegment q m a l x) :
    EqOn (maximalExteriorCurve q m a l x) s.curve (Ioo l s.endpoint) := by
  intro t ht
  have hex : ∃ u : ExteriorSegment q m a l x, t ∈ Ioo l u.endpoint := ⟨s, ht⟩
  rw [maximalExteriorCurve, dif_pos hex]
  exact exterior_segments_agree hq hla hm (Classical.choose hex) s
    ⟨ht.1, lt_min (Classical.choose_spec hex).2 ht.2⟩

/-- The union curve solves the ODE and remains exterior and nonradial on its
entire domain. Existence on this domain is constructed, not a premise. -/
theorem maximalExteriorCurve_solves {q a l : ℝ} {m : ℝ → ℝ} {x : ℝ × ℝ}
    (hq : 1 < q) (hla : l < a)
    (hm : ∀ t ∈ maximalExteriorDomain q m a l x, ContinuousAt m t) :
    ∀ t ∈ maximalExteriorDomain q m a l x,
      HasDerivAt (maximalExteriorCurve q m a l x)
        (directionField q (m t) (maximalExteriorCurve q m a l x t)) t ∧
      0 < m t ∧ 2 * m t < (maximalExteriorCurve q m a l x t).1 ∧
      0 < (maximalExteriorCurve q m a l x t).2 ∧
      (maximalExteriorCurve q m a l x t).2 < Real.pi := by
  intro t ht
  obtain ⟨s, hs⟩ := mem_iUnion.mp ht
  have heq := maximalExteriorCurve_agrees hq hla hm s
  have hgerm : maximalExteriorCurve q m a l x =ᶠ[𝓝 t] s.curve := by
    filter_upwards [isOpen_Ioo.mem_nhds hs] with v hv
    exact heq hv
  have hd := (s.regular t hs).1.congr_of_eventuallyEq hgerm
  rw [← heq hs] at hd
  exact ⟨hd, by simpa only [heq hs] using (s.regular t hs).2⟩

theorem maximalExteriorCurve_preserves_launch {q a l : ℝ} {m : ℝ → ℝ} {x : ℝ × ℝ}
    (hq : 1 < q) (hla : l < a)
    (hm : ∀ t ∈ maximalExteriorDomain q m a l x, ContinuousAt m t)
    (s : ExteriorSegment q m a l x) : maximalExteriorCurve q m a l x a = x :=
  (maximalExteriorCurve_agrees hq hla hm s ⟨hla, s.launch_lt⟩).trans s.initial

/-- Picard-Lindelof supplies a nonempty collection for a regular launch.
The fixed left endpoint is chosen inside the initial local interval. -/
theorem exterior_segment_exists {q a : ℝ} {m : ℝ → ℝ} {x : ℝ × ℝ}
    (hq : 1 < q) (hmpos : 0 < m a) (hext : 2 * m a < x.1)
    (hα : 0 < x.2 ∧ x.2 < Real.pi) (hm : ContDiffAt ℝ 1 m a) :
    ∃ l < a, Nonempty (ExteriorSegment q m a l x) := by
  obtain ⟨f, hbase, ε, hε, hode⟩ :=
    direction_exterior_local_solution_exists hq hmpos hext hα hm
  refine ⟨a - ε, by linarith, ⟨?_⟩⟩
  exact ⟨a + ε, f, by linarith, hbase, hode⟩

/-- For bounded admissible endpoints the constructed domain is exactly the
open interval ending at their supremum. -/
theorem maximalExteriorDomain_eq_Ioo {q a l : ℝ} {m : ℝ → ℝ} {x : ℝ × ℝ}
    (s : ExteriorSegment q m a l x)
    (hbdd : BddAbove (range (ExteriorSegment.endpoint (q := q) (m := m)
      (a := a) (l := l) (x := x)))) :
    maximalExteriorDomain q m a l x =
      Ioo l (sSup (range (ExteriorSegment.endpoint (q := q) (m := m)
        (a := a) (l := l) (x := x)))) := by
  ext t
  constructor
  · intro ht
    obtain ⟨u, hu⟩ := mem_iUnion.mp ht
    exact ⟨hu.1, hu.2.trans_le (le_csSup hbdd (mem_range_self u))⟩
  · intro ht
    obtain ⟨b, hb, htb⟩ := exists_lt_of_lt_csSup ⟨s.endpoint, mem_range_self s⟩ ht.2
    obtain ⟨u, rfl⟩ := hb
    exact mem_iUnion.mpr ⟨u, ht.1, htb⟩

/-- Unbounded admissible endpoints give a genuinely global future domain.
The conclusion does not assume a future global solution. -/
theorem maximalExteriorDomain_eq_Ioi {q a l : ℝ} {m : ℝ → ℝ} {x : ℝ × ℝ}
    (hbdd : ¬ BddAbove (range (ExteriorSegment.endpoint (q := q) (m := m)
      (a := a) (l := l) (x := x)))) :
    maximalExteriorDomain q m a l x = Ioi l := by
  ext t
  constructor
  · intro ht
    obtain ⟨u, hu⟩ := mem_iUnion.mp ht
    exact hu.1
  · intro ht
    have hex : ∃ u : ExteriorSegment q m a l x, t < u.endpoint := by
      by_contra h
      apply hbdd
      refine ⟨t, ?_⟩
      rintro b ⟨u, rfl⟩
      exact le_of_not_gt (fun hut => h ⟨u, hut⟩)
    obtain ⟨u, hu⟩ := hex
    exact mem_iUnion.mpr ⟨u, ht, hu⟩

/-- Every admissible future segment is already included in the union domain.
This is the precise maximality property used to rule out a finite regular
exterior endpoint. -/
theorem exterior_segment_domain_subset_maximal {q a l : ℝ} {m : ℝ → ℝ}
    {x : ℝ × ℝ} (s : ExteriorSegment q m a l x) :
    Ioo l s.endpoint ⊆ maximalExteriorDomain q m a l x :=
  fun _ ht => mem_iUnion.mpr ⟨s, ht⟩

theorem exterior_global_solution_of_unbounded_endpoints {q a l : ℝ} {m : ℝ → ℝ}
    {x : ℝ × ℝ} (hq : 1 < q) (hla : l < a)
    (hm : ∀ t ∈ maximalExteriorDomain q m a l x, ContinuousAt m t)
    (s : ExteriorSegment q m a l x)
    (hbdd : ¬ BddAbove (range (ExteriorSegment.endpoint (q := q) (m := m)
      (a := a) (l := l) (x := x)))) :
    ∃ f : ℝ → ℝ × ℝ, f a = x ∧ ∀ t ∈ Ici a,
      HasDerivAt f (directionField q (m t) (f t)) t ∧
      0 < m t ∧ 2 * m t < (f t).1 ∧ 0 < (f t).2 ∧ (f t).2 < Real.pi := by
  refine ⟨maximalExteriorCurve q m a l x,
    maximalExteriorCurve_preserves_launch hq hla hm s, ?_⟩
  intro t ht
  apply maximalExteriorCurve_solves hq hla hm
  rw [maximalExteriorDomain_eq_Ioi hbdd]
  exact hla.trans_le ht

end ViaB
