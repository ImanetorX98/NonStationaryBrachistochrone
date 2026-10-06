import ViaB.MaximalExteriorConstruction
import ViaB.AngularFiniteInvariance

open Filter Set
open scoped Topology

namespace ViaB

/-- A finite supremum of all exterior continuations is necessarily a capture
endpoint. A strictly exterior endpoint would produce a larger admissible
segment and contradict the supremum. -/
theorem maximal_exterior_finite_endpoint_is_capture {q a l : ℝ} {m : ℝ → ℝ}
    {x : ℝ × ℝ} (hq : 1 < q) (hla : l < a)
    (hmcont : ∀ t ∈ maximalExteriorDomain q m a l x, ContinuousAt m t)
    (hmono : MonotoneOn m (Ici a)) (s : ExteriorSegment q m a l x)
    (hbdd : BddAbove (range (ExteriorSegment.endpoint (q := q) (m := m)
      (a := a) (l := l) (x := x))))
    (hm : ContDiffAt ℝ 1 m
      (sSup (range (ExteriorSegment.endpoint (q := q) (m := m)
        (a := a) (l := l) (x := x))))) :
    let b := sSup (range (ExteriorSegment.endpoint (q := q) (m := m)
      (a := a) (l := l) (x := x)))
    ∃ h : ℝ → ℝ × ℝ, ∃ ε > 0,
      EqOn h (maximalExteriorCurve q m a l x) (Ioo a b) ∧
      (∀ t ∈ Ioo a (b + ε), HasDerivAt h (directionField q (m t) (h t)) t) ∧
      (h b).1 = 2 * m b := by
  let b := sSup (range (ExteriorSegment.endpoint (q := q) (m := m)
    (a := a) (l := l) (x := x)))
  let f := maximalExteriorCurve q m a l x
  have hab : a < b := s.launch_lt.trans_le (le_csSup hbdd (mem_range_self s))
  have hdomain : maximalExteriorDomain q m a l x = Ioo l b :=
    maximalExteriorDomain_eq_Ioo s hbdd
  have hreg : ∀ t ∈ Ioo l b,
      HasDerivAt f (directionField q (m t) (f t)) t ∧
      0 < m t ∧ 2 * m t < (f t).1 ∧ 0 < (f t).2 ∧ (f t).2 < Real.pi := by
    intro t ht
    exact maximalExteriorCurve_solves hq hla hmcont t (hdomain.symm ▸ ht)
  have ha : a ∈ Ioo l b := ⟨hla, hab⟩
  have hsub : Ico a b ⊆ Ioo l b := fun _ ht => ⟨hla.trans_le ht.1, ht.2⟩
  obtain ⟨h, ε, hε, heq, hode, halt⟩ := exterior_finite_endpoint_capture_or_continues
    hq hab (hreg a ha).2.1 (hmono.mono (fun _ ht => ht.1)) hm
    (fun t ht => (hreg t (hsub ht)).1)
    (fun t ht => (hreg t (hsub ht)).2.2.1) (hreg a ha).2.2.2
  refine ⟨h, ε, hε, heq, hode, ?_⟩
  rcases halt with hcapture | ⟨δ, hδ, hδε, hphysical⟩
  · exact hcapture
  · let H : ℝ → ℝ × ℝ := fun t => if a < t then h t else f t
    have heqold : ∀ t < b, H t = f t := by
      intro t htb
      change (if a < t then h t else f t) = f t
      split_ifs with hat
      · exact heq ⟨hat, htb⟩
      · rfl
    have hHinitial : H a = x := by
      simpa only [H, lt_self_iff_false, if_false] using
        maximalExteriorCurve_preserves_launch hq hla hmcont s
    have hHregular : ∀ t ∈ Ioo l (b + δ),
        HasDerivAt H (directionField q (m t) (H t)) t ∧
        0 < m t ∧ 2 * m t < (H t).1 ∧ 0 < (H t).2 ∧ (H t).2 < Real.pi := by
      intro t ht
      by_cases htb : t < b
      · have htold : t ∈ Ioo l b := ⟨ht.1, htb⟩
        have hgerm : H =ᶠ[𝓝 t] f := by
          filter_upwards [Iio_mem_nhds htb] with v hv
          exact heqold v hv
        have hd := (hreg t htold).1.congr_of_eventuallyEq hgerm
        rw [← heqold t htb] at hd
        exact ⟨hd, by simpa only [heqold t htb] using (hreg t htold).2⟩
      · have hat : a < t := hab.trans_le (le_of_not_gt htb)
        have hgerm : H =ᶠ[𝓝 t] h := by
          filter_upwards [Ioi_mem_nhds hat] with v hv
          exact if_pos hv
        have htε : t ∈ Ioo a (b + ε) := ⟨hat, by linarith [ht.2]⟩
        have hd := (hode t htε).congr_of_eventuallyEq hgerm
        have hmt : 0 < m t := (hreg a ha).2.1.trans_le
          (hmono (by simp) hat.le hat.le)
        refine ⟨?_, hmt, ?_⟩
        · simpa only [H, if_pos hat] using hd
        · simpa only [H, if_pos hat] using hphysical t ⟨hat, ht.2⟩
    let u : ExteriorSegment q m a l x :=
      ⟨b + δ, H, by linarith, hHinitial, hHregular⟩
    have hub : u.endpoint ≤ b := le_csSup hbdd (mem_range_self u)
    change b + δ ≤ b at hub
    linarith

end ViaB
