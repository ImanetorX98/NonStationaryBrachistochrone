import ViaB.LocalDirectionExistence
import Mathlib.Topology.Connected.Clopen

open Filter Set
open scoped Topology

namespace ViaB

/-- Local uniqueness propagates over an open connected regular time domain.
No uniform Lipschitz constant on the entire domain is required. -/
theorem direction_solution_unique_on_open_connected {q a : ℝ} {m : ℝ → ℝ}
    {f g : ℝ → ℝ × ℝ} {S : Set ℝ}
    (hq : 1 < q) (hS : IsOpen S) (hconn : IsPreconnected S) (ha : a ∈ S)
    (hm : ContinuousOn m S) (hm0 : ∀ t ∈ S, 0 ≤ m t)
    (hr : ∀ t ∈ S, 0 < (f t).1)
    (hf : ∀ t ∈ S, HasDerivAt f (directionField q (m t) (f t)) t)
    (hg : ∀ t ∈ S, HasDerivAt g (directionField q (m t) (g t)) t)
    (hbase : f a = g a) : EqOn f g S := by
  let E : Set S := {t | f t = g t}
  have hfc : Continuous (fun t : S => f t) :=
    continuousOn_iff_continuous_restrict.mp
      (fun t ht => (hf t ht).continuousAt.continuousWithinAt)
  have hgc : Continuous (fun t : S => g t) :=
    continuousOn_iff_continuous_restrict.mp
      (fun t ht => (hg t ht).continuousAt.continuousWithinAt)
  have hclosed : IsClosed E := isClosed_eq hfc hgc
  have hopen : IsOpen E := by
    apply isOpen_iff_mem_nhds.mpr
    intro t ht
    have hf' : ∀ᶠ s in 𝓝 (t : ℝ),
        HasDerivAt f (directionField q (m s) (f s)) s := by
      filter_upwards [hS.mem_nhds t.property] with s hs
      exact hf s hs
    have hg' : ∀ᶠ s in 𝓝 (t : ℝ),
        HasDerivAt g (directionField q (m s) (g s)) s := by
      filter_upwards [hS.mem_nhds t.property] with s hs
      exact hg s hs
    have heq := direction_local_solution_unique hq (hr t t.property)
      (hm0 t t.property) ((hm t t.property).continuousAt (hS.mem_nhds t.property))
      hf' hg' ht
    exact continuous_subtype_val.continuousAt.eventually heq
  letI : PreconnectedSpace S := Subtype.preconnectedSpace hconn
  have hfull : E = univ := (show IsClopen E from ⟨hclosed, hopen⟩).eq_univ
    ⟨⟨a, ha⟩, hbase⟩
  intro t ht
  have : (⟨t, ht⟩ : S) ∈ E := by rw [hfull]; trivial
  exact this

/-- Two regular solutions with the same launch agree on their whole common
open interval, in both time directions. -/
theorem direction_solution_unique_on_interval {q a l b : ℝ} {m : ℝ → ℝ}
    {f g : ℝ → ℝ × ℝ}
    (hq : 1 < q) (ha : a ∈ Ioo l b)
    (hm : ContinuousOn m (Ioo l b)) (hm0 : ∀ t ∈ Ioo l b, 0 ≤ m t)
    (hr : ∀ t ∈ Ioo l b, 0 < (f t).1)
    (hf : ∀ t ∈ Ioo l b, HasDerivAt f (directionField q (m t) (f t)) t)
    (hg : ∀ t ∈ Ioo l b, HasDerivAt g (directionField q (m t) (g t)) t)
    (hbase : f a = g a) : EqOn f g (Ioo l b) :=
  direction_solution_unique_on_open_connected hq isOpen_Ioo
    (convex_Ioo l b).isPreconnected ha hm hm0 hr hf hg hbase

end ViaB
