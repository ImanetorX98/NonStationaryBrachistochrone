import ViaB.LocalDirectionExistence
import Mathlib.Analysis.Normed.Group.Bounded

open Filter Set
open scoped Topology NNReal

namespace ViaB

/-- A Lipschitz planar path admits a continuous extension and therefore a
finite endpoint limit. The extension is not asserted to solve an ODE. -/
theorem planar_lipschitz_extension {s : Set ℝ} {f : ℝ → ℝ × ℝ} {K : ℝ≥0}
    (hf : LipschitzOnWith K f s) :
    ∃ g : ℝ → ℝ × ℝ, LipschitzWith K g ∧ EqOn f g s := by
  have h₁ : LipschitzOnWith K (fun t => (f t).1) s := by
    simpa only [one_mul] using LipschitzWith.prod_fst.comp_lipschitzOnWith hf
  have h₂ : LipschitzOnWith K (fun t => (f t).2) s := by
    simpa only [one_mul] using LipschitzWith.prod_snd.comp_lipschitzOnWith hf
  obtain ⟨g₁, hg₁, heq₁⟩ := h₁.extend_real
  obtain ⟨g₂, hg₂, heq₂⟩ := h₂.extend_real
  refine ⟨fun t => (g₁ t, g₂ t), ?_, ?_⟩
  · simpa only [max_self] using hg₁.prodMk hg₂
  · intro t ht
    exact Prod.ext (heq₁ ht) (heq₂ ht)

/-- Bounded speed supplies the endpoint limit without assuming one. -/
theorem bounded_speed_endpoint_limit {a b : ℝ} {f df : ℝ → ℝ × ℝ} {K : ℝ≥0}
    (hder : ∀ t ∈ Ioo a b, HasDerivAt f (df t) t)
    (hbound : ∀ t ∈ Ioo a b, ‖df t‖₊ ≤ K) :
    ∃ z : ℝ × ℝ, Tendsto f (𝓝[Ioo a b] b) (𝓝 z) := by
  have hLip := (convex_Ioo a b).lipschitzOnWith_of_nnnorm_hasDerivWithin_le
    (fun t ht => (hder t ht).hasDerivWithinAt) hbound
  obtain ⟨g, hg, heq⟩ := planar_lipschitz_extension hLip
  refine ⟨g b, ?_⟩
  exact tendsto_nhdsWithin_congr (fun t ht => (heq ht).symm)
    hg.continuous.continuousAt.continuousWithinAt

/-- Compactness of a regular state-mass tube gives a uniform speed bound. -/
theorem compact_direction_speed_bound {q : ℝ} {C : Set ((ℝ × ℝ) × ℝ)}
    (hq : 1 < q) (hC : IsCompact C)
    (hr : ∀ z ∈ C, 0 < z.1.1) (hm : ∀ z ∈ C, 0 ≤ z.2) :
    ∃ K : ℝ≥0, ∀ z ∈ C, ‖directionField q z.2 z.1‖₊ ≤ K := by
  have hc : ContinuousOn (fun z : (ℝ × ℝ) × ℝ => directionField q z.2 z.1) C :=
    fun z hz => (directionField_contDiffAt hq (hr z hz) (hm z hz)).continuousAt.continuousWithinAt
  obtain ⟨R, hR⟩ := hC.exists_bound_of_continuousOn hc
  refine ⟨⟨max R 0, le_max_right _ _⟩, ?_⟩
  intro z hz
  change ‖directionField q z.2 z.1‖ ≤ max R 0
  exact (hR z hz).trans (le_max_left _ _)

/-- A trajectory trapped in a regular compact state-mass tube has a regular
finite endpoint and an independently constructed local solution through it.
Agreement with the old path near the endpoint is a separate gluing step. -/
theorem compact_direction_endpoint_restart {q a b : ℝ} {m : ℝ → ℝ}
    {f : ℝ → ℝ × ℝ} {C : Set ((ℝ × ℝ) × ℝ)}
    (hq : 1 < q) (hab : a < b) (hC : IsCompact C)
    (hr : ∀ z ∈ C, 0 < z.1.1) (hm0 : ∀ z ∈ C, 0 ≤ z.2)
    (hm : ContDiffAt ℝ 1 m b)
    (hode : ∀ t ∈ Ioo a b, HasDerivAt f (directionField q (m t) (f t)) t)
    (htube : ∀ t ∈ Ioo a b, (f t, m t) ∈ C) :
    ∃ z : ℝ × ℝ, Tendsto f (𝓝[Ioo a b] b) (𝓝 z) ∧ (z, m b) ∈ C ∧
      ∃ g : ℝ → ℝ × ℝ, g b = z ∧ ∃ ε > 0,
        ∀ t ∈ Ioo (b - ε) (b + ε),
          HasDerivAt g (directionField q (m t) (g t)) t := by
  obtain ⟨K, hK⟩ := compact_direction_speed_bound hq hC hr hm0
  obtain ⟨z, hz⟩ := bounded_speed_endpoint_limit hode (fun t ht => hK (f t, m t) (htube t ht))
  have hb : b ∈ closure (Ioo a b) := by
    rw [closure_Ioo hab.ne]
    exact ⟨hab.le, le_rfl⟩
  haveI : NeBot (𝓝[Ioo a b] b) := mem_closure_iff_nhdsWithin_neBot.mp hb
  have hzm : Tendsto (fun t => (f t, m t)) (𝓝[Ioo a b] b) (𝓝 (z, m b)) :=
    hz.prodMk_nhds hm.continuousAt.continuousWithinAt
  have hzc : (z, m b) ∈ C := hC.isClosed.mem_of_tendsto hzm
    (eventually_nhdsWithin_of_forall htube)
  exact ⟨z, hz, hzc, direction_local_solution_exists hq (hr _ hzc) (hm0 _ hzc) hm⟩

end ViaB
