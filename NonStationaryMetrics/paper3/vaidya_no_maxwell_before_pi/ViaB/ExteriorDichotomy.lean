import ViaB.BoundedLingering

open Filter Topology

namespace ViaB

theorem finite_angle_exterior_escapes {q M a U : ℝ} (flow : DirectionFlow q M a)
    (φ : ℝ → ℝ) (hq : 1 < q) (hM : 0 < M) (hmpositive : 0 < flow.mass a)
    (hexterior : ∀ t ≥ a, 2 * flow.mass t ≤ flow.radius t)
    (hmmono : MonotoneOn flow.mass (Set.Ici a))
    (hφode : ∀ t ≥ a, HasDerivAt φ
      (Real.sqrt (q - 1 + 2 * flow.mass t / flow.radius t) / flow.radius t * Real.sin (flow.angle t)) t)
    (hφupper : ∀ t ≥ a, φ t ≤ U) : Tendsto flow.radius atTop atTop := by
  apply unbounded_radius_escapes flow hq hM
  intro K
  by_contra hn
  have hrupper : ∀ t ≥ a, flow.radius t ≤ K := by
    intro t ht
    by_contra hnot
    exact hn ⟨t, ht, lt_of_not_ge hnot⟩
  have hρ : 0 < 2 * flow.mass a := by positivity
  have hr : ∀ t ≥ a, 2 * flow.mass a ≤ flow.radius t ∧ flow.radius t ≤ K := by
    intro t ht
    have hm := hmmono (by simp : a ∈ Set.Ici a) (Set.mem_Ici.mpr ht) ht
    exact ⟨by linarith [hexterior t ht], hrupper t ht⟩
  exact no_bounded_exterior_finite_angle flow φ hq hρ hr hexterior hmmono hφode hφupper

/-- A global exterior nonradial solution either escapes or accumulates
unbounded swept angle. The alternatives are not asserted to be exclusive. -/
theorem exterior_escape_or_angle_diverges {q M a : ℝ} (flow : DirectionFlow q M a)
    (φ : ℝ → ℝ) (hq : 1 < q) (hM : 0 < M) (hmpositive : 0 < flow.mass a)
    (hexterior : ∀ t ≥ a, 2 * flow.mass t ≤ flow.radius t)
    (hmmono : MonotoneOn flow.mass (Set.Ici a))
    (hφode : ∀ t ≥ a, HasDerivAt φ
      (Real.sqrt (q - 1 + 2 * flow.mass t / flow.radius t) / flow.radius t * Real.sin (flow.angle t)) t) :
    Tendsto flow.radius atTop atTop ∨ Tendsto φ atTop atTop := by
  by_cases hu : ∀ K : ℝ, ∃ t ≥ a, K < flow.radius t
  · exact Or.inl (unbounded_radius_escapes flow hq hM hu)
  · push_neg at hu
    obtain ⟨R, hR⟩ := hu
    exact Or.inr (bounded_exterior_flow_angle_diverges flow φ hq hmpositive hR hexterior hmmono hφode)

end ViaB
