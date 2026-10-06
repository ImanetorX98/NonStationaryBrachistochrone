import ViaB.Algebra
import Mathlib.Topology.Covering
import Mathlib.Topology.LocallyConstant.Basic
import Mathlib.SetTheory.Cardinal.Finite
import Mathlib.Tactic.Ring

namespace ViaB

/-- Linear algebra at a receiver crossing, including radial tangencies.
A,B are the derivatives of r; C,D those of the swept angle. -/
theorem trace_tangent_nondegenerate {A B C D : ℝ} (hdet : A * D - B * C ≠ 0) :
    A * B + B * (-A) = 0 ∧ C * B + D * (-A) ≠ 0 := by
  constructor
  · ring
  · have heq : C * B + D * (-A) = -(A * D - B * C) := by ring
    rw [heq]
    exact neg_ne_zero.mpr hdet

noncomputable def fiberCount {E X : Type*} (f : E → X) (x : X) : ℕ :=
  Nat.card (f ⁻¹' {x})

/-- Covering maps have locally constant fiber cardinality (Nat.card assigns
zero to infinite fibers; the subsequent count-one theorem has no such ambiguity). -/
theorem covering_fiberCount_locallyConstant {E X : Type*}
    [TopologicalSpace E] [TopologicalSpace X] {f : E → X}
    (hf : IsCoveringMap f) : IsLocallyConstant (fiberCount f) := by
  apply (IsLocallyConstant.iff_exists_open _).2
  intro x
  obtain ⟨hdisc, U, hxU, hU, hfU, H, hH⟩ := hf x
  refine ⟨U, hU, hxU, ?_⟩
  intro y hyU
  have hy : IsEvenlyCovered f y (f ⁻¹' {x}) := ⟨hdisc, U, hyU, hU, hfU, H, hH⟩
  exact (Nat.card_congr hy.fiberHomeomorph.toEquiv).symm

/-- The one-sheet step, independently of geometric properness arguments.
It is conditional on an actual covering map and an actual singleton fiber. -/
theorem covering_injective_of_one_fiber {E X : Type*}
    [TopologicalSpace E] [TopologicalSpace X] [PreconnectedSpace X]
    {f : E → X} (hf : IsCoveringMap f) {x₀ : X} (hone : fiberCount f x₀ = 1) :
    Function.Injective f := by
  have hcount : ∀ x, fiberCount f x = 1 := by
    intro x
    rw [(covering_fiberCount_locallyConstant hf).apply_eq_of_isPreconnected
      isPreconnected_univ (Set.mem_univ x) (Set.mem_univ x₀)]
    exact hone
  intro a b hab
  have hsingle : Subsingleton (f ⁻¹' {f a}) :=
    (Nat.card_eq_one_iff_unique.mp (hcount (f a))).1
  have heq : (⟨a, rfl⟩ : f ⁻¹' {f a}) = ⟨b, hab.symm⟩ := hsingle.elim _ _
  exact congrArg Subtype.val heq

end ViaB
