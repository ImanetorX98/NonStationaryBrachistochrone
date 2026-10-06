import ViaB.Wronskian
import Mathlib.Analysis.Calculus.Deriv.Slope
import Mathlib.Tactic.Linarith

open Filter Topology

namespace ViaB

/-- A differentiable function positive immediately to the left of a zero has
a nonpositive derivative at that zero. No monotonicity assumption is needed. -/
theorem derivative_nonpos_at_first_zero {z : ℝ → ℝ} {a b dz : ℝ}
    (hab : a < b) (hz : HasDerivAt z dz b) (hzero : z b = 0)
    (hpos : ∀ t ∈ Set.Ioo a b, 0 < z t) : dz ≤ 0 := by
  have ht : Tendsto (slope z b) (𝓝[<] b) (𝓝 dz) :=
    (hasDerivWithinAt_iff_tendsto_slope' (show b ∉ Set.Iio b by simp)).mp
      (hz.hasDerivWithinAt)
  have hs : ∀ᶠ t in 𝓝[<] b, slope z b t ≤ 0 := by
    filter_upwards [self_mem_nhdsWithin,
      (eventually_gt_nhds hab).filter_mono nhdsWithin_le_nhds] with t ht hta
    simp only [slope_def_field, hzero, sub_zero]
    exact div_nonpos_of_nonneg_of_nonpos (le_of_lt (hpos t ⟨hta, ht⟩))
      (le_of_lt (sub_neg.mpr ht))
  exact le_of_tendsto ht hs

/-- Strict Sturm comparison excludes a candidate first zero before a reference
mode vanishes. Positivity before b is explicit; b is the hypothesized first zero. -/
theorem no_first_zero {u du z dz kperp kpar : ℝ → ℝ} {a b : ℝ}
    (hab : a < b)
    (hu : ∀ x ∈ Set.Icc a b, HasDerivAt u (du x) x)
    (hz : ∀ x ∈ Set.Icc a b, HasDerivAt z (dz x) x)
    (hdu : ∀ x ∈ Set.Icc a b, HasDerivAt du (-kperp x * u x) x)
    (hdz : ∀ x ∈ Set.Icc a b, HasDerivAt dz (-kpar x * z x) x)
    (hgap : ∀ x ∈ Set.Ioo a b, 0 < kperp x - kpar x)
    (hupos : ∀ x ∈ Set.Ioo a b, 0 < u x)
    (hzpos : ∀ x ∈ Set.Ioo a b, 0 < z x)
    (hua : u a = 0) (hza : z a = 0) (hub : 0 < u b) : z b ≠ 0 := by
  intro hzb
  have hmono := wronskian_strictMonoOn hu hz hdu hdz hgap hupos hzpos
  have hwa : wronskian u du z dz a = 0 := by simp [wronskian, hua, hza]
  have hwb : 0 < wronskian u du z dz b := by
    have hlt := hmono (show a ∈ Set.Icc a b from ⟨le_rfl, hab.le⟩)
      (show b ∈ Set.Icc a b from ⟨hab.le, le_rfl⟩) hab
    rwa [hwa] at hlt
  have hdnonpos := derivative_nonpos_at_first_zero hab (hz b ⟨hab.le, le_rfl⟩) hzb hzpos
  have hwbnonpos : wronskian u du z dz b ≤ 0 := by
    simp only [wronskian, hzb, zero_mul, sub_zero]
    exact mul_nonpos_of_nonpos_of_nonneg hdnonpos hub.le
  linarith

end ViaB
