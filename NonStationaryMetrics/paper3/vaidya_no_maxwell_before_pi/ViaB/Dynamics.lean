import ViaB.Sturm
import ViaB.FarField
import Mathlib.Analysis.SpecialFunctions.Log.Deriv
import Mathlib.Tactic.FieldSimp

open Filter Topology

namespace ViaB

theorem increment_ge_of_derivative_lower {f df : ℝ → ℝ} {a b C : ℝ}
    (hab : a ≤ b) (hder : ∀ t ∈ Set.Icc a b, HasDerivAt f (df t) t)
    (hbound : ∀ t ∈ Set.Ioo a b, C ≤ df t) : C * (b - a) ≤ f b - f a := by
  apply (convex_Icc a b).mul_sub_le_image_sub_of_le_deriv
    (fun t ht => (hder t ht).continuousAt.continuousWithinAt)
    (fun t ht => (hder t (interior_subset ht)).differentiableAt.differentiableWithinAt)
    (fun t ht => ?_) a ⟨le_rfl, hab⟩ b ⟨hab, le_rfl⟩ hab
  rw [(hder t (interior_subset ht)).deriv]
  exact hbound t (by simpa only [interior_Icc] using ht)

theorem increment_le_of_derivative_upper {f df : ℝ → ℝ} {a b C : ℝ}
    (hab : a ≤ b) (hder : ∀ t ∈ Set.Icc a b, HasDerivAt f (df t) t)
    (hbound : ∀ t ∈ Set.Ioo a b, df t ≤ C) : f b - f a ≤ C * (b - a) := by
  have hh := increment_ge_of_derivative_lower (f := fun t => -f t) (df := fun t => -df t)
    (C := -C) hab (fun t ht => (hder t ht).neg) (fun t ht => neg_le_neg (hbound t ht))
  linarith

theorem positive_right_of_positive_value {f : ℝ → ℝ} {a : ℝ}
    (hcont : ContinuousAt f a) (hpos : 0 < f a) :
    ∃ c > a, ∀ t ∈ Set.Ioo a c, 0 < f t := by
  have hh : ∀ᶠ t in 𝓝[>] a, 0 < f t :=
    (hcont.eventually (eventually_gt_nhds hpos)).filter_mono nhdsWithin_le_nhds
  exact mem_nhdsGT_iff_exists_Ioo_subset.mp hh

/-- Coupled strict barriers: neither coordinate can have a first zero if its
derivative is strictly positive on its boundary while the other is nonnegative. -/
theorem positive_pair_invariant {f df g dg : ℝ → ℝ} {a b : ℝ}
    (hf : ∀ t ∈ Set.Icc a b, HasDerivAt f (df t) t)
    (hg : ∀ t ∈ Set.Icc a b, HasDerivAt g (dg t) t)
    (hfa : 0 < f a) (hga : 0 < g a)
    (hfb : ∀ t ∈ Set.Icc a b, f t = 0 → 0 ≤ g t → 0 < df t)
    (hgb : ∀ t ∈ Set.Icc a b, g t = 0 → 0 ≤ f t → 0 < dg t) :
    ∀ t ∈ Set.Icc a b, 0 < f t ∧ 0 < g t := by
  intro t ht
  rcases ht.1.eq_or_lt with heq | hat
  · simpa [← heq] using And.intro hfa hga
  have ha : a ∈ Set.Icc a b := ⟨le_rfl, ht.1.trans ht.2⟩
  have hmcont : ContinuousAt (fun x => min (f x) (g x)) a :=
    (hf a ha).continuousAt.min (hg a ha).continuousAt
  obtain ⟨c, hac, hstart⟩ := positive_right_of_positive_value hmcont (lt_min hfa hga)
  by_contra hnpos
  have hmin : min (f t) (g t) ≤ 0 := by
    by_contra hnot
    exact hnpos (lt_min_iff.mp (lt_of_not_ge hnot))
  have hsub : Set.Icc a t ⊆ Set.Icc a b := Set.Icc_subset_Icc le_rfl ht.2
  have hcont : ContinuousOn (fun x => min (f x) (g x)) (Set.Icc a t) := by
    intro x hx
    have hh : ContinuousAt (fun y => min (f y) (g y)) x :=
      (hf x (hsub hx)).continuousAt.min (hg x (hsub hx)).continuousAt
    exact hh.continuousWithinAt
  obtain ⟨d, hd, hdmin, hbefore⟩ := exists_first_zero_of_nonpos hat hac hcont hstart hmin
  have hdb : d ∈ Set.Icc a b := ⟨hd.1.le, hd.2.trans ht.2⟩
  by_cases hfg : f d ≤ g d
  · have hfd : f d = 0 := by rwa [min_eq_left hfg] at hdmin
    have hgd : 0 ≤ g d := by linarith
    have hnonpos := derivative_nonpos_at_first_zero hd.1 (hf d hdb) hfd
      (fun x hx => (lt_min_iff.mp (hbefore x hx)).1)
    exact (not_lt_of_ge hnonpos) (hfb d hdb hfd hgd)
  · have hgf : g d ≤ f d := le_of_not_ge hfg
    have hgd : g d = 0 := by rwa [min_eq_right hgf] at hdmin
    have hfd : 0 ≤ f d := by linarith
    have hnonpos := derivative_nonpos_at_first_zero hd.1 (hg d hdb) hgd
      (fun x hx => (lt_min_iff.mp (hbefore x hx)).2)
    exact (not_lt_of_ge hnonpos) (hgb d hdb hgd hfd)

/-- Integrated excursion estimate proved without an integral: the auxiliary
function alpha + (k/C) log r has nonpositive derivative. -/
theorem excursion_log_bound {r dr α dα : ℝ → ℝ} {a b k C : ℝ}
    (hab : a ≤ b) (hk : 0 < k) (hC : 0 < C)
    (hr : ∀ t ∈ Set.Icc a b, 0 < r t)
    (hdr : ∀ t ∈ Set.Icc a b, HasDerivAt r (dr t) t)
    (hdα : ∀ t ∈ Set.Icc a b, HasDerivAt α (dα t) t)
    (hrbound : ∀ t ∈ Set.Ioo a b, dr t ≤ C)
    (hαbound : ∀ t ∈ Set.Ioo a b, dα t ≤ -k / r t)
    (hαa : α a ≤ Real.pi) (hαb : 0 ≤ α b) :
    k * (Real.log (r b) - Real.log (r a)) ≤ C * Real.pi := by
  let F := fun t => α t + (k / C) * Real.log (r t)
  let dF := fun t => dα t + (k / C) * (dr t / r t)
  have hderF : ∀ t ∈ Set.Icc a b, HasDerivAt F (dF t) t := by
    intro t ht
    exact (hdα t ht).add (((hdr t ht).log (ne_of_gt (hr t ht))).const_mul (k / C))
  have hboundF : ∀ t ∈ Set.Ioo a b, dF t ≤ 0 := by
    intro t ht
    have hrpos := hr t (Set.Ioo_subset_Icc_self ht)
    have hdiv := div_le_div_of_nonneg_right (hrbound t ht) hrpos.le
    have hmul := mul_le_mul_of_nonneg_left hdiv (div_pos hk hC).le
    have heq : -k / r t + (k / C) * (C / r t) = 0 := by
      field_simp
      ring
    dsimp [dF]
    calc
      _ ≤ -k / r t + (k / C) * (C / r t) := add_le_add (hαbound t ht) hmul
      _ = 0 := heq
  have hF := increment_le_of_derivative_upper hab hderF hboundF
  dsimp [F] at hF
  have hmul := mul_le_mul_of_nonneg_left hF hC.le
  have heq : C * (α b + (k / C) * Real.log (r b) - (α a + (k / C) * Real.log (r a))) =
      C * (α b - α a) + k * (Real.log (r b) - Real.log (r a)) := by
    field_simp
    ring
  rw [heq] at hmul
  nlinarith [mul_nonneg hC.le hαb, mul_nonneg hC.le (sub_nonneg.mpr hαa)]

end ViaB
