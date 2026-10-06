import ViaB.RegularEndpointContinuation
import ViaB.NullHamiltonianDynamics

open Set
open scoped Topology

namespace ViaB

/-- The direction system supplies the scalar radial derivative. -/
theorem direction_radial_derivative {q t : ℝ} {m : ℝ → ℝ} {f : ℝ → ℝ × ℝ}
    (hf : HasDerivAt f (directionField q (m t) (f t)) t) :
    HasDerivAt (fun v => (f v).1)
      (radialSpeed q (vaidyaW q (m t) (f t).1) (Real.cos (f t).2)) t := by
  simpa [directionField, vaidyaW] using
    (ContinuousLinearMap.fst ℝ ℝ ℝ).hasFDerivAt.comp_hasDerivAt t hf

/-- The exterior domain itself supplies w∈[q−1,q]. -/
theorem exterior_w_bounds {q m r : ℝ} (hm : 0 ≤ m) (hr : 0 < r)
    (hext : 2 * m ≤ r) : q - 1 ≤ vaidyaW q m r ∧ vaidyaW q m r ≤ q := by
  have hdiv : 2 * m / r ≤ 1 := (div_le_iff₀ hr).mpr (by simpa using hext)
  have hnonneg : 0 ≤ 2 * m / r := by positivity
  unfold vaidyaW
  constructor <;> linarith

/-- Uniform upper speed in the exterior; no radius bound is assumed. -/
theorem exterior_radial_speed_le_charge {q m r α : ℝ}
    (hq : 1 < q) (hm : 0 ≤ m) (hr : 0 < r) (hext : 2 * m ≤ r) :
    radialSpeed q (vaidyaW q m r) (Real.cos α) ≤ q := by
  have hw := exterior_w_bounds (q := q) hm hr hext
  exact radialSpeed_le_charge (by linarith) (by linarith [hw.1]) hw.2 (Real.cos_le_one α)

/-- Integrating the exterior speed bound gives a finite-time radius ceiling. -/
theorem exterior_finite_radius_upper {q a b : ℝ} {m : ℝ → ℝ} {f : ℝ → ℝ × ℝ}
    (hq : 1 < q)
    (hode : ∀ t ∈ Ico a b, HasDerivAt f (directionField q (m t) (f t)) t)
    (hm : ∀ t ∈ Ico a b, 0 ≤ m t)
    (hr : ∀ t ∈ Ico a b, 0 < (f t).1)
    (hext : ∀ t ∈ Ico a b, 2 * m t ≤ (f t).1) :
    ∀ t ∈ Ico a b, (f t).1 ≤ (f a).1 + q * (t - a) := by
  intro t ht
  have hs : Icc a t ⊆ Ico a b := fun u hu => ⟨hu.1, hu.2.trans_lt ht.2⟩
  have hd := fun u hu => direction_radial_derivative (hode u (hs hu))
  have hh := (convex_Icc a t).image_sub_le_mul_sub_of_deriv_le
    (fun u hu => (hd u hu).continuousAt.continuousWithinAt)
    (fun u hu => (hd u (interior_subset hu)).differentiableAt.differentiableWithinAt)
    (fun u hu => by
      rw [(hd u (interior_subset hu)).deriv]
      exact exterior_radial_speed_le_charge hq (hm u (hs (interior_subset hu)))
        (hr u (hs (interior_subset hu))) (hext u (hs (interior_subset hu))))
    a ⟨le_rfl, ht.1⟩ t ⟨ht.1, le_rfl⟩ ht.1
  linarith

/-- Monotone positive mass supplies lower and upper bounds in finite time. -/
theorem exterior_finite_state_bounds {q a b : ℝ} {m : ℝ → ℝ} {f : ℝ → ℝ × ℝ}
    (hq : 1 < q) (hab : a < b) (hmpos : 0 < m a)
    (hmono : MonotoneOn m (Icc a b))
    (hode : ∀ t ∈ Ico a b, HasDerivAt f (directionField q (m t) (f t)) t)
    (hext : ∀ t ∈ Ico a b, 2 * m t < (f t).1) :
    ∀ t ∈ Ico a b,
      2 * m a ≤ (f t).1 ∧ (f t).1 ≤ (f a).1 + q * (b - a) ∧
      0 ≤ m t ∧ m t ≤ m b := by
  have hma : ∀ t ∈ Ico a b, m a ≤ m t :=
    fun t ht => hmono ⟨le_rfl, hab.le⟩ ⟨ht.1, ht.2.le⟩ ht.1
  have hmr : ∀ t ∈ Ico a b, 0 < m t := fun t ht => hmpos.trans_le (hma t ht)
  have hr : ∀ t ∈ Ico a b, 0 < (f t).1 := by
    intro t ht; linarith [hmr t ht, hext t ht]
  have hup := exterior_finite_radius_upper hq hode (fun t ht => (hmr t ht).le)
    hr (fun t ht => (hext t ht).le)
  intro t ht
  have hmb := hmono ⟨ht.1, ht.2.le⟩ ⟨hab.le, le_rfl⟩ ht.2.le
  refine ⟨by linarith [hma t ht, hext t ht], ?_, (hmr t ht).le, hmb⟩
  have htime := mul_le_mul_of_nonneg_left (sub_le_sub_right ht.2.le a)
    (show 0 ≤ q by linarith)
  linarith [hup t ht]

/-- A regular compact tube is derived from exterior evolution and the closed
angular strip. It is not supplied as a separate confinement hypothesis. -/
theorem exterior_finite_regular_compact_tube {q a b : ℝ} {m : ℝ → ℝ}
    {f : ℝ → ℝ × ℝ} (hq : 1 < q) (hab : a < b) (hmpos : 0 < m a)
    (hmono : MonotoneOn m (Icc a b))
    (hode : ∀ t ∈ Ico a b, HasDerivAt f (directionField q (m t) (f t)) t)
    (hext : ∀ t ∈ Ico a b, 2 * m t < (f t).1)
    (hα : ∀ t ∈ Ico a b, 0 ≤ (f t).2 ∧ (f t).2 ≤ Real.pi) :
    ∃ C : Set ((ℝ × ℝ) × ℝ), IsCompact C ∧
      (∀ z ∈ C, 0 < z.1.1) ∧ (∀ z ∈ C, 0 ≤ z.2) ∧
      ∀ t ∈ Ico a b, (f t, m t) ∈ C := by
  let C := (Icc (2 * m a) ((f a).1 + q * (b - a)) ×ˢ Icc 0 Real.pi) ×ˢ Icc 0 (m b)
  have hbnd := exterior_finite_state_bounds hq hab hmpos hmono hode hext
  refine ⟨C, (isCompact_Icc.prod isCompact_Icc).prod isCompact_Icc, ?_, ?_, ?_⟩
  · intro z hz
    have hh : 2 * m a ≤ z.1.1 := hz.1.1.1
    linarith
  · intro z hz; exact hz.2.1
  · intro t ht
    exact ⟨⟨⟨(hbnd t ht).1, (hbnd t ht).2.1⟩, hα t ht⟩,
      ⟨(hbnd t ht).2.2.1, (hbnd t ht).2.2.2⟩⟩

/-- Finite-time exterior evolution cannot lose existence of the coordinate
ODE while remaining in the closed angular strip. Capture is not excluded. -/
theorem exterior_finite_coordinate_continuation {q a b : ℝ} {m : ℝ → ℝ}
    {f : ℝ → ℝ × ℝ} (hq : 1 < q) (hab : a < b) (hmpos : 0 < m a)
    (hmono : MonotoneOn m (Icc a b)) (hm : ContDiffAt ℝ 1 m b)
    (hode : ∀ t ∈ Ico a b, HasDerivAt f (directionField q (m t) (f t)) t)
    (hext : ∀ t ∈ Ico a b, 2 * m t < (f t).1)
    (hα : ∀ t ∈ Ico a b, 0 ≤ (f t).2 ∧ (f t).2 ≤ Real.pi) :
    ∃ h : ℝ → ℝ × ℝ, ∃ ε > 0, EqOn h f (Ioo a b) ∧
      ∀ t ∈ Ioo a (b + ε), HasDerivAt h (directionField q (m t) (h t)) t := by
  obtain ⟨C, hC, hrC, hmC, htube⟩ := exterior_finite_regular_compact_tube
    hq hab hmpos hmono hode hext hα
  exact compact_direction_continues hq hab hC hrC hmC hm
    (fun t ht => hode t ⟨ht.1.le, ht.2⟩)
    (fun t ht => htube t ⟨ht.1.le, ht.2⟩)

end ViaB
