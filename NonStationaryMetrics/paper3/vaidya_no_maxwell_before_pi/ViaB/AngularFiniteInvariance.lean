import ViaB.ExteriorFiniteTube
import ViaB.BoundedAngle

open Set
open scoped Topology NNReal

namespace ViaB

theorem direction_angular_derivative {q t : ℝ} {m : ℝ → ℝ} {f : ℝ → ℝ × ℝ}
    (hf : HasDerivAt f (directionField q (m t) (f t)) t) :
    HasDerivAt (fun v => (f v).2)
      (Real.sin (f t).2 / (f t).1 *
        directionB q (vaidyaW q (m t) (f t).1) (Real.cos (f t).2)) t := by
  simpa [directionField, vaidyaW] using
    (ContinuousLinearMap.snd ℝ ℝ ℝ).hasFDerivAt.comp_hasDerivAt t hf

/-- A scalar linear ODE with bounded coefficient cannot reach zero in finite
time from a nonzero initial value, by uniqueness backwards from the putative zero. -/
theorem linear_ode_initial_nonzero {a b : ℝ} {y k : ℝ → ℝ} {K : ℝ≥0}
    (hode : ∀ t ∈ Icc a b, HasDerivAt y (k t * y t) t)
    (hk : ∀ t ∈ Icc a b, |k t| ≤ (K : ℝ)) (hinit : y a ≠ 0) :
    ∀ t ∈ Icc a b, y t ≠ 0 := by
  intro t ht hzero
  have hs : Icc a t ⊆ Icc a b := fun u hu => ⟨hu.1, hu.2.trans ht.2⟩
  have hLip : ∀ u ∈ Ioc a t,
      LipschitzOnWith K (fun z : ℝ => k u * z) univ := by
    intro u hu
    apply LipschitzOnWith.of_dist_le_mul
    intro z _ w _
    rw [Real.dist_eq, show k u * z - k u * w = k u * (z - w) by ring,
      abs_mul, Real.dist_eq]
    exact mul_le_mul_of_nonneg_right (hk u (hs ⟨hu.1.le, hu.2⟩)) (abs_nonneg _)
  have heq := ODE_solution_unique_of_mem_Icc_left hLip
    (fun u hu => (hode u (hs hu)).continuousAt.continuousWithinAt)
    (fun u hu => (hode u (hs ⟨hu.1.le, hu.2⟩)).hasDerivWithinAt)
    (fun _ _ => mem_univ _)
    (continuousOn_const (c := (0 : ℝ)))
    (fun u _ => by simpa only [mul_zero] using (hasDerivAt_const u (0 : ℝ)).hasDerivWithinAt)
    (fun _ _ => mem_univ _) hzero
  exact hinit (heq ⟨le_rfl, ht.1⟩)

/-- The coefficient in (sin α)'=k sin α is bounded uniformly in α;
no angular-strip hypothesis is used. -/
theorem sine_direction_coefficient_bound {q ρ r w α : ℝ}
    (hq : 1 < q) (hρ : 0 < ρ) (hr : ρ ≤ r)
    (hw : q - 1 ≤ w ∧ w ≤ q) :
    |Real.cos α / r * directionB q w (Real.cos α)| ≤ directionBound q / ρ := by
  have hrpos : 0 < r := hρ.trans_le hr
  have hcos : |Real.cos α| ≤ 1 := abs_le.mpr ⟨Real.neg_one_le_cos α, Real.cos_le_one α⟩
  have hB := directionB_abs_bound hq hw.1 hw.2 ⟨Real.neg_one_le_cos α, Real.cos_le_one α⟩
  have hdiv : |Real.cos α| / r ≤ 1 / ρ :=
    (div_le_div_of_nonneg_right hcos hrpos.le).trans
      (div_le_div_of_nonneg_left (by norm_num) hρ hr)
  have hh := mul_le_mul hdiv hB (abs_nonneg _) (show 0 ≤ 1 / ρ by positivity)
  rw [abs_mul, abs_div, abs_of_pos hrpos]
  convert hh using 1; ring

/-- The nonradial angular strip is invariant on every finite existing
interval. It is a conclusion, not an assumption used to bound the coefficient. -/
theorem finite_direction_angle_stays_nonradial {q a b ρ : ℝ} {r w α : ℝ → ℝ}
    (hq : 1 < q) (hρ : 0 < ρ)
    (hr : ∀ t ∈ Icc a b, ρ ≤ r t)
    (hw : ∀ t ∈ Icc a b, q - 1 ≤ w t ∧ w t ≤ q)
    (hode : ∀ t ∈ Icc a b, HasDerivAt α
      (Real.sin (α t) / r t * directionB q (w t) (Real.cos (α t))) t)
    (hinit : 0 < α a ∧ α a < Real.pi) :
    ∀ t ∈ Icc a b, 0 < α t ∧ α t < Real.pi := by
  let k := fun t => Real.cos (α t) / r t * directionB q (w t) (Real.cos (α t))
  let K : ℝ≥0 := ⟨directionBound q / ρ, (div_pos (directionBound_pos hq) hρ).le⟩
  have hs : ∀ t ∈ Icc a b,
      HasDerivAt (fun u => Real.sin (α u)) (k t * Real.sin (α t)) t := by
    intro t ht
    convert (hode t ht).sin using 1
    dsimp [k]; ring
  have hnonzero := linear_ode_initial_nonzero (K := K) hs
    (fun t ht => sine_direction_coefficient_bound hq hρ (hr t ht) (hw t ht))
    (ne_of_gt (Real.sin_pos_of_pos_of_lt_pi hinit.1 hinit.2))
  intro t ht
  have hsub : Icc a t ⊆ Icc a b := fun u hu => ⟨hu.1, hu.2.trans ht.2⟩
  have hc : ContinuousOn α (Icc a t) := fun u hu => (hode u (hsub hu)).continuousAt.continuousWithinAt
  constructor
  · by_contra hbad
    obtain ⟨u, hu, heq⟩ := intermediate_value_Icc' ht.1 hc ⟨le_of_not_gt hbad, hinit.1.le⟩
    exact hnonzero u (hsub hu) (by simp only [heq, Real.sin_zero])
  · by_contra hbad
    obtain ⟨u, hu, heq⟩ := intermediate_value_Icc ht.1 hc ⟨hinit.2.le, le_of_not_gt hbad⟩
    exact hnonzero u (hsub hu) (by simp only [heq, Real.sin_pi])

/-- Exterior evolution and a nonradial launch derive the angular strip on
the whole finite existing interval, with no separate strip hypothesis. -/
theorem exterior_finite_angle_invariance {q a b : ℝ} {m : ℝ → ℝ} {f : ℝ → ℝ × ℝ}
    (hq : 1 < q) (hab : a < b) (hmpos : 0 < m a)
    (hmono : MonotoneOn m (Icc a b))
    (hode : ∀ t ∈ Ico a b, HasDerivAt f (directionField q (m t) (f t)) t)
    (hext : ∀ t ∈ Ico a b, 2 * m t < (f t).1)
    (hinit : 0 < (f a).2 ∧ (f a).2 < Real.pi) :
    ∀ t ∈ Ico a b, 0 < (f t).2 ∧ (f t).2 < Real.pi := by
  have hbnd := exterior_finite_state_bounds hq hab hmpos hmono hode hext
  intro t ht
  have hs : Icc a t ⊆ Ico a b := fun u hu => ⟨hu.1, hu.2.trans_lt ht.2⟩
  apply finite_direction_angle_stays_nonradial
    (r := fun u => (f u).1) (w := fun u => vaidyaW q (m u) (f u).1)
    hq (show 0 < 2 * m a by positivity)
    (fun u hu => (hbnd u (hs hu)).1)
    (fun u hu => exterior_w_bounds (q := q) (hbnd u (hs hu)).2.2.1
      (by linarith [(hbnd u (hs hu)).1]) (hext u (hs hu)).le)
    (fun u hu => direction_angular_derivative (hode u (hs hu))) hinit t ⟨ht.1, le_rfl⟩

/-- Coordinate continuation now follows from exterior data and the initial
nonradial angle alone: compact confinement and the strip are both derived. -/
theorem exterior_nonradial_finite_coordinate_continuation {q a b : ℝ} {m : ℝ → ℝ}
    {f : ℝ → ℝ × ℝ} (hq : 1 < q) (hab : a < b) (hmpos : 0 < m a)
    (hmono : MonotoneOn m (Icc a b)) (hm : ContDiffAt ℝ 1 m b)
    (hode : ∀ t ∈ Ico a b, HasDerivAt f (directionField q (m t) (f t)) t)
    (hext : ∀ t ∈ Ico a b, 2 * m t < (f t).1)
    (hinit : 0 < (f a).2 ∧ (f a).2 < Real.pi) :
    ∃ h : ℝ → ℝ × ℝ, ∃ ε > 0, EqOn h f (Ioo a b) ∧
      ∀ t ∈ Ioo a (b + ε), HasDerivAt h (directionField q (m t) (h t)) t := by
  have hα := exterior_finite_angle_invariance hq hab hmpos hmono hode hext hinit
  exact exterior_finite_coordinate_continuation hq hab hmpos hmono hm hode hext
    (fun t ht => ⟨(hα t ht).1.le, (hα t ht).2.le⟩)

/-- At a finite exterior endpoint the continued path is still nonradial,
and its radius is at or outside 2m(b). The limit cannot be a radial direction. -/
theorem exterior_nonradial_endpoint_control {q a b : ℝ} {m : ℝ → ℝ}
    {f : ℝ → ℝ × ℝ} (hq : 1 < q) (hab : a < b) (hmpos : 0 < m a)
    (hmono : MonotoneOn m (Icc a b)) (hm : ContDiffAt ℝ 1 m b)
    (hode : ∀ t ∈ Ico a b, HasDerivAt f (directionField q (m t) (f t)) t)
    (hext : ∀ t ∈ Ico a b, 2 * m t < (f t).1)
    (hinit : 0 < (f a).2 ∧ (f a).2 < Real.pi) :
    ∃ h : ℝ → ℝ × ℝ, ∃ ε > 0, EqOn h f (Ioo a b) ∧
      (∀ t ∈ Ioo a (b + ε), HasDerivAt h (directionField q (m t) (h t)) t) ∧
      2 * m b ≤ (h b).1 ∧ 0 < (h b).2 ∧ (h b).2 < Real.pi := by
  obtain ⟨h, ε, hε, heq, hodeh⟩ := exterior_nonradial_finite_coordinate_continuation
    hq hab hmpos hmono hm hode hext hinit
  have hbode := hodeh b ⟨hab, by linarith⟩
  have hbclosure : b ∈ closure (Ioo a b) := by
    rw [closure_Ioo hab.ne]; exact ⟨hab.le, le_rfl⟩
  haveI : Filter.NeBot (𝓝[Ioo a b] b) := mem_closure_iff_nhdsWithin_neBot.mp hbclosure
  have hmc : ContinuousAt (fun t => 2 * m t) b := hm.continuousAt.const_mul 2
  have hborder : 2 * m b ≤ (h b).1 := le_of_tendsto_of_tendsto
    hmc.continuousWithinAt
    hbode.continuousAt.fst.continuousWithinAt (by
      filter_upwards [self_mem_nhdsWithin] with t ht
      rw [heq ht]; exact (hext t ⟨ht.1.le, ht.2⟩).le)
  let c := (a + b) / 2
  have hac : a < c := by dsimp [c]; linarith
  have hcb : c < b := by dsimp [c]; linarith
  have hsub : Icc c b ⊆ Ioo a (b + ε) := fun u hu =>
    ⟨hac.trans_le hu.1, by linarith [hu.2]⟩
  have hma : ∀ u ∈ Icc c b, m a ≤ m u := fun u hu =>
    hmono ⟨le_rfl, hab.le⟩ ⟨hac.le.trans hu.1, hu.2⟩ (hac.le.trans hu.1)
  have hrad : ∀ u ∈ Icc c b, 2 * m u ≤ (h u).1 := by
    intro u hu
    by_cases hub : u = b
    · subst u; exact hborder
    · rw [heq ⟨hac.trans_le hu.1, lt_of_le_of_ne hu.2 hub⟩]
      exact (hext u ⟨hac.le.trans hu.1, lt_of_le_of_ne hu.2 hub⟩).le
  have hρ : 0 < 2 * m a := by positivity
  have hαold := exterior_finite_angle_invariance hq hab hmpos hmono hode hext hinit
  have hαc : 0 < (h c).2 ∧ (h c).2 < Real.pi := by
    rw [heq ⟨hac, hcb⟩]; exact hαold c ⟨hac.le, hcb⟩
  have hαnew := finite_direction_angle_stays_nonradial
    (q := q) (a := c) (b := b) (ρ := 2 * m a)
    (r := fun u => (h u).1) (w := fun u => vaidyaW q (m u) (h u).1)
    hq hρ (fun u hu => by linarith [hma u hu, hrad u hu])
    (fun u hu => exterior_w_bounds (q := q) (hmpos.trans_le (hma u hu)).le
      (by linarith [hma u hu, hrad u hu]) (hrad u hu))
    (fun u hu => direction_angular_derivative (hodeh u (hsub hu))) hαc
    b ⟨hcb.le, le_rfl⟩
  exact ⟨h, ε, hε, heq, hodeh, hborder, hαnew⟩

/-- At the capture surface, every nonradial direction crosses transversally
when the mass-rate is nonnegative. This is a derivative statement at the event. -/
theorem nonradial_capture_gap_derivative {q t μ : ℝ} {m : ℝ → ℝ}
    {f : ℝ → ℝ × ℝ} (hq : 1 < q) (hmpos : 0 < m t) (hμ : 0 ≤ μ)
    (hmder : HasDerivAt m μ t)
    (hf : HasDerivAt f (directionField q (m t) (f t)) t)
    (hborder : (f t).1 = 2 * m t)
    (hα : 0 < (f t).2 ∧ (f t).2 < Real.pi) :
    HasDerivAt (fun v => (f v).1 - 2 * m v)
      (q * (Real.cos (f t).2 - 1) - 2 * μ) t ∧
      q * (Real.cos (f t).2 - 1) - 2 * μ < 0 := by
  have hw : vaidyaW q (m t) (f t).1 = q := by
    unfold vaidyaW
    rw [hborder]
    field_simp
    ring
  have hroot : Real.sqrt (q * q) = q := by
    rw [← pow_two, Real.sqrt_sq (show 0 ≤ q by linarith)]
  have hc : Real.cos (f t).2 < 1 := by
    simpa only [Real.cos_zero] using
      Real.cos_lt_cos_of_nonneg_of_le_pi (by norm_num : (0 : ℝ) ≤ 0) hα.2.le hα.1
  constructor
  · convert (direction_radial_derivative hf).sub (hmder.const_mul 2) using 1
    rw [hw]
    unfold radialSpeed
    rw [hroot]
    ring
  · have hh := mul_neg_of_pos_of_neg (show 0 < q by linarith) (sub_neg.mpr hc)
    linarith

/-- A finite endpoint of an existing exterior nonradial trajectory is either
on the capture surface or admits a longer exterior nonradial interval. -/
theorem exterior_finite_endpoint_capture_or_continues {q a b : ℝ} {m : ℝ → ℝ}
    {f : ℝ → ℝ × ℝ} (hq : 1 < q) (hab : a < b) (hmpos : 0 < m a)
    (hmono : MonotoneOn m (Icc a b)) (hm : ContDiffAt ℝ 1 m b)
    (hode : ∀ t ∈ Ico a b, HasDerivAt f (directionField q (m t) (f t)) t)
    (hext : ∀ t ∈ Ico a b, 2 * m t < (f t).1)
    (hinit : 0 < (f a).2 ∧ (f a).2 < Real.pi) :
    ∃ h : ℝ → ℝ × ℝ, ∃ ε > 0, EqOn h f (Ioo a b) ∧
      (∀ t ∈ Ioo a (b + ε), HasDerivAt h (directionField q (m t) (h t)) t) ∧
      ((h b).1 = 2 * m b ∨ ∃ δ > 0, δ ≤ ε ∧
        ∀ t ∈ Ioo a (b + δ),
          2 * m t < (h t).1 ∧ 0 < (h t).2 ∧ (h t).2 < Real.pi) := by
  obtain ⟨h, ε, hε, heq, hodeh, hbext, hbα₀, hbαπ⟩ :=
    exterior_nonradial_endpoint_control hq hab hmpos hmono hm hode hext hinit
  refine ⟨h, ε, hε, heq, hodeh, ?_⟩
  by_cases hcapture : (h b).1 = 2 * m b
  · exact Or.inl hcapture
  · have hstrict : 2 * m b < (h b).1 := lt_of_le_of_ne hbext (Ne.symm hcapture)
    have hc := (hodeh b ⟨hab, by linarith⟩).continuousAt
    have hmc : ContinuousAt (fun t => 2 * m t) b := hm.continuousAt.const_mul 2
    have hev : {t | 2 * m t < (h t).1 ∧ 0 < (h t).2 ∧ (h t).2 < Real.pi} ∈ 𝓝 b := by
      filter_upwards [hmc.eventually_lt hc.fst hstrict,
        continuousAt_const.eventually_lt hc.snd hbα₀,
        hc.snd.eventually_lt continuousAt_const hbαπ] with t ht₁ ht₂ ht₃
      exact ⟨ht₁, ht₂, ht₃⟩
    obtain ⟨ν, hν, hball⟩ := Metric.mem_nhds_iff.mp hev
    refine Or.inr ⟨min ε ν, lt_min hε hν, min_le_left _ _, ?_⟩
    intro t ht
    by_cases htb : t < b
    · rw [heq ⟨ht.1, htb⟩]
      exact ⟨hext t ⟨ht.1.le, htb⟩,
        exterior_finite_angle_invariance hq hab hmpos hmono hode hext hinit t ⟨ht.1.le, htb⟩⟩
    · apply hball
      rw [Real.ball_eq_Ioo]
      constructor <;> linarith [le_of_not_gt htb, min_le_right ε ν, ht.2]

/-- Restore the launch point and its ODE when an open-interval continuation
has already been constructed. The new curve agrees with f also at a. -/
theorem direction_continuation_preserves_launch {q a b ε : ℝ} {m : ℝ → ℝ}
    {f h : ℝ → ℝ × ℝ} (hab : a < b)
    (haode : HasDerivAt f (directionField q (m a) (f a)) a)
    (heq : EqOn h f (Ioo a b))
    (hode : ∀ t ∈ Ioo a (b + ε), HasDerivAt h (directionField q (m t) (h t)) t) :
    ∃ H : ℝ → ℝ × ℝ, EqOn H f (Ico a b) ∧
      ∀ t ∈ Ico a (b + ε), HasDerivAt H (directionField q (m t) (H t)) t := by
  let H := fun t => if a < t then h t else f t
  have hloca : H =ᶠ[𝓝 a] f := by
    filter_upwards [Iio_mem_nhds hab] with t ht
    change (if a < t then h t else f t) = f t
    split_ifs with hat
    · exact heq ⟨hat, ht⟩
    · rfl
  refine ⟨H, ?_, ?_⟩
  · intro t ht
    change (if a < t then h t else f t) = f t
    split_ifs with hat
    · exact heq ⟨hat, ht.2⟩
    · rfl
  · intro t ht
    by_cases hta : t = a
    · subst t
      simpa only [H, lt_self_iff_false, if_false] using haode.congr_of_eventuallyEq hloca
    · have hat : a < t := lt_of_le_of_ne ht.1 (Ne.symm hta)
      have hloc : H =ᶠ[𝓝 t] h := by
        filter_upwards [Ioi_mem_nhds hat] with u hu
        exact if_pos hu
      simpa only [H, if_pos hat] using (hode t ⟨hat, ht.2⟩).congr_of_eventuallyEq hloc

/-- The finite exterior continuation preserves the full initial-value
solution, including the launch point and its derivative. -/
theorem exterior_nonradial_initial_value_continuation {q a b : ℝ} {m : ℝ → ℝ}
    {f : ℝ → ℝ × ℝ} (hq : 1 < q) (hab : a < b) (hmpos : 0 < m a)
    (hmono : MonotoneOn m (Icc a b)) (hm : ContDiffAt ℝ 1 m b)
    (hode : ∀ t ∈ Ico a b, HasDerivAt f (directionField q (m t) (f t)) t)
    (hext : ∀ t ∈ Ico a b, 2 * m t < (f t).1)
    (hinit : 0 < (f a).2 ∧ (f a).2 < Real.pi) :
    ∃ H : ℝ → ℝ × ℝ, ∃ ε > 0, EqOn H f (Ico a b) ∧
      ∀ t ∈ Ico a (b + ε), HasDerivAt H (directionField q (m t) (H t)) t := by
  obtain ⟨h, ε, hε, heq, hodeh⟩ := exterior_nonradial_finite_coordinate_continuation
    hq hab hmpos hmono hm hode hext hinit
  obtain ⟨H, heqH, hodeH⟩ := direction_continuation_preserves_launch hab
    (hode a ⟨le_rfl, hab⟩) heq hodeh
  exact ⟨H, ε, hε, heqH, hodeH⟩

end ViaB
