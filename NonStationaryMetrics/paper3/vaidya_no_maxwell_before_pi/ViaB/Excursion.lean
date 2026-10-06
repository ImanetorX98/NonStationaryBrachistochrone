import ViaB.ConeEscape

namespace ViaB

theorem outward_direction_lt_pi_div_two {q w α : ℝ} (hw : 0 < w)
    (hα : α ≤ Real.pi) (hout : 0 ≤ radialSpeed q w (Real.cos α)) : α < Real.pi / 2 := by
  by_contra hn
  have hcos := Real.cos_nonpos_of_pi_div_two_le_of_le (le_of_not_gt hn)
    (show α ≤ Real.pi + Real.pi / 2 by linarith [Real.pi_pos])
  have hmul := mul_nonpos_of_nonneg_of_nonpos (Real.sqrt_nonneg (q * w)) hcos
  unfold radialSpeed at hout
  linarith

/-- Every sufficiently large outward excursion staying above R must meet the
angular cone. Choosing such an excursion for an arbitrary unbounded trajectory
is a separate compactness/last-crossing step, not a hypothesis hidden here. -/
theorem large_excursion_enters_angular_cone {q M a b β R κ : ℝ}
    (flow : DirectionFlow q M a) (hq : 1 < q)
    (hab : a ≤ b) (hR : 2 * M < R) (hκ : 0 < κ)
    (hβ : 0 < β ∧ β < Real.pi / 2)
    (hB : ∀ r m c : ℝ, R ≤ r → 0 ≤ m → m ≤ M → c ≤ 1 →
      directionB q (q - 1 + 2 * m / r) c ≤ -κ)
    (hrR : ∀ t ∈ Set.Icc a b, R ≤ flow.radius t)
    (hout : 0 ≤ radialSpeed q (q - 1 + 2 * flow.mass a / flow.radius a)
      (Real.cos (flow.angle a)))
    (hlarge : q * Real.pi < κ * Real.sin β *
      (Real.log (flow.radius b) - Real.log (flow.radius a))) :
    ∃ t ∈ Set.Icc a b, flow.angle t < β := by
  have hqpos : 0 < q := by linarith
  have hsinβ : 0 < Real.sin β := Real.sin_pos_of_pos_of_lt_pi hβ.1 (by linarith [Real.pi_pos])
  have ha : a ∈ Set.Icc a b := ⟨le_rfl, hab⟩
  have hwa : 0 < q - 1 + 2 * flow.mass a / flow.radius a := by
    have hm := (flow.mass_range a le_rfl).1
    have hr := flow.radius_pos a le_rfl
    have hdiv : 0 ≤ 2 * flow.mass a / flow.radius a := by positivity
    linarith
  have hanglea := outward_direction_lt_pi_div_two hwa (flow.angle_range a le_rfl).2.le hout
  have hBpoint : ∀ t ∈ Set.Icc a b,
      directionB q (q - 1 + 2 * flow.mass t / flow.radius t) (Real.cos (flow.angle t)) ≤ -κ := by
    intro t ht
    exact hB _ _ _ (hrR t ht) (flow.mass_range t ht.1).1 (flow.mass_range t ht.1).2
      (Real.cos_le_one _)
  have hnonpos : ∀ t ∈ Set.Icc a b,
      Real.sin (flow.angle t) / flow.radius t *
        directionB q (q - 1 + 2 * flow.mass t / flow.radius t) (Real.cos (flow.angle t)) ≤ 0 := by
    intro t ht
    exact mul_nonpos_of_nonneg_of_nonpos
      (div_pos (Real.sin_pos_of_pos_of_lt_pi (flow.angle_range t ht.1).1
        (flow.angle_range t ht.1).2) (flow.radius_pos t ht.1)).le (by linarith [hBpoint t ht])
  have hangle : ∀ t ∈ Set.Icc a b, flow.angle t ≤ flow.angle a := by
    intro t ht
    have hh := increment_le_of_derivative_upper (C := 0) ht.1
      (fun x hx => flow.angular_ode x hx.1)
      (fun x hx => hnonpos x ⟨hx.1.le, hx.2.le.trans ht.2⟩)
    linarith
  by_contra hnot
  have hangleβ : ∀ t ∈ Set.Icc a b, β ≤ flow.angle t := by
    intro t ht
    by_contra hn
    exact hnot ⟨t, ht, lt_of_not_ge hn⟩
  have hrbound : ∀ t ∈ Set.Ioo a b,
      radialSpeed q (q - 1 + 2 * flow.mass t / flow.radius t) (Real.cos (flow.angle t)) ≤ q := by
    intro t ht
    have hr := flow.radius_pos t ht.1.le
    have hm := flow.mass_range t ht.1.le
    have hm0 := hm.1
    have hrbig := hrR t (Set.Ioo_subset_Icc_self ht)
    have hdiv : 0 ≤ 2 * flow.mass t / flow.radius t := by positivity
    have hdivup : 2 * flow.mass t / flow.radius t ≤ 1 := (div_le_iff₀ hr).2 (by nlinarith)
    exact radialSpeed_le_charge hqpos (by linarith) (by linarith) (Real.cos_le_one _)
  have hαbound : ∀ t ∈ Set.Ioo a b,
      Real.sin (flow.angle t) / flow.radius t *
        directionB q (q - 1 + 2 * flow.mass t / flow.radius t) (Real.cos (flow.angle t)) ≤
      -(κ * Real.sin β) / flow.radius t := by
    intro t ht
    have htcc := Set.Ioo_subset_Icc_self ht
    have hr := flow.radius_pos t ht.1.le
    have hsinle := Real.sin_le_sin_of_le_of_le_pi_div_two
      (show -(Real.pi / 2) ≤ β by linarith [Real.pi_pos])
      ((hangle t htcc).trans hanglea.le) (hangleβ t htcc)
    have hdivle := div_le_div_of_nonneg_right hsinle hr.le
    have hfac := (div_pos (Real.sin_pos_of_pos_of_lt_pi (flow.angle_range t ht.1.le).1
      (flow.angle_range t ht.1.le).2) hr).le
    have hfirst := mul_le_mul_of_nonneg_left (hBpoint t htcc) hfac
    have hsecond := mul_le_mul_of_nonpos_right hdivle (show -κ ≤ 0 by linarith)
    calc
      _ ≤ Real.sin (flow.angle t) / flow.radius t * (-κ) := hfirst
      _ ≤ Real.sin β / flow.radius t * (-κ) := hsecond
      _ = -(κ * Real.sin β) / flow.radius t := by ring
  have hh := excursion_log_bound hab (mul_pos hκ hsinβ) hqpos
    (fun t ht => flow.radius_pos t ht.1) (fun t ht => flow.radial_ode t ht.1)
    (fun t ht => flow.angular_ode t ht.1) hrbound hαbound
    (flow.angle_range a le_rfl).2.le (flow.angle_range b hab).1.le
  exact (not_lt_of_ge hh) hlarge

end ViaB
