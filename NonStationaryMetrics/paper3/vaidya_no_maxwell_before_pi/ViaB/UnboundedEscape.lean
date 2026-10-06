import ViaB.LastExit

open Filter Topology

namespace ViaB

theorem exists_outgoing_cone_angle (q : ℝ) :
    ∃ β : ℝ, (0 < β ∧ β < Real.pi / 2) ∧ q - 1 < q * Real.cos β ^ 2 := by
  have hc : ContinuousAt (fun β => q * Real.cos β ^ 2) 0 := by fun_prop
  have hev : ∀ᶠ β in 𝓝 (0 : ℝ), q - 1 < q * Real.cos β ^ 2 :=
    hc.eventually (eventually_gt_nhds (by simp))
  obtain ⟨ε, heps, hε⟩ := Metric.eventually_nhds_iff.mp hev
  obtain ⟨β, hβ, hβmin⟩ := exists_between (lt_min heps (half_pos Real.pi_pos))
  have hd : dist β (0 : ℝ) < ε := by
    rw [Real.dist_eq, sub_zero, abs_of_pos hβ]
    exact hβmin.trans_le (min_le_left _ _)
  exact ⟨β, ⟨hβ, hβmin.trans_le (min_le_right _ _)⟩, hε hd⟩

/-- Unboundedness is a weak hypothesis: arbitrarily large radii at some future
times. The conclusion is permanent escape, for the declared global ODE solution. -/
theorem unbounded_radius_escapes {q M a : ℝ} (flow : DirectionFlow q M a)
    (hq : 1 < q) (hM : 0 < M)
    (hunbounded : ∀ K : ℝ, ∃ t ≥ a, K < flow.radius t) :
    Tendsto flow.radius atTop atTop := by
  obtain ⟨β, hβ, hcone⟩ := exists_outgoing_cone_angle q
  obtain ⟨R, hR, κ, hκ, δ, hδ, hmargins⟩ := outgoing_cone_margins hq hM hβ hcone
  let S := max R (flow.radius a) + 1
  have hRS : R ≤ S := by dsimp [S]; linarith [le_max_left R (flow.radius a)]
  have hraS : flow.radius a < S := by dsimp [S]; linarith [le_max_right R (flow.radius a)]
  have hS : 2 * M < S := hR.trans_le hRS
  have hSpos : 0 < S := by linarith
  have hmarginsS : ∀ r m c : ℝ, S ≤ r → 0 ≤ m → m ≤ M →
      (c ≤ 1 → directionB q (q - 1 + 2 * m / r) c ≤ -κ) ∧
      (Real.cos β ≤ c → δ ≤ radialSpeed q (q - 1 + 2 * m / r) c) :=
    fun r m c hr hm hmM => hmargins r m c (hRS.trans hr) hm hmM
  have hsinβ : 0 < Real.sin β := Real.sin_pos_of_pos_of_lt_pi hβ.1 (by linarith [Real.pi_pos])
  let k := κ * Real.sin β
  have hk : 0 < k := mul_pos hκ hsinβ
  let L := Real.log S + q * Real.pi / k + 1
  obtain ⟨b, hab, hb⟩ := hunbounded (max S (Real.exp L))
  have hrbS : S < flow.radius b := (le_max_left _ _).trans_lt hb
  have hrexpb : Real.exp L < flow.radius b := (le_max_right _ _).trans_lt hb
  have hLb : L < Real.log (flow.radius b) :=
    (Real.lt_log_iff_exp_lt (flow.radius_pos b hab)).mpr hrexpb
  have hdiff : q * Real.pi / k < Real.log (flow.radius b) - Real.log S := by
    dsimp [L] at hLb
    linarith
  have hlarge : q * Real.pi < κ * Real.sin β * (Real.log (flow.radius b) - Real.log S) := by
    have hh := (div_lt_iff₀ hk).mp hdiff
    simpa only [k, mul_comm] using hh
  have hcont : ContinuousOn flow.radius (Set.Icc a b) :=
    fun t ht => (flow.radial_ode t ht.1).continuousAt.continuousWithinAt
  obtain ⟨d, hd, hrd, hafter⟩ := exists_last_exit hab hcont hraS hrbS
  have houtexit := derivative_nonneg_at_right_exit hd.2 (flow.radial_ode d hd.1) hrd
    (fun t ht => hafter t ⟨ht.1, ht.2.le⟩)
  have hrS : ∀ t ∈ Set.Icc d b, S ≤ flow.radius t := by
    intro t ht
    rcases ht.1.eq_or_lt with heq | hdt
    · rw [← heq, hrd]
    · exact (hafter t ⟨hdt, ht.2⟩).le
  have hlarge' : q * Real.pi < κ * Real.sin β *
      (Real.log (flow.radius b) - Real.log (flow.radius d)) := by rwa [hrd]
  obtain ⟨t, ht, hαt⟩ := large_excursion_enters_angular_cone (flow.restrict d hd.1) hq
    hd.2.le hS hκ hβ (fun r m c hr hm hmM hc => (hmarginsS r m c hr hm hmM).1 hc)
    hrS houtexit hlarge'
  exact cone_escape_from_nonstrict_radius (flow.restrict t (hd.1.trans ht.1)) hβ hκ hδ
    hmarginsS (hrS t ht) hαt

end ViaB
