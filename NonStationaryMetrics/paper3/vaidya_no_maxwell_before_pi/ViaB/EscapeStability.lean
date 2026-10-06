import ViaB.BoundedLingering

open Filter Topology

namespace ViaB

/-- Arbitrarily large future excursions yield strict entry into any fixed
outgoing cone with the declared uniform margins. -/
theorem unbounded_radius_strict_cone_entry {q M a β R κ δ : ℝ}
    (flow : DirectionFlow q M a) (hq : 1 < q) (hRM : 2 * M < R)
    (hM : 0 < M) (hβ : 0 < β ∧ β < Real.pi / 2) (hκ : 0 < κ) (hδ : 0 < δ)
    (hmargins : ∀ r m c : ℝ, R ≤ r → 0 ≤ m → m ≤ M →
      (c ≤ 1 → directionB q (q - 1 + 2 * m / r) c ≤ -κ) ∧
      (Real.cos β ≤ c → δ ≤ radialSpeed q (q - 1 + 2 * m / r) c))
    (hunbounded : ∀ K : ℝ, ∃ t ≥ a, K < flow.radius t) :
    ∃ t ≥ a, R < flow.radius t ∧ flow.angle t < β := by
  let S := max R (flow.radius a) + 1
  have hRS : R ≤ S := by dsimp [S]; linarith [le_max_left R (flow.radius a)]
  have hraS : flow.radius a < S := by dsimp [S]; linarith [le_max_right R (flow.radius a)]
  have hS : 2 * M < S := hRM.trans_le hRS
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
  obtain ⟨s, hts, hrs, hαs⟩ := strict_cone_entry
    (flow.restrict t (hd.1.trans ht.1)) hβ hδ hmarginsS (hrS t ht) hαt
  exact ⟨s, (hd.1.trans ht.1).trans hts.le, hRS.trans_lt hrs, hαs⟩

/-- Finite-time evaluation continuity transfers a strict cone entry to
nearby members of a family of global solutions. -/
theorem cone_entry_stable_in_family {P : Type*} [TopologicalSpace P]
    {q M a β R κ δ s : ℝ} (flows : P → DirectionFlow q M a) (p : P)
    (has : a ≤ s) (hβ : 0 < β ∧ β < Real.pi / 2) (hκ : 0 < κ) (hδ : 0 < δ)
    (hmargins : ∀ r m c : ℝ, R ≤ r → 0 ≤ m → m ≤ M →
      (c ≤ 1 → directionB q (q - 1 + 2 * m / r) c ≤ -κ) ∧
      (Real.cos β ≤ c → δ ≤ radialSpeed q (q - 1 + 2 * m / r) c))
    (hcr : ContinuousAt (fun x => (flows x).radius s) p)
    (hca : ContinuousAt (fun x => (flows x).angle s) p)
    (hrp : R < (flows p).radius s) (hαp : (flows p).angle s < β) :
    ∀ᶠ x in 𝓝 p,
      (∀ t ≥ s, δ * (t - s) ≤ (flows x).radius t - (flows x).radius s) ∧
      Tendsto (flows x).radius atTop atTop := by
  have her := hcr.eventually (eventually_gt_nhds hrp)
  have hea := hca.eventually (eventually_lt_nhds hαp)
  filter_upwards [her, hea] with x hrx hαx
  have hg := outgoing_cone_linear_growth ((flows x).restrict s has)
    hβ hκ hδ hmargins hrx hαx
  exact ⟨hg, radial_growth_tendsto_atTop hδ hg⟩

/-- Escape is open in a family whose finite-time evaluations are continuous.
Existence and evaluation continuity of that family are explicit premises. -/
theorem escaping_set_isOpen {P : Type*} [TopologicalSpace P]
    {q M a : ℝ} (flows : P → DirectionFlow q M a) (hq : 1 < q) (hM : 0 < M)
    (hcr : ∀ s ≥ a, Continuous (fun x => (flows x).radius s))
    (hca : ∀ s ≥ a, Continuous (fun x => (flows x).angle s)) :
    IsOpen {p : P | Tendsto (flows p).radius atTop atTop} := by
  rw [isOpen_iff_mem_nhds]
  intro p hp
  obtain ⟨β, hβ, hcone⟩ := exists_outgoing_cone_angle q
  obtain ⟨R, hRM, κ, hκ, δ, hδ, hmargins⟩ := outgoing_cone_margins hq hM hβ hcone
  have hu : ∀ K : ℝ, ∃ t ≥ a, K < (flows p).radius t := by
    intro K
    have hev := (tendsto_atTop.1 hp) (K + 1)
    obtain ⟨T, hT⟩ := eventually_atTop.1 hev
    refine ⟨max a T, le_max_left _ _, ?_⟩
    have hh := hT (max a T) (le_max_right _ _)
    linarith
  obtain ⟨s, has, hrs, hαs⟩ := unbounded_radius_strict_cone_entry
    (flows p) hq hRM hM hβ hκ hδ hmargins hu
  have hev := cone_entry_stable_in_family flows p has hβ hκ hδ hmargins
    (hcr s has).continuousAt (hca s has).continuousAt hrs hαs
  filter_upwards [hev] with x hx
  exact hx.2

end ViaB
