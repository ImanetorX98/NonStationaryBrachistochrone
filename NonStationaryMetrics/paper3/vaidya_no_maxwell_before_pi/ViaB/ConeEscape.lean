import ViaB.Dynamics

open Filter Topology

namespace ViaB

/-- The real ODE system to be connected later to geometric null geodesics.
The future domain and the nonradial angle range are explicit assumptions. -/
structure DirectionFlow (q M a : ℝ) where
  radius : ℝ → ℝ
  angle : ℝ → ℝ
  mass : ℝ → ℝ
  radius_pos : ∀ t ≥ a, 0 < radius t
  angle_range : ∀ t ≥ a, 0 < angle t ∧ angle t < Real.pi
  mass_range : ∀ t ≥ a, 0 ≤ mass t ∧ mass t ≤ M
  radial_ode : ∀ t ≥ a, HasDerivAt radius
    (radialSpeed q (q - 1 + 2 * mass t / radius t) (Real.cos (angle t))) t
  angular_ode : ∀ t ≥ a, HasDerivAt angle
    (Real.sin (angle t) / radius t * directionB q (q - 1 + 2 * mass t / radius t)
      (Real.cos (angle t))) t

theorem outgoing_cone_margins {q M β : ℝ} (hq : 1 < q) (hM : 0 < M)
    (hβ : 0 < β ∧ β < Real.pi / 2) (hcone : q - 1 < q * Real.cos β ^ 2) :
    ∃ R > 2 * M, ∃ κ > 0, ∃ δ > 0, ∀ r m c : ℝ,
      R ≤ r → 0 ≤ m → m ≤ M →
      (c ≤ 1 → directionB q (q - 1 + 2 * m / r) c ≤ -κ) ∧
      (Real.cos β ≤ c → δ ≤ radialSpeed q (q - 1 + 2 * m / r) c) := by
  have hcβ : 0 < Real.cos β := Real.cos_pos_of_mem_Ioo ⟨by linarith [Real.pi_pos], hβ.2⟩
  obtain ⟨εB, hεB, κ, hk, hB⟩ := directionB_uniform_far hq
  obtain ⟨εS, hεS, δ, hδ, hS⟩ := radialSpeed_uniform_far hq hcβ hcone
  let ε := min εB εS
  have heps : 0 < ε := lt_min hεB hεS
  obtain ⟨R, hR⟩ := exists_gt (max (2 * M) (2 * M / ε))
  have hRM : 2 * M < R := (le_max_left _ _).trans_lt hR
  have hReps : 2 * M / ε < R := (le_max_right _ _).trans_lt hR
  have hmass : 2 * M < R * ε := (div_lt_iff₀ heps).mp hReps
  refine ⟨R, hRM, κ, hk, δ, hδ, ?_⟩
  intro r m c hr hm hmM
  have hrpos : 0 < r := by linarith
  have hdiv : 0 ≤ 2 * m / r := by positivity
  have hdivlt : 2 * m / r < ε := (div_lt_iff₀ hrpos).2 (by nlinarith)
  have hw : q - 1 ≤ q - 1 + 2 * m / r := by linarith
  constructor
  · intro hc
    exact hB _ _ hw (by have hh := min_le_left εB εS; dsimp [ε] at hdivlt; linarith) hc
  · intro hc
    exact hS _ _ hw (by have hh := min_le_right εB εS; dsimp [ε] at hdivlt; linarith) hc

theorem outgoing_cone_invariant {q M a β R κ δ : ℝ} (flow : DirectionFlow q M a)
    (hβ : 0 < β ∧ β < Real.pi / 2) (hκ : 0 < κ) (hδ : 0 < δ)
    (hmargins : ∀ r m c : ℝ, R ≤ r → 0 ≤ m → m ≤ M →
      (c ≤ 1 → directionB q (q - 1 + 2 * m / r) c ≤ -κ) ∧
      (Real.cos β ≤ c → δ ≤ radialSpeed q (q - 1 + 2 * m / r) c))
    (hra : R < flow.radius a) (hαa : flow.angle a < β) :
    ∀ t ≥ a, R < flow.radius t ∧ flow.angle t < β := by
  intro t hat
  have hb := positive_pair_invariant
    (a := a) (b := t)
    (f := fun x => flow.radius x - R) (g := fun x => β - flow.angle x)
    (fun x hx => (flow.radial_ode x hx.1).sub_const R)
    (fun x hx => (flow.angular_ode x hx.1).const_sub β)
    (sub_pos.mpr hra) (sub_pos.mpr hαa)
    (fun x hx hzero hother => ?_) (fun x hx hzero hother => ?_)
  · obtain ⟨hr, hα⟩ := hb t ⟨hat, le_rfl⟩
    exact ⟨by linarith, by linarith⟩
  · have hxr : flow.radius x = R := by linarith
    have hangle : flow.angle x ≤ β := by linarith
    have hc := Real.cos_le_cos_of_nonneg_of_le_pi (flow.angle_range x hx.1).1.le
      (show β ≤ Real.pi by linarith [Real.pi_pos]) hangle
    have hmass := flow.mass_range x hx.1
    have hmargin := (hmargins (flow.radius x) (flow.mass x) (Real.cos (flow.angle x))
      (by linarith) hmass.1 hmass.2).2 hc
    exact hδ.trans_le hmargin
  · have hangle : flow.angle x = β := by linarith
    have hrR : R ≤ flow.radius x := by linarith
    have hmass := flow.mass_range x hx.1
    have hB := (hmargins (flow.radius x) (flow.mass x) (Real.cos (flow.angle x))
      hrR hmass.1 hmass.2).1 (Real.cos_le_one (flow.angle x))
    have hBneg : directionB q (q - 1 + 2 * flow.mass x / flow.radius x)
        (Real.cos (flow.angle x)) < 0 := by linarith
    have hsin : 0 < Real.sin (flow.angle x) :=
      Real.sin_pos_of_pos_of_lt_pi (flow.angle_range x hx.1).1 (flow.angle_range x hx.1).2
    exact neg_pos.mpr (mul_neg_of_pos_of_neg (div_pos hsin (flow.radius_pos x hx.1)) hBneg)

theorem outgoing_cone_linear_growth {q M a β R κ δ : ℝ} (flow : DirectionFlow q M a)
    (hβ : 0 < β ∧ β < Real.pi / 2) (hκ : 0 < κ) (hδ : 0 < δ)
    (hmargins : ∀ r m c : ℝ, R ≤ r → 0 ≤ m → m ≤ M →
      (c ≤ 1 → directionB q (q - 1 + 2 * m / r) c ≤ -κ) ∧
      (Real.cos β ≤ c → δ ≤ radialSpeed q (q - 1 + 2 * m / r) c))
    (hra : R < flow.radius a) (hαa : flow.angle a < β) :
    ∀ t ≥ a, δ * (t - a) ≤ flow.radius t - flow.radius a := by
  have hinv := outgoing_cone_invariant flow hβ hκ hδ hmargins hra hαa
  intro t hat
  apply increment_ge_of_derivative_lower hat (fun x hx => flow.radial_ode x hx.1)
  intro x hx
  have hcone := hinv x hx.1.le
  have hmass := flow.mass_range x hx.1.le
  have hc := Real.cos_le_cos_of_nonneg_of_le_pi (flow.angle_range x hx.1.le).1.le
    (show β ≤ Real.pi by linarith [Real.pi_pos]) hcone.2.le
  exact (hmargins _ _ _ hcone.1.le hmass.1 hmass.2).2 hc

theorem radial_growth_tendsto_atTop {r : ℝ → ℝ} {a δ : ℝ} (hδ : 0 < δ)
    (hgrowth : ∀ t ≥ a, δ * (t - a) ≤ r t - r a) : Tendsto r atTop atTop := by
  apply tendsto_atTop.2
  intro K
  filter_upwards [eventually_ge_atTop (max a (a + (K - r a) / δ))] with t ht
  have hat : a ≤ t := (le_max_left _ _).trans ht
  have hdiv : (K - r a) / δ ≤ t - a := by have hh := (le_max_right _ _).trans ht; linarith
  have hK : K - r a ≤ (t - a) * δ := (div_le_iff₀ hδ).mp hdiv
  have hg := hgrowth t hat
  nlinarith

/-- Permanent escape after strict cone entry, for the declared ODE system. -/
theorem outgoing_cone_escape {q M a β : ℝ} (flow : DirectionFlow q M a)
    (hq : 1 < q) (hM : 0 < M) (hβ : 0 < β ∧ β < Real.pi / 2)
    (hcone : q - 1 < q * Real.cos β ^ 2) :
    ∃ R > 2 * M, R < flow.radius a → flow.angle a < β →
      (∀ t ≥ a, R < flow.radius t ∧ flow.angle t < β) ∧ Tendsto flow.radius atTop atTop := by
  obtain ⟨R, hR, κ, hk, δ, hδ, hmargins⟩ := outgoing_cone_margins hq hM hβ hcone
  refine ⟨R, hR, ?_⟩
  intro hra hαa
  exact ⟨outgoing_cone_invariant flow hβ hk hδ hmargins hra hαa,
    radial_growth_tendsto_atTop hδ (outgoing_cone_linear_growth flow hβ hk hδ hmargins hra hαa)⟩

end ViaB
