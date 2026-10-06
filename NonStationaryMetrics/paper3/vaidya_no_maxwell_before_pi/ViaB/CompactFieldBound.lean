import ViaB.TubeBootstrap

open Filter Topology
open scoped NNReal

namespace ViaB

theorem finite_positive_lower_bound {ι : Type*} (s : Finset ι) (ρ : ι → ℝ)
    (hρ : ∀ i ∈ s, 0 < ρ i) : ∃ ε > 0, ∀ i ∈ s, ε ≤ ρ i := by
  classical
  induction s using Finset.induction_on with
  | empty => exact ⟨1, by norm_num, by simp⟩
  | @insert i s hi ih =>
    obtain ⟨ε, hε, hbound⟩ := ih (fun j hj => hρ j (Finset.mem_insert_of_mem hj))
    refine ⟨min ε (ρ i), lt_min hε (hρ i (Finset.mem_insert_self _ _)), ?_⟩
    intro j hj
    rcases Finset.mem_insert.mp hj with rfl | hj
    · exact min_le_right _ _
    · exact (min_le_left _ _).trans (hbound j hj)

/-- Compactness turns pointwise local Lipschitz neighborhoods into one
constant and one positive tube width around the whole compact set. -/
theorem compact_uniform_local_lipschitz {X Y : Type*}
    [PseudoMetricSpace X] [PseudoMetricSpace Y] {C : Set X} (f : X → Y)
    (hC : IsCompact C)
    (hloc : ∀ z ∈ C, ∃ K : ℝ≥0, ∃ U ∈ 𝓝 z, LipschitzOnWith K f U) :
    ∃ K : ℝ≥0, ∃ ε > 0, ∀ z ∈ C, ∀ y,
      dist y z ≤ ε → dist (f y) (f z) ≤ (K : ℝ) * dist y z := by
  classical
  have hballs : ∀ z : C, ∃ K : ℝ≥0, ∃ ρ > 0,
      LipschitzOnWith K f (Metric.ball (z : X) ρ) := by
    intro z
    obtain ⟨K, U, hU, hLip⟩ := hloc z z.property
    obtain ⟨ρ, hρ, hball⟩ := Metric.mem_nhds_iff.mp hU
    exact ⟨K, ρ, hρ, hLip.mono hball⟩
  choose L ρ hρ hLip using hballs
  obtain ⟨s, hcover⟩ := hC.elim_finite_subcover
    (fun z : C => Metric.ball (z : X) (ρ z / 3))
    (fun _ => Metric.isOpen_ball) (fun z hz => by
      exact Set.mem_iUnion.mpr ⟨⟨z, hz⟩,
        Metric.mem_ball_self (div_pos (hρ ⟨z, hz⟩) (by norm_num))⟩)
  obtain ⟨ε, hε, hwidth⟩ := finite_positive_lower_bound s (fun i => ρ i / 3)
    (fun i _ => div_pos (hρ i) (by norm_num))
  let K := s.sup L
  refine ⟨K, ε, hε, ?_⟩
  intro z hz y hy
  obtain ⟨i, hi, hzi⟩ := Set.mem_iUnion₂.mp (hcover hz)
  have hzi' : dist z (i : X) < ρ i / 3 := hzi
  have hwi := hwidth i hi
  have hyi : y ∈ Metric.ball (i : X) (ρ i) := by
    have ht := dist_triangle y z (i : X)
    change dist y (i : X) < ρ i
    linarith [hρ i]
  have hzir : z ∈ Metric.ball (i : X) (ρ i) := by
    change dist z (i : X) < ρ i
    linarith [hρ i]
  have hL : L i ≤ K := Finset.le_sup hi
  exact ((hLip i).dist_le_mul y hyi z hzir).trans
    (mul_le_mul_of_nonneg_right (show (L i : ℝ) ≤ (K : ℝ) from hL) dist_nonneg)

/-- Uniform local bound along any continuous regular finite-time reference.
No ODE, neighboring trajectories, or future existence is assumed. -/
theorem finite_reference_uniform_field_bound {q a b : ℝ}
    (x : ℝ → ℝ × ℝ) (m : ℝ → ℝ) (hq : 1 < q)
    (hxcont : ContinuousOn x (Set.Icc a b)) (hmcont : ContinuousOn m (Set.Icc a b))
    (hr : ∀ t ∈ Set.Icc a b, 0 < (x t).1)
    (hm : ∀ t ∈ Set.Icc a b, 0 ≤ m t) :
    ∃ K : ℝ≥0, ∃ ε > 0, ∀ t ∈ Set.Icc a b, ∀ y : ℝ × ℝ,
      dist y (x t) ≤ ε →
      dist (directionField q (m t) y) (directionField q (m t) (x t)) ≤
        (K : ℝ) * dist y (x t) := by
  let z := fun t => (x t, m t)
  have hz : ContinuousOn z (Set.Icc a b) := hxcont.prodMk hmcont
  have hcompact : IsCompact (z '' Set.Icc a b) := isCompact_Icc.image_of_continuousOn hz
  obtain ⟨K, ε, hε, hbound⟩ := compact_uniform_local_lipschitz
    (fun x : (ℝ × ℝ) × ℝ => directionField q x.2 x.1) hcompact (by
      intro x hx
      obtain ⟨t, ht, rfl⟩ := hx
      exact directionField_locally_lipschitz hq (hr t ht) (hm t ht))
  refine ⟨K, ε, hε, ?_⟩
  intro t ht y hy
  have heq : dist (y, m t) (z t) = dist y (x t) := by
    simp [z, Prod.dist_eq, dist_self]
  have hh := hbound (z t) ⟨t, ht, rfl⟩ (y, m t) (by rwa [heq])
  simpa only [heq] using hh

/-- Uniform local state bound along one regular finite-time solution.
Only the reference trajectory and its continuous mass history enter here. -/
theorem directionFlow_uniform_field_bound {q M a b : ℝ} (flow : DirectionFlow q M a)
    (hq : 1 < q) (hmcont : ContinuousOn flow.mass (Set.Icc a b)) :
    ∃ K : ℝ≥0, ∃ ε > 0, ∀ t ∈ Set.Icc a b, ∀ y : ℝ × ℝ,
      dist y (flow.radius t, flow.angle t) ≤ ε →
      dist (directionField q (flow.mass t) y)
        (directionField q (flow.mass t) (flow.radius t, flow.angle t)) ≤
        (K : ℝ) * dist y (flow.radius t, flow.angle t) := by
  exact finite_reference_uniform_field_bound
    (fun t => (flow.radius t, flow.angle t)) flow.mass hq
    ((HasDerivAt.continuousOn (fun t ht => flow.radial_ode t ht.1)).prodMk
      (HasDerivAt.continuousOn (fun t ht => flow.angular_ode t ht.1))) hmcont
    (fun t ht => flow.radius_pos t ht.1) (fun t ht => (flow.mass_range t ht.1).1)

/-- Both the tube width and the finite-time dependence constant are derived
from the reference flow. Tube membership of the second flow is a conclusion. -/
theorem directionFlow_finite_time_dependence {q M a b : ℝ} (flow : DirectionFlow q M a)
    (hq : 1 < q) (hab : a ≤ b) (hmcont : ContinuousOn flow.mass (Set.Icc a b)) :
    ∃ K : ℝ≥0, ∃ ε > 0, ∀ other : DirectionFlow q M a,
      (∀ t ∈ Set.Icc a b, other.mass t = flow.mass t) →
      dist (other.radius a, other.angle a) (flow.radius a, flow.angle a) *
        Real.exp ((K : ℝ) * (b - a)) < ε →
      ∀ t ∈ Set.Icc a b,
        dist (other.radius t, other.angle t) (flow.radius t, flow.angle t) < ε ∧
        dist (other.radius t, other.angle t) (flow.radius t, flow.angle t) ≤
          dist (other.radius a, other.angle a) (flow.radius a, flow.angle a) *
            Real.exp ((K : ℝ) * (t - a)) := by
  obtain ⟨K, ε, hε, hfield⟩ := directionFlow_uniform_field_bound flow hq hmcont
  refine ⟨K, ε, hε, ?_⟩
  intro other hshared hsmall
  have hf : ∀ t ∈ Set.Icc a b,
      HasDerivAt (fun v => (other.radius v, other.angle v))
        (directionField q (flow.mass t) (other.radius t, other.angle t)) t := by
    intro t ht
    have hh := directionFlow_state_derivative other ht.1
    rwa [hshared t ht] at hh
  have hg := fun t (ht : t ∈ Set.Icc a b) => directionFlow_state_derivative flow ht.1
  have hgap : ∀ t ∈ Set.Ico a b,
      dist (other.radius t, other.angle t) (flow.radius t, flow.angle t) ≤ ε →
      ‖directionField q (flow.mass t) (other.radius t, other.angle t) -
        directionField q (flow.mass t) (flow.radius t, flow.angle t)‖ ≤
        (K : ℝ) * ‖(other.radius t, other.angle t) - (flow.radius t, flow.angle t)‖ := by
    intro t ht hclose
    simpa only [dist_eq_norm] using hfield t ⟨ht.1, ht.2.le⟩ _ hclose
  have htube := trajectory_stays_in_tube hab K.coe_nonneg hf hg hgap hsmall
  have hbound := trajectory_gap_bound hf hg
    (fun t ht => hgap t ht (htube t ⟨ht.1, ht.2.le⟩).le)
  exact fun t ht => ⟨htube t ht, hbound t ht⟩

/-- Finite-time evaluation continuity follows from initial continuity and
the continuous reference mass; no uniform field bound is assumed. -/
theorem directionFlow_eval_continuousAt_from_initial {P : Type*} [TopologicalSpace P]
    {q M a b : ℝ} (flows : P → DirectionFlow q M a) (p : P)
    (hq : 1 < q) (hab : a ≤ b)
    (hmcont : ContinuousOn (flows p).mass (Set.Icc a b))
    (hshared : ∀ x t, t ∈ Set.Icc a b → (flows x).mass t = (flows p).mass t)
    (hinitial : ContinuousAt (fun x => ((flows x).radius a, (flows x).angle a)) p) :
    ContinuousAt (fun x => ((flows x).radius b, (flows x).angle b)) p := by
  obtain ⟨K, ε, hε, hfield⟩ := directionFlow_uniform_field_bound (flows p) hq hmcont
  apply trajectory_eval_continuousAt_of_local_gap
    (f := fun x t => ((flows x).radius t, (flows x).angle t))
    (df := fun x t => directionField q ((flows p).mass t) ((flows x).radius t, (flows x).angle t))
    p hab K.coe_nonneg hε
  · intro x t ht
    have hh := directionFlow_state_derivative (flows x) ht.1
    rwa [hshared x t ht] at hh
  · exact Filter.Eventually.of_forall (fun x t ht hclose => by
      simpa only [dist_eq_norm] using hfield t ⟨ht.1, ht.2.le⟩ _ hclose)
  · exact hinitial

/-- Strict cone entry of the reference flow is stable from initial data.
The uniform local field bound and the tube width are derived internally. -/
theorem cone_escape_stable_from_regular_reference {P : Type*} [TopologicalSpace P]
    {q M a s β R κ δ : ℝ} (flows : P → DirectionFlow q M a) (p : P)
    (hq : 1 < q) (has : a ≤ s)
    (hβ : 0 < β ∧ β < Real.pi / 2) (hκ : 0 < κ) (hδ : 0 < δ)
    (hmargins : ∀ r m c : ℝ, R ≤ r → 0 ≤ m → m ≤ M →
      (c ≤ 1 → directionB q (q - 1 + 2 * m / r) c ≤ -κ) ∧
      (Real.cos β ≤ c → δ ≤ radialSpeed q (q - 1 + 2 * m / r) c))
    (hrp : R < (flows p).radius s) (hαp : (flows p).angle s < β)
    (hmcont : ContinuousOn (flows p).mass (Set.Icc a s))
    (hshared : ∀ x t, t ∈ Set.Icc a s → (flows x).mass t = (flows p).mass t)
    (hinitial : ContinuousAt (fun x => ((flows x).radius a, (flows x).angle a)) p) :
    ∀ᶠ x in 𝓝 p,
      (∀ t ≥ s, δ * (t - s) ≤ (flows x).radius t - (flows x).radius s) ∧
      Tendsto (flows x).radius atTop atTop := by
  have hc := directionFlow_eval_continuousAt_from_initial flows p hq has hmcont hshared hinitial
  exact cone_entry_stable_in_family flows p has hβ hκ hδ hmargins hc.fst hc.snd hrp hαp

/-- Openness of escape in a family of existing global flows, derived from
initial continuity rather than assuming continuity at every finite time. -/
theorem escaping_set_isOpen_from_initial {P : Type*} [TopologicalSpace P]
    {q M a : ℝ} (flows : P → DirectionFlow q M a) (hq : 1 < q) (hM : 0 < M)
    (hmcont : ∀ p b, a ≤ b → ContinuousOn (flows p).mass (Set.Icc a b))
    (hshared : ∀ x p t, a ≤ t → (flows x).mass t = (flows p).mass t)
    (hinitial : Continuous (fun x => ((flows x).radius a, (flows x).angle a))) :
    IsOpen {p : P | Tendsto (flows p).radius atTop atTop} := by
  have hc : ∀ s ≥ a, Continuous (fun x => ((flows x).radius s, (flows x).angle s)) := by
    intro s has
    rw [continuous_iff_continuousAt]
    intro p
    exact directionFlow_eval_continuousAt_from_initial flows p hq has (hmcont p s has)
      (fun x t ht => hshared x p t ht.1) hinitial.continuousAt
  exact escaping_set_isOpen flows hq hM
    (fun s has => (hc s has).fst) (fun s has => (hc s has).snd)

end ViaB
