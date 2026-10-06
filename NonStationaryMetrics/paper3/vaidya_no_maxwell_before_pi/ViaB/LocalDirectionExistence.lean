import ViaB.FlowDependence
import ViaB.TimeDependentHamiltonian
import Mathlib.Analysis.ODE.PicardLindelof

open Filter Set
open scoped Topology

namespace ViaB

/-- Add the time coordinate to make the nonautonomous field autonomous. -/
noncomputable def liftedDirectionField (q : ℝ) (m : ℝ → ℝ)
    (z : (ℝ × ℝ) × ℝ) : (ℝ × ℝ) × ℝ :=
  (directionField q (m z.2) z.1, 1)

theorem liftedDirectionField_contDiffAt {q a : ℝ} {m : ℝ → ℝ} {x : ℝ × ℝ}
    (hq : 1 < q) (hr : 0 < x.1) (hm0 : 0 ≤ m a)
    (hm : ContDiffAt ℝ 1 m a) :
    ContDiffAt ℝ 1 (liftedDirectionField q m) (x, a) := by
  have hmap : ContDiffAt ℝ 1 (fun z : (ℝ × ℝ) × ℝ => (z.1, m z.2)) (x, a) :=
    contDiffAt_fst.prodMk (hm.comp (x, a) contDiffAt_snd)
  exact ((directionField_contDiffAt hq hr hm0).comp (x, a) hmap).prodMk contDiffAt_const

/-- Local existence for the nonstationary direction system, obtained from
Picard-Lindelof. No pre-existing trajectory is a premise. -/
theorem direction_local_solution_exists {q a : ℝ} {m : ℝ → ℝ} {x : ℝ × ℝ}
    (hq : 1 < q) (hr : 0 < x.1) (hm0 : 0 ≤ m a)
    (hm : ContDiffAt ℝ 1 m a) :
    ∃ f : ℝ → ℝ × ℝ, f a = x ∧ ∃ ε > 0,
      ∀ t ∈ Ioo (a - ε) (a + ε),
        HasDerivAt f (directionField q (m t) (f t)) t := by
  obtain ⟨α, hbase, ε, hε, hode⟩ :=
    ContDiffAt.exists_forall_mem_closedBall_exists_eq_forall_mem_Ioo_hasDerivAt₀
      (liftedDirectionField_contDiffAt hq hr hm0 hm) a
  have ha : a ∈ Ioo (a - ε) (a + ε) := by constructor <;> linarith
  have htime : ∀ t ∈ Ioo (a - ε) (a + ε), (α t).2 = t := by
    intro t ht
    have hz := hamiltonian_null_level_propagates
      (H := fun s => (α s).2 - s) (convex_Ioo (a - ε) (a + ε)) ha ht
      (fun s hs => by
        have hd : HasDerivAt (fun u => (α u).2) 1 s := by
          simpa [liftedDirectionField] using
            (ContinuousLinearMap.snd ℝ (ℝ × ℝ) ℝ).hasFDerivAt.comp_hasDerivAt s (hode s hs)
        simpa only [sub_self] using hd.sub (hasDerivAt_id s))
      (by simp only [hbase, sub_self])
    exact sub_eq_zero.mp hz
  refine ⟨fun t => (α t).1, congrArg Prod.fst hbase, ε, hε, ?_⟩
  intro t ht
  have hd := (ContinuousLinearMap.fst ℝ (ℝ × ℝ) ℝ).hasFDerivAt.comp_hasDerivAt t (hode t ht)
  simpa [liftedDirectionField, htime t ht] using hd


/-- Local uniqueness for a prescribed continuous mass history. The spatial
Lipschitz bound is derived at the regular launch, not assumed separately. -/
theorem direction_local_solution_unique {q a : ℝ} {m : ℝ → ℝ}
    {f g : ℝ → ℝ × ℝ}
    (hq : 1 < q) (hr : 0 < (f a).1) (hm0 : 0 ≤ m a)
    (hm : ContinuousAt m a)
    (hf : ∀ᶠ t in 𝓝 a, HasDerivAt f (directionField q (m t) (f t)) t)
    (hg : ∀ᶠ t in 𝓝 a, HasDerivAt g (directionField q (m t) (g t)) t)
    (hbase : f a = g a) : f =ᶠ[𝓝 a] g := by
  obtain ⟨K, U, hU, hLip⟩ := directionField_locally_lipschitz (z := (f a, m a)) hq hr hm0
  have hfc : ContinuousAt f a := hf.self_of_nhds.continuousAt
  have hgc : ContinuousAt g a := hg.self_of_nhds.continuousAt
  have hfe : ∀ᶠ t in 𝓝 a, (f t, m t) ∈ U := (hfc.prodMk hm).eventually hU
  have hge : ∀ᶠ t in 𝓝 a, (g t, m t) ∈ U := by
    have hUg : U ∈ 𝓝 (g a, m a) := by rwa [← hbase]
    exact (hgc.prodMk hm).eventually hUg
  have hl : ∀ t : ℝ, LipschitzOnWith K (directionField q (m t))
      {y : ℝ × ℝ | (y, m t) ∈ U} := by
    intro t
    apply LipschitzOnWith.of_dist_le_mul
    intro y hy z hz
    have hdist : dist (y, m t) (z, m t) = dist y z := by
      rw [Prod.dist_eq, dist_self, max_eq_left dist_nonneg]
    simpa only [hdist] using hLip.dist_le_mul (y, m t) hy (z, m t) hz
  exact ODE_solution_unique_of_eventually (Eventually.of_forall hl)
    (hf.and hfe) (hg.and hge) hbase

/-- A strictly exterior nonradial launch has a local solution staying in the
physical open domain. Positivity on the interval is a conclusion. -/
theorem direction_exterior_local_solution_exists {q a : ℝ} {m : ℝ → ℝ}
    {x : ℝ × ℝ} (hq : 1 < q) (hmpos : 0 < m a)
    (hext : 2 * m a < x.1) (hα : 0 < x.2 ∧ x.2 < Real.pi)
    (hm : ContDiffAt ℝ 1 m a) :
    ∃ f : ℝ → ℝ × ℝ, f a = x ∧ ∃ ε > 0,
      ∀ t ∈ Ioo (a - ε) (a + ε),
        HasDerivAt f (directionField q (m t) (f t)) t ∧
        0 < m t ∧ 2 * m t < (f t).1 ∧ 0 < (f t).2 ∧ (f t).2 < Real.pi := by
  have hr : 0 < x.1 := by linarith
  obtain ⟨f, hbase, ε, hε, hode⟩ := direction_local_solution_exists hq hr hmpos.le hm
  have ha : a ∈ Ioo (a - ε) (a + ε) := by constructor <;> linarith
  have hc := (hode a ha).continuousAt
  have hmp : ∀ᶠ t in 𝓝 a, 0 < m t :=
    continuousAt_const.eventually_lt hm.continuousAt hmpos
  have hrp : ∀ᶠ t in 𝓝 a, 2 * m t < (f t).1 :=
    (hm.continuousAt.const_mul 2).eventually_lt hc.fst (by simpa only [hbase] using hext)
  have hαp : ∀ᶠ t in 𝓝 a, 0 < (f t).2 :=
    continuousAt_const.eventually_lt hc.snd (by simpa only [hbase] using hα.1)
  have hαπ : ∀ᶠ t in 𝓝 a, (f t).2 < Real.pi :=
    hc.snd.eventually_lt continuousAt_const (by simpa only [hbase] using hα.2)
  have hprops : {t | 0 < m t ∧ 2 * m t < (f t).1 ∧
      0 < (f t).2 ∧ (f t).2 < Real.pi} ∈ 𝓝 a := by
    filter_upwards [hmp, hrp, hαp, hαπ] with t ht₁ ht₂ ht₃ ht₄
    exact ⟨ht₁, ht₂, ht₃, ht₄⟩
  obtain ⟨δ, hδ, hball⟩ := Metric.mem_nhds_iff.mp hprops
  refine ⟨f, hbase, min ε δ, lt_min hε hδ, ?_⟩
  intro t ht
  have htε : t ∈ Ioo (a - ε) (a + ε) := by
    constructor <;> linarith [min_le_left ε δ, ht.1, ht.2]
  have htδ : t ∈ Metric.ball a δ := by
    rw [Real.ball_eq_Ioo]
    constructor <;> linarith [min_le_right ε δ, ht.1, ht.2]
  exact ⟨hode t htε, hball htδ⟩

/-- Existence and uniqueness of the local solution germ for regular data.
The solution and its Lipschitz neighborhood are both constructed. -/
theorem direction_local_initial_value_problem {q a : ℝ} {m : ℝ → ℝ} {x : ℝ × ℝ}
    (hq : 1 < q) (hr : 0 < x.1) (hm0 : 0 ≤ m a)
    (hm : ContDiffAt ℝ 1 m a) :
    ∃ f : ℝ → ℝ × ℝ, f a = x ∧
      (∀ᶠ t in 𝓝 a, HasDerivAt f (directionField q (m t) (f t)) t) ∧
      ∀ g : ℝ → ℝ × ℝ, g a = x →
        (∀ᶠ t in 𝓝 a, HasDerivAt g (directionField q (m t) (g t)) t) →
        f =ᶠ[𝓝 a] g := by
  obtain ⟨f, hbase, ε, hε, hode⟩ := direction_local_solution_exists hq hr hm0 hm
  have hf : ∀ᶠ t in 𝓝 a, HasDerivAt f (directionField q (m t) (f t)) t := by
    filter_upwards [Ioo_mem_nhds (show a - ε < a by linarith)
      (show a < a + ε by linarith)] with t ht
    exact hode t ht
  refine ⟨f, hbase, hf, ?_⟩
  intro g hg hgode
  exact direction_local_solution_unique hq (by simpa only [hbase] using hr) hm0
    hm.continuousAt hf hgode (hbase.trans hg.symm)

end ViaB
