import ViaB.ConstructedAffineClock
import ViaB.MaximalConeEscape

open Filter
open scoped Topology

namespace ViaB

/-- The complete coordinate Hamiltonian equations, including the time force,
the null constraint and the future orientation, at one affine event. -/
structure HamiltonianLiftAt (q L : ℝ) (m μ v R Φ pv p : ℝ → ℝ) (s : ℝ) : Prop where
  time : HasDerivAt v (-pv s / vaidyaW q (m (v s)) (R s) + p s) s
  radial : HasDerivAt R (pv s + (q - vaidyaW q (m (v s)) (R s)) * p s) s
  angular : HasDerivAt Φ (L / (R s) ^ 2) s
  time_momentum : HasDerivAt pv
    (-μ (v s) * massHamiltonianGradient q (m (v s)) (R s) (pv s) (p s)) s
  radial_momentum : HasDerivAt p
    (radialMomentumAffine q (m (v s)) (R s) (pv s) (p s) L) s
  angular_momentum : HasDerivAt (fun _ : ℝ => L) 0 s
  null : cometricQuad q (vaidyaW q (m (v s)) (R s)) (R s) (pv s) (p s) L = 0
  future : 0 < -pv s / vaidyaW q (m (v s)) (R s) + p s

theorem direction_angular_primitive_exists {q a : ℝ} {m r α : ℝ → ℝ}
    (hr : 0 < r a) (hα : 0 < α a ∧ α a < Real.pi)
    (hw : 0 < vaidyaW q (m a) (r a))
    (hm : ∀ᶠ t in 𝓝 a, ContinuousAt m t)
    (hrc : ∀ᶠ t in 𝓝 a, ContinuousAt r t)
    (hαc : ∀ᶠ t in 𝓝 a, ContinuousAt α t) :
    ∃ φ : ℝ → ℝ, φ a = 0 ∧ ∀ᶠ t in 𝓝 a,
      HasDerivAt φ (Real.sqrt (vaidyaW q (m t) (r t)) / r t * Real.sin (α t)) t := by
  let F : ℝ → ℝ := fun t => Real.sqrt (vaidyaW q (m t) (r t)) / r t * Real.sin (α t)
  have hF : ∀ᶠ t in 𝓝 a, ContinuousAt F t := by
    filter_upwards [hm, hrc, hαc,
      continuousAt_const.eventually_lt hrc.self_of_nhds hr] with t hmt hrt hαt hrpos
    have hcw : ContinuousAt (fun v => vaidyaW q (m v) (r v)) t :=
      continuousAt_const.add ((hmt.const_mul 2).div hrt (ne_of_gt hrpos))
    exact (hcw.sqrt.div hrt (ne_of_gt hrpos)).mul
      (Real.continuous_sin.continuousAt.comp hαt)
  have hs := Real.sin_pos_of_pos_of_lt_pi hα.1 hα.2
  have hpos : 0 < F a := mul_pos (div_pos (Real.sqrt_pos.2 hw) hr) hs
  obtain ⟨φ, _, hφ, _, hder, _⟩ := positive_clock_with_inverse_exists hF hpos
  exact ⟨φ, hφ, hder⟩

/-- Chain rule converts the reconstructed equations to the complete affine
Hamiltonian equations. The supplied time derivative is discharged by the
constructed inverse-clock theorem in the following germ result. -/
theorem direction_to_affine_hamiltonian_at {q L s : ℝ} {m μ r α φ v : ℝ → ℝ}
    (hq : 0 < q) (hL : 0 < L) (hr : 0 < r (v s))
    (hα : 0 < α (v s) ∧ α (v s) < Real.pi) (hw : 0 < vaidyaW q (m (v s)) (r (v s)))
    (hm : HasDerivAt m (μ (v s)) (v s))
    (hrder : HasDerivAt r (radialSpeed q (vaidyaW q (m (v s)) (r (v s)))
      (Real.cos (α (v s)))) (v s))
    (hαder : HasDerivAt α (Real.sin (α (v s)) / r (v s) *
      directionB q (vaidyaW q (m (v s)) (r (v s))) (Real.cos (α (v s)))) (v s))
    (hφ : HasDerivAt φ (Real.sqrt (vaidyaW q (m (v s)) (r (v s))) /
      r (v s) * Real.sin (α (v s))) (v s))
    (hv : HasDerivAt v
      (reconstructedTimeSpeed L (vaidyaW q (m (v s)) (r (v s))) (r (v s)) (α (v s))) s) :
    HamiltonianLiftAt q L m μ v (fun u => r (v u)) (fun u => φ (v u))
      (fun u => reconstructedPv q L (m (v u)) (r (v u)) (α (v u)))
      (fun u => reconstructedPr q L (r (v u)) (α (v u))) s := by
  have hnull := reconstructed_future_null hq hL hr hα hw
  have hs := Real.sin_pos_of_pos_of_lt_pi hα.1 hα.2
  have htime := reconstructed_time_speed_identity (q := q) (L := L)
    (ne_of_gt hr) (ne_of_gt hs) hw
  have hV : reconstructedTimeSpeed L (vaidyaW q (m (v s)) (r (v s)))
      (r (v s)) (α (v s)) ≠ 0 := ne_of_gt (htime ▸ hnull.2)
  have hids := reconstructed_normalized_hamiltonian_identities hq hL hr hα hw
  refine ⟨?_, ?_, ?_, ?_, ?_, hasDerivAt_const s L, hnull.1, hnull.2⟩
  · simpa only [htime] using hv
  · exact (hrder.comp s hv).congr_deriv ((div_eq_iff hV).mp hids.1).symm
  · exact (hφ.comp s hv).congr_deriv ((div_eq_iff hV).mp hids.2.1).symm
  · have hpv := reconstructed_time_momentum_dynamics hq hL hr hα hw hm hrder hαder
    simpa only [div_mul_cancel₀ _ hV] using hpv.comp s hv
  · have hp := reconstructed_radial_momentum_dynamics hq hL hr hα hw hrder hαder
    exact (hp.comp s hv).congr_deriv ((div_eq_iff hV).mp hids.2.2).symm

/-- Construct the angle primitive, affine clock and its one chosen inverse,
then prove the complete Hamiltonian equations throughout an affine
neighborhood. No Hamiltonian trajectory or reparametrization is a premise. -/
theorem direction_affine_hamiltonian_germ_exists {q L a : ℝ} {m μ r α : ℝ → ℝ}
    (hq : 0 < q) (hL : 0 < L) (hr : 0 < r a)
    (hα : 0 < α a ∧ α a < Real.pi) (hw : 0 < vaidyaW q (m a) (r a))
    (hm : ∀ᶠ t in 𝓝 a, HasDerivAt m (μ t) t)
    (hrder : ∀ᶠ t in 𝓝 a, HasDerivAt r
      (radialSpeed q (vaidyaW q (m t) (r t)) (Real.cos (α t))) t)
    (hαder : ∀ᶠ t in 𝓝 a, HasDerivAt α
      (Real.sin (α t) / r t * directionB q (vaidyaW q (m t) (r t)) (Real.cos (α t))) t) :
    ∃ ℓ v φ : ℝ → ℝ, ℓ a = 0 ∧ v 0 = a ∧ φ a = 0 ∧
      (∀ᶠ t in 𝓝 a, v (ℓ t) = t) ∧
      (∀ᶠ s in 𝓝 (0 : ℝ), ℓ (v s) = s) ∧
      (∀ᶠ s in 𝓝 (0 : ℝ), HamiltonianLiftAt q L m μ v
        (fun u => r (v u)) (fun u => φ (v u))
        (fun u => reconstructedPv q L (m (v u)) (r (v u)) (α (v u)))
        (fun u => reconstructedPr q L (r (v u)) (α (v u))) s) := by
  have hmc := hm.mono (fun _ ht => ht.continuousAt)
  have hrc := hrder.mono (fun _ ht => ht.continuousAt)
  have hαc := hαder.mono (fun _ ht => ht.continuousAt)
  obtain ⟨ℓ, v, hℓ, hva, _, hleft, hright, hv, hvnear⟩ :=
    direction_affine_clock_exists hL hr hα hw hmc hrc hαc
  obtain ⟨φ, hφa, hφ⟩ := direction_angular_primitive_exists hr hα hw hmc hrc hαc
  have hcw : ContinuousAt (fun t => vaidyaW q (m t) (r t)) a :=
    continuousAt_const.add ((hmc.self_of_nhds.const_mul 2).div
      hrc.self_of_nhds (ne_of_gt hr))
  have hnear : ∀ᶠ t in 𝓝 a,
      0 < r t ∧ (0 < α t ∧ α t < Real.pi) ∧ 0 < vaidyaW q (m t) (r t) ∧
      HasDerivAt m (μ t) t ∧
      HasDerivAt r (radialSpeed q (vaidyaW q (m t) (r t)) (Real.cos (α t))) t ∧
      HasDerivAt α (Real.sin (α t) / r t *
        directionB q (vaidyaW q (m t) (r t)) (Real.cos (α t))) t ∧
      HasDerivAt φ (Real.sqrt (vaidyaW q (m t) (r t)) / r t * Real.sin (α t)) t := by
    filter_upwards [hm, hrder, hαder, hφ,
      continuousAt_const.eventually_lt hrc.self_of_nhds hr,
      continuousAt_const.eventually_lt hαc.self_of_nhds hα.1,
      hαc.self_of_nhds.eventually_lt continuousAt_const hα.2,
      continuousAt_const.eventually_lt hcw hw] with t hmt hrt hαt hφt hrp hαp hαπ hwp
    exact ⟨hrp, ⟨hαp, hαπ⟩, hwp, hmt, hrt, hαt, hφt⟩
  have hvt : Tendsto v (𝓝 (0 : ℝ)) (𝓝 a) := by
    simpa only [hva] using hv.continuousAt.tendsto
  refine ⟨ℓ, v, φ, hℓ, hva, hφa, hleft, hright, ?_⟩
  filter_upwards [hvt.eventually hnear, hvnear] with s hs hvs
  exact direction_to_affine_hamiltonian_at hq hL hs.1 hs.2.1 hs.2.2.1
    hs.2.2.2.1 hs.2.2.2.2.1 hs.2.2.2.2.2.1 hs.2.2.2.2.2.2 hvs

/-- Apply the lift to the maximal exterior curve already constructed from
local solutions. Radius positivity, angle range and the direction ODE are
conclusions of that construction, not new premises here. -/
theorem maximal_exterior_affine_hamiltonian_germ_exists {q L a l c : ℝ}
    {m : ℝ → ℝ} {x : ℝ × ℝ} (hq : 1 < q) (hL : 0 < L) (hla : l < a)
    (hmC1 : ∀ t ∈ maximalExteriorDomain q m a l x, ContDiffAt ℝ 1 m t)
    (hc : c ∈ maximalExteriorDomain q m a l x) :
    let f := maximalExteriorCurve q m a l x
    ∃ ℓ v φ : ℝ → ℝ, ℓ c = 0 ∧ v 0 = c ∧ φ c = 0 ∧
      (∀ᶠ t in 𝓝 c, v (ℓ t) = t) ∧
      (∀ᶠ s in 𝓝 (0 : ℝ), ℓ (v s) = s) ∧
      (∀ᶠ s in 𝓝 (0 : ℝ), HamiltonianLiftAt q L m (deriv m) v
        (fun u => (f (v u)).1) (fun u => φ (v u))
        (fun u => reconstructedPv q L (m (v u)) (f (v u)).1 (f (v u)).2)
        (fun u => reconstructedPr q L (f (v u)).1 (f (v u)).2) s) := by
  let f := maximalExteriorCurve q m a l x
  have hmc : ∀ t ∈ maximalExteriorDomain q m a l x, ContinuousAt m t :=
    fun t ht => (hmC1 t ht).continuousAt
  have hreg := maximalExteriorCurve_solves hq hla hmc c hc
  have hr : 0 < (f c).1 := by linarith [hreg.2.1, hreg.2.2.1]
  have hα : 0 < (f c).2 ∧ (f c).2 < Real.pi := hreg.2.2.2
  have hw : 0 < vaidyaW q (m c) (f c).1 := vaidyaW_pos hq hreg.2.1.le hr
  have hnear : ∀ᶠ t in 𝓝 c, t ∈ maximalExteriorDomain q m a l x :=
    (maximalExteriorDomain_isOpen q m a l x).mem_nhds hc
  have hm : ∀ᶠ t in 𝓝 c, HasDerivAt m (deriv m t) t := by
    filter_upwards [hnear] with t ht
    exact ((hmC1 t ht).differentiableAt le_rfl).hasDerivAt
  have hrder : ∀ᶠ t in 𝓝 c, HasDerivAt (fun u => (f u).1)
      (radialSpeed q (vaidyaW q (m t) (f t).1) (Real.cos (f t).2)) t := by
    filter_upwards [hnear] with t ht
    exact direction_radial_derivative (maximalExteriorCurve_solves hq hla hmc t ht).1
  have hαder : ∀ᶠ t in 𝓝 c, HasDerivAt (fun u => (f u).2)
      (Real.sin (f t).2 / (f t).1 *
        directionB q (vaidyaW q (m t) (f t).1) (Real.cos (f t).2)) t := by
    filter_upwards [hnear] with t ht
    exact direction_angular_derivative (maximalExteriorCurve_solves hq hla hmc t ht).1
  exact direction_affine_hamiltonian_germ_exists (by linarith) hL hr hα hw hm hrder hαder

end ViaB
