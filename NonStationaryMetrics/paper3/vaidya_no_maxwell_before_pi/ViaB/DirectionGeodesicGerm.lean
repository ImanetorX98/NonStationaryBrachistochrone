import ViaB.HamiltonianGeodesic

open Filter
open scoped Topology

namespace ViaB

/-- A future-null affine geodesic equation in the fixed Fermat chart.
The connection coefficients are the independently verified Koszul ones. -/
structure CoordinateNullGeodesicAt (q : ℝ) (m μ v R Φ : ℝ → ℝ) (s : ℝ) : Prop where
  equation : ∀ i : CoordinateIndex,
    HasDerivAt (fun t => ![deriv v t, deriv R t, deriv Φ t] i)
      (-fermatChristoffelContraction q (m (v s)) (μ (v s)) (R s)
        (deriv v s) (deriv R s) (deriv Φ s) i) s
  null : metricQuad q (vaidyaW q (m (v s)) (R s)) (R s)
    (deriv v s) (deriv R s) (deriv Φ s) = 0
  future : 0 < deriv v s

/-- Actual second derivatives and the null geodesic equation hold throughout
the same Hamiltonian neighborhood. No C2 coordinate curve is assumed. -/
theorem hamiltonian_germ_coordinate_null_geodesic {q L s : ℝ}
    {m μ v R Φ pv p : ℝ → ℝ} (hq : q ≠ 0) (hr : R s ≠ 0)
    (hw : vaidyaW q (m (v s)) (R s) ≠ 0)
    (hm : ∀ᶠ t in 𝓝 (v s), HasDerivAt m (μ t) t)
    (h : ∀ᶠ t in 𝓝 s, HamiltonianLiftAt q L m μ v R Φ pv p t) :
    ∀ᶠ t in 𝓝 s, CoordinateNullGeodesicAt q m μ v R Φ t := by
  have hs := h.self_of_nhds
  have hmc : ContinuousAt (fun t => m (v t)) s :=
    hm.self_of_nhds.continuousAt.comp hs.time.continuousAt
  have hcw : ContinuousAt (fun t => vaidyaW q (m (v t)) (R t)) s :=
    continuousAt_const.add ((hmc.const_mul 2).div hs.radial.continuousAt hr)
  filter_upwards [eventually_eventually_nhds.mpr h,
    hs.time.continuousAt.eventually hm,
    hs.radial.continuousAt.eventually (eventually_ne_nhds hr),
    hcw.eventually (eventually_ne_nhds hw)] with t ht hmt hrt hwt
  have hpoint := ht.self_of_nhds
  exact ⟨hamiltonian_germ_coordinate_geodesic hq hrt hwt hmt ht,
    hamiltonian_lift_null_coordinate_velocity hq hrt hwt hpoint,
    by simpa only [hpoint.time.deriv] using hpoint.future⟩

/-- Constructed direction trajectories are genuine affine null geodesics
in the coordinate Levi-Civita sense, throughout one local affine interval.
Clock, inverse, momenta and second derivatives are all conclusions. -/
theorem direction_coordinate_null_geodesic_germ_exists {q L a : ℝ}
    {m μ r α : ℝ → ℝ} (hq : 0 < q) (hL : 0 < L) (hr : 0 < r a)
    (hα : 0 < α a ∧ α a < Real.pi) (hw : 0 < vaidyaW q (m a) (r a))
    (hm : ∀ᶠ t in 𝓝 a, HasDerivAt m (μ t) t)
    (hrder : ∀ᶠ t in 𝓝 a, HasDerivAt r
      (radialSpeed q (vaidyaW q (m t) (r t)) (Real.cos (α t))) t)
    (hαder : ∀ᶠ t in 𝓝 a, HasDerivAt α
      (Real.sin (α t) / r t * directionB q (vaidyaW q (m t) (r t)) (Real.cos (α t))) t) :
    ∃ ℓ v φ : ℝ → ℝ, ℓ a = 0 ∧ v 0 = a ∧ φ a = 0 ∧
      (∀ᶠ t in 𝓝 a, v (ℓ t) = t) ∧
      (∀ᶠ s in 𝓝 (0 : ℝ), ℓ (v s) = s) ∧
      (∀ᶠ s in 𝓝 (0 : ℝ), CoordinateNullGeodesicAt q m μ v
        (fun u => r (v u)) (fun u => φ (v u)) s) := by
  obtain ⟨ℓ, v, φ, hℓ, hva, hφ, hleft, hright, hHam⟩ :=
    direction_affine_hamiltonian_germ_exists hq hL hr hα hw hm hrder hαder
  refine ⟨ℓ, v, φ, hℓ, hva, hφ, hleft, hright, ?_⟩
  apply hamiltonian_germ_coordinate_null_geodesic (ne_of_gt hq) (by simpa only [hva] using ne_of_gt hr)
    (by simpa only [hva] using ne_of_gt hw) (by simpa only [hva] using hm) hHam

/-- The previously constructed maximal exterior curve has affine null
geodesic germs at every event in its domain, for prescribed C1 mass. -/
theorem maximal_exterior_coordinate_null_geodesic_germ_exists {q L a l c : ℝ}
    {m : ℝ → ℝ} {x : ℝ × ℝ} (hq : 1 < q) (hL : 0 < L) (hla : l < a)
    (hmC1 : ∀ t ∈ maximalExteriorDomain q m a l x, ContDiffAt ℝ 1 m t)
    (hc : c ∈ maximalExteriorDomain q m a l x) :
    let f := maximalExteriorCurve q m a l x
    ∃ ℓ v φ : ℝ → ℝ, ℓ c = 0 ∧ v 0 = c ∧ φ c = 0 ∧
      (∀ᶠ t in 𝓝 c, v (ℓ t) = t) ∧
      (∀ᶠ s in 𝓝 (0 : ℝ), ℓ (v s) = s) ∧
      (∀ᶠ s in 𝓝 (0 : ℝ), CoordinateNullGeodesicAt q m (deriv m) v
        (fun u => (f (v u)).1) (fun u => φ (v u)) s) := by
  let f := maximalExteriorCurve q m a l x
  obtain ⟨ℓ, v, φ, hℓ, hva, hφ, hleft, hright, hHam⟩ :=
    maximal_exterior_affine_hamiltonian_germ_exists hq hL hla hmC1 hc
  have hmc : ∀ t ∈ maximalExteriorDomain q m a l x, ContinuousAt m t :=
    fun t ht => (hmC1 t ht).continuousAt
  have hreg := maximalExteriorCurve_solves hq hla hmc c hc
  have hr : 0 < (f c).1 := by linarith [hreg.2.1, hreg.2.2.1]
  have hw : 0 < vaidyaW q (m c) (f c).1 := vaidyaW_pos hq hreg.2.1.le hr
  have hm : ∀ᶠ t in 𝓝 c, HasDerivAt m (deriv m t) t := by
    filter_upwards [(maximalExteriorDomain_isOpen q m a l x).mem_nhds hc] with t ht
    exact ((hmC1 t ht).differentiableAt le_rfl).hasDerivAt
  refine ⟨ℓ, v, φ, hℓ, hva, hφ, hleft, hright, ?_⟩
  apply hamiltonian_germ_coordinate_null_geodesic (by linarith : q ≠ 0)
    (by simpa only [hva] using ne_of_gt hr) (by simpa only [hva] using ne_of_gt hw)
    (by simpa only [hva] using hm) hHam

end ViaB
