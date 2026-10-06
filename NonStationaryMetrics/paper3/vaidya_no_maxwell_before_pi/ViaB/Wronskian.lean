import Mathlib.Analysis.Calculus.Deriv.Inv
import Mathlib.Analysis.Calculus.Deriv.MeanValue
import Mathlib.Tactic.Ring

namespace ViaB

def wronskian (u du z dz : ℝ → ℝ) (x : ℝ) : ℝ := dz x * u x - z x * du x

/-- Exact differential identity for two scalar Jacobi equations. The ODE
assumptions are explicit; this does not derive them from the Vaidya metric. -/
theorem wronskian_hasDerivAt {u du z dz : ℝ → ℝ} {kperp kpar x : ℝ}
    (hu : HasDerivAt u (du x) x) (hz : HasDerivAt z (dz x) x)
    (hdu : HasDerivAt du (-kperp * u x) x)
    (hdz : HasDerivAt dz (-kpar * z x) x) :
    HasDerivAt (wronskian u du z dz) ((kperp - kpar) * u x * z x) x := by
  convert (hdz.fun_mul hu).sub (hz.fun_mul hdu) using 1
  ring

theorem ratio_hasDerivAt {u du z dz : ℝ → ℝ} {x : ℝ}
    (hu : HasDerivAt u (du x) x) (hz : HasDerivAt z (dz x) x)
    (hne : u x ≠ 0) :
    HasDerivAt (fun t => z t / u t) (wronskian u du z dz x / u x ^ 2) x := by
  simpa only [wronskian] using hz.div hu hne

/-- On any positive comparison interval the Wronskian is strictly increasing.
Positivity of z is a hypothesis here, not a conclusion of the full Sturm lemma. -/
theorem wronskian_strictMonoOn {u du z dz kperp kpar : ℝ → ℝ} {a b : ℝ}
    (hu : ∀ x ∈ Set.Icc a b, HasDerivAt u (du x) x)
    (hz : ∀ x ∈ Set.Icc a b, HasDerivAt z (dz x) x)
    (hdu : ∀ x ∈ Set.Icc a b, HasDerivAt du (-kperp x * u x) x)
    (hdz : ∀ x ∈ Set.Icc a b, HasDerivAt dz (-kpar x * z x) x)
    (hgap : ∀ x ∈ Set.Ioo a b, 0 < kperp x - kpar x)
    (hupos : ∀ x ∈ Set.Ioo a b, 0 < u x)
    (hzpos : ∀ x ∈ Set.Ioo a b, 0 < z x) :
    StrictMonoOn (wronskian u du z dz) (Set.Icc a b) := by
  have hw : ∀ x ∈ Set.Icc a b,
      HasDerivAt (wronskian u du z dz) ((kperp x - kpar x) * u x * z x) x :=
    fun x hx => wronskian_hasDerivAt (hu x hx) (hz x hx) (hdu x hx) (hdz x hx)
  apply strictMonoOn_of_deriv_pos (convex_Icc a b)
  · exact fun x hx => (hw x hx).continuousAt.continuousWithinAt
  · intro x hx
    have hi : x ∈ Set.Ioo a b := by simpa only [interior_Icc] using hx
    rw [(hw x (Set.Ioo_subset_Icc_self hi)).deriv]
    exact mul_pos (mul_pos (hgap x hi) (hupos x hi)) (hzpos x hi)

/-- Positive Wronskian and a nonzero reference mode force a positive ratio derivative. -/
theorem ratio_deriv_pos {u du z dz : ℝ → ℝ} {x : ℝ}
    (hu : HasDerivAt u (du x) x) (hz : HasDerivAt z (dz x) x)
    (hupos : 0 < u x) (hwpos : 0 < wronskian u du z dz x) :
    0 < deriv (fun t => z t / u t) x := by
  rw [(ratio_hasDerivAt hu hz (ne_of_gt hupos)).deriv]
  exact div_pos hwpos (sq_pos_of_pos hupos)

end ViaB
