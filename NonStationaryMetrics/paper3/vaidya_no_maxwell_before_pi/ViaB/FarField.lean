import ViaB.Algebra
import Mathlib.Analysis.SpecialFunctions.Trigonometric.Basic
import Mathlib.Tactic.FunProp

open Filter Topology

namespace ViaB

noncomputable def radialSpeed (q w c : ℝ) : ℝ := -w + Real.sqrt (q * w) * c

theorem directionB_continuousAt {q w c : ℝ} (hw : w ≠ 0) :
    ContinuousAt (fun p : ℝ × ℝ => directionB q p.1 p.2) (w, c) := by
  unfold directionB
  fun_prop (disch := simpa using hw)

theorem directionB_le_at_one {q w c : ℝ} (hq : 1 < q)
    (hw : q - 1 ≤ w) (hc : c ≤ 1) : directionB q w c ≤ directionB q w 1 := by
  have hcoef : 0 ≤ 2 * w - q + 1 := by linarith
  have hp := mul_nonneg (sub_nonneg.mpr hc) hcoef
  unfold directionB
  nlinarith

/-- A uniform negative margin, valid for all c ≤ 1 at once. -/
theorem directionB_uniform_far {q : ℝ} (hq : 1 < q) :
    ∃ ε > 0, ∃ κ > 0, ∀ w c : ℝ,
      q - 1 ≤ w → w < q - 1 + ε → c ≤ 1 → directionB q w c ≤ -κ := by
  have hneg := directionB_at_infinity_negative hq (le_refl (1 : ℝ))
  let κ := -directionB q (q - 1) 1 / 2
  have hk : 0 < κ := by dsimp [κ]; linarith
  have hcont : ContinuousAt (fun w => directionB q w 1) (q - 1) := by
    unfold directionB
    fun_prop (disch := linarith)
  have hevent : ∀ᶠ w in 𝓝 (q - 1), directionB q w 1 < -κ :=
    hcont.eventually (eventually_lt_nhds (show directionB q (q - 1) 1 < -κ by
      dsimp [κ]; linarith))
  obtain ⟨ε, heps, hε⟩ := Metric.eventually_nhds_iff.mp hevent
  refine ⟨ε, heps, κ, hk, ?_⟩
  intro w c hwlow hwup hc
  have hd : dist w (q - 1) < ε := by
    rw [Real.dist_eq, abs_of_nonneg (sub_nonneg.mpr hwlow)]
    linarith
  exact (directionB_le_at_one hq hwlow hc).trans (hε hd).le

/-- Bounded mass gives a single radius valid for every time and direction.
No hypothesis on the mass derivative occurs in this statement. -/
theorem directionB_large_radius {q M : ℝ} (hq : 1 < q) (hM : 0 < M) :
    ∃ R > 2 * M, ∃ κ > 0, ∀ r m c : ℝ,
      R ≤ r → 0 ≤ m → m ≤ M → c ≤ 1 →
      directionB q (q - 1 + 2 * m / r) c ≤ -κ := by
  obtain ⟨ε, heps, κ, hk, hfar⟩ := directionB_uniform_far hq
  obtain ⟨R, hR⟩ := exists_gt (max (2 * M) (2 * M / ε))
  have hRM : 2 * M < R := (le_max_left _ _).trans_lt hR
  have hReps : 2 * M / ε < R := (le_max_right _ _).trans_lt hR
  have hRpos : 0 < R := by linarith
  have hmass : 2 * M < R * ε := (div_lt_iff₀ heps).mp hReps
  refine ⟨R, hRM, κ, hk, ?_⟩
  intro r m c hr hm hmM hc
  have hrpos : 0 < r := hRpos.trans_le hr
  have hdiv : 0 ≤ 2 * m / r := by positivity
  have hdivlt : 2 * m / r < ε := (div_lt_iff₀ hrpos).2 (by nlinarith)
  exact hfar _ _ (by linarith) (by linarith) hc

theorem radialSpeed_at_infinity_pos {q c : ℝ} (hq : 1 < q)
    (hc : 0 < c) (hcone : q - 1 < q * c ^ 2) :
    0 < radialSpeed q (q - 1) c := by
  have hw : 0 < q - 1 := by linarith
  have hqpos : 0 < q := by linarith
  have hprod : 0 ≤ q * (q - 1) := (mul_pos hqpos hw).le
  have hs := Real.sq_sqrt hprod
  have hsnonneg := Real.sqrt_nonneg (q * (q - 1))
  have hh := mul_pos hw (sub_pos.mpr hcone)
  have hsqc : (Real.sqrt (q * (q - 1)) * c) ^ 2 = q * (q - 1) * c ^ 2 := by
    rw [mul_pow, hs]
  have hlt : q - 1 < Real.sqrt (q * (q - 1)) * c :=
    (sq_lt_sq₀ hw.le (mul_nonneg hsnonneg hc.le)).mp (by nlinarith [hh, hsqc])
  unfold radialSpeed
  linarith

theorem radialSpeed_le_charge {q w c : ℝ} (hq : 0 < q)
    (hw : 0 ≤ w) (hwq : w ≤ q) (hc : c ≤ 1) : radialSpeed q w c ≤ q := by
  have hsnonneg := Real.sqrt_nonneg (q * w)
  have hssq := Real.sq_sqrt (mul_nonneg hq.le hw)
  have hsq : Real.sqrt (q * w) ≤ q := by nlinarith
  have hmul := mul_nonneg hsnonneg (sub_nonneg.mpr hc)
  unfold radialSpeed
  nlinarith

/-- The outgoing cone admits a strictly positive radial speed margin. -/
theorem radialSpeed_uniform_far {q c₀ : ℝ} (hq : 1 < q)
    (hc₀ : 0 < c₀) (hcone : q - 1 < q * c₀ ^ 2) :
    ∃ ε > 0, ∃ δ > 0, ∀ w c : ℝ,
      q - 1 ≤ w → w < q - 1 + ε → c₀ ≤ c → δ ≤ radialSpeed q w c := by
  have hpos := radialSpeed_at_infinity_pos hq hc₀ hcone
  let δ := radialSpeed q (q - 1) c₀ / 2
  have hδ : 0 < δ := by dsimp [δ]; linarith
  have hcont : ContinuousAt (fun w => radialSpeed q w c₀) (q - 1) := by
    unfold radialSpeed
    fun_prop
  have hevent : ∀ᶠ w in 𝓝 (q - 1), δ < radialSpeed q w c₀ :=
    hcont.eventually (eventually_gt_nhds (show δ < radialSpeed q (q - 1) c₀ by
      dsimp [δ]; linarith))
  obtain ⟨ε, heps, hε⟩ := Metric.eventually_nhds_iff.mp hevent
  refine ⟨ε, heps, δ, hδ, ?_⟩
  intro w c hwlow hwup hc
  have hd : dist w (q - 1) < ε := by
    rw [Real.dist_eq, abs_of_nonneg (sub_nonneg.mpr hwlow)]
    linarith
  have hs := Real.sqrt_nonneg (q * w)
  have hp := mul_nonneg hs (sub_nonneg.mpr hc)
  have hm : radialSpeed q w c₀ ≤ radialSpeed q w c := by unfold radialSpeed; nlinarith
  exact (hε hd).le.trans hm

end ViaB
