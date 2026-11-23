import Mathlib.Probability.BrownianMotion.Basic
import Mathlib.Probability.Martingale.Basic
import Mathlib.MeasureTheory.Integral.IntervalIntegral.IntegrationByParts
import Mathlib.Analysis.SpecialFunctions.Integrals.Basic
import Mathlib.Tactic

/-!
# Brownian compensator: verified components of the recorded problem

This module contains genuine analytic, measure-theoretic and algebraic proofs.
The full stochastic-integral construction, Ito formula and terminal compensator
martingale are open obligations; no theorem here assumes these conclusions and
then advertises them as proved. See VERIFICATION.md for the coverage ledger.
-/

open MeasureTheory ProbabilityTheory Set
open scoped BigOperators NNReal ENNReal Topology Interval

namespace BrownianCompensator

noncomputable section

def digital (c x : ℝ) : ℝ := if x < c then 1 else 0

theorem digital_bounds (c x : ℝ) : 0 ≤ digital c x ∧ digital c x ≤ 1 := by
  unfold digital
  split_ifs <;> constructor <;> norm_num

theorem digital_measurable (c : ℝ) : Measurable (digital c) := by
  unfold digital
  exact measurable_const.ite (measurableSet_lt measurable_id measurable_const) measurable_const

def atom (a b x : ℝ) : ℝ :=
  -(1 / 4 : ℝ) * max (a + x) 0 + (1 / 2 : ℝ) * max (b - a) 0 -
    (1 / 4 : ℝ) * digital ((1 - a - b) / 2) x

def atomEnvelope (a b : ℝ) : ℝ :=
  (1 / 4 : ℝ) * |a| + (1 / 2 : ℝ) * max (b - a) 0 + 1 / 4

theorem atom_measurable (a b : ℝ) : Measurable (atom a b) := by
  unfold atom
  exact (((measurable_const.add measurable_id).max measurable_const).const_mul _).add
    measurable_const |>.sub ((digital_measurable _).const_mul _)

theorem atomEnvelope_nonneg (a b : ℝ) : 0 ≤ atomEnvelope a b := by
  unfold atomEnvelope
  positivity

theorem atom_growth (a b x : ℝ) : |atom a b x| ≤ atomEnvelope a b + |x| / 4 := by
  have hp : 0 ≤ max (a + x) 0 := le_max_right _ _
  have hq : 0 ≤ max (b - a) 0 := le_max_right _ _
  have hmax : max (a + x) 0 ≤ |a| + |x| := by
    exact max_le (add_le_add (le_abs_self _) (le_abs_self _)) (by positivity)
  obtain ⟨hd0, hd1⟩ := digital_bounds ((1 - a - b) / 2) x
  unfold atom atomEnvelope
  apply abs_le.mpr
  constructor <;> nlinarith [le_abs_self a, le_abs_self x]

def payoff (n : ℕ) (a : Fin n → ℝ) (x : ℝ) : ℝ :=
  ((n : ℝ) ^ 2)⁻¹ * ∑ i : Fin n, ∑ j : Fin n, atom (a i) (a j) x

def envelope (n : ℕ) (a : Fin n → ℝ) : ℝ :=
  ((n : ℝ) ^ 2)⁻¹ * ∑ i : Fin n, ∑ j : Fin n, atomEnvelope (a i) (a j)

theorem payoff_measurable (n : ℕ) (a : Fin n → ℝ) : Measurable (payoff n a) := by
  unfold payoff
  exact Measurable.const_mul
    (Finset.measurable_sum _ fun i _ => Finset.measurable_sum _ fun j _ =>
      atom_measurable _ _) _

theorem envelope_nonneg (n : ℕ) (a : Fin n → ℝ) : 0 ≤ envelope n a := by
  unfold envelope
  apply mul_nonneg (by positivity)
  exact Finset.sum_nonneg fun i _ => Finset.sum_nonneg fun j _ => atomEnvelope_nonneg _ _

theorem payoff_linear_growth (n : ℕ) (hn : 0 < n) (a : Fin n → ℝ) (x : ℝ) :
    |payoff n a x| ≤ envelope n a + |x| / 4 := by
  have hnR : (n : ℝ) ≠ 0 := Nat.cast_ne_zero.mpr (Nat.ne_of_gt hn)
  calc
    |payoff n a x|
        ≤ ((n : ℝ) ^ 2)⁻¹ * ∑ i : Fin n, ∑ j : Fin n, |atom (a i) (a j) x| := by
          unfold payoff
          rw [abs_mul, abs_of_nonneg (inv_nonneg.mpr (sq_nonneg _))]
          apply mul_le_mul_of_nonneg_left _ (by positivity)
          exact (Finset.abs_sum_le_sum_abs _ _).trans
            (Finset.sum_le_sum fun i _ => Finset.abs_sum_le_sum_abs _ _)
    _ ≤ ((n : ℝ) ^ 2)⁻¹ * ∑ i : Fin n, ∑ j : Fin n,
          (atomEnvelope (a i) (a j) + |x| / 4) := by
          apply mul_le_mul_of_nonneg_left _ (by positivity)
          exact Finset.sum_le_sum fun i _ => Finset.sum_le_sum fun j _ => atom_growth _ _ _
    _ = envelope n a + |x| / 4 := by
          simp only [Finset.sum_add_distrib, Finset.sum_const, Finset.card_univ,
            Fintype.card_fin, nsmul_eq_mul, envelope]
          field_simp

theorem atom_original_indicator (a b x : ℝ) :
    digital ((1 - a - b) / 2) x = (if a + b + 2 * x < 1 then 1 else 0) := by
  have h : x < (1 - a - b) / 2 ↔ a + b + 2 * x < 1 := by constructor <;> intro h <;> linarith
  simp only [digital, h]

def primitive (S : ℝ → ℝ) (x : ℝ) : ℝ := ∫ r in (0 : ℝ)..x, r * S r

theorem primitive_zero (S : ℝ → ℝ) : primitive S 0 = 0 := by simp [primitive]

theorem primitive_hasDerivAt (S : ℝ → ℝ) (hS : Continuous S) (x : ℝ) :
    HasDerivAt (primitive S) (x * S x) x := by
  have hf : Continuous (fun r : ℝ => r * S r) := continuous_id.mul hS
  exact intervalIntegral.integral_hasDerivAt_right (hf.intervalIntegrable 0 x)
    hf.stronglyMeasurable.stronglyMeasurableAtFilter hf.continuousAt

theorem primitive_deriv (S : ℝ → ℝ) (hS : Continuous S) (x : ℝ) :
    deriv (primitive S) x = x * S x := (primitive_hasDerivAt S hS x).deriv

theorem primitive_second_deriv (S : ℝ → ℝ) (hS : Continuous S)
    (x sx : ℝ) (h : HasDerivAt S sx x) :
    deriv (fun y => deriv (primitive S) y) x = S x + x * sx := by
  simp_rw [primitive_deriv S hS]
  have he : (fun y : ℝ => y * S y) = id * S := by funext y; rfl
  rw [he]
  simpa only [id_eq, one_mul] using ((hasDerivAt_id x).mul h).deriv

theorem weighted_integration_by_parts (S Sx Sxx : ℝ → ℝ)
    (hS : ∀ x, HasDerivAt S (Sx x) x)
    (hSx : ∀ x, HasDerivAt Sx (Sxx x) x)
    (hSxx : Continuous Sxx) (x : ℝ) :
    (∫ r in (0 : ℝ)..x, r * Sxx r) = x * Sx x - S x + S 0 := by
  have cS : Continuous S := continuous_iff_continuousAt.mpr fun y => (hS y).continuousAt
  have cSx : Continuous Sx := continuous_iff_continuousAt.mpr fun y => (hSx y).continuousAt
  have h1 := intervalIntegral.integral_mul_deriv_eq_deriv_mul_of_hasDerivAt
    continuous_id.continuousOn cSx.continuousOn
    (fun y _ => hasDerivAt_id y) (fun y _ => hSx y)
    (continuous_const.intervalIntegrable 0 x) (hSxx.intervalIntegrable 0 x)
  have h2 : (∫ r in (0 : ℝ)..x, Sx r) = S x - S 0 :=
    intervalIntegral.integral_eq_sub_of_hasDerivAt (fun y _ => hS y)
      (cSx.intervalIntegrable 0 x)
  simp only [id_eq, one_mul, zero_mul, sub_zero] at h1
  rw [h2] at h1
  linarith

theorem generator_cancellation (s sx s0 x : ℝ) :
    -(1 / 2 : ℝ) * (x * sx - s + s0) + (1 / 2 : ℝ) * (s + x * sx)
      = s - (1 / 2 : ℝ) * s0 := by ring

theorem integrand_square_bound (C x y : ℝ) (hC : 0 ≤ C)
    (hy : |y| ≤ C * (1 + |x|)) :
    (x * y)^2 ≤ 2 * C^2 * (x^2 + x^4) := by
  have hx : 0 ≤ |x| := abs_nonneg x
  have hbase : 0 ≤ C * (1 + |x|) := mul_nonneg hC (by positivity)
  have hs : y^2 ≤ (C * (1 + |x|))^2 := by
    nlinarith [mul_nonneg (sub_nonneg.mpr hy)
      (add_nonneg hbase (abs_nonneg y)), sq_abs y]
  have hpoly : (1 + |x|)^2 ≤ 2 * (1 + x^2) := by
    nlinarith [sq_nonneg (|x| - 1), sq_abs x]
  have hxx : 0 ≤ x^2 := sq_nonneg x
  nlinarith [mul_nonneg hxx (sub_nonneg.mpr hs),
    mul_nonneg (mul_nonneg (sq_nonneg C) hxx) (sub_nonneg.mpr hpoly)]

theorem integrated_moment_majorant (C T : ℝ) :
    (∫ s in (0 : ℝ)..T, 2 * C^2 * (s + 3*s^2)) = C^2 * (T^2 + 2*T^3) := by
  have hf (s : ℝ) : HasDerivAt (fun u : ℝ => C^2*(u^2+2*u^3))
      (2*C^2*(s+3*s^2)) s := by
    convert (((hasDerivAt_pow 2 s).add ((hasDerivAt_pow 3 s).const_mul 2)).const_mul (C^2)) using 1; ring
  have hc : Continuous (fun s : ℝ => 2*C^2*(s+3*s^2)) := by fun_prop
  have h := intervalIntegral.integral_eq_sub_of_hasDerivAt (fun s _ => hf s)
    (hc.intervalIntegrable 0 T)
  simpa using h

theorem compensated_path_identity (P K I : ℝ) (h : P = K + I) : P - K = I := by linarith

variable {Ω : Type*} [MeasurableSpace Ω]

theorem brownian_finite_moments (B : ℝ≥0 → Ω → ℝ) (μ : Measure Ω)
    (hB : IsBrownianReal B μ) (t : ℝ≥0) (p : ℝ≥0∞) (hp : p ≠ ∞) :
    MemLp (B t) p μ :=
  (hB.toIsPreBrownianReal.hasLaw_eval t).memLp (memLp_id_gaussianReal' p hp)

theorem terminal_payoff_memLp (n : ℕ) (hn : 0 < n) (a : Fin n → ℝ)
    (B : ℝ≥0 → Ω → ℝ) (μ : Measure Ω) [IsProbabilityMeasure μ]
    (hB : IsBrownianReal B μ) (T : ℝ≥0) (p : ℝ≥0∞) (hp : p ≠ ∞) :
    MemLp (fun ω => payoff n a (B T ω)) p μ := by
  have hBmem := brownian_finite_moments B μ hB T p hp
  have hmajor : MemLp (fun ω => envelope n a + |B T ω| / 4) p μ := by
    have hh := (memLp_const (envelope n a)).add (hBmem.norm.mul_const (4 : ℝ)⁻¹)
    convert hh using 1
    funext ω
    simp only [Pi.add_apply, Real.norm_eq_abs, div_eq_mul_inv]
  apply hmajor.mono'
    ((payoff_measurable n a).comp_aemeasurable (hB.toIsPreBrownianReal.aemeasurable T)).aestronglyMeasurable
  exact ae_of_all _ fun ω => by
    simpa only [Real.norm_eq_abs, Function.comp_apply] using payoff_linear_growth n hn a (B T ω)

/-- The actual terminal payoff has a genuine conditional-expectation martingale.
Identifying this with the Gaussian convolution is a separate open obligation. -/
theorem terminal_payoff_condExp_martingale (n : ℕ) (a : Fin n → ℝ)
    (B : ℝ≥0 → Ω → ℝ) (μ : Measure Ω)
    (ℱ : Filtration ℝ≥0 ‹MeasurableSpace Ω›) [SigmaFiniteFiltration μ ℱ] (T : ℝ≥0) :
    Martingale (fun t => μ[(fun ω => payoff n a (B T ω)) | ℱ t]) ℱ μ :=
  martingale_condExp _ ℱ μ

end
end BrownianCompensator
