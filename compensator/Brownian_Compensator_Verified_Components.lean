import Mathlib.Probability.BrownianMotion.Basic
import Mathlib.Probability.Martingale.Basic
import Mathlib.MeasureTheory.Integral.IntervalIntegral.IntegrationByParts
import Mathlib.Analysis.SpecialFunctions.Integrals.Basic
import Mathlib.Tactic
import Mathlib.Analysis.SpecialFunctions.Log.Deriv
import Mathlib.Analysis.SpecialFunctions.ExpDeriv

/-! Compiled components of the recorded Brownian compensator problem.
Full analytic solution: paper.pdf. Exact formal coverage: VERIFICATION.md.
The complete stochastic-integral compensator theorem remains open in Lean. -/

-- Source module: CompensatorCore.lean

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


-- Source module: HeatKernel.lean

/-! Actual derivatives of the logarithmic Gaussian heat kernel.
For v > 0 and c = (sqrt (2*pi))^-1 this is the usual heat kernel.
This module verifies the kernel PDE, not interchange of derivatives with the
convolution integral for the original payoff. -/

namespace BrownianCompensator
noncomputable section

def logKernel (c v x : ℝ) : ℝ :=
  c * Real.exp (-Real.log v / 2 - x^2 / (2*v))

theorem kernel_exponent_time (v x : ℝ) (hv : 0 < v) :
    HasDerivAt (fun u : ℝ => -Real.log u / 2 - x^2 / (2*u))
      ((x^2-v)/(2*v^2)) v := by
  have h1 := (Real.hasDerivAt_log hv.ne').neg.div_const 2
  have h2 := (hasDerivAt_const v (x^2)).div
    ((hasDerivAt_id v).const_mul 2) (by positivity : 2*v ≠ 0)
  have hh := h1.sub h2
  convert hh using 1
  · funext u; rfl
  · simp only [id_eq]
    field_simp
    ring

theorem kernel_exponent_space (v x : ℝ) (hv : 0 < v) :
    HasDerivAt (fun y : ℝ => -Real.log v / 2 - y^2 / (2*v)) (-x/v) x := by
  convert (hasDerivAt_const x (-Real.log v / 2)).sub
    ((hasDerivAt_pow 2 x).div_const (2*v)) using 1; field_simp; ring

theorem kernel_time_hasDerivAt (c v x : ℝ) (hv : 0 < v) :
    HasDerivAt (fun u => logKernel c u x)
      (logKernel c v x * ((x^2-v)/(2*v^2))) v := by
  simpa [logKernel, mul_assoc] using ((kernel_exponent_time v x hv).exp.const_mul c)

theorem kernel_space_hasDerivAt (c v x : ℝ) (hv : 0 < v) :
    HasDerivAt (fun y => logKernel c v y) (logKernel c v x * (-x/v)) x := by
  simpa [logKernel, mul_assoc] using ((kernel_exponent_space v x hv).exp.const_mul c)

theorem kernel_space_deriv (c v x : ℝ) (hv : 0 < v) :
    deriv (fun y => logKernel c v y) x = logKernel c v x * (-x/v) :=
  (kernel_space_hasDerivAt c v x hv).deriv

theorem kernel_second_deriv (c v x : ℝ) (hv : 0 < v) :
    deriv (fun y => deriv (fun z => logKernel c v z) y) x =
      logKernel c v x * (x^2/v^2 - 1/v) := by
  simp_rw [kernel_space_deriv c v _ hv]
  have h1 := kernel_space_hasDerivAt c v x hv
  have h2 := (hasDerivAt_id x).neg.div_const v
  have hh := (h1.mul h2).deriv
  have he : (fun y : ℝ => logKernel c v y * (-y/v)) =
      (fun y => logKernel c v y) * (fun y : ℝ => (-id) y / v) := by funext y; rfl
  simp only [Pi.neg_apply, id_eq] at hh
  rw [he]
  simp only [Pi.neg_apply, id_eq]
  rw [hh]
  field_simp
  ring

theorem kernel_heat_equation (c v x : ℝ) (hv : 0 < v) :
    deriv (fun u => logKernel c u x) v =
      (1/2 : ℝ) * deriv (fun y => deriv (fun z => logKernel c v z) y) x := by
  rw [(kernel_time_hasDerivAt c v x hv).deriv, kernel_second_deriv c v x hv]
  field_simp

theorem kernel_backward_heat (c T t x : ℝ) (ht : t < T) :
    deriv (fun u => logKernel c (T-u) x) t +
      (1/2 : ℝ) * deriv (fun y => deriv (fun z => logKernel c (T-t) z) y) x = 0 := by
  have hv : 0 < T-t := sub_pos.mpr ht
  have h := (kernel_time_hasDerivAt c (T-t) x hv).comp t
    ((hasDerivAt_const t T).sub (hasDerivAt_id t))
  have he : (fun u : ℝ => logKernel c (T-u) x) =
      (fun u => logKernel c u x) ∘ HSub.hSub T := by funext u; rfl
  rw [he, h.deriv, kernel_second_deriv c (T-t) x hv]
  field_simp
  ring

end
end BrownianCompensator


-- Source module: QuadraticTest.lean
open MeasureTheory ProbabilityTheory
open scoped NNReal ENNReal
namespace BrownianQuadraticTest

theorem gaussian_second_moment (v : ℝ≥0) :
    (∫ z : ℝ, z ^ 2 ∂gaussianReal 0 v) = (v : ℝ) := by
  have h : Var[(fun z : ℝ => z); gaussianReal 0 v] = v :=
    variance_fun_id_gaussianReal
  rw [variance_eq_integral measurable_id'.aemeasurable] at h
  simpa using h

theorem gaussian_convolution_square (v : ℝ≥0) (x : ℝ) :
    (∫ z : ℝ, (x + z)^2 ∂gaussianReal 0 v) = x^2 + (v : ℝ) := by
  let ν : Measure ℝ := gaussianReal 0 v
  have hz : Integrable (fun z : ℝ => z) ν := by
    exact memLp_one_iff_integrable.mp (memLp_id_gaussianReal' 1 (by norm_num))
  have hz2 : Integrable (fun z : ℝ => z ^ 2) ν := by
    exact (memLp_id_gaussianReal' 2 (by norm_num)).integrable_sq
  have hzlin : Integrable (fun z : ℝ => (2 * x) * z) ν := hz.const_mul _
  have hc : Integrable (fun _ : ℝ => x ^ 2) ν := integrable_const _
  calc
    (∫ z : ℝ, (x + z)^2 ∂ν)
        = ∫ z : ℝ, ((z ^ 2 + (2*x)*z) + x^2) ∂ν := by
            congr 1
            funext z
            ring
    _ = (∫ z : ℝ, z^2 ∂ν) + (2*x) * (∫ z : ℝ, z ∂ν) + x^2 := by
          have hi1 : (∫ z : ℝ, (z^2 + (2*x)*z) + x^2 ∂ν) =
              (∫ z : ℝ, z^2 + (2*x)*z ∂ν) + (∫ z : ℝ, x^2 ∂ν) := by
            simpa only [Pi.add_apply] using (integral_add (hz2.add hzlin) hc)
          have hi2 : (∫ z : ℝ, z^2 + (2*x)*z ∂ν) =
              (∫ z : ℝ, z^2 ∂ν) + (∫ z : ℝ, (2*x)*z ∂ν) := by
            simpa only [Pi.add_apply] using (integral_add hz2 hzlin)
          rw [hi1, hi2, integral_const_mul, integral_const]
          simp [ν]
    _ = x^2 + (v : ℝ) := by
          simp [ν, gaussian_second_moment, integral_id_gaussianReal]
          ring

/-- Quadratic Gaussian-semigroup test payoff; not the original nonsmooth payoff. -/
def quadraticS (T t x : ℝ) : ℝ := x ^ 2 + T - t

theorem quadraticS_time_derivative (T t x : ℝ) :
    deriv (fun u : ℝ => quadraticS T u x) t = -1 := by
  have hf : (fun u : ℝ => quadraticS T u x) =
      (fun u : ℝ => (x ^ 2 + T) - u) := by
    funext u
    simp only [quadraticS]
  rw [hf]
  convert ((hasDerivAt_const t (x ^ 2 + T)).sub (hasDerivAt_id t)).deriv using 1 <;>
    norm_num

theorem quadraticS_space_derivative (T t x : ℝ) :
    deriv (fun y : ℝ => quadraticS T t y) x = 2 * x := by
  have hf : (fun y : ℝ => quadraticS T t y) =
      (fun y : ℝ => y ^ 2 + (T - t)) := by
    funext y
    dsimp [quadraticS]
    ring
  rw [hf]
  convert ((hasDerivAt_pow 2 x).add_const (T - t)).deriv using 1; norm_num

theorem quadraticS_second_derivative (T t x : ℝ) :
    deriv (fun y : ℝ => deriv (fun z : ℝ => quadraticS T t z) y) x = 2 := by
  have h : (fun y : ℝ => deriv (fun z : ℝ => quadraticS T t z) y) =
      (fun y : ℝ => y * 2) := by
    funext y
    rw [quadraticS_space_derivative]
    ring
  rw [h]
  simpa only [id_eq, one_mul] using ((hasDerivAt_id x).mul_const 2).deriv

theorem quadratic_backward_heat (T t x : ℝ) :
    deriv (fun u : ℝ => quadraticS T u x) t +
      (1 / 2 : ℝ) * deriv
        (fun y : ℝ => deriv (fun z : ℝ => quadraticS T t z) y) x = 0 := by
  rw [quadraticS_time_derivative, quadraticS_second_derivative]
  ring

theorem integrand_polynomial_bound (x s C : ℝ)
    (h : s ^ 2 ≤ C ^ 2 * (1 + x ^ 2)) :
    (x * s) ^ 2 ≤ C ^ 2 * (x ^ 2 + x ^ 4) := by
  have hx : 0 ≤ x ^ 2 := sq_nonneg x
  have hmul := mul_le_mul_of_nonneg_left h hx
  nlinarith [sq_nonneg C]

theorem brownian_finite_moments
    {Ω : Type*} [MeasurableSpace Ω]
    (B : ℝ≥0 → Ω → ℝ) (μ : Measure Ω)
    (hB : IsBrownianReal B μ) (t : ℝ≥0)
    (p : ℝ≥0∞) (hp : p ≠ ∞) :
    MemLp (B t) p μ :=
  (hB.toIsPreBrownianReal.hasLaw_eval t).memLp
    (memLp_id_gaussianReal' p hp)

theorem brownian_L4
    {Ω : Type*} [MeasurableSpace Ω]
    (B : ℝ≥0 → Ω → ℝ) (μ : Measure Ω)
    (hB : IsBrownianReal B μ) (t : ℝ≥0) :
    MemLp (B t) 4 μ :=
  brownian_finite_moments B μ hB t 4 (by norm_num)

/-- Gaussian convolution agrees with the explicit quadratic heat solution. -/
theorem quadratic_gaussian_representation (T t x : ℝ) (ht : t ≤ T) :
    (∫ z : ℝ, (x + z)^2 ∂gaussianReal 0 (Real.toNNReal (T - t))) =
       quadraticS T t x := by
  have h := gaussian_convolution_square (Real.toNNReal (T - t)) x
  convert h using 1
  simp [quadraticS, max_eq_left (sub_nonneg.mpr ht)]
  ring

end BrownianQuadraticTest


-- Dependency audit of all included theorem statements.
#print axioms BrownianCompensator.digital_bounds
#print axioms BrownianCompensator.digital_measurable
#print axioms BrownianCompensator.atom_measurable
#print axioms BrownianCompensator.atomEnvelope_nonneg
#print axioms BrownianCompensator.atom_growth
#print axioms BrownianCompensator.payoff_measurable
#print axioms BrownianCompensator.envelope_nonneg
#print axioms BrownianCompensator.payoff_linear_growth
#print axioms BrownianCompensator.atom_original_indicator
#print axioms BrownianCompensator.primitive_zero
#print axioms BrownianCompensator.primitive_hasDerivAt
#print axioms BrownianCompensator.primitive_deriv
#print axioms BrownianCompensator.primitive_second_deriv
#print axioms BrownianCompensator.weighted_integration_by_parts
#print axioms BrownianCompensator.generator_cancellation
#print axioms BrownianCompensator.integrand_square_bound
#print axioms BrownianCompensator.integrated_moment_majorant
#print axioms BrownianCompensator.compensated_path_identity
#print axioms BrownianCompensator.brownian_finite_moments
#print axioms BrownianCompensator.terminal_payoff_memLp
#print axioms BrownianCompensator.terminal_payoff_condExp_martingale
#print axioms BrownianCompensator.kernel_exponent_time
#print axioms BrownianCompensator.kernel_exponent_space
#print axioms BrownianCompensator.kernel_time_hasDerivAt
#print axioms BrownianCompensator.kernel_space_hasDerivAt
#print axioms BrownianCompensator.kernel_space_deriv
#print axioms BrownianCompensator.kernel_second_deriv
#print axioms BrownianCompensator.kernel_heat_equation
#print axioms BrownianCompensator.kernel_backward_heat
#print axioms BrownianQuadraticTest.gaussian_second_moment
#print axioms BrownianQuadraticTest.gaussian_convolution_square
#print axioms BrownianQuadraticTest.quadraticS_time_derivative
#print axioms BrownianQuadraticTest.quadraticS_space_derivative
#print axioms BrownianQuadraticTest.quadraticS_second_derivative
#print axioms BrownianQuadraticTest.quadratic_backward_heat
#print axioms BrownianQuadraticTest.integrand_polynomial_bound
#print axioms BrownianQuadraticTest.brownian_finite_moments
#print axioms BrownianQuadraticTest.brownian_L4
#print axioms BrownianQuadraticTest.quadratic_gaussian_representation
