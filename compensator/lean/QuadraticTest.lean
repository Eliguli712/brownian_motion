import Mathlib.Probability.BrownianMotion.Basic
import Mathlib.Tactic
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
