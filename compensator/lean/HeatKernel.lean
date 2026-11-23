import Mathlib.Analysis.SpecialFunctions.Log.Deriv
import Mathlib.Analysis.SpecialFunctions.ExpDeriv
import Mathlib.Tactic

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
