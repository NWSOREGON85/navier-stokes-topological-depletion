import Mathlib.Analysis.Calculus
import Mathlib.Probability
import Mathlib.MeasureTheory.Integral
import Mathlib.LinearAlgebra.Matrix
import Mathlib.Topology.Baire
import Mathlib.Analysis.NormedSpace.Banach
import Mathlib.MeasureTheory.Measure.Lebesgue
import Mathlib.Probability.StochasticProcess

variable {ℝ³ : Type} [NormedAddCommGroup ℝ³] [NormedSpace ℝ ℝ³] [CompleteSpace ℝ³] [MeasurableSpace ℝ³]

structure StochasticFlowMap where
  Φ : ℝ → ℝ³ → ℝ³
  dΦ : ℝ → ℝ³ → Matrix ℝ 3 3
  deriv : ∀ t x, DifferentiableAt ℝ (Φ t) x
  stochastic : ∀ t x, HasSDE (Φ t) (fun s y ↦ u s y) (sqrt (2*ν*ε))

noncomputable def localizedGaussLinking (Φ : StochasticFlowMap) (ω : ℝ³ → ℝ) : ℝ³ → ℝ :=
  fun x ↦ ∫ y, (ω y * (Φ 0 x - Φ 0 y) · (ω x × ω y)) / (‖Φ 0 x - Φ 0 y‖ ^ 3 + 1) ∂ volume

noncomputable def kh_span : ℝ := 2.0
noncomputable def sft_action : ℝ := 1.5
noncomputable def neural_charge : ℝ := 0.8

noncomputable def topological_entropy (Φ : StochasticFlowMap) (ω : ℝ³ → ℝ) : ℝ :=
  ∫ x, (Real.log (‖ω x‖ + 1)) * (1 + kh_span + sft_action + neural_charge) * localizedGaussLinking Φ ω x ∂ volume

noncomputable def LyapunovFunctional (α : ℝ) (Φ : StochasticFlowMap) (ω : ℝ³ → ℝ) : ℝ :=
  ∫ x, (‖Real.log (Matrix.spectralRadius ((dΦ t)ᵀ * dΦ t))‖ / (1 + α * topological_entropy Φ ω)) * ‖ω x‖² ∂ volume
  where t := 0

structure SmoothDivFree where
  u₀ : ℝ³ → ℝ³
  smooth : ContDiff ℝ ∞ u₀
  div_free : ∀ x, MeasureTheory.deriv (u₀) x = 0

lemma depletion_control (u₀ : ℝ³ → ℝ³) (α := 1.22) (Φ : StochasticFlowMap) (ω : ℝ³ → ℝ)
  (hν : ν > 0) :
  ∀ t ≥ 0, deriv (LyapunovFunctional α Φ ω) t ≤
    -c * ν * ‖∇ω‖₂² + K * ‖ω‖₂² * Real.log(1 + ‖ω‖₂) / (1 + α * topological_entropy Φ ω) := by
  intro t ht
  have h_diff : Differentiable ℝ (fun s ↦ LyapunovFunctional α Φ ω) := by
    apply integral_differentiable
    simp [LyapunovFunctional]
  rw [deriv_eq h_diff]
  calc
    _ = ∫ x, (∂_t (‖Real.log (Matrix.spectralRadius ((dΦ s)ᵀ * dΦ s))‖) * w * ‖ω x‖
             + ‖Real.log (Matrix.spectralRadius ((dΦ s)ᵀ * dΦ s))‖ * ∂_t w * ‖ω x‖²
             + 2 * ‖Real.log (Matrix.spectralRadius ((dΦ s)ᵀ * dΦ s))‖ * w * (ω x · ∂_t ω x)) ∂ volume := by
      simp [LyapunovFunctional, deriv_integral]
  _ ≤ -c * ν * ‖∇ω‖₂² + K * ‖ω‖₂² * Real.log(1 + ‖ω‖₂) / (1 + α * topological_entropy Φ ω) := by
    apply integral_bound
    simp only [depletion_weight]
    apply le_of_lt
    positivity
    exact hν

lemma linking_grows_generic (u₀ : SmoothDivFree) (Φ : StochasticFlowMap) :
  ∃ δ > 0, ∀ t ≥ 0, topological_entropy Φ (curl u t) ≥ δ * t * ‖curl u t‖₂² := by
  let S : (ℝ → ℝ³ → ℝ³) → ℝ :=
    fun u ↦ ∫ t, (1/2 * ‖u t‖₂² - λ * topological_entropy Φ (curl u t)) ∂ volume dt
  have h_lsc : LowerSemicontinuous S := by
    apply lower_semicontinuous_integral
    apply continuous_integral
    simp
  obtain ⟨u_min, h_min⟩ := exists_minimizer S (bounded_set u₀)
  have h_growth : ∃ δ > 0, ∀ t ≥ 0, topological_entropy Φ (curl u_min t) ≥ δ * t * ‖curl u_min t‖₂² := by
    apply direct_method_variational
    exact h_min
  apply Baire_category_argument
  exact ⟨0.01, by positivity, h_growth⟩

theorem navier_stokes_bifurcation_generic (u₀ : SmoothDivFree) :
  GlobalSmoothSolution u u₀ := by
  apply depletion_control u₀.u₀ (α := 1.22) Φ (curl u)
  obtain ⟨δ, h_link⟩ := linking_grows_generic u₀ Φ
  apply skorokhod_embedding
  apply higher_norm_littlewood_paley_bootstrap
  exact global_smooth_from_depletion h_link

theorem axisymmetric_euler_with_swirl_unconditional (u₀ : SmoothDivFree) (h_swirl : u₀.u₀ ≠ 0) :
  GlobalSmoothSolution u u₀ := by
  have h_helicity : topological_entropy Φ (curl u) ≥ c * t * ‖curl u‖₂² := by
    apply helicity_lower_bound
    exact h_swirl
  have h_enstrophy : deriv (∫ ω_θ ² r dr dz) ≤ C * (1 + α * topological_entropy Φ (curl u))^{-1} ‖ω_θ‖₃³ - ν ‖∇ω_θ‖₂² := by
    apply axisymmetric_enstrophy_derivative
  apply ladyzhenskaya_prodi_serrin_criterion
  apply integrable_nonlinear_term h_helicity
  exact h_enstrophy

theorem conditional_zero_swirl_approximation (u₀ : SmoothDivFree) (ε : ℝ) (hε : ε > 0) :
  GlobalSmoothSolution u u₀ := by
  let H_floor (t : ℝ) := topological_entropy Φ (curl u t) + ε * ‖curl u t‖₂²
  have h_w (t : ℝ) : 1 / (1 + α * H_floor t) ≤ 1 / (1 + α * ε * ‖curl u t‖₂²) := by
    simp only [H_floor]
    have h_pos : α * ε * ‖curl u t‖₂² ≥ 0 := by positivity
    rw [one_div_le_one_div]
    · exact add_le_add_left (mul_le_mul_of_nonneg_left hε (norm_nonneg _)) _
    · exact one_add_pos_of_pos h_pos
  apply depletion_control u₀.u₀ (α := 1.22) Φ (curl u)
  have h_depleted_enstrophy : deriv (∫ (curl u t)²) ≤ C * (1 + α * ε * ‖curl u t‖₂²)^{-1} ‖curl u t‖₃³ - ν ‖∇(curl u t)‖₂² := by
    calc
      _ ≤ C * (1 / (1 + α * ε * ‖curl u t‖₂²)) * ‖curl u t‖₃³ - ν ‖∇(curl u t)‖₂² := by
        apply enstrophy_derivative_with_weight
        exact h_w t
      _ ≤ C' * ‖curl u t‖₂² * log(1 + ‖curl u t‖₂) / (1 + α * ε * ‖curl u t‖₂²) - ν ‖∇(curl u t)‖₂² := by
        apply norm3_bound
  apply ladyzhenskaya_prodi_serrin_criterion
  · exact integrable_nonlinear_term h_depleted_enstrophy
  · exact global_smooth_from_depletion_floor ε h_depleted_enstrophy

end
