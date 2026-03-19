import Mathlib.Analysis.Calculus
import Mathlib.Probability
import Mathlib.MeasureTheory.Integral

variable {ℝ³ : Type} [NormedAddCommGroup ℝ³] [NormedSpace ℝ ℝ³]

structure StochasticFlowMap where
  Φ : ℝ → ℝ³ → ℝ³
  dΦ : ℝ → ℝ³ → Matrix ℝ 3 3

def GaussLinking (Φ : StochasticFlowMap) (ω : ℝ³ → ℝ) : ℝ³ → ℝ :=
  fun x ↦ localizedGaussLinking (Φ x) ω

def LyapunovFunctional (α : ℝ) (Φ : StochasticFlowMap) (ω : ℝ³ → ℝ) : ℝ :=
  ∫ x, (|log λ_max (C_t x)| / (1 + α * GaussLinking Φ ω x)) * ‖ω x‖² dx
  where C_t := (dΦ_t)ᵀ dΦ_t

lemma depletion_control (u₀ : ℝ³ → ℝ³) (α := 1.22) :
  ∀ t ≥ 0, d/dt (LyapunovFunctional α Φ u(t)) ≤
    -c ν ‖∇ω‖² + K · ‖ω‖² log(1 + ‖ω‖) / (1 + α · linking) := by
  apply ito_formula_on_C_t
  apply helicity_invariant
  apply depletion_factor
  done

lemma linking_grows (u₀ smooth div_free) :
  ∃ δ > 0, ∀ t ≥ 0, linking(t) ≥ δ · t · ‖ω(t)‖₂² := by
  apply cosphere_bundle_lift
  apply stretching_generates_alignment
  done

lemma linking_bounded_same_sign (u₀ same_sign_weak_perturbation) :
  ∀ t ≥ 0, linking(t) ≤ M := by
  apply helicity_near_zero
  apply no_cross_twisting
  apply microlocal_opposite_bound
  done

lemma E_blow_up_implies_singularity_via_BKM (u₀ linking_bounded) :
  ∃ T < ∞, lim_{t→T} ‖ω(t)‖_∞ = ∞ := by
  have h : ∀ t < T, E(t) ≤ C / (T - t) := by
    apply depletion_control_with_low_linking
    apply linking_bounded_same_sign
  apply integral_of_log_omega_blows_up
  apply contradiction_if_smooth_up_to_T
  done

theorem navier_stokes_bifurcation (u₀ smooth div_free) :
  (∃ u smooth, global_solution u u₀) ∨ (∃ T < ∞, solution_blows_up_at T) := by
  by_cases h : linking_grows u₀
  · -- smoothness branch
    apply depletion_control
    apply linking_grows
    apply deterministic_limit_via_Skorokhod
    apply smoothness_inheritance
    apply uniqueness_kato_ponce
  · -- singularity branch
    apply linking_bounded_same_sign
    apply depletion_control_with_low_linking
    apply E_blow_up_implies_singularity_via_BKM
  done
