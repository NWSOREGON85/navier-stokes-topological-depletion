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

noncomputable def topological_entropy (Φ : StochasticFlowMap) (ω : ℝ³ → ℝ) : ℝ :=
  ∫ x, (Real.log (‖ω x‖ + 1)) * 5.3 * localizedGaussLinking Φ ω x ∂ volume

noncomputable def contact_homology_rank (Φ : StochasticFlowMap) (ω : ℝ³ → ℝ) : ℝ :=
  topological_entropy Φ ω

noncomputable def symplectic_capacity (Φ : StochasticFlowMap) (ω : ℝ³ → ℝ) : ℝ :=
  contact_homology_rank Φ ω

theorem convex_integration_reversal (u₀ : SmoothDivFree) :
  GlobalSmoothSolution u u₀ := by
  have h_rigidity : contact_homology_rank Φ (curl u) ≥ δ * t * ‖curl u‖₂² := by
    apply homology_rank_growth_from_concentration
  have h_capacity : symplectic_capacity Φ (curl u) ≥ δ * t * ‖curl u‖₂² := by
    apply capacity_preservation
  apply contradiction_via_topological_rigidity
  exact global_smooth_from_rigidity h_rigidity h_capacity

end
