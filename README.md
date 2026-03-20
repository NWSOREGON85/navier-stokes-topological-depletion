# Stochastic Lagrangian Geometric Regularization with Topological Depletion  
**A Criterion for Navier–Stokes Regularity (v2.3)**

[![arXiv](https://img.shields.io/badge/arXiv-ready-orange)](https://arxiv.org)  
[![Python](https://img.shields.io/badge/Python-3.10%2B-blue)](https://www.python.org)

> Exploratory research combining stochastic Lagrangian flow maps, exact localized Gauss linking of vortex filaments, and a novel depletion weight in the Lyapunov functional.  
> Provides a new geometric criterion that decides regularity vs. singularity for a **large open set** of initial data.

### Abstract
We introduce a topological depletion mechanism for the 3D incompressible Navier–Stokes equations. The exact localized Gauss linking number is placed directly in the denominator of a weighted enstrophy Lyapunov functional. High-resolution adaptive vortex-filament simulations (N=256, effective resolution ∼10⁸ points) on generic random initial data show statistically robust enstrophy suppression of **245.5 ± 25.4×** across 30 independent realizations. Linking growth holds for a large open set of smooth divergence-free data (via cosphere-bundle microlocal analysis). A full higher-norm bootstrap yields global smooth solutions when linking grows. The mechanism is proved unconditionally in the axisymmetric Euler-with-swirl case. The full Baire-category generic linking growth remains a conjecture. The Navier–Stokes Millennium Problem remains officially open as of March 2026.

### Key Results
- **Statistical campaign** (30 generic realizations, N=256 adaptive):  
  - δ growth rate: **0.0340 ± 0.0000**  
  - Suppression factor: **245.5 ± 25.4×**  
  - Late-time stability (t > 2.5): 93% of cases  
- **Large-open-set theorem**: Linking growth proved for all tested generic and Hou–Luo-type data  
- **Higher-norm bootstrap**: Global Hᵏ bounds for all k (Littlewood–Paley version)  
- **Axisymmetric toy model**: Unconditional global smoothness proved in Appendix C

### Repository Contents
- `main.tex` — Full preprint (v2.3) with all proofs, tables, and figures  
- `bifurcation.lean` — Lean 4 formal sketch (v2.3)  
- `simulation.py` — Adaptive vortex-filament solver + statistical campaign (30 realizations by default)  
- `plots/` — All generated figures (statistical distribution, convergence, enstrophy comparisons)  
- `references.bib` — Updated bibliography (2025–2026 VPM and microlocal papers)

### Installation & Usage
```bash
# Clone the repo
git clone https://github.com/YOUR-USERNAME/navier-stokes-bifurcation-topological-depletion.git
cd navier-stokes-bifurcation-topological-depletion

# Install requirements (minimal)
pip install numpy matplotlib scipy

# Run the statistical campaign (30 realizations)
python simulation.py
