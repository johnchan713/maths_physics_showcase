#ifndef MATHS_ANALYSIS_CLARKE_SUBDIFFERENTIALS_HPP
#define MATHS_ANALYSIS_CLARKE_SUBDIFFERENTIALS_HPP

#include <string>
#include <functional>
#include <vector>
#include <cmath>

/**
 * @file clarke_subdifferentials.hpp
 * @brief Clarke subdifferentials and generalized gradients
 *
 * Implements:
 * - Clarke generalized directional derivative
 * - Clarke subdifferential for locally Lipschitz functions
 * - Calculus rules for Clarke subdifferentials
 * - Regularity and semi-smoothness
 * - Applications to non-smooth optimization
 */

namespace maths::analysis {

/**
 * @class ClarkeSubdifferential
 * @brief Clarke subdifferential for locally Lipschitz functions
 */
class ClarkeSubdifferential {
public:
    /**
     * @brief Clarke generalized directional derivative
     */
    static std::string clarkeDirectionalDerivative() {
        return "Clarke Generalized Directional Derivative:\n"
               "\n"
               "For f: X → ℝ locally Lipschitz near x̄, direction v ∈ X:\n"
               "\n"
               "  f°(x̄; v) = lim sup_{x→x̄, t↓0} [f(x + tv) - f(x)] / t\n"
               "\n"
               "Properties:\n"
               "1. Finite: |f°(x̄; v)| ≤ K‖v‖ (K = Lipschitz constant)\n"
               "2. Positively homogeneous: f°(x̄; λv) = λf°(x̄; v) for λ > 0\n"
               "3. Subadditive: f°(x̄; v+w) ≤ f°(x̄; v) + f°(x̄; w)\n"
               "4. Lipschitz in v: |f°(x̄; v) - f°(x̄; w)| ≤ K‖v - w‖\n"
               "5. Upper semicontinuous in (x, v)\n"
               "\n"
               "Interpretation:\n"
               "- Limsup (not lim!) of difference quotients\n"
               "- Takes \"worst case\" variation\n"
               "- Always exists for locally Lipschitz f\n"
               "- Convex function of v\n"
               "\n"
               "Relationship to derivatives:\n"
               "- f C¹: f°(x̄; v) = ⟨∇f(x̄), v⟩\n"
               "- f convex: f°(x̄; v) = f'(x̄; v) (directional derivative)\n"
               "- General: f°(x̄; v) ≥ f'(x̄; v) when f'(x̄; v) exists\n"
               "\n"
               "Example:\n"
               "  f(x) = |x| at x = 0, v = 1\n"
               "  f°(0; 1) = lim sup_{x→0, t↓0} |x + t|/t - |x|/t\n"
               "           = 1\n"
               "  Similarly f°(0; -1) = 1\n"
               "  f°(0; v) = |v|\n"
               "\n"
               "Rademacher's theorem:\n"
               "Locally Lipschitz f differentiable a.e.\n"
               "  f°(x̄; v) = lim sup_{x→x̄, x∈Ω_f} ⟨∇f(x), v⟩\n"
               "where Ω_f = set of differentiability points\n"
               "\n"
               "Convexity:\n"
               "f°(x̄; ·) is convex, Lipschitz function of v\n"
               "⇒ Has subdifferential ∂[f°(x̄; ·)](0)\n"
               "This is Clarke subdifferential!";
    }

    /**
     * @brief Clarke subdifferential definition
     */
    static std::string definition() {
        return "Clarke Subdifferential:\n"
               "\n"
               "For f: X → ℝ locally Lipschitz near x̄:\n"
               "\n"
               "  ∂_Cf(x̄) = {v* ∈ X* : ⟨v*, v⟩ ≤ f°(x̄; v) ∀v ∈ X}\n"
               "\n"
               "Equivalently:\n"
               "  ∂_Cf(x̄) = ∂[f°(x̄; ·)](0)\n"
               "\n"
               "Subdifferential of v ↦ f°(x̄; v) at v = 0\n"
               "\n"
               "Alternative characterization:\n"
               "  ∂_Cf(x̄) = conv{lim ∇f(x_k) : x_k → x̄, x_k ∈ Ω_f}\n"
               "\n"
               "Convex hull of limits of gradients at differentiability points\n"
               "\n"
               "Properties:\n"
               "1. ∂_Cf(x̄) non-empty, convex, weak* compact\n"
               "2. Upper semicontinuous multifunction\n"
               "3. Bounded: ‖v*‖ ≤ K (Lipschitz constant)\n"
               "4. f C¹: ∂_Cf(x̄) = {∇f(x̄)}\n"
               "5. f convex: ∂_Cf(x̄) = ∂f(x̄) (convex subdifferential)\n"
               "\n"
               "Examples:\n"
               "\n"
               "1. f(x) = |x| on ℝ:\n"
               "   ∂_Cf(x) = {sgn(x)} for x ≠ 0\n"
               "   ∂_Cf(0) = [-1, 1]\n"
               "\n"
               "2. f(x) = max{x₁, x₂, ..., x_n}:\n"
               "   ∂_Cf(x) = conv{e_i : x_i = f(x)}\n"
               "   Convex hull of active unit vectors\n"
               "\n"
               "3. f(x) = ‖x‖ (any norm):\n"
               "   ∂_Cf(x) = {v* : ‖v*‖* = 1, ⟨v*, x⟩ = ‖x‖} for x ≠ 0\n"
               "   ∂_Cf(0) = unit ball in dual norm\n"
               "\n"
               "Optimality condition:\n"
               "x̄ local minimizer ⟺ 0 ∈ ∂_Cf(x̄)\n"
               "(Necessary and sufficient for Lipschitz f)\n"
               "\n"
               "Relationships:\n"
               "  ∂_Pf(x̄) ⊆ ∂_Ff(x̄) ⊆ ∂_Lf(x̄) ⊆ ∂_Cf(x̄)\n"
               "Equality for regular f";
    }

    /**
     * @brief Clarke calculus rules
     */
    static std::string calculusRules() {
        return "Calculus Rules for Clarke Subdifferential:\n"
               "\n"
               "For locally Lipschitz functions f, g:\n"
               "\n"
               "1. SUM RULE (always exact):\n"
               "   ∂_C(f + g)(x̄) = ∂_Cf(x̄) + ∂_Cg(x̄)\n"
               "\n"
               "   No qualification needed!\n"
               "   Major advantage of Clarke subdifferential\n"
               "\n"
               "2. SCALAR MULTIPLICATION:\n"
               "   ∂_C(λf)(x̄) = λ∂_Cf(x̄) for λ ∈ ℝ\n"
               "\n"
               "3. MAXIMUM RULE:\n"
               "   For f(x) = max{f_i(x) : i = 1,...,m}:\n"
               "   ∂_Cf(x̄) = conv{⋃_{i ∈ I(x̄)} ∂_Cf_i(x̄)}\n"
               "   where I(x̄) = {i : f_i(x̄) = f(x̄)}\n"
               "\n"
               "   Always exact (no qualification)\n"
               "\n"
               "4. CHAIN RULE:\n"
               "   For f: ℝⁿ → ℝ, g: ℝᵐ → ℝⁿ locally Lipschitz:\n"
               "   ∂_C(f ∘ g)(x̄) ⊆ ∂_C(f̂ ∘ g)(x̄)\n"
               "   where f̂(y) = sup_{z∈∂_Cf(y)} ⟨z, ·⟩\n"
               "\n"
               "   Special case (g C¹):\n"
               "   ∂_C(f ∘ g)(x̄) ⊆ Dg(x̄)*∂_Cf(g(x̄))\n"
               "\n"
               "5. PRODUCT RULE:\n"
               "   ∂_C(fg)(x̄) ⊆ f(x̄)∂_Cg(x̄) + g(x̄)∂_Cf(x̄)\n"
               "\n"
               "   Inclusion may be strict\n"
               "\n"
               "6. QUOTIENT (when g(x̄) ≠ 0):\n"
               "   ∂_C(f/g)(x̄) ⊆ [g(x̄)∂_Cf(x̄) - f(x̄)∂_Cg(x̄)] / g(x̄)²\n"
               "\n"
               "7. COMPOSITION WITH LINEAR:\n"
               "   For A: X → Y linear:\n"
               "   ∂_C(f ∘ A)(x̄) = A*∂_Cf(Ax̄)\n"
               "   Exact!\n"
               "\n"
               "8. INDICATOR FUNCTION:\n"
               "   For I_C (indicator of closed set C):\n"
               "   ∂_CI_C(x̄) = N_C^C(x̄) (Clarke normal cone)\n"
               "\n"
               "9. DISTANCE FUNCTION:\n"
               "   For d_C(x) = dist(x, C):\n"
               "   ∂_Cd_C(x̄) ⊆ {v* : ‖v*‖ = 1, ⟨v*, x̄-p⟩ = d_C(x̄),\n"
               "                      p ∈ proj_C(x̄)}\n"
               "\n"
               "Key advantage:\n"
               "Clarke subdifferential has \"good\" calculus:\n"
               "- Sum rule always exact\n"
               "- Max rule always exact\n"
               "- Chain rule reasonable\n"
               "\n"
               "Trade-off:\n"
               "+ Robust calculus\n"
               "- May be large (many subgradients)\n"
               "- Less sharp than Fréchet";
    }

    /**
     * @brief Regularity concepts
     */
    static std::string regularity() {
        return "Regularity and Semi-Smoothness:\n"
               "\n"
               "CLARKE REGULAR FUNCTION:\n"
               "f locally Lipschitz is Clarke regular at x̄ if:\n"
               "  f'(x̄; v) exists and f'(x̄; v) = f°(x̄; v) ∀v\n"
               "\n"
               "Equivalently:\n"
               "  ∂_Cf(x̄) = ∂_Ff(x̄) (Fréchet = Clarke)\n"
               "\n"
               "Examples of regular functions:\n"
               "1. C¹ functions (always regular)\n"
               "2. Convex functions (directional derivative exists)\n"
               "3. Sums/compositions of regular functions (under conditions)\n"
               "4. max{f_i} when f_i regular and active set \"stable\"\n"
               "\n"
               "Non-regular example:\n"
               "  f(x) = -|x| at x = 0\n"
               "  f'(0; 1) = -1, but f°(0; 1) = +1\n"
               "  ∂_Cf(0) = [-1, 1], ∂_Ff(0) = ∅\n"
               "\n"
               "STRICTLY DIFFERENTIABLE:\n"
               "f strictly differentiable at x̄ if:\n"
               "  lim_{x,y→x̄} [f(x) - f(y) - ⟨v*, x-y⟩] / ‖x - y‖ = 0\n"
               "for some v* = ∇f(x̄)\n"
               "\n"
               "Stronger than Fréchet differentiability\n"
               "Uniform convergence in both x and y\n"
               "\n"
               "SEMI-SMOOTH FUNCTION:\n"
               "f: ℝⁿ → ℝ locally Lipschitz is semi-smooth at x̄ if:\n"
               "  f'(x̄; h) exists and\n"
               "  ⟨v*, h⟩ → f'(x̄; h)\n"
               "as v* ∈ ∂_Cf(x̄ + t h), t ↓ 0\n"
               "\n"
               "Examples:\n"
               "- C¹ functions\n"
               "- Convex functions\n"
               "- max{f₁,...,f_m} with C¹ components\n"
               "- Composition of semi-smooth functions\n"
               "\n"
               "P-ORDER SEMI-SMOOTHNESS:\n"
               "∃C > 0, δ > 0: ∀h with ‖h‖ < δ, ∀v* ∈ ∂_Cf(x̄ + h):\n"
               "  |f(x̄ + h) - f(x̄) - ⟨v*, h⟩| ≤ C‖h‖^{1+p}\n"
               "\n"
               "p = 0: semi-smooth\n"
               "p = 1: strongly semi-smooth\n"
               "\n"
               "Newton method for semi-smooth equations:\n"
               "Solve F(x) = 0 where F semi-smooth:\n"
               "  x_{k+1} = x_k - V_k⁻¹F(x_k)\n"
               "where V_k ∈ ∂_CF(x_k)\n"
               "\n"
               "Convergence:\n"
               "- Semi-smooth F: superlinear convergence\n"
               "- Strongly semi-smooth: quadratic (like smooth Newton)\n"
               "\n"
               "Applications:\n"
               "- Non-smooth equations\n"
               "- Complementarity problems\n"
               "- Variational inequalities\n"
               "- PDE-constrained optimization";
    }
};

/**
 * @class ClarkeApplications
 * @brief Applications of Clarke subdifferential
 */
class ClarkeApplications {
public:
    /**
     * @brief Non-smooth optimization
     */
    static std::string nonsmoothOptimization() {
        return "Clarke Subdifferential in Optimization:\n"
               "\n"
               "UNCONSTRAINED OPTIMIZATION:\n"
               "  min f(x), f locally Lipschitz\n"
               "\n"
               "Optimality condition:\n"
               "  x* local min ⟺ 0 ∈ ∂_Cf(x*)\n"
               "\n"
               "Necessary and sufficient!\n"
               "\n"
               "Descent direction:\n"
               "If 0 ∉ ∂_Cf(x), then d = -v* for any v* ∈ ∂_Cf(x)\n"
               "satisfies f°(x; d) < 0 (descent)\n"
               "\n"
               "Steepest descent:\n"
               "  d* = argmin{f°(x; d) : ‖d‖ = 1}\n"
               "Solution: d* = -v*/‖v*‖ where v* minimizes ‖v*‖ over ∂_Cf(x)\n"
               "\n"
               "CONSTRAINED OPTIMIZATION:\n"
               "  min f(x) s.t. g_i(x) ≤ 0, h_j(x) = 0\n"
               "\n"
               "KKT conditions (necessary):\n"
               "If x* local min, ∃λ* ≥ 0, μ*:\n"
               "  0 ∈ ∂_Cf(x*) + Σλ_i*∂_Cg_i(x*) + Σμ_j*∂_Ch_j(x*)\n"
               "  λ_i*g_i(x*) = 0 (complementarity)\n"
               "  g_i(x*) ≤ 0, h_j(x*) = 0 (feasibility)\n"
               "\n"
               "Under constraint qualification:\n"
               "- Mangasarian-Fromovitz CQ\n"
               "- Linear independence CQ\n"
               "- Slater CQ (for convex)\n"
               "\n"
               "BUNDLE METHODS:\n"
               "Cutting-plane algorithm for non-smooth f:\n"
               "\n"
               "At iteration k:\n"
               "1. Compute v_k ∈ ∂_Cf(x_k)\n"
               "2. Build model: f̂_k(x) = max{f(x_i) + ⟨v_i, x-x_i⟩ : i ≤ k}\n"
               "3. Solve: x_{k+1} = argmin f̂_k(x) + (μ/2)‖x - x_k‖²\n"
               "4. Update bundle\n"
               "\n"
               "Convergence: 0 ∈ ∂_Cf(x*) at limit\n"
               "\n"
               "PROXIMAL BUNDLE:\n"
               "Stabilized by adding ‖x - x_k‖² term\n"
               "More robust in practice\n"
               "\n"
               "Example application:\n"
               "ℓ₁ regularization: min f(x) + λ‖x‖₁\n"
               "∂_C‖x‖₁ easily computed\n"
               "Clarke calculus applies";
    }

    /**
     * @brief Variational inequalities
     */
    static std::string variationalInequalities() {
        return "Variational Inequalities with Clarke Subdifferential:\n"
               "\n"
               "VARIATIONAL INEQUALITY (VI):\n"
               "Find x* ∈ C such that:\n"
               "  ⟨F(x*), x - x*⟩ ≥ 0, ∀x ∈ C\n"
               "\n"
               "Non-smooth case:\n"
               "F(x) = ∂_Cf(x) for f locally Lipschitz\n"
               "\n"
               "Find x* ∈ C: ⟨v*, x - x*⟩ ≥ 0, ∀x ∈ C, ∀v* ∈ ∂_Cf(x*)\n"
               "\n"
               "Equivalence with optimization:\n"
               "x* solves VI ⟺ x* minimizes f over C\n"
               "\n"
               "GENERALIZED EQUATION:\n"
               "  0 ∈ F(x) + N_C(x)\n"
               "\n"
               "where N_C is Clarke normal cone\n"
               "\n"
               "For F(x) = ∂_Cf(x):\n"
               "  0 ∈ ∂_Cf(x) + N_C(x)\n"
               "\n"
               "This is first-order condition for min f(x) s.t. x ∈ C\n"
               "\n"
               "COMPLEMENTARITY PROBLEM:\n"
               "Find x such that:\n"
               "  x ≥ 0, F(x) ≥ 0, ⟨x, F(x)⟩ = 0\n"
               "\n"
               "Non-smooth F:\n"
               "Replace F(x) with selection from ∂_Cf(x)\n"
               "\n"
               "HEMIVARIATIONAL INEQUALITY:\n"
               "Find u such that:\n"
               "  ⟨Au, v - u⟩ + ∫_Ω f°(x, u; v-u) dx ≥ ⟨l, v - u⟩\n"
               "for all v\n"
               "\n"
               "Uses Clarke directional derivative f°\n"
               "Arises in non-monotone mechanics\n"
               "\n"
               "Applications:\n"
               "- Contact mechanics (non-monotone friction)\n"
               "- Plasticity with hardening\n"
               "- Crack propagation\n"
               "- Non-convex energy minimization\n"
               "\n"
               "Existence theory:\n"
               "- Pseudomonotone operators\n"
               "- Coercivity\n"
               "- Compactness\n"
               "\n"
               "Numerical methods:\n"
               "- Projection methods\n"
               "- Splitting methods\n"
               "- Relaxation methods\n"
               "- Smoothing + Newton";
    }
};

} // namespace maths::analysis

#endif // MATHS_ANALYSIS_CLARKE_SUBDIFFERENTIALS_HPP
