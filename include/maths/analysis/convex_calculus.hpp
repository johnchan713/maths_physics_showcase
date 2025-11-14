#ifndef MATHS_ANALYSIS_CONVEX_CALCULUS_HPP
#define MATHS_ANALYSIS_CONVEX_CALCULUS_HPP

#include <string>
#include <functional>
#include <vector>
#include <cmath>

/**
 * @file convex_calculus.hpp
 * @brief Advanced convex calculus: fuzzy rules, exact rules, mean value theorems
 *
 * Implements:
 * - Fuzzy calculus rules in convex analysis
 * - Exact calculus rules
 * - Mean value theorems for subdifferentials
 * - Smoothness of norms
 * - Favorable classes of Banach spaces
 */

namespace maths::analysis {

/**
 * @class FuzzyCalculus
 * @brief Fuzzy (approximate) calculus rules for subdifferentials
 */
class FuzzyCalculus {
public:
    /**
     * @brief Fuzzy sum rule
     */
    static std::string fuzzySum() {
        return "Fuzzy Sum Rule:\n"
               "\n"
               "For f, g: X → ℝ ∪ {+∞} lsc, x̄ ∈ dom(f+g):\n"
               "\n"
               "∀ε > 0, ∀v* ∈ ∂_ε(f+g)(x̄), ∃x_f, x_g ∈ B_ε(x̄):\n"
               "  v* ∈ ∂_{2ε}f(x_f) + ∂_{2ε}g(x_g)\n"
               "\n"
               "where ∂_εf(x) = ε-subdifferential:\n"
               "  {v* : f(y) ≥ f(x) + ⟨v*, y-x⟩ - ε, ∀y}\n"
               "\n"
               "Interpretation:\n"
               "- Approximate subdifferential of sum\n"
               "- Can be decomposed into approximate subdifferentials of summands\n"
               "- Points x_f, x_g may differ from x̄\n"
               "\n"
               "Key insight:\n"
               "Fuzzy rules trade exactness for generality:\n"
               "- No constraint qualifications needed\n"
               "- Always valid for lsc functions\n"
               "- Useful when exact rules fail\n"
               "\n"
               "Applications:\n"
               "- Non-smooth optimization algorithms\n"
               "- Approximate optimality conditions\n"
               "- Bundle methods";
    }

    /**
     * @brief Fuzzy chain rule
     */
    static std::string fuzzyChain() {
        return "Fuzzy Chain Rule:\n"
               "\n"
               "For f: Y → ℝ lsc, g: X → Y, x̄ ∈ X:\n"
               "\n"
               "Let h = f ∘ g. For ε > 0, v* ∈ ∂_εh(x̄):\n"
               "\n"
               "∃y ∈ B_ε(g(x̄)), w* ∈ ∂_{2ε}f(y):\n"
               "  v* ∈ g'*(x̄; w*) + εB*\n"
               "\n"
               "where g'*(x̄; w*) is coderivative.\n"
               "\n"
               "Smooth case:\n"
               "If g Fréchet differentiable at x̄:\n"
               "  v* ∈ Dg(x̄)*w* + εB*\n"
               "\n"
               "Special case (linear g = A):\n"
               "  ∂_ε(f ∘ A)(x) ⊆ A*∂_{2ε}f(Ax + δ) + εB*\n"
               "  for some δ with ‖δ‖ ≤ ε\n"
               "\n"
               "Error accumulation:\n"
               "- Input ε-subdifferential\n"
               "- Output 2ε (or worse) subdifferential\n"
               "- Need careful error management\n"
               "\n"
               "Applications:\n"
               "- Composite optimization\n"
               "- Regularized problems\n"
               "- Penalty methods";
    }

    /**
     * @brief Fuzzy maximum rule
     */
    static std::string fuzzyMaximum() {
        return "Fuzzy Maximum Rule:\n"
               "\n"
               "For f(x) = max{f_i(x) : i ∈ I} with |I| < ∞:\n"
               "\n"
               "Let I_ε(x̄) = {i : f_i(x̄) ≥ f(x̄) - ε} (ε-active indices)\n"
               "\n"
               "Then: ∂_εf(x̄) ⊆ conv{⋃_{i ∈ I_ε(x̄)} ∂_{2ε}f_i(x̄)}\n"
               "\n"
               "Interpretation:\n"
               "- Nearly active functions contribute\n"
               "- Exact active set I(x̄) replaced by ε-active set I_ε(x̄)\n"
               "- Convex hull of approximate subdifferentials\n"
               "\n"
               "Example (two functions):\n"
               "  f(x) = max{x², 1-x}\n"
               "  At x = 0.9: f₁(0.9) = 0.81, f₂(0.9) = 0.1\n"
               "  For small ε, only f₁ is ε-active\n"
               "  ∂_εf(0.9) ≈ {2·0.9} = {1.8}\n"
               "\n"
               "Infinite case:\n"
               "For f(x) = sup_{i∈I} f_i(x) with compact I:\n"
               "Similar rule with limiting process\n"
               "\n"
               "Applications:\n"
               "- Minimax problems\n"
               "- Support vector machines\n"
               "- Robust optimization";
    }

    /**
     * @brief Approximate minimization rule
     */
    static std::string approximateMinimization() {
        return "Approximate Minimization Rule (Ekeland-type):\n"
               "\n"
               "For f: X → ℝ ∪ {+∞} lsc, bounded below:\n"
               "\n"
               "∀ε > 0, ∀x̄ with f(x̄) ≤ inf f + ε:\n"
               "∃x_ε with ‖x_ε - x̄‖ ≤ √ε:\n"
               "  0 ∈ ∂_{√ε}f(x_ε)\n"
               "\n"
               "Interpretation:\n"
               "- Approximate minimizers have small subdifferentials\n"
               "- Foundation of Ekeland's variational principle\n"
               "- Quantifies \"nearly critical\"\n"
               "\n"
               "Stronger version (Ekeland):\n"
               "∃x_ε with:\n"
               "  f(x_ε) ≤ f(x̄)\n"
               "  ‖x_ε - x̄‖ ≤ √ε\n"
               "  f(y) ≥ f(x_ε) - √ε‖y - x_ε‖, ∀y\n"
               "\n"
               "Consequences:\n"
               "1. ∂f(x_ε) ∩ √εB* ≠ ∅\n"
               "2. Existence of approximate critical points\n"
               "3. Basis for descent methods\n"
               "\n"
               "Algorithm implications:\n"
               "- Stopping criteria: ‖∂f(x)‖ ≤ ε\n"
               "- Convergence: f(x_k) ↓ inf f ⇒ dist(0, ∂f(x_k)) → 0\n"
               "- Practical implementable conditions";
    }
};

/**
 * @class ExactCalculus
 * @brief Exact calculus rules with constraint qualifications
 */
class ExactCalculus {
public:
    /**
     * @brief Exact sum rule with qualification
     */
    static std::string exactSum() {
        return "Exact Sum Rule:\n"
               "\n"
               "For f, g: X → ℝ ∪ {+∞} lsc proper, x̄ ∈ dom f ∩ dom g:\n"
               "\n"
               "Constraint qualification:\n"
               "  0 ∈ int(dom f - dom g)\n"
               "or equivalently:\n"
               "  ∃r > 0: rB ⊆ dom f - dom g\n"
               "\n"
               "Then: ∂(f + g)(x̄) = ∂f(x̄) + ∂g(x̄)\n"
               "\n"
               "Weaker qualifications (one suffices):\n"
               "1. Attouch-Brézis: One function continuous at x̄\n"
               "2. Relative interior: ri(dom f) ∩ ri(dom g) ≠ ∅\n"
               "3. Finite-dimensional: X = ℝⁿ, both closed\n"
               "\n"
               "Counterexample (equality fails):\n"
               "  f(x) = I_{(-∞,0]}(x), g(x) = I_{[0,∞)}(x) on ℝ\n"
               "  dom f ∩ dom g = {0}\n"
               "  ∂f(0) = (-∞,0], ∂g(0) = [0,∞)\n"
               "  ∂(f+g)(0) = ℝ but ∂f(0) + ∂g(0) = ℝ\n"
               "  (Actually equal here, but pathological)\n"
               "\n"
               "Practical check:\n"
               "If one function is Lipschitz, sum rule exact\n"
               "\n"
               "Applications:\n"
               "- Optimality conditions: 0 ∈ ∂(f+g)(x*)\n"
               "- Proximal splitting: separate f and g\n"
               "- Duality theory";
    }

    /**
     * @brief Exact chain rule
     */
    static std::string exactChain() {
        return "Exact Chain Rule:\n"
               "\n"
               "Case 1: Smooth outer function\n"
               "If f: ℝᵐ → ℝ differentiable, g: ℝⁿ → ℝᵐ:\n"
               "  ∂(f ∘ g)(x) ⊇ (∂g(x))*∇f(g(x))\n"
               "\n"
               "Exact when:\n"
               "- g subdifferentially regular at x\n"
               "- f continuously differentiable\n"
               "\n"
               "Case 2: Smooth inner function\n"
               "If f: ℝᵐ → ℝ convex, g: ℝⁿ → ℝᵐ C¹:\n"
               "  ∂(f ∘ g)(x) = Dg(x)*∂f(g(x))\n"
               "\n"
               "Qualification for equality:\n"
               "  g(x) ∈ int(dom f)\n"
               "or\n"
               "  Dg(x) surjective onto aff(dom f - g(x))\n"
               "\n"
               "Case 3: Linear inner function\n"
               "For linear A: X → Y, f: Y → ℝ lsc:\n"
               "  ∂(f ∘ A)(x) = A*∂f(Ax)\n"
               "\n"
               "Always exact! No qualification needed.\n"
               "\n"
               "Example:\n"
               "  f(y) = ‖y‖₁, A matrix, h(x) = ‖Ax‖₁\n"
               "  ∂h(x) = A*∂(‖·‖₁)(Ax)\n"
               "  Exact for all A, x\n"
               "\n"
               "Applications:\n"
               "- Least squares: f(x) = ‖Ax - b‖²\n"
               "- Compressed sensing: min ‖x‖₁ s.t. Ax = b\n"
               "- Regularization: f(x) + λ‖Dx‖";
    }

    /**
     * @brief Exact maximum rule
     */
    static std::string exactMaximum() {
        return "Exact Maximum Rule:\n"
               "\n"
               "For f(x) = max{f₁(x), ..., f_m(x)} where f_i lsc:\n"
               "\n"
               "Let I(x̄) = {i : f_i(x̄) = f(x̄)} (active index set)\n"
               "\n"
               "Then: ∂f(x̄) = conv{⋃_{i ∈ I(x̄)} ∂f_i(x̄)}\n"
               "\n"
               "Qualification:\n"
               "Functions f_i are subdifferentially regular at x̄\n"
               "\n"
               "Proof idea:\n"
               "1. Upper bound: ∂f ⊆ conv{⋃ ∂f_i}\n"
               "2. Lower bound: Each ∂f_i(x̄) ⊆ ∂f(x̄) for i ∈ I(x̄)\n"
               "3. Convexity: conv preserved\n"
               "\n"
               "Infinite index set:\n"
               "For f(x) = sup_{t∈T} f_t(x) with T compact:\n"
               "\n"
               "  ∂f(x̄) ⊇ conv{⋃_{t: f_t(x̄)=f(x̄)} ∂f_t(x̄)}\n"
               "\n"
               "Equality requires:\n"
               "- Uniform continuity in t\n"
               "- Compactness of active set\n"
               "\n"
               "Danskin's theorem (differentiable case):\n"
               "If f_t(x) differentiable in x, argmax unique:\n"
               "  ∇f(x̄) = ∇_x f_t(x̄) where t = argmax\n"
               "\n"
               "Example:\n"
               "  f(x) = max{x, -x, 0}\n"
               "  At x = 0: all three active\n"
               "  ∂f(0) = conv{1, -1, 0} = [-1, 1]\n"
               "\n"
               "Applications:\n"
               "- ℓ∞ norm: ‖x‖∞ = max_i |x_i|\n"
               "- ReLU: max{0, x}\n"
               "- Minimax optimization";
    }
};

/**
 * @class SubdifferentialMVT
 * @brief Mean value theorems for subdifferentials
 */
class SubdifferentialMVT {
public:
    /**
     * @brief Convex mean value theorem
     */
    static std::string convexMVT() {
        return "Mean Value Theorem for Convex Functions:\n"
               "\n"
               "Let f: X → ℝ convex, x, y ∈ X.\n"
               "\n"
               "Then ∃λ ∈ (0, 1), v* ∈ ∂f(z) where z = λx + (1-λ)y:\n"
               "\n"
               "  f(y) - f(x) = ⟨v*, y - x⟩\n"
               "\n"
               "Interpretation:\n"
               "- Subgradient exactly captures slope\n"
               "- Point z on line segment [x, y]\n"
               "- Generalizes smooth MVT: f(y) - f(x) = f'(z)(y - x)\n"
               "\n"
               "Proof sketch:\n"
               "Consider h(t) = f(x + t(y-x)) for t ∈ [0,1]\n"
               "By convexity: h'₊(t) exists\n"
               "h(1) - h(0) = ∫₀¹ h'₊(t)dt\n"
               "Intermediate value gives λ with h'₊(λ) = h(1) - h(0)\n"
               "\n"
               "Consequence (subdifferential monotonicity):\n"
               "For x ≠ y, v* ∈ ∂f(x), w* ∈ ∂f(y):\n"
               "  ⟨v* - w*, x - y⟩ ≥ 0\n"
               "\n"
               "Strong convexity:\n"
               "If f strongly convex with modulus μ:\n"
               "  ⟨v* - w*, x - y⟩ ≥ μ‖x - y‖²\n"
               "\n"
               "Applications:\n"
               "- Convergence of gradient descent\n"
               "- Monotone operator theory\n"
               "- Variational inequalities";
    }

    /**
     * @brief Non-smooth mean value theorem
     */
    static std::string nonsmoothMVT() {
        return "Mean Value Theorem for Non-Smooth Functions:\n"
               "\n"
               "Lebourg's MVT:\n"
               "Let f: ℝⁿ → ℝ locally Lipschitz. For x, y ∈ ℝⁿ:\n"
               "\n"
               "∃λ ∈ (0, 1), z = x + λ(y - x), v* ∈ ∂_Cf(z):\n"
               "  f(y) - f(x) = ⟨v*, y - x⟩\n"
               "\n"
               "where ∂_Cf is Clarke subdifferential.\n"
               "\n"
               "Key properties:\n"
               "1. No convexity assumption\n"
               "2. Requires local Lipschitz continuity\n"
               "3. Uses Clarke (not Fréchet) subdifferential\n"
               "4. Point z in interior of segment\n"
               "\n"
               "Proof uses:\n"
               "- Rademacher's theorem (a.e. differentiability)\n"
               "- Generalized directional derivatives\n"
               "- Clarke's calculus\n"
               "\n"
               "Example:\n"
               "  f(x) = |x| on ℝ\n"
               "  f(1) - f(-1) = 2\n"
               "  At z = 0: ∂_Cf(0) = [-1, 1]\n"
               "  Can choose v* = 1, gives ⟨1, 2⟩ = 2 ✓\n"
               "\n"
               "Fuzzy version:\n"
               "For lsc f (not necessarily Lipschitz):\n"
               "∀ε > 0, ∃z ∈ [x, y], v* ∈ ∂f(z):\n"
               "  |f(y) - f(x) - ⟨v*, y - x⟩| ≤ ε‖y - x‖\n"
               "\n"
               "Applications:\n"
               "- Error bounds\n"
               "- Descent lemmas\n"
               "- Convergence analysis";
    }

    /**
     * @brief Approximate mean value theorem
     */
    static std::string approximateMVT() {
        return "Approximate Mean Value Inequality:\n"
               "\n"
               "For f: X → ℝ lsc, x, y ∈ X, ε > 0:\n"
               "\n"
               "∃z ∈ [x, y], v* ∈ ∂_εf(z):\n"
               "  f(y) - f(x) ≤ ⟨v*, y - x⟩ + 2ε‖y - x‖\n"
               "\n"
               "Interpretation:\n"
               "- One-sided inequality (upper bound)\n"
               "- Uses ε-subdifferential\n"
               "- Error term 2ε‖y - x‖\n"
               "- Valid for any lsc function\n"
               "\n"
               "Proof idea:\n"
               "Apply Ekeland's principle to:\n"
               "  φ(w) = f(w) - ⟨v₀*, w⟩\n"
               "along segment [x, y]\n"
               "\n"
               "Refinement (both bounds):\n"
               "∃z₁, z₂ ∈ [x, y], v₁* ∈ ∂_εf(z₁), v₂* ∈ ∂_εf(z₂):\n"
               "  ⟨v₁*, y - x⟩ - ε‖y - x‖ ≤ f(y) - f(x)\n"
               "                           ≤ ⟨v₂*, y - x⟩ + ε‖y - x‖\n"
               "\n"
               "Letting ε → 0:\n"
               "Recovers classical MVT when limits exist\n"
               "\n"
               "Application to descent:\n"
               "If ‖∂_εf(x_k)‖ → 0 and f(x_k) bounded:\n"
               "  f(x_{k+1}) ≈ f(x_k) (by MVT with small ∂_εf)\n"
               "  ⇒ sequence converges to stationary point\n"
               "\n"
               "Algorithmic use:\n"
               "Stopping criterion: dist(0, ∂_εf(x)) < tol\n"
               "Justified by approximate MVT";
    }
};

/**
 * @class SmoothnessNorms
 * @brief Smoothness properties of norms in Banach spaces
 */
class SmoothnessNorms {
public:
    /**
     * @brief Smooth and strictly convex norms
     */
    static std::string smoothNorms() {
        return "Smoothness of Norms:\n"
               "\n"
               "A norm ‖·‖ on Banach space X is:\n"
               "\n"
               "1. Gâteaux differentiable at x ≠ 0 if:\n"
               "   ∃!v* ∈ ∂‖x‖ = {f ∈ X* : ‖f‖ ≤ 1, f(x) = ‖x‖}\n"
               "   Equivalently: Support functional unique\n"
               "\n"
               "2. Fréchet differentiable at x ≠ 0 if:\n"
               "   ‖x + h‖ = ‖x‖ + ⟨D‖·‖(x), h⟩ + o(‖h‖)\n"
               "   Uniformly in directions\n"
               "\n"
               "3. Uniformly smooth if:\n"
               "   ρ(τ)/τ → 0 as τ → 0\n"
               "   where ρ(τ) = sup{(‖x+y‖ + ‖x-y‖)/2 - 1 : ‖x‖=1, ‖y‖=τ}\n"
               "\n"
               "Examples:\n"
               "\n"
               "ℓᵖ norms (1 < p < ∞):\n"
               "  ‖x‖_p = (Σ|x_i|ᵖ)^(1/p)\n"
               "  - Uniformly smooth and strictly convex\n"
               "  - Fréchet differentiable away from 0\n"
               "  - ∂‖x‖_p = {v : v_i = |x_i|^(p-1)sgn(x_i)/‖x‖_p^(p-1)}\n"
               "\n"
               "ℓ² (Hilbert):\n"
               "  ‖x‖₂ = √(Σx_i²)\n"
               "  - Optimally smooth\n"
               "  - ∇‖x‖₂ = x/‖x‖₂\n"
               "  - Modulus ρ(τ) = √(1 + τ²) - 1 ≈ τ²/2\n"
               "\n"
               "ℓ¹ and ℓ∞:\n"
               "  - NOT differentiable (many subdifferentials)\n"
               "  - ∂‖x‖₁ = {v : ‖v‖∞ ≤ 1, v_i = sgn(x_i) if x_i ≠ 0}\n"
               "  - ∂‖x‖∞ multi-valued at non-smooth points\n"
               "\n"
               "Duality:\n"
               "  X smooth ⟺ X* strictly convex\n"
               "  X uniformly smooth ⟺ X* uniformly convex\n"
               "\n"
               "Consequences:\n"
               "- Smooth norms: unique subdifferential\n"
               "- Projection onto balls well-defined\n"
               "- Optimization algorithms simpler";
    }

    /**
     * @brief Favorable classes of Banach spaces
     */
    static std::string favorableSpaces() {
        return "Favorable Classes of Banach Spaces:\n"
               "\n"
               "1. ASPLUND SPACES\n"
               "   X is Asplund if every convex continuous f on open set\n"
               "   is Fréchet differentiable on dense G_δ set.\n"
               "\n"
               "   Characterizations:\n"
               "   - X* separable ⟹ X Asplund\n"
               "   - X reflexive ⟹ X Asplund\n"
               "   - All ℓᵖ, Lᵖ (1 < p < ∞) are Asplund\n"
               "\n"
               "   Calculus advantage:\n"
               "   Fuzzy calculus rules become exact on dense sets\n"
               "\n"
               "2. HILBERT SPACES\n"
               "   Inner product structure: ⟨·,·⟩\n"
               "\n"
               "   Special properties:\n"
               "   - Norm: ‖x‖² = ⟨x, x⟩\n"
               "   - Riesz representation: X ≅ X*\n"
               "   - Projection theorem: unique nearest point\n"
               "   - Orthogonality: ⟨x, y⟩ = 0\n"
               "\n"
               "   Subdifferential simplification:\n"
               "   ∂f(x) identified with gradient ∇f(x)\n"
               "\n"
               "3. FINITE-DIMENSIONAL SPACES\n"
               "   X = ℝⁿ\n"
               "\n"
               "   Advantages:\n"
               "   - All norms equivalent\n"
               "   - Closed = compact (Heine-Borel)\n"
               "   - No constraint qualifications needed\n"
               "   - Exact sum rule always\n"
               "   - Clarke = Fréchet subdifferential\n"
               "\n"
               "4. SMOOTH BANACH SPACES\n"
               "   Norm Fréchet differentiable away from 0\n"
               "\n"
               "   Examples:\n"
               "   - ℓᵖ, Lᵖ for 1 < p < ∞\n"
               "   - All Hilbert spaces\n"
               "   - Sobolev spaces W^{k,p} (1 < p < ∞)\n"
               "\n"
               "   Benefit:\n"
               "   Duality mapping single-valued\n"
               "   Facilitates optimization algorithms\n"
               "\n"
               "5. REFLEXIVE SPACES\n"
               "   X = X** (canonical isomorphism)\n"
               "\n"
               "   Examples:\n"
               "   - ℓᵖ, Lᵖ for 1 < p < ∞\n"
               "   - All Hilbert spaces\n"
               "\n"
               "   Advantages:\n"
               "   - Weak compactness: bounded ⇒ weakly compact\n"
               "   - Existence of minimizers\n"
               "   - Strong duality in optimization\n"
               "\n"
               "Hierarchy:\n"
               "  Hilbert ⊂ Uniformly smooth & convex\n"
               "         ⊂ Reflexive\n"
               "         ⊂ Asplund\n"
               "\n"
               "Non-favorable:\n"
               "- ℓ¹, ℓ∞: Not reflexive\n"
               "- C[0,1]: Not separable dual\n"
               "- L¹: Not reflexive";
    }
};

} // namespace maths::analysis

#endif // MATHS_ANALYSIS_CONVEX_CALCULUS_HPP
