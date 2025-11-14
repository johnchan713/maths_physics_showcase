#ifndef MATHS_ANALYSIS_ELEMENTARY_SUBDIFFERENTIALS_HPP
#define MATHS_ANALYSIS_ELEMENTARY_SUBDIFFERENTIALS_HPP

#include <string>
#include <functional>
#include <vector>
#include <cmath>

/**
 * @file elementary_subdifferentials.hpp
 * @brief Elementary, proximal, and viscosity subdifferentials
 *
 * Implements:
 * - Elementary subderivatives and subdifferentials
 * - Proximal and Fréchet subdifferentials
 * - Coderivatives of set-valued mappings
 * - Viscosity subdifferentials for Hamilton-Jacobi equations
 * - Applications to optimization and control
 */

namespace maths::analysis {

/**
 * @class ElementarySubdifferentials
 * @brief Elementary (Fréchet) subdifferentials and subderivatives
 */
class ElementarySubdifferentials {
public:
    /**
     * @brief Fréchet subdifferential definition
     */
    static std::string frechetSubdifferential() {
        return "Fréchet (Elementary) Subdifferential:\n"
               "\n"
               "For f: X → ℝ ∪ {+∞}, x̄ ∈ dom f:\n"
               "\n"
               "v* ∈ ∂_Ff(x̄) (Fréchet subdifferential) if:\n"
               "\n"
               "  lim inf_{x → x̄} [f(x) - f(x̄) - ⟨v*, x - x̄⟩] / ‖x - x̄‖ ≥ 0\n"
               "\n"
               "Equivalently:\n"
               "  f(x) ≥ f(x̄) + ⟨v*, x - x̄⟩ + o(‖x - x̄‖)\n"
               "\n"
               "as x → x̄.\n"
               "\n"
               "Interpretation:\n"
               "- v* is \"tangent\" functional at x̄\n"
               "- Linear underestimate with o(‖x - x̄‖) error\n"
               "- Strict (not just lim sup)\n"
               "- Localized to neighborhood of x̄\n"
               "\n"
               "Properties:\n"
               "1. ∂_Ff(x̄) is weak* closed convex (possibly empty)\n"
               "2. f convex ⇒ ∂_Ff = ∂f (classical subdifferential)\n"
               "3. f C¹ at x̄ ⇒ ∂_Ff(x̄) = {∇f(x̄)}\n"
               "4. Smaller than limiting/Clarke subdifferentials\n"
               "\n"
               "Example:\n"
               "  f(x) = -√|x| at x = 0 on ℝ\n"
               "  ∂_Ff(0) = ∅ (no Fréchet subgradients)\n"
               "  But ∂_Cf(0) = [-∞, +∞] (Clarke)\n"
               "\n"
               "Optimality:\n"
               "  x̄ local minimizer ⇒ 0 ∈ ∂_Ff(x̄)\n"
               "  (First-order necessary condition)\n"
               "\n"
               "Advantages:\n"
               "- Good calculus in smooth spaces\n"
               "- Respects local geometry\n"
               "- Natural for differentiable approximations\n"
               "\n"
               "Disadvantages:\n"
               "- May be empty (non-robust)\n"
               "- Calculus requires qualifications\n"
               "- Not closed under limits";
    }

    /**
     * @brief Proximal subdifferential
     */
    static std::string proximalSubdifferential() {
        return "Proximal Subdifferential:\n"
               "\n"
               "For f: X → ℝ ∪ {+∞}, x̄ ∈ dom f:\n"
               "\n"
               "v* ∈ ∂_Pf(x̄) (proximal subdifferential) if:\n"
               "\n"
               "∃σ > 0, η > 0:\n"
               "  f(x) ≥ f(x̄) + ⟨v*, x - x̄⟩ - σ‖x - x̄‖²\n"
               "\n"
               "for all x with ‖x - x̄‖ < η.\n"
               "\n"
               "Interpretation:\n"
               "- Parabolic (quadratic) lower bound\n"
               "- v* is \"approximately gradient\" with curvature -σ\n"
               "- Stronger than Fréchet: σ‖x - x̄‖² vs o(‖x - x̄‖)\n"
               "\n"
               "Relationship:\n"
               "  ∂_Pf(x̄) ⊆ ∂_Ff(x̄) ⊆ ∂_Cf(x̄)\n"
               "\n"
               "Geometric meaning:\n"
               "v* ∈ ∂_Pf(x̄) ⟺ x̄ local minimizer of\n"
               "  f(x) - ⟨v*, x⟩ + (σ/2)‖x - x̄‖²\n"
               "\n"
               "Connection to proximal mapping:\n"
               "  prox_{σf}(x̄ + σv*) = x̄\n"
               "\n"
               "Properties:\n"
               "1. Convex function: ∂_Pf = ∂f\n"
               "2. Hilbert space: simpler characterization\n"
               "3. Robust under perturbations\n"
               "4. Useful for perturbed optimization\n"
               "\n"
               "Example:\n"
               "  f(x) = |x| at x = 0 on ℝ\n"
               "  ∂_Pf(0) = {0} (parabola touches from below)\n"
               "  ∂_Ff(0) = ∅ (no o(|x|) bound)\n"
               "\n"
               "Applications:\n"
               "- Hamilton-Jacobi equations\n"
               "- Viscosity solutions\n"
               "- Perturbed problems\n"
               "- Proximal point algorithms\n"
               "\n"
               "Hilbert space formula:\n"
               "v* ∈ ∂_Pf(x̄) ⟺ ∃σ, r > 0:\n"
               "  x̄ = argmin_{‖x-x̄‖<r} [f(x) - ⟨v*, x⟩ + (σ/2)‖x - x̄‖²]";
    }

    /**
     * @brief Limiting subdifferential
     */
    static std::string limitingSubdifferential() {
        return "Limiting (Mordukhovich) Subdifferential:\n"
               "\n"
               "For f: X → ℝ ∪ {+∞}, x̄ ∈ dom f:\n"
               "\n"
               "v* ∈ ∂_Lf(x̄) if ∃sequences:\n"
               "  x_k → x̄ with f(x_k) → f(x̄)\n"
               "  v_k* → v* with v_k* ∈ ∂_Ff(x_k)\n"
               "\n"
               "Notation: Also denoted ∂f, ∂̂f, ∂̄f in literature\n"
               "\n"
               "Interpretation:\n"
               "- Sequential closure of Fréchet subdifferential\n"
               "- Captures limiting behavior\n"
               "- Robust under approximation\n"
               "- Never empty for lsc functions (in Asplund spaces)\n"
               "\n"
               "Properties:\n"
               "1. ∂_Ff(x̄) ⊆ ∂_Lf(x̄) always\n"
               "2. Convex f: ∂_Lf = ∂f (convex subdifferential)\n"
               "3. Closed graph: graph(∂_Lf) closed in X × X*\n"
               "4. Calculus rules (sum, chain) available\n"
               "5. Non-empty for lsc f in Asplund spaces\n"
               "\n"
               "Example:\n"
               "  f(x) = -|x|^(1/2) at x = 0\n"
               "  ∂_Ff(0) = ∅\n"
               "  ∂_Lf(0) = ℝ (all limits of gradients as x → 0)\n"
               "\n"
               "Relationship to Clarke:\n"
               "  ∂_Lf(x̄) ⊆ ∂_Cf(x̄)\n"
               "Equality for locally Lipschitz f\n"
               "\n"
               "Optimality conditions:\n"
               "x̄ local min ⇒ 0 ∈ ∂_Lf(x̄)\n"
               "Converse true under regularity\n"
               "\n"
               "Calculus (fuzzy sum rule):\n"
               "For lsc f, g in Asplund space:\n"
               "  ∂_L(f+g)(x̄) ⊆ lim sup [∂_Lf(x) + ∂_Lg(y)]\n"
               "                x→x̄, y→x̄\n"
               "Exact under qualification\n"
               "\n"
               "Applications:\n"
               "- Generalized equations\n"
               "- Variational analysis\n"
               "- Non-smooth mechanics\n"
               "- Optimal control\n"
               "\n"
               "Advantages:\n"
               "+ Robust (closed, never empty)\n"
               "+ Good calculus\n"
               "+ Optimality conditions\n"
               "- May be large (many elements)";
    }
};

/**
 * @class Coderivatives
 * @brief Coderivatives of set-valued mappings
 */
class Coderivatives {
public:
    /**
     * @brief Fréchet coderivative
     */
    static std::string frechetCoderivative() {
        return "Fréchet Coderivative:\n"
               "\n"
               "For set-valued map F: X ⇉ Y, (x̄, ȳ) ∈ graph(F):\n"
               "\n"
               "The Fréchet coderivative D*F(x̄, ȳ): Y* → X* is:\n"
               "\n"
               "  x* ∈ D*F(x̄, ȳ)(y*) ⟺\n"
               "  (x*, -y*) ∈ N̂_F(x̄, ȳ)\n"
               "\n"
               "where N̂_F is Fréchet normal cone to graph(F):\n"
               "\n"
               "  (x*, y*) ∈ N̂((x̄, ȳ), graph(F)) if:\n"
               "  lim sup_{(x,y)→(x̄,ȳ), (x,y)∈graph(F)}\n"
               "    [⟨x*, x-x̄⟩ + ⟨y*, y-ȳ⟩] / ‖(x-x̄, y-ȳ)‖ ≤ 0\n"
               "\n"
               "Single-valued case:\n"
               "If F(x) = {f(x)} (function), then:\n"
               "  D*f(x̄)(y*) = {Df(x̄)*y*}\n"
               "Adjoint of derivative!\n"
               "\n"
               "Properties:\n"
               "1. D*F(x̄, ȳ) is positively homogeneous\n"
               "2. Closed graph\n"
               "3. D*F(x̄, ȳ)(0) = N̂_F(x̄, ȳ) ∩ (X* × {0})\n"
               "\n"
               "Chain rule (composition):\n"
               "For F: X ⇉ Y, G: Y ⇉ Z, (x̄, ȳ, z̄):\n"
               "  D*(G ∘ F)(x̄, z̄)(z*) ⊇ D*F(x̄, ȳ)(D*G(ȳ, z̄)(z*))\n"
               "\n"
               "Sum rule:\n"
               "  D*(F + G)(x̄, ȳ+z̄)(w*) ⊇\n"
               "    D*F(x̄, ȳ)(w*) + D*G(x̄, z̄)(w*)\n"
               "\n"
               "Example (normal cone map):\n"
               "F(x) = N_C(x) for convex set C\n"
               "D*F(x̄, v̄)(v*) = ∂σ_C(v̄) where σ_C is support function\n"
               "\n"
               "Applications:\n"
               "- Sensitivity analysis\n"
               "- Generalized derivatives\n"
               "- Implicit function theorems\n"
               "- Optimality conditions";
    }

    /**
     * @brief Limiting coderivative
     */
    static std::string limitingCoderivative() {
        return "Limiting (Mordukhovich) Coderivative:\n"
               "\n"
               "For F: X ⇉ Y, (x̄, ȳ) ∈ graph(F):\n"
               "\n"
               "  x* ∈ D*_LF(x̄, ȳ)(y*) if:\n"
               "  ∃(x_k, y_k) → (x̄, ȳ), y_k ∈ F(x_k)\n"
               "  ∃x_k* → x*, y_k* → y*\n"
               "  with x_k* ∈ D*_FF(x_k, y_k)(y_k*)\n"
               "\n"
               "Sequential closure of Fréchet coderivative\n"
               "\n"
               "Properties:\n"
               "1. Closed graph\n"
               "2. Positively homogeneous\n"
               "3. Upper semicontinuous (outer limit)\n"
               "4. D*_FF ⊆ D*_LF\n"
               "\n"
               "Calculus (chain rule):\n"
               "Under metric regularity:\n"
               "  D*_L(G ∘ F)(x̄, z̄) = D*_LF(x̄, ȳ) ∘ D*_LG(ȳ, z̄)\n"
               "\n"
               "Criterion for metric regularity:\n"
               "F metrically regular at (x̄, ȳ) ⟺\n"
               "  D*_LF(x̄, ȳ)(0) = {0}\n"
               "\n"
               "Criterion for pseudo-Lipschitz:\n"
               "F pseudo-Lipschitz at (x̄, ȳ) ⟺\n"
               "  y* = 0 ⇒ D*_LF⁻¹(ȳ, x̄)(y*) = {0}\n"
               "\n"
               "Example (subdifferential map):\n"
               "F(x) = ∂f(x) for lsc f\n"
               "D*_LF(x̄, v̄) relates to second-order subdifferential\n"
               "\n"
               "Optimality conditions:\n"
               "For constrained problem min f(x) s.t. g(x) ∈ C:\n"
               "  0 ∈ ∂_Lf(x̄) + D*_Lg(x̄, ḡ)(N_C(ḡ))\n"
               "\n"
               "Applications:\n"
               "- Stability and sensitivity\n"
               "- Implicit functions\n"
               "- Constraint qualifications\n"
               "- Variational inequalities";
    }

    /**
     * @brief Normal cone and tangent cone
     */
    static std::string normalTangentCones() {
        return "Normal and Tangent Cones:\n"
               "\n"
               "For set C ⊆ X, x̄ ∈ C:\n"
               "\n"
               "FRÉCHET NORMAL CONE:\n"
               "  N̂_C(x̄) = {v* : lim sup_{x→x̄, x∈C} ⟨v*, x-x̄⟩/‖x-x̄‖ ≤ 0}\n"
               "\n"
               "Interpretation: v* \"outward normal\" to C at x̄\n"
               "\n"
               "LIMITING NORMAL CONE:\n"
               "  N_C(x̄) = lim sup_{x→x̄, x∈C} N̂_C(x)\n"
               "\n"
               "Sequential closure of Fréchet normals\n"
               "\n"
               "TANGENT CONE (Bouligand):\n"
               "  T_C(x̄) = {v : lim inf_{t↓0} d(x̄+tv, C)/t = 0}\n"
               "\n"
               "Interpretation: feasible directions from x̄\n"
               "\n"
               "POLAR RELATIONSHIP:\n"
               "  N̂_C(x̄) = (T_C(x̄))⁻ (negative polar)\n"
               "  N_C(x̄) ⊆ (T_C(x̄))⁻\n"
               "\n"
               "Properties:\n"
               "1. Convex C: N_C(x̄) = normal cone (convex analysis)\n"
               "2. Smooth boundary: N_C(x̄) = ℝ₊·∇g(x̄) if C = {g ≤ 0}\n"
               "3. Product: N_{C₁×C₂}(x̄₁, x̄₂) = N_{C₁}(x̄₁) × N_{C₂}(x̄₂)\n"
               "\n"
               "Example (half-space):\n"
               "  C = {x : ⟨a, x⟩ ≤ b}, x̄ on boundary\n"
               "  N_C(x̄) = ℝ₊·a (ray in direction a)\n"
               "  T_C(x̄) = {v : ⟨a, v⟩ ≤ 0}\n"
               "\n"
               "Subdifferential as normal cone:\n"
               "  ∂f(x̄) = {v* : (v*, -1) ∈ N_{epi f}(x̄, f(x̄))}\n"
               "\n"
               "Calculus:\n"
               "  N_{C₁∩C₂}(x̄) ⊇ N_{C₁}(x̄) + N_{C₂}(x̄)\n"
               "Equality under qualification\n"
               "\n"
               "Applications:\n"
               "- Constrained optimization\n"
               "- Variational inequalities\n"
               "- Equilibrium problems\n"
               "- Differential inclusions";
    }
};

/**
 * @class ViscositySubdifferentials
 * @brief Viscosity subdifferentials for PDEs
 */
class ViscositySubdifferentials {
public:
    /**
     * @brief Viscosity subdifferential definition
     */
    static std::string definition() {
        return "Viscosity Subdifferential:\n"
               "\n"
               "For f: X → ℝ, x̄ ∈ X (typically X = ℝⁿ):\n"
               "\n"
               "VISCOSITY SUBDIFFERENTIAL ∂⁻f(x̄):\n"
               "v* ∈ ∂⁻f(x̄) if ∃φ ∈ C¹(X) with:\n"
               "  1. φ(x̄) = f(x̄)\n"
               "  2. φ(x) ≤ f(x) for x near x̄\n"
               "  3. ∇φ(x̄) = v*\n"
               "\n"
               "φ is called \"test function\" or \"smooth minorant\"\n"
               "\n"
               "VISCOSITY SUPERDIFFERENTIAL ∂⁺f(x̄):\n"
               "v* ∈ ∂⁺f(x̄) if -v* ∈ ∂⁻(-f)(x̄)\n"
               "Equivalently: smooth majorant with gradient v*\n"
               "\n"
               "Interpretation:\n"
               "- v* = gradient of smooth function touching from below/above\n"
               "- \"Touched by smooth function\"\n"
               "- Generalizes classical derivative\n"
               "\n"
               "Properties:\n"
               "1. f C¹ at x̄ ⇒ ∂⁻f(x̄) = ∂⁺f(x̄) = {∇f(x̄)}\n"
               "2. ∂⁻f(x̄) = ∂_Pf(x̄) (proximal subdifferential)\n"
               "3. Non-empty for lsc f at local minima\n"
               "4. Closed under + λI for λ > 0\n"
               "\n"
               "Example:\n"
               "  f(x) = |x| at x = 0\n"
               "  ∂⁻f(0) = {0} (parabola φ(x) = x²/2ε touches)\n"
               "  ∂⁺f(0) = [-1, 1] (all tangent lines)\n"
               "\n"
               "Relationship to other concepts:\n"
               "  ∂⁻f ⊆ ∂_Ff ⊆ ∂_Lf ⊆ ∂_Cf\n"
               "  ∂⁺f ⊇ ∂_Cf (Clarke)\n"
               "\n"
               "Viscosity solutions:\n"
               "u is viscosity subsolution of H(x, u, Du) = 0 if:\n"
               "  ∀x, ∀p ∈ ∂⁺u(x): H(x, u(x), p) ≤ 0\n"
               "\n"
               "Applications:\n"
               "- Hamilton-Jacobi-Bellman equations\n"
               "- Optimal control (value functions)\n"
               "- Differential games\n"
               "- Level set methods";
    }

    /**
     * @brief Viscosity solutions of PDEs
     */
    static std::string viscositySolutions() {
        return "Viscosity Solutions of PDEs:\n"
               "\n"
               "Consider first-order PDE:\n"
               "  H(x, u(x), Du(x)) = 0 in Ω ⊂ ℝⁿ\n"
               "\n"
               "Classical solution requires u ∈ C¹.\n"
               "Viscosity solution: generalize to continuous u.\n"
               "\n"
               "DEFINITION:\n"
               "u: Ω → ℝ continuous is:\n"
               "\n"
               "Viscosity subsolution if ∀φ ∈ C¹:\n"
               "  u - φ has local max at x̄ ⇒ H(x̄, u(x̄), ∇φ(x̄)) ≤ 0\n"
               "\n"
               "Viscosity supersolution if ∀φ ∈ C¹:\n"
               "  u - φ has local min at x̄ ⇒ H(x̄, u(x̄), ∇φ(x̄)) ≥ 0\n"
               "\n"
               "Viscosity solution: both sub and supersolution\n"
               "\n"
               "Equivalent characterization:\n"
               "- Subsolution: H(x, u(x), p) ≤ 0 ∀p ∈ ∂⁺u(x)\n"
               "- Supersolution: H(x, u(x), p) ≥ 0 ∀p ∈ ∂⁻u(x)\n"
               "\n"
               "Key properties:\n"
               "1. Stability: limits of viscosity solutions are solutions\n"
               "2. Uniqueness: under comparison principle\n"
               "3. Existence: Perron's method\n"
               "4. C¹ solutions are viscosity solutions\n"
               "\n"
               "Example (eikonal equation):\n"
               "  |Du| = 1 in Ω, u = 0 on ∂Ω\n"
               "  u(x) = dist(x, ∂Ω) is viscosity solution\n"
               "  Not C¹ at medial axis!\n"
               "\n"
               "Hamilton-Jacobi-Bellman:\n"
               "  ∂_t u + H(x, Du) = 0\n"
               "  u(0, x) = u₀(x)\n"
               "\n"
               "Arises in optimal control:\n"
               "  u(t, x) = value function\n"
               "  H = Hamiltonian\n"
               "\n"
               "Comparison principle:\n"
               "If H proper, u subsolution, v supersolution:\n"
               "  u ≤ v on ∂Ω ⇒ u ≤ v in Ω\n"
               "\n"
               "Perron's method (existence):\n"
               "  u*(x) = sup{v(x) : v viscosity subsolution, v ≤ g on ∂Ω}\n"
               "  u* is viscosity solution\n"
               "\n"
               "Numerical methods:\n"
               "- Monotone schemes converge to viscosity solution\n"
               "- Semi-Lagrangian schemes\n"
               "- Level set methods\n"
               "\n"
               "Applications:\n"
               "- Optimal control\n"
               "- Differential games\n"
               "- Image processing\n"
               "- Front propagation";
    }
};

} // namespace maths::analysis

#endif // MATHS_ANALYSIS_ELEMENTARY_SUBDIFFERENTIALS_HPP
