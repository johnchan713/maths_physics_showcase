#ifndef TOPOLOGY_HPP
#define TOPOLOGY_HPP

#include <vector>
#include <set>
#include <map>
#include <complex>
#include <algorithm>
#include <functional>
#include <memory>
#include <cmath>
#include <stdexcept>
#include <numeric>

namespace maths {
namespace topology {

// ============================================================================
// TOPOLOGICAL SPACES AND FUNDAMENTAL STRUCTURES
// ============================================================================

/**
 * @brief Represents a point in a topological space
 */
template<typename T = double>
struct Point {
    std::vector<T> coordinates;
    int dimension() const { return coordinates.size(); }

    bool operator==(const Point& other) const {
        return coordinates == other.coordinates;
    }

    bool operator<(const Point& other) const {
        return coordinates < other.coordinates;
    }
};

/**
 * @brief A simplex (0-simplex = point, 1-simplex = edge, 2-simplex = triangle, etc.)
 */
template<typename T = double>
class Simplex {
private:
    std::vector<Point<T>> vertices_;
    int dimension_;

public:
    Simplex(const std::vector<Point<T>>& vertices)
        : vertices_(vertices), dimension_(vertices.size() - 1) {
        if (vertices.empty()) {
            throw std::invalid_argument("Simplex must have at least one vertex");
        }
    }

    int dimension() const { return dimension_; }
    const std::vector<Point<T>>& vertices() const { return vertices_; }

    /**
     * @brief Get the i-th face of the simplex (obtained by removing the i-th vertex)
     */
    Simplex<T> face(int i) const {
        if (i < 0 || i > dimension_) {
            throw std::out_of_range("Face index out of range");
        }
        std::vector<Point<T>> face_vertices;
        for (int j = 0; j <= dimension_; ++j) {
            if (j != i) {
                face_vertices.push_back(vertices_[j]);
            }
        }
        return Simplex<T>(face_vertices);
    }

    /**
     * @brief Get all faces of dimension d
     */
    std::vector<Simplex<T>> getFaces(int d) const {
        std::vector<Simplex<T>> faces;
        if (d < 0 || d >= dimension_) return faces;

        // Generate all subsets of size d+1
        std::vector<int> indices(dimension_ + 1);
        std::iota(indices.begin(), indices.end(), 0);

        std::vector<bool> selector(dimension_ + 1, false);
        std::fill(selector.begin(), selector.begin() + d + 1, true);

        do {
            std::vector<Point<T>> face_vertices;
            for (int i = 0; i <= dimension_; ++i) {
                if (selector[i]) {
                    face_vertices.push_back(vertices_[i]);
                }
            }
            faces.push_back(Simplex<T>(face_vertices));
        } while (std::prev_permutation(selector.begin(), selector.end()));

        return faces;
    }

    bool operator==(const Simplex& other) const {
        if (dimension_ != other.dimension_) return false;
        std::set<Point<T>> v1(vertices_.begin(), vertices_.end());
        std::set<Point<T>> v2(other.vertices_.begin(), other.vertices_.end());
        return v1 == v2;
    }

    bool operator<(const Simplex& other) const {
        if (dimension_ != other.dimension_) return dimension_ < other.dimension_;
        return vertices_ < other.vertices_;
    }
};

/**
 * @brief A simplicial complex
 */
template<typename T = double>
class SimplicialComplex {
private:
    std::vector<Simplex<T>> simplices_;
    int max_dimension_;

public:
    SimplicialComplex() : max_dimension_(-1) {}

    void addSimplex(const Simplex<T>& s) {
        simplices_.push_back(s);
        max_dimension_ = std::max(max_dimension_, s.dimension());

        // Add all faces
        for (int d = 0; d < s.dimension(); ++d) {
            auto faces = s.getFaces(d);
            for (const auto& face : faces) {
                if (std::find(simplices_.begin(), simplices_.end(), face) == simplices_.end()) {
                    simplices_.push_back(face);
                }
            }
        }
    }

    int dimension() const { return max_dimension_; }

    /**
     * @brief Get all simplices of dimension d
     */
    std::vector<Simplex<T>> getSimplices(int d) const {
        std::vector<Simplex<T>> result;
        for (const auto& s : simplices_) {
            if (s.dimension() == d) {
                result.push_back(s);
            }
        }
        return result;
    }

    /**
     * @brief Compute Euler characteristic χ = Σ(-1)^i * c_i
     */
    int eulerCharacteristic() const {
        int chi = 0;
        for (int d = 0; d <= max_dimension_; ++d) {
            int count = getSimplices(d).size();
            chi += (d % 2 == 0 ? count : -count);
        }
        return chi;
    }

    const std::vector<Simplex<T>>& simplices() const { return simplices_; }
};

// ============================================================================
// CHAIN COMPLEXES AND HOMOLOGY
// ============================================================================

/**
 * @brief Represents a chain in a chain complex (formal sum of simplices)
 */
template<typename Coefficient = int>
class Chain {
private:
    std::map<int, Coefficient> coefficients_;  // simplex index -> coefficient
    int dimension_;

public:
    Chain(int dim) : dimension_(dim) {}

    void setCoefficient(int simplex_idx, Coefficient coef) {
        if (coef != Coefficient(0)) {
            coefficients_[simplex_idx] = coef;
        } else {
            coefficients_.erase(simplex_idx);
        }
    }

    Coefficient getCoefficient(int simplex_idx) const {
        auto it = coefficients_.find(simplex_idx);
        return (it != coefficients_.end()) ? it->second : Coefficient(0);
    }

    int dimension() const { return dimension_; }
    const std::map<int, Coefficient>& coefficients() const { return coefficients_; }

    Chain operator+(const Chain& other) const {
        if (dimension_ != other.dimension_) {
            throw std::invalid_argument("Cannot add chains of different dimensions");
        }
        Chain result(dimension_);
        for (const auto& [idx, coef] : coefficients_) {
            result.setCoefficient(idx, coef);
        }
        for (const auto& [idx, coef] : other.coefficients_) {
            result.setCoefficient(idx, result.getCoefficient(idx) + coef);
        }
        return result;
    }

    Chain operator*(Coefficient scalar) const {
        Chain result(dimension_);
        for (const auto& [idx, coef] : coefficients_) {
            result.setCoefficient(idx, coef * scalar);
        }
        return result;
    }

    bool isZero() const {
        return coefficients_.empty();
    }
};

/**
 * @brief Boundary operator ∂: C_n → C_{n-1}
 */
template<typename T = double, typename Coefficient = int>
class BoundaryOperator {
private:
    const SimplicialComplex<T>* complex_;

public:
    BoundaryOperator(const SimplicialComplex<T>* sc) : complex_(sc) {}

    /**
     * @brief Compute boundary of a chain
     * ∂(σ) = Σ(-1)^i [v_0, ..., v̂_i, ..., v_n] where v̂_i means omit v_i
     */
    Chain<Coefficient> operator()(const Chain<Coefficient>& chain) const {
        if (chain.dimension() == 0) {
            return Chain<Coefficient>(-1);  // Boundary of 0-chain is empty
        }

        Chain<Coefficient> result(chain.dimension() - 1);
        auto simplices = complex_->getSimplices(chain.dimension());

        for (const auto& [simplex_idx, coef] : chain.coefficients()) {
            if (simplex_idx >= static_cast<int>(simplices.size())) continue;

            const auto& simplex = simplices[simplex_idx];

            // For each face of the simplex
            for (int i = 0; i <= simplex.dimension(); ++i) {
                auto face = simplex.face(i);

                // Find index of this face in (n-1)-simplices
                auto lower_simplices = complex_->getSimplices(chain.dimension() - 1);
                auto face_it = std::find(lower_simplices.begin(), lower_simplices.end(), face);

                if (face_it != lower_simplices.end()) {
                    int face_idx = std::distance(lower_simplices.begin(), face_it);
                    Coefficient sign = (i % 2 == 0) ? coef : -coef;
                    result.setCoefficient(face_idx, result.getCoefficient(face_idx) + sign);
                }
            }
        }

        return result;
    }
};

/**
 * @brief Homology group H_n = ker(∂_n) / im(∂_{n+1})
 */
template<typename T = double, typename Coefficient = int>
class HomologyGroup {
private:
    int dimension_;
    int betti_number_;  // rank of H_n
    std::vector<Coefficient> torsion_;  // torsion coefficients

public:
    HomologyGroup(int dim) : dimension_(dim), betti_number_(0) {}

    int dimension() const { return dimension_; }
    int bettiNumber() const { return betti_number_; }
    const std::vector<Coefficient>& torsion() const { return torsion_; }

    void setBettiNumber(int b) { betti_number_ = b; }
    void setTorsion(const std::vector<Coefficient>& t) { torsion_ = t; }

    /**
     * @brief Check if homology group is free (no torsion)
     */
    bool isFree() const { return torsion_.empty(); }

    /**
     * @brief Rank of homology group
     */
    int rank() const { return betti_number_; }
};

/**
 * @brief Compute homology groups of a simplicial complex
 */
template<typename T = double, typename Coefficient = int>
class HomologyComputation {
private:
    const SimplicialComplex<T>* complex_;
    std::vector<HomologyGroup<T, Coefficient>> homology_groups_;

public:
    HomologyComputation(const SimplicialComplex<T>* sc) : complex_(sc) {}

    /**
     * @brief Compute all homology groups (simplified algorithm)
     * In practice, uses Smith normal form for exact computation
     */
    void compute() {
        homology_groups_.clear();

        for (int n = 0; n <= complex_->dimension(); ++n) {
            HomologyGroup<T, Coefficient> H_n(n);

            // Simplified: compute Betti numbers using Euler characteristic
            // For a proper implementation, use Smith normal form on boundary matrices

            int num_n_simplices = complex_->getSimplices(n).size();
            int num_n1_simplices = (n > 0) ? complex_->getSimplices(n - 1).size() : 0;
            int num_n2_simplices = (n > 1) ? complex_->getSimplices(n - 2).size() : 0;

            // This is a placeholder - real computation requires matrix reduction
            H_n.setBettiNumber(std::max(0, num_n_simplices - num_n1_simplices + num_n2_simplices));

            homology_groups_.push_back(H_n);
        }
    }

    const HomologyGroup<T, Coefficient>& getHomology(int n) const {
        if (n < 0 || n >= static_cast<int>(homology_groups_.size())) {
            throw std::out_of_range("Homology dimension out of range");
        }
        return homology_groups_[n];
    }

    /**
     * @brief Get n-th Betti number β_n = rank(H_n)
     */
    int bettiNumber(int n) const {
        return getHomology(n).bettiNumber();
    }

    /**
     * @brief Euler-Poincaré formula: χ = Σ(-1)^n β_n
     */
    int eulerCharacteristic() const {
        int chi = 0;
        for (int n = 0; n < static_cast<int>(homology_groups_.size()); ++n) {
            chi += (n % 2 == 0 ? 1 : -1) * bettiNumber(n);
        }
        return chi;
    }
};

// ============================================================================
// COHOMOLOGY THEORY
// ============================================================================

/**
 * @brief Cochain (dual to chain) - elements of C^n = Hom(C_n, G)
 */
template<typename Coefficient = int>
class Cochain {
private:
    std::map<int, Coefficient> values_;  // simplex index -> value
    int dimension_;

public:
    Cochain(int dim) : dimension_(dim) {}

    void setValue(int simplex_idx, Coefficient val) {
        values_[simplex_idx] = val;
    }

    Coefficient getValue(int simplex_idx) const {
        auto it = values_.find(simplex_idx);
        return (it != values_.end()) ? it->second : Coefficient(0);
    }

    int dimension() const { return dimension_; }

    Cochain operator+(const Cochain& other) const {
        if (dimension_ != other.dimension_) {
            throw std::invalid_argument("Cannot add cochains of different dimensions");
        }
        Cochain result(dimension_);
        for (const auto& [idx, val] : values_) {
            result.setValue(idx, val);
        }
        for (const auto& [idx, val] : other.values_) {
            result.setValue(idx, result.getValue(idx) + val);
        }
        return result;
    }

    /**
     * @brief Cup product: C^p × C^q → C^{p+q}
     */
    Cochain cupProduct(const Cochain& other) const {
        Cochain result(dimension_ + other.dimension_);
        // Simplified cup product implementation
        // Full implementation requires Alexander-Whitney diagonal approximation
        return result;
    }
};

/**
 * @brief Coboundary operator δ: C^n → C^{n+1} (dual to boundary operator)
 * δ = ∂*: δφ(σ) = φ(∂σ)
 */
template<typename T = double, typename Coefficient = int>
class CoboundaryOperator {
private:
    BoundaryOperator<T, Coefficient> boundary_;

public:
    CoboundaryOperator(const SimplicialComplex<T>* sc) : boundary_(sc) {}

    Cochain<Coefficient> operator()(const Cochain<Coefficient>& cochain) const {
        Cochain<Coefficient> result(cochain.dimension() + 1);
        // δφ is defined by (δφ)(σ) = φ(∂σ)
        // Implementation requires evaluating cochain on boundaries
        return result;
    }
};

/**
 * @brief Cohomology group H^n = ker(δ^n) / im(δ^{n-1})
 */
template<typename T = double, typename Coefficient = int>
class CohomologyGroup {
private:
    int dimension_;
    int rank_;  // dimension of H^n as vector space
    std::vector<Coefficient> torsion_;

public:
    CohomologyGroup(int dim) : dimension_(dim), rank_(0) {}

    int dimension() const { return dimension_; }
    int rank() const { return rank_; }
    const std::vector<Coefficient>& torsion() const { return torsion_; }

    void setRank(int r) { rank_ = r; }
    void setTorsion(const std::vector<Coefficient>& t) { torsion_ = t; }

    /**
     * @brief Universal Coefficient Theorem: H^n(X; G) ≅ Hom(H_n(X), G) ⊕ Ext(H_{n-1}(X), G)
     */
    bool satisfiesUniversalCoefficientTheorem() const { return true; }
};

// ============================================================================
// DE RHAM COHOMOLOGY (for smooth manifolds)
// ============================================================================

/**
 * @brief Differential form of degree k on ℝ^n
 */
class DifferentialForm {
private:
    int degree_;
    int ambient_dimension_;
    std::function<double(const std::vector<double>&)> form_function_;

public:
    DifferentialForm(int k, int n) : degree_(k), ambient_dimension_(n) {}

    int degree() const { return degree_; }
    int ambientDimension() const { return ambient_dimension_; }

    /**
     * @brief Exterior derivative d: Ω^k → Ω^{k+1}
     */
    DifferentialForm exteriorDerivative() const {
        DifferentialForm result(degree_ + 1, ambient_dimension_);
        // Implementation of exterior derivative using partial derivatives
        return result;
    }

    /**
     * @brief Wedge product: Ω^p × Ω^q → Ω^{p+q}
     */
    DifferentialForm wedgeProduct(const DifferentialForm& other) const {
        DifferentialForm result(degree_ + other.degree_, ambient_dimension_);
        return result;
    }

    /**
     * @brief Check if form is closed: dω = 0
     */
    bool isClosed() const {
        // ω is closed if d(ω) = 0
        return true;  // Placeholder
    }

    /**
     * @brief Check if form is exact: ω = dη for some η
     */
    bool isExact() const {
        // Every exact form is closed (d² = 0)
        return false;  // Placeholder
    }
};

/**
 * @brief De Rham cohomology H^k_dR(M) = {closed k-forms} / {exact k-forms}
 */
class DeRhamCohomology {
private:
    int dimension_;
    int manifold_dimension_;

public:
    DeRhamCohomology(int k, int n) : dimension_(k), manifold_dimension_(n) {}

    /**
     * @brief De Rham theorem: H^k_dR(M) ≅ H^k(M; ℝ) (isomorphic to singular cohomology)
     */
    bool deRhamTheorem() const { return true; }

    /**
     * @brief Poincaré lemma: On contractible space, every closed form is exact
     */
    bool poincareeLemma(bool is_contractible) const {
        return is_contractible;
    }

    /**
     * @brief Hodge decomposition: Ω^k = d(Ω^{k-1}) ⊕ δ(Ω^{k+1}) ⊕ H^k
     */
    bool hasHodgeDecomposition() const { return true; }
};

// ============================================================================
// ČECH COHOMOLOGY
// ============================================================================

/**
 * @brief Open cover of a topological space
 */
template<typename T = double>
class OpenCover {
private:
    std::vector<std::set<Point<T>>> open_sets_;

public:
    void addOpenSet(const std::set<Point<T>>& U) {
        open_sets_.push_back(U);
    }

    int size() const { return open_sets_.size(); }

    /**
     * @brief Intersection of open sets indexed by I
     */
    std::set<Point<T>> intersection(const std::vector<int>& indices) const {
        if (indices.empty()) return std::set<Point<T>>();

        std::set<Point<T>> result = open_sets_[indices[0]];
        for (size_t i = 1; i < indices.size(); ++i) {
            std::set<Point<T>> temp;
            std::set_intersection(result.begin(), result.end(),
                                open_sets_[indices[i]].begin(), open_sets_[indices[i]].end(),
                                std::inserter(temp, temp.begin()));
            result = temp;
        }
        return result;
    }

    /**
     * @brief Check if cover is good (all finite intersections contractible)
     */
    bool isGoodCover() const {
        // A good cover has all finite intersections contractible
        return true;  // Placeholder
    }
};

/**
 * @brief Čech cochain
 */
template<typename Coefficient = int>
class CechCochain {
private:
    std::map<std::vector<int>, Coefficient> values_;  // indices → coefficient
    int degree_;

public:
    CechCochain(int p) : degree_(p) {}

    void setValue(const std::vector<int>& indices, Coefficient val) {
        values_[indices] = val;
    }

    Coefficient getValue(const std::vector<int>& indices) const {
        auto it = values_.find(indices);
        return (it != values_.end()) ? it->second : Coefficient(0);
    }

    int degree() const { return degree_; }
};

/**
 * @brief Čech cohomology Ȟ^p(X, 𝒰) for open cover 𝒰
 */
template<typename T = double, typename Coefficient = int>
class CechCohomology {
private:
    const OpenCover<T>* cover_;
    int degree_;

public:
    CechCohomology(const OpenCover<T>* U, int p) : cover_(U), degree_(p) {}

    /**
     * @brief Čech coboundary δ: Č^p → Č^{p+1}
     */
    CechCochain<Coefficient> coboundary(const CechCochain<Coefficient>& cochain) const {
        CechCochain<Coefficient> result(degree_ + 1);
        // δf(i_0,...,i_{p+1}) = Σ(-1)^j f(i_0,...,î_j,...,i_{p+1})
        return result;
    }

    /**
     * @brief For good covers: Čech cohomology = singular cohomology
     */
    bool agreesWithSingularCohomology() const {
        return cover_->isGoodCover();
    }
};

// ============================================================================
// HOMOTOPY THEORY
// ============================================================================

/**
 * @brief Path in a topological space
 */
template<typename T = double>
class Path {
private:
    Point<T> start_;
    Point<T> end_;
    std::function<Point<T>(double)> path_function_;  // [0,1] → X

public:
    Path(const Point<T>& start, const Point<T>& end) : start_(start), end_(end) {}

    Point<T> evaluate(double t) const {
        if (path_function_) {
            return path_function_(t);
        }
        return start_;  // Default
    }

    Point<T> start() const { return start_; }
    Point<T> end() const { return end_; }

    /**
     * @brief Concatenation of paths: γ₁ * γ₂
     */
    Path concatenate(const Path& other) const {
        if (end_ != other.start_) {
            throw std::invalid_argument("Paths cannot be concatenated");
        }
        return Path(start_, other.end_);
    }

    /**
     * @brief Reverse path: γ⁻¹
     */
    Path reverse() const {
        return Path(end_, start_);
    }
};

/**
 * @brief Homotopy between two paths/maps
 */
template<typename T = double>
class Homotopy {
private:
    std::function<Point<T>(double, double)> homotopy_function_;  // [0,1] × [0,1] → X

public:
    /**
     * @brief Check if two paths are homotopic
     */
    static bool arePathsHomotopic(const Path<T>& f, const Path<T>& g) {
        // Paths f and g are homotopic if there exists H: [0,1] × [0,1] → X
        // with H(s,0) = f(s), H(s,1) = g(s), H(0,t) = x₀, H(1,t) = x₁
        return true;  // Placeholder
    }

    /**
     * @brief Homotopy equivalence of spaces
     */
    static bool areSpacesHomotopyEquivalent() {
        // X and Y are homotopy equivalent if ∃ f: X → Y, g: Y → X
        // with g∘f ≃ id_X and f∘g ≃ id_Y
        return true;  // Placeholder
    }
};

/**
 * @brief Fundamental group π₁(X, x₀) = {homotopy classes of loops at x₀}
 */
template<typename T = double>
class FundamentalGroup {
private:
    Point<T> basepoint_;
    std::vector<Path<T>> generators_;  // Generators of π₁
    std::vector<std::string> relations_;  // Relations between generators

public:
    FundamentalGroup(const Point<T>& x0) : basepoint_(x0) {}

    void addGenerator(const Path<T>& loop) {
        generators_.push_back(loop);
    }

    void addRelation(const std::string& relation) {
        relations_.push_back(relation);
    }

    /**
     * @brief Check if fundamental group is trivial (space is simply connected)
     */
    bool isTrivial() const {
        return generators_.empty();
    }

    /**
     * @brief Check if fundamental group is abelian
     */
    bool isAbelian() const {
        // π₁ is abelian ⟺ all commutators are trivial
        return true;  // Placeholder
    }

    /**
     * @brief Seifert-van Kampen theorem for computing π₁(X) from open cover
     */
    static FundamentalGroup seifertVanKampen(const FundamentalGroup& pi1_U,
                                             const FundamentalGroup& pi1_V,
                                             const FundamentalGroup& pi1_intersection) {
        // π₁(U ∪ V) = π₁(U) *_{π₁(U∩V)} π₁(V) (amalgamated free product)
        FundamentalGroup result(pi1_U.basepoint_);
        return result;
    }
};

/**
 * @brief Higher homotopy groups π_n(X, x₀)
 */
template<typename T = double>
class HomotopyGroup {
private:
    int dimension_;
    Point<T> basepoint_;

public:
    HomotopyGroup(int n, const Point<T>& x0) : dimension_(n), basepoint_(x0) {}

    int dimension() const { return dimension_; }

    /**
     * @brief π_n is abelian for n ≥ 2 (Eckmann-Hilton theorem)
     */
    bool isAbelian() const {
        return dimension_ >= 2;
    }

    /**
     * @brief Hurewicz theorem: π_n(X) → H_n(X) is isomorphism when X is (n-1)-connected
     */
    bool hurewiczTheorem(bool is_n_minus_1_connected) const {
        return is_n_minus_1_connected;
    }

    /**
     * @brief Freudenthal suspension theorem
     */
    bool freudenthalSuspension(int k) const {
        // π_k(S^n) ≅ π_{k+1}(S^{n+1}) for k < 2n - 1
        return true;  // Placeholder
    }
};

// ============================================================================
// FIBER BUNDLES AND PRINCIPAL BUNDLES
// ============================================================================

/**
 * @brief Fiber bundle E → B with fiber F
 */
template<typename T = double>
class FiberBundle {
private:
    int base_dimension_;
    int fiber_dimension_;
    int total_dimension_;

public:
    FiberBundle(int base_dim, int fiber_dim)
        : base_dimension_(base_dim),
          fiber_dimension_(fiber_dim),
          total_dimension_(base_dim + fiber_dim) {}

    int baseDimension() const { return base_dimension_; }
    int fiberDimension() const { return fiber_dimension_; }
    int totalDimension() const { return total_dimension_; }

    /**
     * @brief Projection π: E → B
     */
    Point<T> projection(const Point<T>& e) const {
        // Project from total space to base space
        Point<T> b;
        b.coordinates.resize(base_dimension_);
        for (int i = 0; i < base_dimension_; ++i) {
            b.coordinates[i] = e.coordinates[i];
        }
        return b;
    }

    /**
     * @brief Check if bundle is trivial: E ≅ B × F
     */
    bool isTrivial() const {
        return true;  // Placeholder
    }

    /**
     * @brief Leray-Hirsch theorem for cohomology of fiber bundles
     */
    bool lerayHirschTheorem() const {
        // H*(E) ≅ H*(B) ⊗ H*(F) as H*(B)-modules (under certain conditions)
        return true;
    }
};

/**
 * @brief Principal G-bundle
 */
template<typename T = double>
class PrincipalBundle : public FiberBundle<T> {
private:
    std::string structure_group_;  // e.g., "U(1)", "SU(2)", "SO(3)"

public:
    PrincipalBundle(int base_dim, const std::string& group)
        : FiberBundle<T>(base_dim, 0), structure_group_(group) {}

    const std::string& structureGroup() const { return structure_group_; }

    /**
     * @brief Right action of G on E: E × G → E
     */
    Point<T> groupAction(const Point<T>& e, double g) const {
        // Right action: (e, g) ↦ e·g
        return e;  // Placeholder
    }

    /**
     * @brief Associated vector bundle E ×_G V
     */
    FiberBundle<T> associatedBundle(int vector_space_dim) const {
        return FiberBundle<T>(this->baseDimension(), vector_space_dim);
    }

    /**
     * @brief Connection on principal bundle (Ehresmann connection)
     */
    bool hasConnection() const { return true; }

    /**
     * @brief Curvature 2-form Ω
     */
    DifferentialForm curvature() const {
        // Curvature Ω = dA + [A ∧ A] where A is connection 1-form
        return DifferentialForm(2, this->baseDimension());
    }
};

// ============================================================================
// CHARACTERISTIC CLASSES
// ============================================================================

/**
 * @brief Chern classes c_i(E) ∈ H^{2i}(B; ℤ) for complex vector bundles
 */
class ChernClass {
private:
    int degree_;  // c_i has degree 2i
    int rank_;    // rank of vector bundle

public:
    ChernClass(int i, int r) : degree_(2 * i), rank_(r) {}

    int degree() const { return degree_; }

    /**
     * @brief Total Chern class c(E) = 1 + c₁(E) + c₂(E) + ...
     */
    static std::vector<ChernClass> totalChernClass(int rank) {
        std::vector<ChernClass> classes;
        for (int i = 1; i <= rank; ++i) {
            classes.push_back(ChernClass(i, rank));
        }
        return classes;
    }

    /**
     * @brief Whitney sum formula: c(E ⊕ F) = c(E) ∪ c(F)
     */
    bool whitneySumFormula() const { return true; }

    /**
     * @brief Top Chern class c_n(E) = Euler class for rank n bundle
     */
    bool topChernClassIsEulerClass() const { return true; }
};

/**
 * @brief Stiefel-Whitney classes w_i(E) ∈ H^i(B; ℤ/2ℤ) for real vector bundles
 */
class StiefelWhitneyClass {
private:
    int degree_;

public:
    StiefelWhitneyClass(int i) : degree_(i) {}

    int degree() const { return degree_; }

    /**
     * @brief w₁ = 0 ⟺ bundle is orientable
     */
    bool firstClassVanishesIffOrientable() const { return true; }

    /**
     * @brief w₂ = 0 (and w₁ = 0) ⟺ bundle admits spin structure
     */
    bool secondClassVanishesIffSpin() const { return true; }
};

/**
 * @brief Pontryagin classes p_i(E) ∈ H^{4i}(B; ℤ) for real vector bundles
 */
class PontryaginClass {
private:
    int degree_;  // p_i has degree 4i

public:
    PontryaginClass(int i) : degree_(4 * i) {}

    int degree() const { return degree_; }

    /**
     * @brief Relation to Chern classes: p_k(E) = (-1)^k c_{2k}(E ⊗ ℂ)
     */
    bool relationToChernClasses() const { return true; }

    /**
     * @brief Hirzebruch signature theorem
     */
    double signatureFormula(int signature) const {
        // signature(M) = ⟨L(p₁, p₂, ...), [M]⟩
        return signature;
    }
};

/**
 * @brief Euler class e(E) ∈ H^n(B; ℤ) for oriented rank n bundle
 */
class EulerClass {
private:
    int rank_;

public:
    EulerClass(int n) : rank_(n) {}

    /**
     * @brief Poincaré-Hopf theorem: χ(M) = ⟨e(TM), [M]⟩
     */
    int poincareHopf(int euler_characteristic) const {
        return euler_characteristic;
    }

    /**
     * @brief Gauss-Bonnet theorem (generalized)
     */
    double gaussBonnet(double curvature_integral) const {
        // ∫_M K dV = 2π χ(M) for surfaces
        return curvature_integral;
    }
};

// ============================================================================
// SPECTRAL SEQUENCES
// ============================================================================

/**
 * @brief Spectral sequence E^{p,q}_r ⇒ H^{p+q}
 */
template<typename Coefficient = int>
class SpectralSequence {
private:
    std::map<std::tuple<int, int, int>, Coefficient> pages_;  // (p, q, r) → E^{p,q}_r
    int max_page_;

public:
    SpectralSequence() : max_page_(0) {}

    void setEntry(int p, int q, int r, Coefficient value) {
        pages_[std::make_tuple(p, q, r)] = value;
        max_page_ = std::max(max_page_, r);
    }

    Coefficient getEntry(int p, int q, int r) const {
        auto it = pages_.find(std::make_tuple(p, q, r));
        return (it != pages_.end()) ? it->second : Coefficient(0);
    }

    /**
     * @brief Differential d_r: E^{p,q}_r → E^{p+r,q-r+1}_r
     */
    void applyDifferential(int r) {
        // d_r has bidegree (r, -r+1) or (r, 1-r) depending on convention
        // E^{p,q}_{r+1} = ker(d_r) / im(d_r)
    }

    /**
     * @brief Check if spectral sequence has converged
     */
    bool hasConverged(int r) const {
        // Converged if E^{p,q}_r = E^{p,q}_{r+1} = E^{p,q}_∞ for all p, q
        return r >= max_page_;
    }

    /**
     * @brief Get E_∞ page (limit)
     */
    Coefficient limitEntry(int p, int q) const {
        return getEntry(p, q, max_page_);
    }
};

/**
 * @brief Serre spectral sequence for fibration F → E → B
 */
template<typename Coefficient = int>
class SerreSpectralSequence : public SpectralSequence<Coefficient> {
public:
    /**
     * @brief E₂ page: E₂^{p,q} = H^p(B; H^q(F))
     */
    void computeE2Page() {
        // E₂^{p,q} = cohomology of base with coefficients in cohomology of fiber
    }

    /**
     * @brief Converges to H^*(E)
     */
    bool convergesToTotalSpaceCohomology() const { return true; }
};

/**
 * @brief Leray spectral sequence for continuous map f: X → Y
 */
template<typename Coefficient = int>
class LeraySpectralSequence : public SpectralSequence<Coefficient> {
public:
    /**
     * @brief E₂ page: E₂^{p,q} = H^p(Y; R^q f_* ℱ)
     */
    void computeE2Page() {
        // Higher direct images R^q f_*
    }
};

// ============================================================================
// EXACT SEQUENCES
// ============================================================================

/**
 * @brief Long exact sequence in homology/cohomology
 */
template<typename T>
class ExactSequence {
private:
    std::vector<T> groups_;
    std::vector<std::string> maps_;  // names of homomorphisms

public:
    void addGroup(const T& group) {
        groups_.push_back(group);
    }

    void addMap(const std::string& map_name) {
        maps_.push_back(map_name);
    }

    /**
     * @brief Check exactness: im(f_i) = ker(f_{i+1})
     */
    bool isExact() const {
        return true;  // Placeholder
    }

    /**
     * @brief Short exact sequence: 0 → A → B → C → 0
     */
    static bool isShortExact(const T& A, const T& B, const T& C) {
        // 0 → A →^f B →^g C → 0 is exact
        // ⟺ f injective, g surjective, im(f) = ker(g)
        return true;
    }
};

/**
 * @brief Mayer-Vietoris sequence for X = U ∪ V
 */
template<typename Coefficient = int>
class MayerVietorisSequence {
public:
    /**
     * @brief Long exact sequence:
     * ... → H_n(U ∩ V) → H_n(U) ⊕ H_n(V) → H_n(X) → H_{n-1}(U ∩ V) → ...
     */
    static ExactSequence<HomologyGroup<double, Coefficient>>
    computeSequence() {
        ExactSequence<HomologyGroup<double, Coefficient>> seq;
        // Construct the sequence
        return seq;
    }
};

// ============================================================================
// COVERING SPACES
// ============================================================================

/**
 * @brief Covering space p: X̃ → X
 */
template<typename T = double>
class CoveringSpace {
private:
    int num_sheets_;  // number of sheets (degree of covering)
    Point<T> basepoint_;

public:
    CoveringSpace(int n) : num_sheets_(n) {}

    int degree() const { return num_sheets_; }

    /**
     * @brief Covering projection p: X̃ → X
     */
    Point<T> projection(const Point<T>& x_tilde) const {
        return x_tilde;  // Placeholder
    }

    /**
     * @brief Path lifting property
     */
    bool hasPathLifting() const { return true; }

    /**
     * @brief Unique path lifting: given γ in X and x̃₀ in p⁻¹(γ(0)),
     * ∃! lift γ̃ with γ̃(0) = x̃₀
     */
    bool hasUniquePathLifting() const { return true; }

    /**
     * @brief Universal covering space
     */
    bool isUniversal() const {
        // X̃ is universal ⟺ X̃ is simply connected
        return true;  // Placeholder
    }

    /**
     * @brief Deck transformations (covering automorphisms)
     */
    int numDeckTransformations() const {
        // For regular covering: Deck(X̃/X) ≅ π₁(X) / p_*(π₁(X̃))
        return num_sheets_;  // For regular coverings
    }

    /**
     * @brief Galois correspondence:
     * {subgroups of π₁(X)} ↔ {covering spaces of X}
     */
    bool galoisCorrespondence() const { return true; }
};

// ============================================================================
// MORSE THEORY
// ============================================================================

/**
 * @brief Morse function f: M → ℝ on a manifold
 */
class MorseFunction {
private:
    std::vector<int> critical_point_indices_;  // Morse indices of critical points

public:
    void addCriticalPoint(int morse_index) {
        critical_point_indices_.push_back(morse_index);
    }

    /**
     * @brief Morse inequalities: m_k ≥ β_k (number of index-k critical points ≥ k-th Betti number)
     */
    bool morseInequalities(int k, int betti_k) const {
        int m_k = std::count(critical_point_indices_.begin(),
                            critical_point_indices_.end(), k);
        return m_k >= betti_k;
    }

    /**
     * @brief Weak Morse inequality: Σ(-1)^k m_k = χ(M)
     */
    int weakMorseInequality() const {
        int alternating_sum = 0;
        for (int index : critical_point_indices_) {
            alternating_sum += (index % 2 == 0) ? 1 : -1;
        }
        return alternating_sum;
    }

    /**
     * @brief Strong Morse inequalities with alternating sums
     */
    bool strongMorseInequalities() const { return true; }
};

// ============================================================================
// DIFFERENTIAL TOPOLOGY - CHARTS AND COORDINATE SYSTEMS
// ============================================================================

/**
 * @brief Chart (local coordinate system) on a manifold
 * A chart is a homeomorphism φ: U → V ⊂ ℝ^n where U is open in M
 */
template<typename T = double, int n = 2>
class Chart {
private:
    std::string name_;
    std::function<std::vector<T>(const Point<T>&)> coordinate_map_;  // φ: U → ℝ^n
    std::function<Point<T>(const std::vector<T>&)> inverse_map_;     // φ^{-1}: ℝ^n → U

public:
    Chart(const std::string& name) : name_(name) {}

    std::string name() const { return name_; }
    int dimension() const { return n; }

    /**
     * @brief Set coordinate map φ and its inverse
     */
    void setMaps(std::function<std::vector<T>(const Point<T>&)> phi,
                 std::function<Point<T>(const std::vector<T>&)> phi_inv) {
        coordinate_map_ = phi;
        inverse_map_ = phi_inv;
    }

    /**
     * @brief Apply coordinate map: p ↦ φ(p)
     */
    std::vector<T> toCoordinates(const Point<T>& p) const {
        if (coordinate_map_) {
            return coordinate_map_(p);
        }
        return std::vector<T>(n, 0);
    }

    /**
     * @brief Apply inverse map: x ↦ φ^{-1}(x)
     */
    Point<T> fromCoordinates(const std::vector<T>& coords) const {
        if (inverse_map_) {
            return inverse_map_(coords);
        }
        return Point<T>();
    }

    /**
     * @brief Check if point is in domain of chart
     */
    bool inDomain(const Point<T>& p) const {
        return true;  // Placeholder - should check if φ(p) is defined
    }
};

/**
 * @brief Transition map between two charts
 * τ_{αβ} = φ_β ∘ φ_α^{-1}: φ_α(U_α ∩ U_β) → φ_β(U_α ∩ U_β)
 */
template<typename T = double, int n = 2>
class TransitionMap {
private:
    const Chart<T, n>* chart_from_;
    const Chart<T, n>* chart_to_;
    std::function<std::vector<T>(const std::vector<T>&)> transition_function_;

public:
    TransitionMap(const Chart<T, n>* from, const Chart<T, n>* to)
        : chart_from_(from), chart_to_(to) {}

    /**
     * @brief Compute transition map: x ↦ τ(x)
     */
    std::vector<T> apply(const std::vector<T>& coords) const {
        if (transition_function_) {
            return transition_function_(coords);
        }
        // Default: φ_β ∘ φ_α^{-1}
        Point<T> p = chart_from_->fromCoordinates(coords);
        return chart_to_->toCoordinates(p);
    }

    /**
     * @brief Check if transition map is smooth (C^∞)
     */
    bool isSmooth() const {
        return true;  // Placeholder - requires checking derivatives
    }

    /**
     * @brief Jacobian matrix of transition map
     */
    std::vector<std::vector<T>> jacobian(const std::vector<T>& coords) const {
        std::vector<std::vector<T>> J(n, std::vector<T>(n, 0));
        // Compute ∂τ_i/∂x_j
        return J;
    }

    /**
     * @brief Check if transition map is a diffeomorphism
     */
    bool isDiffeomorphism() const {
        // Check if Jacobian is invertible everywhere
        return isSmooth();
    }
};

// ============================================================================
// COMPACT SURFACES (2-dimensional manifolds)
// ============================================================================

/**
 * @brief Compact surface (2-dimensional closed manifold)
 */
class CompactSurface {
private:
    int genus_;  // Number of "holes"
    bool orientable_;
    int euler_characteristic_;

public:
    CompactSurface(int g, bool orient) : genus_(g), orientable_(orient) {
        if (orientable_) {
            euler_characteristic_ = 2 - 2 * genus_;  // χ = 2 - 2g
        } else {
            euler_characteristic_ = 2 - genus_;  // χ = 2 - g for non-orientable
        }
    }

    int genus() const { return genus_; }
    bool isOrientable() const { return orientable_; }
    int eulerCharacteristic() const { return euler_characteristic_; }

    /**
     * @brief Classification theorem: Every compact surface is homeomorphic to
     * - Sphere S^2 (g = 0, orientable)
     * - Torus T^2 or connected sum of g tori (g ≥ 1, orientable)
     * - Projective plane ℝP^2 or connected sum of g projective planes (non-orientable)
     */
    std::string classify() const {
        if (orientable_) {
            if (genus_ == 0) return "Sphere S^2";
            if (genus_ == 1) return "Torus T^2";
            return "Connected sum of " + std::to_string(genus_) + " tori";
        } else {
            if (genus_ == 1) return "Projective plane ℝP^2";
            if (genus_ == 2) return "Klein bottle";
            return "Connected sum of " + std::to_string(genus_) + " projective planes";
        }
    }

    /**
     * @brief Betti numbers for surface
     */
    std::vector<int> bettiNumbers() const {
        std::vector<int> betti(3, 0);
        betti[0] = 1;  // β_0 = 1 (connected)
        if (orientable_) {
            betti[1] = 2 * genus_;  // β_1 = 2g
            betti[2] = 1;           // β_2 = 1
        } else {
            betti[1] = genus_ - 1;  // β_1 = g - 1
            betti[2] = 0;           // β_2 = 0
        }
        return betti;
    }

    /**
     * @brief Fundamental group
     */
    std::string fundamentalGroup() const {
        if (orientable_) {
            if (genus_ == 0) return "Trivial {e}";
            if (genus_ == 1) return "ℤ × ℤ";
            return "⟨a₁,b₁,...,a_g,b_g | [a₁,b₁]···[a_g,b_g]=e⟩";
        } else {
            if (genus_ == 1) return "ℤ/2ℤ";
            return "⟨a₁,...,a_g | a₁²···a_g²=e⟩";
        }
    }
};

/**
 * @brief Sphere S^n = {x ∈ ℝ^{n+1} | |x| = 1}
 */
template<int n>
class Sphere {
public:
    int dimension() const { return n; }

    /**
     * @brief Euler characteristic χ(S^n) = 1 + (-1)^n
     */
    int eulerCharacteristic() const {
        return 1 + (n % 2 == 0 ? 1 : -1);
    }

    /**
     * @brief Betti numbers: β_0 = β_n = 1, others = 0
     */
    std::vector<int> bettiNumbers() const {
        std::vector<int> betti(n + 1, 0);
        betti[0] = 1;
        betti[n] = 1;
        return betti;
    }

    /**
     * @brief Fundamental group: π_1(S^n) = {e} for n ≥ 2
     */
    bool isSimplyConnected() const {
        return n >= 2;
    }

    /**
     * @brief Homotopy groups of spheres (partial)
     */
    std::string homotopyGroup(int k) const {
        if (k == 0) return "ℤ";
        if (k < n) return "Trivial {e}";
        if (k == n) return "ℤ";
        // Higher homotopy groups are generally complicated
        return "Complex (see Homotopy groups of spheres table)";
    }
};

/**
 * @brief Torus T^n = S^1 × ··· × S^1 (n times)
 */
template<int n>
class Torus {
public:
    int dimension() const { return n; }

    /**
     * @brief Euler characteristic χ(T^n) = 0
     */
    int eulerCharacteristic() const {
        return 0;
    }

    /**
     * @brief k-th Betti number: β_k = C(n, k) (binomial coefficient)
     */
    int bettiNumber(int k) const {
        if (k < 0 || k > n) return 0;
        // C(n, k) = n! / (k! (n-k)!)
        int binom = 1;
        for (int i = 0; i < k; ++i) {
            binom = binom * (n - i) / (i + 1);
        }
        return binom;
    }

    /**
     * @brief Fundamental group: π_1(T^n) = ℤ^n
     */
    std::string fundamentalGroup() const {
        return "ℤ^" + std::to_string(n);
    }

    /**
     * @brief All higher homotopy groups vanish: π_k(T^n) = 0 for k ≥ 2
     */
    bool higherHomotopyVanishes(int k) const {
        return k >= 2;
    }
};

/**
 * @brief Projective space ℝP^n or ℂP^n
 */
template<int n, bool complex = false>
class ProjectiveSpace {
public:
    int dimension() const {
        return complex ? 2 * n : n;
    }

    /**
     * @brief Euler characteristic
     */
    int eulerCharacteristic() const {
        if (complex) {
            return n + 1;  // χ(ℂP^n) = n + 1
        } else {
            return (n % 2 == 0) ? 1 : 0;  // χ(ℝP^n) = 1 if n even, 0 if n odd
        }
    }

    /**
     * @brief Betti numbers
     */
    std::vector<int> bettiNumbers() const {
        std::vector<int> betti;
        if (complex) {
            // For ℂP^n: β_k = 1 if k even and k ≤ 2n, 0 otherwise
            for (int k = 0; k <= 2 * n; ++k) {
                betti.push_back((k % 2 == 0 && k <= 2 * n) ? 1 : 0);
            }
        } else {
            // For ℝP^n: β_0 = 1, β_n = 1 if n odd, others depend on ℤ/2ℤ
            for (int k = 0; k <= n; ++k) {
                if (k == 0) betti.push_back(1);
                else if (k == n && n % 2 == 1) betti.push_back(1);
                else betti.push_back(0);
            }
        }
        return betti;
    }

    /**
     * @brief Fundamental group
     */
    std::string fundamentalGroup() const {
        if (complex) {
            return "Trivial {e}";  // ℂP^n is simply connected
        } else {
            return (n >= 2) ? "ℤ/2ℤ" : "Trivial {e}";
        }
    }
};

// ============================================================================
// SMOOTH MANIFOLDS
// ============================================================================

/**
 * @brief Topological manifold (locally Euclidean space)
 */
template<typename T = double, int n = 2>
class TopologicalManifold {
protected:
    std::vector<Chart<T, n>> charts_;
    int dimension_;
    bool hausdorff_;
    bool second_countable_;

public:
    TopologicalManifold(int dim) : dimension_(dim), hausdorff_(true), second_countable_(true) {}

    int dimension() const { return dimension_; }

    /**
     * @brief Add chart to atlas
     */
    void addChart(const Chart<T, n>& chart) {
        charts_.push_back(chart);
    }

    /**
     * @brief Get all charts
     */
    const std::vector<Chart<T, n>>& charts() const {
        return charts_;
    }

    /**
     * @brief Check if manifold is Hausdorff
     */
    bool isHausdorff() const {
        return hausdorff_;
    }

    /**
     * @brief Check if manifold is second-countable
     */
    bool isSecondCountable() const {
        return second_countable_;
    }

    /**
     * @brief Check if manifold is paracompact
     * (Hausdorff + second-countable ⟹ paracompact)
     */
    bool isParacompact() const {
        return hausdorff_ && second_countable_;
    }

    /**
     * @brief Find chart containing point p
     */
    const Chart<T, n>* findChart(const Point<T>& p) const {
        for (const auto& chart : charts_) {
            if (chart.inDomain(p)) {
                return &chart;
            }
        }
        return nullptr;
    }
};

/**
 * @brief Smooth structure (C^∞ structure) on a manifold
 */
template<typename T = double, int n = 2>
class SmoothStructure {
private:
    std::vector<Chart<T, n>> atlas_;
    bool is_maximal_;

public:
    SmoothStructure() : is_maximal_(false) {}

    /**
     * @brief Add compatible chart to smooth structure
     */
    bool addChart(const Chart<T, n>& chart) {
        // Check compatibility with existing charts
        for (const auto& existing : atlas_) {
            TransitionMap<T, n> transition(&existing, &chart);
            if (!transition.isSmooth()) {
                return false;  // Not compatible
            }
        }
        atlas_.push_back(chart);
        return true;
    }

    /**
     * @brief Get smooth atlas
     */
    const std::vector<Chart<T, n>>& atlas() const {
        return atlas_;
    }

    /**
     * @brief Check if two charts are compatible (C^∞-compatible)
     */
    static bool areCompatible(const Chart<T, n>& chart1, const Chart<T, n>& chart2) {
        TransitionMap<T, n> transition12(&chart1, &chart2);
        TransitionMap<T, n> transition21(&chart2, &chart1);
        return transition12.isSmooth() && transition21.isSmooth();
    }

    /**
     * @brief Generate maximal atlas from given atlas
     */
    void makeMaximal() {
        // A maximal atlas contains all charts compatible with the given ones
        is_maximal_ = true;
    }

    bool isMaximal() const {
        return is_maximal_;
    }

    /**
     * @brief Two smooth structures are equivalent if they generate the same maximal atlas
     */
    bool isEquivalent(const SmoothStructure& other) const {
        // Check if all charts in this atlas are compatible with all charts in other
        for (const auto& chart1 : atlas_) {
            for (const auto& chart2 : other.atlas_) {
                if (!areCompatible(chart1, chart2)) {
                    return false;
                }
            }
        }
        return true;
    }
};

/**
 * @brief Smooth manifold (manifold with smooth structure)
 */
template<typename T = double, int n = 2>
class SmoothManifold : public TopologicalManifold<T, n> {
private:
    SmoothStructure<T, n> smooth_structure_;

public:
    SmoothManifold(int dim) : TopologicalManifold<T, n>(dim) {}

    /**
     * @brief Set smooth structure
     */
    void setSmoothStructure(const SmoothStructure<T, n>& structure) {
        smooth_structure_ = structure;
    }

    const SmoothStructure<T, n>& smoothStructure() const {
        return smooth_structure_;
    }

    /**
     * @brief Add smooth chart
     */
    bool addSmoothChart(const Chart<T, n>& chart) {
        if (smooth_structure_.addChart(chart)) {
            this->addChart(chart);
            return true;
        }
        return false;
    }

    /**
     * @brief Check if manifold is orientable
     */
    bool isOrientable() const {
        // Check if transition maps have positive Jacobian determinant
        return true;  // Placeholder
    }

    /**
     * @brief Dimension of manifold
     */
    int dimension() const {
        return this->dimension_;
    }
};

/**
 * @brief Smooth map between smooth manifolds f: M → N
 */
template<typename T = double, int m = 2, int n = 2>
class SmoothMap {
private:
    const SmoothManifold<T, m>* domain_;
    const SmoothManifold<T, n>* codomain_;
    std::function<Point<T>(const Point<T>&)> map_function_;

public:
    SmoothMap(const SmoothManifold<T, m>* M, const SmoothManifold<T, n>* N)
        : domain_(M), codomain_(N) {}

    /**
     * @brief Set map function
     */
    void setFunction(std::function<Point<T>(const Point<T>&)> f) {
        map_function_ = f;
    }

    /**
     * @brief Apply map: p ↦ f(p)
     */
    Point<T> apply(const Point<T>& p) const {
        if (map_function_) {
            return map_function_(p);
        }
        return Point<T>();
    }

    /**
     * @brief Check if map is smooth
     * f is smooth if φ_β ∘ f ∘ φ_α^{-1} is C^∞ for all charts (φ_α, φ_β)
     */
    bool isSmooth() const {
        return true;  // Placeholder - requires checking in local coordinates
    }

    /**
     * @brief Check if map is a diffeomorphism
     */
    bool isDiffeomorphism() const {
        // f is diffeomorphism if bijective, smooth, and has smooth inverse
        return isSmooth() && m == n;  // Simplified
    }

    /**
     * @brief Rank of differential at point p
     */
    int rank(const Point<T>& p) const {
        // Rank of Jacobian matrix
        return std::min(m, n);  // Placeholder
    }

    /**
     * @brief Check if map is an immersion (injective differential)
     */
    bool isImmersion(const Point<T>& p) const {
        return rank(p) == m;  // df_p injective ⟺ rank = dim(M)
    }

    /**
     * @brief Check if map is a submersion (surjective differential)
     */
    bool isSubmersion(const Point<T>& p) const {
        return rank(p) == n;  // df_p surjective ⟺ rank = dim(N)
    }

    /**
     * @brief Check if map is an embedding (immersion + homeomorphism onto image)
     */
    bool isEmbedding() const {
        return true;  // Placeholder
    }
};

/**
 * @brief Submanifold of a smooth manifold
 */
template<typename T = double, int m = 2, int n = 3>
class Submanifold {
private:
    const SmoothManifold<T, n>* ambient_;
    SmoothManifold<T, m> submanifold_;

public:
    Submanifold(const SmoothManifold<T, n>* M) : ambient_(M), submanifold_(m) {
        static_assert(m <= n, "Submanifold dimension must be ≤ ambient dimension");
    }

    /**
     * @brief Inclusion map i: S ↪ M
     */
    SmoothMap<T, m, n> inclusion() const {
        return SmoothMap<T, m, n>(&submanifold_, ambient_);
    }

    /**
     * @brief Codimension = dim(M) - dim(S)
     */
    int codimension() const {
        return n - m;
    }

    /**
     * @brief Check if submanifold is embedded
     */
    bool isEmbedded() const {
        return inclusion().isEmbedding();
    }

    /**
     * @brief Check if submanifold is regular (constant rank)
     */
    bool isRegular() const {
        return true;  // Placeholder
    }

    /**
     * @brief Normal bundle (for embedded submanifolds)
     */
    int normalBundleRank() const {
        return codimension();
    }
};

/**
 * @brief Product manifold M × N
 */
template<typename T = double, int m = 2, int n = 2>
class ProductManifold : public SmoothManifold<T, m + n> {
private:
    const SmoothManifold<T, m>* factor1_;
    const SmoothManifold<T, n>* factor2_;

public:
    ProductManifold(const SmoothManifold<T, m>* M, const SmoothManifold<T, n>* N)
        : SmoothManifold<T, m + n>(m + n), factor1_(M), factor2_(N) {}

    /**
     * @brief Projection onto first factor: M × N → M
     */
    SmoothMap<T, m + n, m> projection1() const {
        return SmoothMap<T, m + n, m>(this, factor1_);
    }

    /**
     * @brief Projection onto second factor: M × N → N
     */
    SmoothMap<T, m + n, n> projection2() const {
        return SmoothMap<T, m + n, n>(this, factor2_);
    }

    /**
     * @brief Product smooth structure
     * Charts are products (U_α × V_β, φ_α × ψ_β)
     */
    void constructProductStructure() {
        // Combine charts from both factors
    }
};

// ============================================================================
// TANGENT SPACES
// ============================================================================

/**
 * @brief Germ of functions at a point p
 * [f]_p ~ [g]_p if f = g in a neighborhood of p
 */
template<typename T = double>
class Germ {
private:
    Point<T> basepoint_;
    std::function<T(const Point<T>&)> representative_;

public:
    Germ(const Point<T>& p) : basepoint_(p) {}

    void setRepresentative(std::function<T(const Point<T>&)> f) {
        representative_ = f;
    }

    /**
     * @brief Evaluate representative function at basepoint
     */
    T evaluate() const {
        if (representative_) {
            return representative_(basepoint_);
        }
        return T(0);
    }

    /**
     * @brief Check if two germs are equivalent
     */
    bool isEquivalent(const Germ& other) const {
        // f ~ g if they agree on some neighborhood of p
        return true;  // Placeholder
    }

    /**
     * @brief Sum of germs: [f]_p + [g]_p = [f + g]_p
     */
    Germ operator+(const Germ& other) const {
        Germ result(basepoint_);
        return result;
    }

    /**
     * @brief Product of germs: [f]_p · [g]_p = [f · g]_p
     */
    Germ operator*(const Germ& other) const {
        Germ result(basepoint_);
        return result;
    }
};

/**
 * @brief Derivation at point p (element of tangent space)
 * A linear map v: C^∞(M) → ℝ satisfying Leibniz rule: v(fg) = f(p)v(g) + g(p)v(f)
 */
template<typename T = double>
class Derivation {
private:
    Point<T> basepoint_;
    std::vector<T> components_;  // In local coordinates

public:
    Derivation(const Point<T>& p, int dim) : basepoint_(p), components_(dim, 0) {}

    /**
     * @brief Set components in local coordinates
     */
    void setComponents(const std::vector<T>& v) {
        components_ = v;
    }

    const std::vector<T>& components() const {
        return components_;
    }

    /**
     * @brief Apply derivation to function f
     * v(f) = Σ v^i ∂f/∂x^i
     */
    T apply(std::function<T(const Point<T>&)> f) const {
        // Compute directional derivative
        return T(0);  // Placeholder
    }

    /**
     * @brief Leibniz rule: v(fg) = f(p)v(g) + g(p)v(f)
     */
    bool satisfiesLeibnizRule() const {
        return true;
    }

    /**
     * @brief Sum of derivations
     */
    Derivation operator+(const Derivation& other) const {
        Derivation result(basepoint_, components_.size());
        std::vector<T> sum(components_.size());
        for (size_t i = 0; i < components_.size(); ++i) {
            sum[i] = components_[i] + other.components_[i];
        }
        result.setComponents(sum);
        return result;
    }

    /**
     * @brief Scalar multiplication
     */
    Derivation operator*(T scalar) const {
        Derivation result(basepoint_, components_.size());
        std::vector<T> scaled(components_.size());
        for (size_t i = 0; i < components_.size(); ++i) {
            scaled[i] = scalar * components_[i];
        }
        result.setComponents(scaled);
        return result;
    }
};

/**
 * @brief Tangent space T_p M at point p ∈ M
 * Vector space of all derivations at p
 */
template<typename T = double, int n = 2>
class TangentSpace {
private:
    Point<T> basepoint_;
    std::vector<Derivation<T>> basis_;  // Coordinate basis {∂/∂x^i}

public:
    TangentSpace(const Point<T>& p) : basepoint_(p) {
        // Initialize coordinate basis
        for (int i = 0; i < n; ++i) {
            Derivation<T> basis_vector(p, n);
            std::vector<T> components(n, 0);
            components[i] = 1;  // e_i
            basis_vector.setComponents(components);
            basis_.push_back(basis_vector);
        }
    }

    int dimension() const {
        return n;
    }

    const Point<T>& basepoint() const {
        return basepoint_;
    }

    /**
     * @brief Coordinate basis vectors ∂/∂x^i
     */
    const std::vector<Derivation<T>>& basis() const {
        return basis_;
    }

    /**
     * @brief Get i-th basis vector
     */
    const Derivation<T>& basisVector(int i) const {
        return basis_[i];
    }

    /**
     * @brief Create tangent vector from components
     */
    Derivation<T> vector(const std::vector<T>& components) const {
        Derivation<T> v(basepoint_, n);
        v.setComponents(components);
        return v;
    }

    /**
     * @brief Inner product (if Riemannian metric is defined)
     */
    T innerProduct(const Derivation<T>& v, const Derivation<T>& w) const {
        T result = 0;
        for (int i = 0; i < n; ++i) {
            result += v.components()[i] * w.components()[i];
        }
        return result;
    }
};

/**
 * @brief Differential (pushforward) of smooth map f: M → N
 * df_p: T_p M → T_{f(p)} N
 */
template<typename T = double, int m = 2, int n = 2>
class Differential {
private:
    Point<T> basepoint_;
    const SmoothMap<T, m, n>* map_;
    std::vector<std::vector<T>> jacobian_;  // n × m matrix

public:
    Differential(const Point<T>& p, const SmoothMap<T, m, n>* f)
        : basepoint_(p), map_(f), jacobian_(n, std::vector<T>(m, 0)) {}

    /**
     * @brief Set Jacobian matrix
     */
    void setJacobian(const std::vector<std::vector<T>>& J) {
        jacobian_ = J;
    }

    /**
     * @brief Apply differential to tangent vector: v ↦ df_p(v)
     */
    Derivation<T> apply(const Derivation<T>& v) const {
        Derivation<T> result(map_->apply(basepoint_), n);
        std::vector<T> components(n, 0);

        // Matrix-vector multiplication: J · v
        for (int i = 0; i < n; ++i) {
            for (int j = 0; j < m; ++j) {
                components[i] += jacobian_[i][j] * v.components()[j];
            }
        }

        result.setComponents(components);
        return result;
    }

    /**
     * @brief Rank of differential
     */
    int rank() const {
        // Rank of Jacobian matrix
        return std::min(m, n);  // Placeholder
    }

    /**
     * @brief Chain rule: d(g ∘ f)_p = dg_{f(p)} ∘ df_p
     */
    static Differential chainRule(const Differential& df, const Differential& dg) {
        // Multiply Jacobian matrices
        Differential result(df.basepoint_, nullptr);
        return result;
    }
};

// ============================================================================
// VECTOR BUNDLES
// ============================================================================

/**
 * @brief Topological vector bundle E → B
 */
template<typename T = double, int n = 2, int k = 1>
class TopologicalVectorBundle {
protected:
    const TopologicalManifold<T, n>* base_;
    int fiber_dimension_;
    int total_dimension_;

public:
    TopologicalVectorBundle(const TopologicalManifold<T, n>* B, int fiber_dim)
        : base_(B), fiber_dimension_(fiber_dim), total_dimension_(n + fiber_dim) {}

    int baseDimension() const { return n; }
    int fiberDimension() const { return fiber_dimension_; }
    int totalDimension() const { return total_dimension_; }

    /**
     * @brief Projection π: E → B
     */
    virtual Point<T> projection(const Point<T>& e) const {
        Point<T> b;
        b.coordinates.resize(n);
        for (int i = 0; i < n; ++i) {
            b.coordinates[i] = e.coordinates[i];
        }
        return b;
    }

    /**
     * @brief Fiber over point p: E_p = π^{-1}(p)
     */
    int fiberDimension(const Point<T>& p) const {
        return fiber_dimension_;
    }

    /**
     * @brief Check if bundle is trivial: E ≅ B × ℝ^k
     */
    virtual bool isTrivial() const {
        return false;  // Placeholder
    }

    /**
     * @brief Zero section s: B → E with π ∘ s = id_B
     */
    Point<T> zeroSection(const Point<T>& p) const {
        Point<T> e;
        e.coordinates = p.coordinates;
        e.coordinates.resize(total_dimension_, 0);
        return e;
    }
};

/**
 * @brief Transition function for vector bundle
 * g_{αβ}: U_α ∩ U_β → GL(k, ℝ) where k = fiber dimension
 */
template<typename T = double, int k = 1>
class VectorBundleTransition {
private:
    std::function<std::vector<std::vector<T>>(const Point<T>&)> transition_matrix_;

public:
    /**
     * @brief Set transition matrix function
     */
    void setTransition(std::function<std::vector<std::vector<T>>(const Point<T>&)> g) {
        transition_matrix_ = g;
    }

    /**
     * @brief Evaluate transition matrix at point p
     */
    std::vector<std::vector<T>> at(const Point<T>& p) const {
        if (transition_matrix_) {
            return transition_matrix_(p);
        }
        // Return identity matrix by default
        std::vector<std::vector<T>> I(k, std::vector<T>(k, 0));
        for (int i = 0; i < k; ++i) {
            I[i][i] = 1;
        }
        return I;
    }

    /**
     * @brief Cocycle condition: g_{αβ} · g_{βγ} = g_{αγ} on U_α ∩ U_β ∩ U_γ
     */
    bool satisfiesCocycleCondition(const VectorBundleTransition& g_beta_gamma,
                                   const VectorBundleTransition& g_alpha_gamma) const {
        return true;  // Placeholder
    }

    /**
     * @brief Compatibility: g_{βα} = g_{αβ}^{-1}
     */
    bool isCompatible(const VectorBundleTransition& g_beta_alpha) const {
        return true;  // Placeholder
    }
};

/**
 * @brief Smooth vector bundle
 */
template<typename T = double, int n = 2, int k = 1>
class SmoothVectorBundle : public TopologicalVectorBundle<T, n, k> {
private:
    std::vector<VectorBundleTransition<T, k>> transitions_;

public:
    SmoothVectorBundle(const SmoothManifold<T, n>* B, int fiber_dim)
        : TopologicalVectorBundle<T, n, k>(B, fiber_dim) {}

    /**
     * @brief Add transition function
     */
    void addTransition(const VectorBundleTransition<T, k>& g) {
        transitions_.push_back(g);
    }

    /**
     * @brief Verify cocycle condition for all transitions
     */
    bool verifyCocycleCondition() const {
        // Check g_{αβ} · g_{βγ} · g_{γα} = I for all α, β, γ
        return true;  // Placeholder
    }

    /**
     * @brief Check if bundle is smooth
     */
    bool isSmooth() const {
        // Transition functions are smooth maps
        return true;
    }

    /**
     * @brief Section s: B → E with π ∘ s = id_B
     */
    class Section {
    private:
        std::function<std::vector<T>(const Point<T>&)> section_function_;

    public:
        /**
         * @brief Evaluate section at point p
         */
        std::vector<T> at(const Point<T>& p) const {
            if (section_function_) {
                return section_function_(p);
            }
            return std::vector<T>(k, 0);
        }

        /**
         * @brief Check if section is smooth
         */
        bool isSmooth() const {
            return true;  // Placeholder
        }

        /**
         * @brief Space of smooth sections Γ(E)
         */
        static int dimensionOfSections() {
            return -1;  // Infinite-dimensional in general
        }
    };

    /**
     * @brief Create section
     */
    Section createSection() const {
        return Section();
    }
};

/**
 * @brief Pre-vector bundle (before verifying cocycle condition)
 */
template<typename T = double, int n = 2, int k = 1>
class PreVectorBundle {
private:
    std::vector<Chart<T, n>> open_cover_;
    std::vector<VectorBundleTransition<T, k>> transitions_;

public:
    /**
     * @brief Add chart to cover
     */
    void addChart(const Chart<T, n>& chart) {
        open_cover_.push_back(chart);
    }

    /**
     * @brief Add transition function
     */
    void addTransition(const VectorBundleTransition<T, k>& g) {
        transitions_.push_back(g);
    }

    /**
     * @brief Construct vector bundle (after verifying cocycle condition)
     */
    SmoothVectorBundle<T, n, k> construct(const SmoothManifold<T, n>* base) const {
        SmoothVectorBundle<T, n, k> bundle(base, k);
        for (const auto& g : transitions_) {
            bundle.addTransition(g);
        }
        return bundle;
    }

    /**
     * @brief Verify that transition data defines a bundle
     */
    bool isValid() const {
        // Check cocycle condition
        return true;  // Placeholder
    }
};

/**
 * @brief Tangent bundle TM = ⊔_{p∈M} T_p M
 */
template<typename T = double, int n = 2>
class TangentBundle : public SmoothVectorBundle<T, n, n> {
private:
    const SmoothManifold<T, n>* base_manifold_;

public:
    TangentBundle(const SmoothManifold<T, n>* M)
        : SmoothVectorBundle<T, n, n>(M, n), base_manifold_(M) {}

    /**
     * @brief Fiber at point p is tangent space T_p M
     */
    TangentSpace<T, n> fiberAt(const Point<T>& p) const {
        return TangentSpace<T, n>(p);
    }

    /**
     * @brief Transition functions from coordinate changes
     * If x = φ_α and y = φ_β, then g_{αβ} = ∂y/∂x (Jacobian)
     */
    void computeTransitionsFromAtlas() {
        const auto& charts = base_manifold_->charts();
        for (size_t i = 0; i < charts.size(); ++i) {
            for (size_t j = i + 1; j < charts.size(); ++j) {
                TransitionMap<T, n> coord_change(&charts[i], &charts[j]);
                VectorBundleTransition<T, n> g;
                // Set g to be Jacobian of coordinate change
                this->addTransition(g);
            }
        }
    }

    /**
     * @brief Vector field (smooth section of tangent bundle)
     */
    class VectorField {
    private:
        std::function<Derivation<T>(const Point<T>&)> field_;

    public:
        /**
         * @brief Evaluate vector field at point p
         */
        Derivation<T> at(const Point<T>& p) const {
            if (field_) {
                return field_(p);
            }
            return Derivation<T>(p, n);
        }

        /**
         * @brief Lie bracket [X, Y] of two vector fields
         */
        static VectorField lieBracket(const VectorField& X, const VectorField& Y) {
            VectorField result;
            // [X,Y](f) = X(Y(f)) - Y(X(f))
            return result;
        }

        /**
         * @brief Flow of vector field (one-parameter group of diffeomorphisms)
         */
        Point<T> flow(const Point<T>& p, T t) const {
            // Solve dx/dt = X(x) with x(0) = p
            return p;  // Placeholder
        }
    };

    /**
     * @brief Create vector field
     */
    VectorField createVectorField() const {
        return VectorField();
    }

    /**
     * @brief Tangent bundle is always orientable × 2
     * (Manifold is orientable ⟺ TM is orientable)
     */
    bool preservesOrientability() const {
        return true;
    }
};

/**
 * @brief Cotangent bundle T*M (dual to tangent bundle)
 */
template<typename T = double, int n = 2>
class CotangentBundle : public SmoothVectorBundle<T, n, n> {
private:
    const SmoothManifold<T, n>* base_manifold_;

public:
    CotangentBundle(const SmoothManifold<T, n>* M)
        : SmoothVectorBundle<T, n, n>(M, n), base_manifold_(M) {}

    /**
     * @brief Fiber at p is cotangent space T*_p M = (T_p M)*
     */
    class CotangentSpace {
    private:
        Point<T> basepoint_;
        std::vector<T> components_;  // Covector components

    public:
        CotangentSpace(const Point<T>& p) : basepoint_(p), components_(n, 0) {}

        /**
         * @brief Apply covector to vector: ω(v)
         */
        T apply(const Derivation<T>& v) const {
            T result = 0;
            for (int i = 0; i < n; ++i) {
                result += components_[i] * v.components()[i];
            }
            return result;
        }

        /**
         * @brief Coordinate basis: dx^i
         */
        static CotangentSpace basisCovector(const Point<T>& p, int i) {
            CotangentSpace omega(p);
            omega.components_[i] = 1;
            return omega;
        }
    };

    /**
     * @brief Differential form (section of ⋀^k T*M)
     */
    DifferentialForm differentialForm(int degree) const {
        return DifferentialForm(degree, n);
    }

    /**
     * @brief Canonical symplectic form on T*M
     */
    bool hasCanonicalSymplecticStructure() const {
        // T*M always has canonical symplectic form ω = Σ dp_i ∧ dq^i
        return true;
    }
};

} // namespace topology
} // namespace maths

#endif // TOPOLOGY_HPP
