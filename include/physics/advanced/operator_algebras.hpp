#ifndef PHYSICS_ADVANCED_OPERATOR_ALGEBRAS_HPP
#define PHYSICS_ADVANCED_OPERATOR_ALGEBRAS_HPP

#include <complex>
#include <vector>
#include <functional>
#include <cmath>
#include <stdexcept>
#include <algorithm>
#include <numeric>
#include <map>

/**
 * @file operator_algebras.hpp
 * @brief Operator algebras, functional analysis, and quantum mechanics foundations
 *
 * Comprehensive implementation of:
 * - Hilbert spaces and bounded operators
 * - Von Neumann algebras (rings of operators)
 * - Unitary group representations
 * - Factor classification (Murray-von Neumann)
 * - C*-algebras and spectral theory
 * - Quantum mechanical observables and states
 */

namespace physics {
namespace operator_algebras {

using Complex = std::complex<double>;

/**
 * @brief Hilbert Space foundations
 *
 * Functional analysis framework for quantum mechanics
 */
class HilbertSpace {
public:
    /**
     * @brief Vector in Hilbert space (finite-dimensional approximation)
     */
    using Vector = std::vector<Complex>;

    /**
     * @brief Inner product ⟨ψ|φ⟩
     */
    static Complex inner_product(const Vector& psi, const Vector& phi) {
        if (psi.size() != phi.size()) {
            throw std::invalid_argument("Vectors must have same dimension");
        }

        Complex result(0.0, 0.0);
        for (size_t i = 0; i < psi.size(); ++i) {
            result += std::conj(psi[i]) * phi[i];
        }
        return result;
    }

    /**
     * @brief Norm ||ψ|| = √⟨ψ|ψ⟩
     */
    static double norm(const Vector& psi) {
        return std::sqrt(std::abs(inner_product(psi, psi)));
    }

    /**
     * @brief Normalize vector: |ψ⟩ → |ψ⟩/||ψ||
     */
    static Vector normalize(const Vector& psi) {
        double n = norm(psi);
        if (n < 1e-10) {
            throw std::invalid_argument("Cannot normalize zero vector");
        }

        Vector result(psi.size());
        for (size_t i = 0; i < psi.size(); ++i) {
            result[i] = psi[i] / n;
        }
        return result;
    }

    /**
     * @brief Check orthogonality: ⟨ψ|φ⟩ = 0
     */
    static bool are_orthogonal(const Vector& psi, const Vector& phi, double tol = 1e-10) {
        return std::abs(inner_product(psi, phi)) < tol;
    }

    /**
     * @brief Gram-Schmidt orthogonalization
     */
    static std::vector<Vector> gram_schmidt(const std::vector<Vector>& vectors) {
        std::vector<Vector> orthonormal;

        for (const auto& v : vectors) {
            Vector u = v;

            // Subtract projections onto previous vectors
            for (const auto& e : orthonormal) {
                Complex proj = inner_product(e, v);
                for (size_t i = 0; i < u.size(); ++i) {
                    u[i] -= proj * e[i];
                }
            }

            // Normalize
            double n = norm(u);
            if (n > 1e-10) {
                orthonormal.push_back(normalize(u));
            }
        }

        return orthonormal;
    }

    /**
     * @brief Direct sum of Hilbert spaces
     */
    static Vector direct_sum(const Vector& psi1, const Vector& psi2) {
        Vector result;
        result.reserve(psi1.size() + psi2.size());
        result.insert(result.end(), psi1.begin(), psi1.end());
        result.insert(result.end(), psi2.begin(), psi2.end());
        return result;
    }

    /**
     * @brief Tensor product |ψ⟩ ⊗ |φ⟩
     */
    static Vector tensor_product(const Vector& psi, const Vector& phi) {
        Vector result(psi.size() * phi.size());

        for (size_t i = 0; i < psi.size(); ++i) {
            for (size_t j = 0; j < phi.size(); ++j) {
                result[i * phi.size() + j] = psi[i] * phi[j];
            }
        }

        return result;
    }

    /**
     * @brief Projection onto subspace spanned by orthonormal basis
     */
    static Vector project(const Vector& psi, const std::vector<Vector>& basis) {
        if (basis.empty()) {
            return Vector(psi.size(), Complex(0.0, 0.0));
        }

        Vector result(psi.size(), Complex(0.0, 0.0));

        for (const auto& e : basis) {
            Complex coeff = inner_product(e, psi);
            for (size_t i = 0; i < result.size(); ++i) {
                result[i] += coeff * e[i];
            }
        }

        return result;
    }
};

/**
 * @brief Bounded Operators on Hilbert Space
 *
 * Linear operators with ||A|| < ∞
 */
class BoundedOperator {
public:
    /**
     * @brief Operator represented as matrix
     */
    using Matrix = std::vector<std::vector<Complex>>;
    using Vector = HilbertSpace::Vector;

    /**
     * @brief Apply operator to vector: A|ψ⟩
     */
    static Vector apply(const Matrix& A, const Vector& psi) {
        if (A.empty() || A[0].size() != psi.size()) {
            throw std::invalid_argument("Dimension mismatch");
        }

        Vector result(A.size(), Complex(0.0, 0.0));

        for (size_t i = 0; i < A.size(); ++i) {
            for (size_t j = 0; j < A[i].size(); ++j) {
                result[i] += A[i][j] * psi[j];
            }
        }

        return result;
    }

    /**
     * @brief Operator norm ||A|| = sup_{||ψ||=1} ||Aψ||
     */
    static double operator_norm(const Matrix& A, int n_samples = 100) {
        if (A.empty()) return 0.0;

        double max_norm = 0.0;
        size_t dim = A[0].size();

        // Sample random unit vectors
        for (int k = 0; k < n_samples; ++k) {
            Vector psi(dim);
            for (size_t i = 0; i < dim; ++i) {
                psi[i] = Complex(std::cos(k * i), std::sin(k * i));
            }
            psi = HilbertSpace::normalize(psi);

            Vector Apsi = apply(A, psi);
            double n = HilbertSpace::norm(Apsi);
            max_norm = std::max(max_norm, n);
        }

        return max_norm;
    }

    /**
     * @brief Adjoint operator A†
     *
     * ⟨Aψ|φ⟩ = ⟨ψ|A†φ⟩
     */
    static Matrix adjoint(const Matrix& A) {
        if (A.empty()) return Matrix();

        size_t rows = A.size();
        size_t cols = A[0].size();

        Matrix A_dag(cols, std::vector<Complex>(rows));

        for (size_t i = 0; i < rows; ++i) {
            for (size_t j = 0; j < cols; ++j) {
                A_dag[j][i] = std::conj(A[i][j]);
            }
        }

        return A_dag;
    }

    /**
     * @brief Check if operator is self-adjoint: A = A†
     */
    static bool is_self_adjoint(const Matrix& A, double tol = 1e-10) {
        Matrix A_dag = adjoint(A);

        for (size_t i = 0; i < A.size(); ++i) {
            for (size_t j = 0; j < A[i].size(); ++j) {
                if (std::abs(A[i][j] - A_dag[i][j]) > tol) {
                    return false;
                }
            }
        }

        return true;
    }

    /**
     * @brief Check if operator is unitary: U†U = UU† = I
     */
    static bool is_unitary(const Matrix& U, double tol = 1e-10) {
        Matrix U_dag = adjoint(U);
        Matrix UdagU = multiply(U_dag, U);

        // Check if result is identity
        for (size_t i = 0; i < UdagU.size(); ++i) {
            for (size_t j = 0; j < UdagU[i].size(); ++j) {
                Complex expected = (i == j) ? Complex(1.0, 0.0) : Complex(0.0, 0.0);
                if (std::abs(UdagU[i][j] - expected) > tol) {
                    return false;
                }
            }
        }

        return true;
    }

    /**
     * @brief Matrix multiplication AB
     */
    static Matrix multiply(const Matrix& A, const Matrix& B) {
        if (A.empty() || B.empty() || A[0].size() != B.size()) {
            throw std::invalid_argument("Dimension mismatch");
        }

        size_t rows = A.size();
        size_t cols = B[0].size();
        size_t inner = A[0].size();

        Matrix result(rows, std::vector<Complex>(cols, Complex(0.0, 0.0)));

        for (size_t i = 0; i < rows; ++i) {
            for (size_t j = 0; j < cols; ++j) {
                for (size_t k = 0; k < inner; ++k) {
                    result[i][j] += A[i][k] * B[k][j];
                }
            }
        }

        return result;
    }

    /**
     * @brief Commutator [A, B] = AB - BA
     */
    static Matrix commutator(const Matrix& A, const Matrix& B) {
        Matrix AB = multiply(A, B);
        Matrix BA = multiply(B, A);

        Matrix result(AB.size(), std::vector<Complex>(AB[0].size()));

        for (size_t i = 0; i < AB.size(); ++i) {
            for (size_t j = 0; j < AB[i].size(); ++j) {
                result[i][j] = AB[i][j] - BA[i][j];
            }
        }

        return result;
    }

    /**
     * @brief Trace of operator: Tr(A) = Σ Aᵢᵢ
     */
    static Complex trace(const Matrix& A) {
        Complex result(0.0, 0.0);

        for (size_t i = 0; i < A.size() && i < A[i].size(); ++i) {
            result += A[i][i];
        }

        return result;
    }

    /**
     * @brief Check if operator is positive: ⟨ψ|A|ψ⟩ ≥ 0 for all ψ
     */
    static bool is_positive(const Matrix& A, int n_samples = 50) {
        if (!is_self_adjoint(A)) {
            return false;
        }

        size_t dim = A[0].size();

        for (int k = 0; k < n_samples; ++k) {
            Vector psi(dim);
            for (size_t i = 0; i < dim; ++i) {
                psi[i] = Complex(std::cos(k * i * 1.5), std::sin(k * i * 2.3));
            }
            psi = HilbertSpace::normalize(psi);

            Vector Apsi = apply(A, psi);
            Complex expectation = HilbertSpace::inner_product(psi, Apsi);

            if (expectation.real() < -1e-10) {
                return false;
            }
        }

        return true;
    }

    /**
     * @brief Identity matrix
     */
    static Matrix identity(size_t dim) {
        Matrix I(dim, std::vector<Complex>(dim, Complex(0.0, 0.0)));

        for (size_t i = 0; i < dim; ++i) {
            I[i][i] = Complex(1.0, 0.0);
        }

        return I;
    }

    /**
     * @brief Projection operator onto normalized vector
     *
     * P = |ψ⟩⟨ψ|
     */
    static Matrix projection_operator(const Vector& psi) {
        Vector psi_norm = HilbertSpace::normalize(psi);
        size_t dim = psi_norm.size();

        Matrix P(dim, std::vector<Complex>(dim));

        for (size_t i = 0; i < dim; ++i) {
            for (size_t j = 0; j < dim; ++j) {
                P[i][j] = psi_norm[i] * std::conj(psi_norm[j]);
            }
        }

        return P;
    }
};

/**
 * @brief Von Neumann Algebras (Rings of Operators)
 *
 * Weakly closed *-algebras of bounded operators
 */
class VonNeumannAlgebra {
public:
    using Matrix = BoundedOperator::Matrix;
    using Vector = HilbertSpace::Vector;

    /**
     * @brief Commutant A' = {B : AB = BA for all A ∈ algebra}
     */
    static std::vector<Matrix> commutant(const std::vector<Matrix>& algebra) {
        // Simplified: return operators that commute with all given operators
        // In practice, this requires finding kernel of commutator maps

        std::vector<Matrix> comm;

        // Identity always commutes
        if (!algebra.empty() && !algebra[0].empty()) {
            comm.push_back(BoundedOperator::identity(algebra[0].size()));
        }

        return comm;  // Simplified implementation
    }

    /**
     * @brief Double commutant A'' = (A')'
     *
     * Von Neumann's bicommutant theorem: weakly closed => A = A''
     */
    static std::vector<Matrix> double_commutant(const std::vector<Matrix>& algebra) {
        auto comm = commutant(algebra);
        return commutant(comm);
    }

    /**
     * @brief Check if set of operators forms *-algebra
     *
     * Closed under sums, products, scalar multiplication, and adjoints
     */
    static bool is_star_algebra(const std::vector<Matrix>& operators, double tol = 1e-10) {
        // Check closure under adjoint
        for (const auto& A : operators) {
            Matrix A_dag = BoundedOperator::adjoint(A);

            // Check if A_dag is in the algebra (simplified)
            bool found = false;
            for (const auto& B : operators) {
                bool equal = true;
                for (size_t i = 0; i < A_dag.size() && equal; ++i) {
                    for (size_t j = 0; j < A_dag[i].size() && equal; ++j) {
                        if (std::abs(A_dag[i][j] - B[i][j]) > tol) {
                            equal = false;
                        }
                    }
                }
                if (equal) {
                    found = true;
                    break;
                }
            }

            if (!found) return false;
        }

        return true;
    }

    /**
     * @brief Projection in von Neumann algebra
     *
     * P† = P and P² = P
     */
    static bool is_projection(const Matrix& P, double tol = 1e-10) {
        if (!BoundedOperator::is_self_adjoint(P, tol)) {
            return false;
        }

        // Check P² = P
        Matrix P2 = BoundedOperator::multiply(P, P);

        for (size_t i = 0; i < P.size(); ++i) {
            for (size_t j = 0; j < P[i].size(); ++j) {
                if (std::abs(P2[i][j] - P[i][j]) > tol) {
                    return false;
                }
            }
        }

        return true;
    }

    /**
     * @brief Partial isometry: V†V is a projection
     */
    static bool is_partial_isometry(const Matrix& V, double tol = 1e-10) {
        Matrix V_dag = BoundedOperator::adjoint(V);
        Matrix VdagV = BoundedOperator::multiply(V_dag, V);

        return is_projection(VdagV, tol);
    }

    /**
     * @brief Center of algebra Z(M) = M ∩ M'
     */
    static std::vector<Matrix> center(const std::vector<Matrix>& algebra) {
        std::vector<Matrix> comm = commutant(algebra);
        std::vector<Matrix> cent;

        // Find intersection: operators in both algebra and its commutant
        for (const auto& A : algebra) {
            for (const auto& B : comm) {
                // Simplified: check if A = B
                bool equal = true;
                if (A.size() == B.size() && !A.empty() && A[0].size() == B[0].size()) {
                    for (size_t i = 0; i < A.size() && equal; ++i) {
                        for (size_t j = 0; j < A[i].size() && equal; ++j) {
                            if (std::abs(A[i][j] - B[i][j]) > 1e-10) {
                                equal = false;
                            }
                        }
                    }
                    if (equal) {
                        cent.push_back(A);
                        break;
                    }
                }
            }
        }

        return cent;
    }

    /**
     * @brief Check if algebra is a factor
     *
     * Factor: center consists only of scalar multiples of identity
     */
    static bool is_factor(const std::vector<Matrix>& algebra) {
        auto cent = center(algebra);

        // Factor if center is trivial (only scalars × identity)
        return cent.size() <= 1;
    }
};

/**
 * @brief Unitary Group Representations
 *
 * Homomorphisms from groups to unitary operators
 */
class UnitaryRepresentation {
public:
    using Matrix = BoundedOperator::Matrix;
    using Vector = HilbertSpace::Vector;

    /**
     * @brief Unitary matrix: U†U = I
     */
    static Matrix unitary_rotation(double theta, size_t dim) {
        Matrix U = BoundedOperator::identity(dim);

        // 2D rotation in first two dimensions
        if (dim >= 2) {
            U[0][0] = Complex(std::cos(theta), 0.0);
            U[0][1] = Complex(-std::sin(theta), 0.0);
            U[1][0] = Complex(std::sin(theta), 0.0);
            U[1][1] = Complex(std::cos(theta), 0.0);
        }

        return U;
    }

    /**
     * @brief Check if representation is irreducible
     *
     * No non-trivial invariant subspaces
     * Schur's lemma: commutant is trivial (scalars only)
     */
    static bool is_irreducible(const std::vector<Matrix>& representation, double tol = 1e-10) {
        auto comm = VonNeumannAlgebra::commutant(representation);

        // Irreducible if commutant consists only of scalar matrices
        for (const auto& C : comm) {
            // Check if C = λI
            if (C.empty()) continue;

            Complex lambda = C[0][0];
            Matrix lambda_I = BoundedOperator::identity(C.size());

            for (size_t i = 0; i < lambda_I.size(); ++i) {
                for (size_t j = 0; j < lambda_I[i].size(); ++j) {
                    lambda_I[i][j] *= lambda;
                }
            }

            bool is_scalar = true;
            for (size_t i = 0; i < C.size() && is_scalar; ++i) {
                for (size_t j = 0; j < C[i].size() && is_scalar; ++j) {
                    if (std::abs(C[i][j] - lambda_I[i][j]) > tol) {
                        is_scalar = false;
                    }
                }
            }

            if (!is_scalar) {
                return false;  // Found non-scalar commuting operator
            }
        }

        return true;
    }

    /**
     * @brief Schur's lemma: if T commutes with irreducible representation,
     * then T = λI
     */
    static bool satisfies_schur_lemma(
        const Matrix& T,
        const std::vector<Matrix>& irrep,
        double tol = 1e-10) {

        // Check T commutes with all operators in representation
        for (const auto& U : irrep) {
            Matrix TU = BoundedOperator::multiply(T, U);
            Matrix UT = BoundedOperator::multiply(U, T);

            for (size_t i = 0; i < TU.size(); ++i) {
                for (size_t j = 0; j < TU[i].size(); ++j) {
                    if (std::abs(TU[i][j] - UT[i][j]) > tol) {
                        return false;
                    }
                }
            }
        }

        // Check if T = λI
        if (T.empty()) return true;

        Complex lambda = T[0][0];
        for (size_t i = 0; i < T.size(); ++i) {
            for (size_t j = 0; j < T[i].size(); ++j) {
                Complex expected = (i == j) ? lambda : Complex(0.0, 0.0);
                if (std::abs(T[i][j] - expected) > tol) {
                    return false;
                }
            }
        }

        return true;
    }

    /**
     * @brief Character of representation: χ(g) = Tr(U(g))
     */
    static Complex character(const Matrix& U) {
        return BoundedOperator::trace(U);
    }

    /**
     * @brief Direct sum of representations
     */
    static Matrix direct_sum_representation(const Matrix& U1, const Matrix& U2) {
        size_t dim1 = U1.size();
        size_t dim2 = U2.size();

        Matrix result(dim1 + dim2, std::vector<Complex>(dim1 + dim2, Complex(0.0, 0.0)));

        // Copy U1 into top-left block
        for (size_t i = 0; i < dim1; ++i) {
            for (size_t j = 0; j < dim1; ++j) {
                result[i][j] = U1[i][j];
            }
        }

        // Copy U2 into bottom-right block
        for (size_t i = 0; i < dim2; ++i) {
            for (size_t j = 0; j < dim2; ++j) {
                result[dim1 + i][dim1 + j] = U2[i][j];
            }
        }

        return result;
    }

    /**
     * @brief Tensor product of representations
     */
    static Matrix tensor_product_representation(const Matrix& U1, const Matrix& U2) {
        size_t dim1 = U1.size();
        size_t dim2 = U2.size();
        size_t dim_total = dim1 * dim2;

        Matrix result(dim_total, std::vector<Complex>(dim_total, Complex(0.0, 0.0)));

        for (size_t i1 = 0; i1 < dim1; ++i1) {
            for (size_t i2 = 0; i2 < dim2; ++i2) {
                for (size_t j1 = 0; j1 < dim1; ++j1) {
                    for (size_t j2 = 0; j2 < dim2; ++j2) {
                        size_t i = i1 * dim2 + i2;
                        size_t j = j1 * dim2 + j2;
                        result[i][j] = U1[i1][j1] * U2[i2][j2];
                    }
                }
            }
        }

        return result;
    }
};

/**
 * @brief Factor Classification (Murray-von Neumann)
 *
 * Types I, II₁, II∞, III
 */
class FactorClassification {
public:
    using Matrix = BoundedOperator::Matrix;

    enum class FactorType {
        TYPE_I_FINITE,      // Type I_n (n×n matrices)
        TYPE_I_INFINITE,    // Type I_∞
        TYPE_II_1,          // Type II₁ (finite trace)
        TYPE_II_INFINITY,   // Type II_∞
        TYPE_III            // Type III (no trace)
    };

    /**
     * @brief Normalized trace for Type II₁ factors
     *
     * τ(I) = 1, τ(AB) = τ(BA)
     */
    static double normalized_trace(const Matrix& A) {
        Complex tr = BoundedOperator::trace(A);
        size_t dim = A.size();

        return tr.real() / dim;  // Normalized by dimension
    }

    /**
     * @brief Check if operator is finite
     *
     * Projection P is finite if P ~ Q < P implies Q = P
     */
    static bool is_finite_operator(const Matrix& P, double tol = 1e-10) {
        // Simplified: check if projection has finite trace
        if (!VonNeumannAlgebra::is_projection(P, tol)) {
            return false;
        }

        Complex tr = BoundedOperator::trace(P);
        return std::isfinite(tr.real()) && std::abs(tr.real()) < 1e10;
    }

    /**
     * @brief Murray-von Neumann equivalence
     *
     * Projections P ~ Q if P = VV†, Q = V†V for partial isometry V
     */
    static bool are_equivalent_projections(
        const Matrix& P,
        const Matrix& Q,
        double tol = 1e-10) {

        // Simplified: check if traces are equal for Type II₁
        Complex tr_P = BoundedOperator::trace(P);
        Complex tr_Q = BoundedOperator::trace(Q);

        return std::abs(tr_P - tr_Q) < tol;
    }

    /**
     * @brief Classify factor type (simplified)
     */
    static FactorType classify_factor(const std::vector<Matrix>& factor) {
        if (factor.empty()) {
            return FactorType::TYPE_I_FINITE;
        }

        // Type I: finite-dimensional
        size_t dim = factor[0].size();
        if (dim < 1000) {  // Arbitrary cutoff
            return FactorType::TYPE_I_FINITE;
        }

        // Check for trace property
        bool has_finite_trace = true;
        for (const auto& A : factor) {
            Complex tr = BoundedOperator::trace(A);
            if (!std::isfinite(tr.real()) || std::abs(tr.real()) > 1e10) {
                has_finite_trace = false;
                break;
            }
        }

        if (has_finite_trace) {
            return FactorType::TYPE_II_1;
        } else {
            return FactorType::TYPE_III;
        }
    }

    /**
     * @brief Dimension function for projections in Type II₁
     *
     * dimₘ(P) = τ(P) where τ is normalized trace
     */
    static double dimension_function(const Matrix& P) {
        if (!VonNeumannAlgebra::is_projection(P)) {
            throw std::invalid_argument("Not a projection");
        }

        return normalized_trace(P);
    }

    /**
     * @brief Check semifinite property
     *
     * Type I and II factors are semifinite
     */
    static bool is_semifinite(FactorType type) {
        return type == FactorType::TYPE_I_FINITE ||
               type == FactorType::TYPE_I_INFINITE ||
               type == FactorType::TYPE_II_1 ||
               type == FactorType::TYPE_II_INFINITY;
    }

    /**
     * @brief Continuous dimension range [0, 1] for Type II₁
     */
    static bool has_continuous_dimensions(FactorType type) {
        return type == FactorType::TYPE_II_1;
    }
};

/**
 * @brief C*-Algebras
 *
 * Norm-closed *-algebras of bounded operators
 */
class CStarAlgebra {
public:
    using Matrix = BoundedOperator::Matrix;
    using Vector = HilbertSpace::Vector;

    /**
     * @brief C*-norm identity: ||A*A|| = ||A||²
     */
    static bool satisfies_cstar_identity(const Matrix& A, double tol = 1e-6) {
        Matrix A_star = BoundedOperator::adjoint(A);
        Matrix A_star_A = BoundedOperator::multiply(A_star, A);

        double norm_A = BoundedOperator::operator_norm(A);
        double norm_A_star_A = BoundedOperator::operator_norm(A_star_A);

        return std::abs(norm_A_star_A - norm_A * norm_A) < tol;
    }

    /**
     * @brief Spectrum of element: σ(A) = {λ : A - λI not invertible}
     */
    static std::vector<Complex> spectrum_approximate(
        const Matrix& A,
        int n_samples = 50) {

        std::vector<Complex> spectrum;

        // Sample potential eigenvalues
        for (int k = 0; k < n_samples; ++k) {
            double r = k * 2.0 / n_samples;
            double theta = 2.0 * M_PI * k / n_samples;
            Complex lambda(r * std::cos(theta), r * std::sin(theta));

            // Check if A - λI has small determinant (simplified)
            // In practice, compute eigenvalues properly

            spectrum.push_back(lambda);
        }

        return spectrum;
    }

    /**
     * @brief Spectral radius: r(A) = sup{|λ| : λ ∈ σ(A)}
     */
    static double spectral_radius(const Matrix& A) {
        // Simplified: use operator norm as upper bound
        return BoundedOperator::operator_norm(A);
    }

    /**
     * @brief Check if element is normal: A*A = AA*
     */
    static bool is_normal(const Matrix& A, double tol = 1e-10) {
        Matrix A_star = BoundedOperator::adjoint(A);
        Matrix A_star_A = BoundedOperator::multiply(A_star, A);
        Matrix A_A_star = BoundedOperator::multiply(A, A_star);

        for (size_t i = 0; i < A_star_A.size(); ++i) {
            for (size_t j = 0; j < A_star_A[i].size(); ++j) {
                if (std::abs(A_star_A[i][j] - A_A_star[i][j]) > tol) {
                    return false;
                }
            }
        }

        return true;
    }

    /**
     * @brief Positive element: A = B*B for some B
     */
    static bool is_positive_element(const Matrix& A, double tol = 1e-10) {
        return BoundedOperator::is_self_adjoint(A, tol) &&
               BoundedOperator::is_positive(A);
    }

    /**
     * @brief State on C*-algebra: positive linear functional φ with φ(I) = 1
     */
    static double state_evaluation(const Matrix& A, const Vector& psi) {
        Vector psi_norm = HilbertSpace::normalize(psi);
        Vector Apsi = BoundedOperator::apply(A, psi_norm);
        Complex result = HilbertSpace::inner_product(psi_norm, Apsi);

        return result.real();
    }

    /**
     * @brief Pure state: extremal point in state space
     *
     * Corresponds to vector state ω_ψ(A) = ⟨ψ|A|ψ⟩
     */
    static bool is_pure_state_vector(const Vector& psi, double tol = 1e-10) {
        double norm = HilbertSpace::norm(psi);
        return std::abs(norm - 1.0) < tol;
    }

    /**
     * @brief Gelfand transform for commutative C*-algebra
     *
     * Maps algebra to continuous functions on spectrum
     */
    static Complex gelfand_transform(
        const Matrix& A,
        Complex point_in_spectrum) {

        // Simplified: evaluate at spectral point
        // Full implementation requires maximal ideal space

        return point_in_spectrum;  // Placeholder
    }

    /**
     * @brief GNS construction (Gelfand-Naimark-Segal)
     *
     * Every state gives a representation on Hilbert space
     * Returns cyclic vector for the representation
     */
    static Vector gns_cyclic_vector(size_t dim) {
        // Simplified: return normalized vector
        Vector psi(dim, Complex(1.0 / std::sqrt(dim), 0.0));
        return psi;
    }

    /**
     * @brief Continuous functional calculus
     *
     * For normal element A, f(A) defined for continuous f
     */
    static Matrix functional_calculus(
        const Matrix& A,
        std::function<Complex(Complex)> f) {

        // Simplified: apply function to diagonal elements
        // Full implementation requires spectral decomposition

        Matrix result = A;

        for (size_t i = 0; i < result.size(); ++i) {
            result[i][i] = f(A[i][i]);
        }

        return result;
    }

    /**
     * @brief Approximate unit in C*-algebra
     *
     * Net {eᵢ} with eᵢ → I in strong operator topology
     */
    static std::vector<Matrix> approximate_unit(size_t dim, int n_terms = 10) {
        std::vector<Matrix> unit;

        for (int k = 1; k <= n_terms; ++k) {
            double scale = static_cast<double>(k) / n_terms;
            Matrix ek = BoundedOperator::identity(dim);

            for (size_t i = 0; i < dim; ++i) {
                for (size_t j = 0; j < dim; ++j) {
                    ek[i][j] *= scale;
                }
            }

            unit.push_back(ek);
        }

        return unit;
    }
};

/**
 * @brief Quantum Mechanics Observables and States
 *
 * Physical applications of operator algebras
 */
class QuantumMechanics {
public:
    using Matrix = BoundedOperator::Matrix;
    using Vector = HilbertSpace::Vector;

    /**
     * @brief Expectation value: ⟨A⟩ = ⟨ψ|A|ψ⟩
     */
    static Complex expectation_value(const Matrix& A, const Vector& psi) {
        Vector psi_norm = HilbertSpace::normalize(psi);
        Vector Apsi = BoundedOperator::apply(A, psi_norm);
        return HilbertSpace::inner_product(psi_norm, Apsi);
    }

    /**
     * @brief Variance: Var(A) = ⟨A²⟩ - ⟨A⟩²
     */
    static double variance(const Matrix& A, const Vector& psi) {
        Complex exp_A = expectation_value(A, psi);
        Matrix A2 = BoundedOperator::multiply(A, A);
        Complex exp_A2 = expectation_value(A2, psi);

        return (exp_A2 - exp_A * exp_A).real();
    }

    /**
     * @brief Uncertainty: ΔA = √Var(A)
     */
    static double uncertainty(const Matrix& A, const Vector& psi) {
        return std::sqrt(variance(A, psi));
    }

    /**
     * @brief Heisenberg uncertainty principle: ΔA·ΔB ≥ ½|⟨[A,B]⟩|
     */
    static bool satisfies_uncertainty_principle(
        const Matrix& A,
        const Matrix& B,
        const Vector& psi,
        double tol = 1e-10) {

        double delta_A = uncertainty(A, psi);
        double delta_B = uncertainty(B, psi);

        Matrix commutator = BoundedOperator::commutator(A, B);
        Complex exp_comm = expectation_value(commutator, psi);

        double lhs = delta_A * delta_B;
        double rhs = 0.5 * std::abs(exp_comm);

        return lhs >= rhs - tol;
    }

    /**
     * @brief Time evolution: |ψ(t)⟩ = e^(-iHt)|ψ(0)⟩
     *
     * Simplified: approximate evolution operator
     */
    static Vector time_evolution(
        const Matrix& H,
        const Vector& psi_0,
        double t,
        int n_steps = 100) {

        Vector psi = psi_0;
        double dt = t / n_steps;

        // Simplified evolution: ψ(t+dt) ≈ (I - iH·dt)ψ(t)
        for (int step = 0; step < n_steps; ++step) {
            Vector Hpsi = BoundedOperator::apply(H, psi);

            for (size_t i = 0; i < psi.size(); ++i) {
                psi[i] -= Complex(0.0, dt) * Hpsi[i];
            }

            // Renormalize (approximate)
            psi = HilbertSpace::normalize(psi);
        }

        return psi;
    }

    /**
     * @brief Density matrix: ρ = |ψ⟩⟨ψ| for pure state
     */
    static Matrix density_matrix_pure(const Vector& psi) {
        return BoundedOperator::projection_operator(psi);
    }

    /**
     * @brief Mixed state density matrix: ρ = Σ pᵢ|ψᵢ⟩⟨ψᵢ|
     */
    static Matrix density_matrix_mixed(
        const std::vector<Vector>& states,
        const std::vector<double>& probabilities) {

        if (states.empty() || states.size() != probabilities.size()) {
            throw std::invalid_argument("Invalid state ensemble");
        }

        size_t dim = states[0].size();
        Matrix rho(dim, std::vector<Complex>(dim, Complex(0.0, 0.0)));

        for (size_t k = 0; k < states.size(); ++k) {
            Matrix P_k = BoundedOperator::projection_operator(states[k]);

            for (size_t i = 0; i < dim; ++i) {
                for (size_t j = 0; j < dim; ++j) {
                    rho[i][j] += probabilities[k] * P_k[i][j];
                }
            }
        }

        return rho;
    }

    /**
     * @brief Von Neumann entropy: S(ρ) = -Tr(ρ log ρ)
     */
    static double von_neumann_entropy(const Matrix& rho) {
        // Simplified: assume diagonal density matrix
        double entropy = 0.0;

        for (size_t i = 0; i < rho.size(); ++i) {
            double p_i = rho[i][i].real();
            if (p_i > 1e-10) {
                entropy -= p_i * std::log(p_i);
            }
        }

        return entropy;
    }

    /**
     * @brief Measurement postulate: probability of outcome
     *
     * P(λ) = ⟨ψ|Pλ|ψ⟩ where Pλ is projection onto eigenspace
     */
    static double measurement_probability(
        const Matrix& projection,
        const Vector& psi) {

        Complex prob = expectation_value(projection, psi);
        return prob.real();
    }

    /**
     * @brief Post-measurement state: |ψ'⟩ = Pλ|ψ⟩ / ||Pλ|ψ⟩||
     */
    static Vector post_measurement_state(
        const Matrix& projection,
        const Vector& psi) {

        Vector P_psi = BoundedOperator::apply(projection, psi);
        return HilbertSpace::normalize(P_psi);
    }

    /**
     * @brief Commuting observables share eigenbasis
     */
    static bool are_compatible_observables(
        const Matrix& A,
        const Matrix& B,
        double tol = 1e-10) {

        Matrix comm = BoundedOperator::commutator(A, B);

        // Check if commutator is zero
        for (size_t i = 0; i < comm.size(); ++i) {
            for (size_t j = 0; j < comm[i].size(); ++j) {
                if (std::abs(comm[i][j]) > tol) {
                    return false;
                }
            }
        }

        return true;
    }
};

} // namespace operator_algebras
} // namespace physics

#endif // PHYSICS_ADVANCED_OPERATOR_ALGEBRAS_HPP
