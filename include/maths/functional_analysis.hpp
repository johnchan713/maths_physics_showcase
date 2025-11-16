#ifndef MATHS_FUNCTIONAL_ANALYSIS_HPP
#define MATHS_FUNCTIONAL_ANALYSIS_HPP

#include <vector>
#include <functional>
#include <cmath>
#include <complex>
#include <algorithm>
#include <numeric>
#include <limits>
#include <stdexcept>

namespace maths {
namespace functional_analysis {

using Complex = std::complex<double>;

// ============================================================================
// NORMED SPACES
// ============================================================================

template<typename T>
class NormedSpace {
public:
    virtual double norm(const std::vector<T>& x) const = 0;
    virtual ~NormedSpace() = default;

    double distance(const std::vector<T>& x, const std::vector<T>& y) const {
        if (x.size() != y.size()) {
            throw std::invalid_argument("Vectors must have same dimension");
        }
        std::vector<T> diff(x.size());
        for (size_t i = 0; i < x.size(); ++i) {
            diff[i] = x[i] - y[i];
        }
        return norm(diff);
    }

    bool isCauchy(const std::vector<std::vector<T>>& sequence,
                  double epsilon = 1e-6) const {
        int n = sequence.size();
        for (int i = n - 10; i < n; ++i) {
            for (int j = i + 1; j < n; ++j) {
                if (distance(sequence[i], sequence[j]) > epsilon) {
                    return false;
                }
            }
        }
        return true;
    }
};

template<typename T>
class LpNorm : public NormedSpace<T> {
private:
    double p_;

public:
    LpNorm(double p) : p_(p) {
        if (p < 1.0) throw std::invalid_argument("p must be >= 1");
    }

    double norm(const std::vector<T>& x) const override {
        if (std::isinf(p_)) {
            double max_val = 0.0;
            for (const auto& val : x) {
                max_val = std::max(max_val, std::abs(static_cast<double>(val)));
            }
            return max_val;
        }

        double sum = 0.0;
        for (const auto& val : x) {
            sum += std::pow(std::abs(static_cast<double>(val)), p_);
        }
        return std::pow(sum, 1.0 / p_);
    }
};

// ============================================================================
// BANACH SPACES
// ============================================================================

template<typename T>
class BanachSpace : public NormedSpace<T> {
public:
    std::vector<T> limit(const std::vector<std::vector<T>>& sequence) const {
        if (sequence.empty()) {
            throw std::invalid_argument("Empty sequence");
        }

        int n = sequence.size();
        size_t dim = sequence[0].size();

        if (n < 2) return sequence[0];

        if (!this->isCauchy(sequence)) {
            return sequence.back();
        }

        std::vector<T> lim(dim);
        for (size_t i = 0; i < dim; ++i) {
            double sum = 0.0;
            int count = std::min(10, n);
            for (int j = n - count; j < n; ++j) {
                sum += static_cast<double>(sequence[j][i]);
            }
            lim[i] = static_cast<T>(sum / count);
        }
        return lim;
    }

    bool isComplete(const std::vector<std::vector<T>>& test_sequences,
                   double tol = 1e-6) const {
        for (const auto& seq : test_sequences) {
            if (this->isCauchy({seq}) && this->norm(limit({seq})) < 1e-10) {
                return false;
            }
        }
        return true;
    }

    std::vector<T> fixedPoint(
        std::function<std::vector<T>(const std::vector<T>&)> T_map,
        const std::vector<T>& x0,
        double contraction_constant,
        int max_iter = 1000,
        double tol = 1e-8) const {

        if (contraction_constant >= 1.0) {
            throw std::invalid_argument("Not a contraction");
        }

        std::vector<T> x = x0;
        for (int iter = 0; iter < max_iter; ++iter) {
            std::vector<T> x_new = T_map(x);
            double dist = this->distance(x, x_new);
            if (dist < tol) {
                return x_new;
            }
            x = x_new;
        }
        return x;
    }
};

// ============================================================================
// HILBERT SPACES
// ============================================================================

template<typename T = double>
class HilbertSpace : public BanachSpace<T> {
public:
    virtual T innerProduct(const std::vector<T>& x,
                          const std::vector<T>& y) const = 0;

    double norm(const std::vector<T>& x) const override {
        T ip = innerProduct(x, x);
        return std::sqrt(std::abs(static_cast<double>(ip)));
    }

    bool areOrthogonal(const std::vector<T>& x, const std::vector<T>& y,
                      double tol = 1e-10) const {
        T ip = innerProduct(x, y);
        return std::abs(static_cast<double>(ip)) < tol;
    }

    std::vector<T> project(const std::vector<T>& x,
                          const std::vector<T>& v) const {
        T ip_xv = innerProduct(x, v);
        T ip_vv = innerProduct(v, v);
        double coeff = static_cast<double>(ip_xv) / static_cast<double>(ip_vv);

        std::vector<T> proj(v.size());
        for (size_t i = 0; i < v.size(); ++i) {
            proj[i] = static_cast<T>(coeff * static_cast<double>(v[i]));
        }
        return proj;
    }

    std::vector<std::vector<T>> gramSchmidt(
        const std::vector<std::vector<T>>& basis) const {

        std::vector<std::vector<T>> orthonormal;
        for (const auto& v : basis) {
            std::vector<T> u = v;

            for (const auto& e : orthonormal) {
                auto proj = project(v, e);
                for (size_t i = 0; i < u.size(); ++i) {
                    u[i] -= proj[i];
                }
            }

            double n = norm(u);
            if (n > 1e-10) {
                for (auto& component : u) {
                    component = static_cast<T>(static_cast<double>(component) / n);
                }
                orthonormal.push_back(u);
            }
        }
        return orthonormal;
    }
};

class RealHilbertSpace : public HilbertSpace<double> {
public:
    double innerProduct(const std::vector<double>& x,
                       const std::vector<double>& y) const override {
        if (x.size() != y.size()) {
            throw std::invalid_argument("Dimension mismatch");
        }
        return std::inner_product(x.begin(), x.end(), y.begin(), 0.0);
    }
};

class ComplexHilbertSpace : public HilbertSpace<Complex> {
public:
    Complex innerProduct(const std::vector<Complex>& x,
                        const std::vector<Complex>& y) const override {
        if (x.size() != y.size()) {
            throw std::invalid_argument("Dimension mismatch");
        }
        Complex sum = 0.0;
        for (size_t i = 0; i < x.size(); ++i) {
            sum += std::conj(x[i]) * y[i];
        }
        return sum;
    }
};

// ============================================================================
// LINEAR OPERATORS
// ============================================================================

template<typename T>
class LinearOperator {
protected:
    size_t dim_;

public:
    LinearOperator(size_t dim) : dim_(dim) {}
    virtual ~LinearOperator() = default;

    virtual std::vector<T> apply(const std::vector<T>& x) const = 0;

    size_t dimension() const { return dim_; }

    double operatorNorm(const NormedSpace<T>& space, int samples = 1000) const {
        double max_ratio = 0.0;
        for (int i = 0; i < samples; ++i) {
            std::vector<T> x(dim_);
            for (size_t j = 0; j < dim_; ++j) {
                x[j] = static_cast<T>(static_cast<double>(rand()) / RAND_MAX - 0.5);
            }
            double norm_x = space.norm(x);
            if (norm_x > 1e-10) {
                auto Tx = apply(x);
                double norm_Tx = space.norm(Tx);
                max_ratio = std::max(max_ratio, norm_Tx / norm_x);
            }
        }
        return max_ratio;
    }

    bool isBounded(const NormedSpace<T>& space, double tol = 1e6) const {
        return operatorNorm(space) < tol;
    }
};

template<typename T>
class MatrixOperator : public LinearOperator<T> {
private:
    std::vector<std::vector<T>> matrix_;

public:
    MatrixOperator(const std::vector<std::vector<T>>& A)
        : LinearOperator<T>(A.size()), matrix_(A) {}

    std::vector<T> apply(const std::vector<T>& x) const override {
        if (x.size() != this->dim_) {
            throw std::invalid_argument("Dimension mismatch");
        }
        std::vector<T> result(this->dim_, 0);
        for (size_t i = 0; i < this->dim_; ++i) {
            for (size_t j = 0; j < this->dim_; ++j) {
                result[i] += matrix_[i][j] * x[j];
            }
        }
        return result;
    }

    const std::vector<std::vector<T>>& getMatrix() const { return matrix_; }

    bool isSelfAdjoint(double tol = 1e-10) const {
        for (size_t i = 0; i < this->dim_; ++i) {
            for (size_t j = 0; j < this->dim_; ++j) {
                T diff = matrix_[i][j] - matrix_[j][i];
                if (std::abs(static_cast<double>(diff)) > tol) {
                    return false;
                }
            }
        }
        return true;
    }

    bool isCompact() const {
        return true;
    }
};

// ============================================================================
// SPECTRAL THEORY
// ============================================================================

template<typename T>
class SpectralTheory {
public:
    struct Eigenvalue {
        T value;
        std::vector<T> eigenvector;
        int multiplicity;
    };

    static std::vector<Eigenvalue> powerIteration(
        const MatrixOperator<T>& A,
        int max_eigenvalues = 5,
        int max_iter = 1000,
        double tol = 1e-8) {

        std::vector<Eigenvalue> eigenvalues;
        size_t n = A.dimension();
        std::vector<T> v(n);

        for (int k = 0; k < max_eigenvalues; ++k) {
            for (size_t i = 0; i < n; ++i) {
                v[i] = static_cast<T>(static_cast<double>(rand()) / RAND_MAX);
            }

            for (const auto& prev_ev : eigenvalues) {
                T ip = 0;
                for (size_t i = 0; i < n; ++i) {
                    ip += v[i] * prev_ev.eigenvector[i];
                }
                for (size_t i = 0; i < n; ++i) {
                    v[i] -= ip * prev_ev.eigenvector[i];
                }
            }

            T lambda = 0;
            T lambda_old;
            for (int iter = 0; iter < max_iter; ++iter) {
                auto Av = A.apply(v);

                T norm_sq = 0;
                for (const auto& val : Av) {
                    norm_sq += val * val;
                }
                T norm_Av = std::sqrt(norm_sq);

                for (auto& val : Av) {
                    val /= norm_Av;
                }

                lambda_old = lambda;
                lambda = 0;
                for (size_t i = 0; i < n; ++i) {
                    lambda += v[i] * A.apply(v)[i];
                }

                if (std::abs(static_cast<double>(lambda - lambda_old)) < tol) {
                    v = Av;
                    break;
                }
                v = Av;
            }

            if (std::abs(static_cast<double>(lambda)) > tol) {
                eigenvalues.push_back({lambda, v, 1});
            } else {
                break;
            }
        }

        return eigenvalues;
    }

    static double spectralRadius(const MatrixOperator<T>& A) {
        auto eigenvalues = powerIteration(A, 1);
        if (eigenvalues.empty()) return 0.0;
        return std::abs(static_cast<double>(eigenvalues[0].value));
    }

    static std::vector<T> resolvent(
        const MatrixOperator<T>& A,
        T lambda,
        const std::vector<T>& x) {

        size_t n = A.dimension();
        auto A_matrix = A.getMatrix();

        std::vector<std::vector<T>> A_lambda(n, std::vector<T>(n));
        for (size_t i = 0; i < n; ++i) {
            for (size_t j = 0; j < n; ++j) {
                A_lambda[i][j] = A_matrix[i][j];
                if (i == j) {
                    A_lambda[i][j] -= lambda;
                }
            }
        }

        std::vector<T> b = x;
        std::vector<T> y(n, 0);

        for (int iter = 0; iter < 100; ++iter) {
            for (size_t i = 0; i < n; ++i) {
                T sum = b[i];
                for (size_t j = 0; j < n; ++j) {
                    if (i != j) {
                        sum -= A_lambda[i][j] * y[j];
                    }
                }
                y[i] = sum / A_lambda[i][i];
            }
        }

        return y;
    }
};

// ============================================================================
// DUAL SPACES
// ============================================================================

template<typename T>
class LinearFunctional {
protected:
    size_t dim_;

public:
    LinearFunctional(size_t dim) : dim_(dim) {}
    virtual ~LinearFunctional() = default;

    virtual T apply(const std::vector<T>& x) const = 0;

    T operator()(const std::vector<T>& x) const {
        return apply(x);
    }

    double norm(const NormedSpace<T>& space, int samples = 1000) const {
        double max_val = 0.0;
        for (int i = 0; i < samples; ++i) {
            std::vector<T> x(dim_);
            for (size_t j = 0; j < dim_; ++j) {
                x[j] = static_cast<T>(static_cast<double>(rand()) / RAND_MAX);
            }
            double norm_x = space.norm(x);
            if (norm_x > 1e-10) {
                T val = apply(x);
                max_val = std::max(max_val, std::abs(static_cast<double>(val)) / norm_x);
            }
        }
        return max_val;
    }
};

template<typename T>
class RieszRepresentation : public LinearFunctional<T> {
private:
    std::vector<T> representative_;
    const HilbertSpace<T>& space_;

public:
    RieszRepresentation(const std::vector<T>& y, const HilbertSpace<T>& H)
        : LinearFunctional<T>(y.size()), representative_(y), space_(H) {}

    T apply(const std::vector<T>& x) const override {
        return space_.innerProduct(representative_, x);
    }

    const std::vector<T>& getRepresentative() const { return representative_; }
};

// ============================================================================
// WEAK CONVERGENCE
// ============================================================================

template<typename T>
class WeakConvergence {
private:
    const HilbertSpace<T>& space_;

public:
    WeakConvergence(const HilbertSpace<T>& H) : space_(H) {}

    bool convergesWeakly(
        const std::vector<std::vector<T>>& sequence,
        const std::vector<T>& limit,
        const std::vector<LinearFunctional<T>*>& test_functionals,
        double tol = 1e-6) const {

        for (const auto& f : test_functionals) {
            T f_limit = f->apply(limit);
            int n = sequence.size();
            T f_xn = f->apply(sequence[n-1]);

            if (std::abs(static_cast<double>(f_xn - f_limit)) > tol) {
                return false;
            }
        }
        return true;
    }

    bool isWeaklyCauchy(
        const std::vector<std::vector<T>>& sequence,
        const std::vector<LinearFunctional<T>*>& test_functionals,
        double tol = 1e-6) const {

        int n = sequence.size();
        for (const auto& f : test_functionals) {
            for (int i = n - 10; i < n; ++i) {
                for (int j = i + 1; j < n; ++j) {
                    T f_xi = f->apply(sequence[i]);
                    T f_xj = f->apply(sequence[j]);
                    if (std::abs(static_cast<double>(f_xi - f_xj)) > tol) {
                        return false;
                    }
                }
            }
        }
        return true;
    }
};

// ============================================================================
// SOBOLEV SPACES
// ============================================================================

class SobolevSpace {
private:
    int k_;
    double p_;

public:
    SobolevSpace(int k, double p) : k_(k), p_(p) {
        if (k < 0) throw std::invalid_argument("k must be non-negative");
        if (p < 1.0) throw std::invalid_argument("p must be >= 1");
    }

    double norm(std::function<double(double)> u,
                std::function<std::vector<double>(double)> derivatives,
                double a, double b,
                int n = 1000) const {

        double sum = 0.0;
        double dx = (b - a) / n;

        for (int j = 0; j <= k_; ++j) {
            for (int i = 0; i < n; ++i) {
                double x = a + i * dx;
                double val;
                if (j == 0) {
                    val = u(x);
                } else {
                    auto derivs = derivatives(x);
                    val = derivs[j-1];
                }
                sum += std::pow(std::abs(val), p_) * dx;
            }
        }

        return std::pow(sum, 1.0 / p_);
    }

    double seminorm(std::function<double(double)> u_k_th_derivative,
                   double a, double b,
                   int n = 1000) const {

        double sum = 0.0;
        double dx = (b - a) / n;

        for (int i = 0; i < n; ++i) {
            double x = a + i * dx;
            double val = u_k_th_derivative(x);
            sum += std::pow(std::abs(val), p_) * dx;
        }

        return std::pow(sum, 1.0 / p_);
    }

    bool satisfiesSobolevEmbedding(double domain_dim, double target_p) const {
        if (k_ * p_ > domain_dim) {
            return target_p == std::numeric_limits<double>::infinity();
        } else if (k_ * p_ == domain_dim) {
            return target_p < std::numeric_limits<double>::infinity();
        } else {
            double critical_p = domain_dim * p_ / (domain_dim - k_ * p_);
            return target_p <= critical_p;
        }
    }
};

// ============================================================================
// THEOREM VERIFICATION
// ============================================================================

class TheoremVerification {
public:
    template<typename T>
    static bool verifyOpenMappingTheorem(
        const LinearOperator<T>& A,
        const BanachSpace<T>& X,
        const BanachSpace<T>& Y,
        int samples = 100) {

        std::vector<std::vector<T>> test_open_sets;
        for (int i = 0; i < samples; ++i) {
            std::vector<T> center(A.dimension());
            for (size_t j = 0; j < A.dimension(); ++j) {
                center[j] = static_cast<T>(static_cast<double>(rand()) / RAND_MAX);
            }
            test_open_sets.push_back(center);
        }

        for (const auto& x : test_open_sets) {
            auto Ax = A.apply(x);
            if (Y.norm(Ax) < 1e-10) {
                return false;
            }
        }
        return true;
    }

    template<typename T>
    static bool verifyClosedGraphTheorem(
        const LinearOperator<T>& A,
        const std::vector<std::pair<std::vector<T>, std::vector<T>>>& graph_points,
        const BanachSpace<T>& X,
        const BanachSpace<T>& Y,
        double tol = 1e-6) {

        for (const auto& [x, y] : graph_points) {
            auto Ax = A.apply(x);
            if (Y.distance(Ax, y) > tol) {
                return false;
            }
        }
        return true;
    }

    static bool verifyHahnBanach(
        const std::vector<double>& subspace_functional,
        const std::vector<double>& extension_functional,
        const std::vector<double>& subspace_point,
        double tol = 1e-10) {

        double sub_val = std::inner_product(
            subspace_functional.begin(), subspace_functional.end(),
            subspace_point.begin(), 0.0);

        double ext_val = std::inner_product(
            extension_functional.begin(), extension_functional.end(),
            subspace_point.begin(), 0.0);

        return std::abs(sub_val - ext_val) < tol;
    }
};

} // namespace functional_analysis
} // namespace maths

#endif // MATHS_FUNCTIONAL_ANALYSIS_HPP
