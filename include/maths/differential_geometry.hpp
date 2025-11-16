#ifndef MATHS_DIFFERENTIAL_GEOMETRY_HPP
#define MATHS_DIFFERENTIAL_GEOMETRY_HPP

#include <vector>
#include <array>
#include <functional>
#include <cmath>
#include <algorithm>
#include <numeric>
#include <stdexcept>

namespace maths {
namespace differential_geometry {

template<int N>
using Point = std::array<double, N>;

template<int N>
using Vector = std::array<double, N>;

template<int N>
using Covector = std::array<double, N>;

template<int N>
using Tensor2 = std::array<std::array<double, N>, N>;

template<int N>
using Tensor3 = std::array<std::array<std::array<double, N>, N>, N>;

template<int N>
using Tensor4 = std::array<std::array<std::array<std::array<double, N>, N>, N>, N>;

// ============================================================================
// RIEMANNIAN METRIC
// ============================================================================

template<int N>
class RiemannianMetric {
public:
    virtual Tensor2<N> at(const Point<N>& p) const = 0;
    virtual ~RiemannianMetric() = default;

    Tensor2<N> inverse(const Point<N>& p) const {
        auto g = at(p);
        Tensor2<N> g_inv;

        if constexpr (N == 2) {
            double det = g[0][0] * g[1][1] - g[0][1] * g[1][0];
            if (std::abs(det) < 1e-15) throw std::runtime_error("Singular metric");
            g_inv[0][0] = g[1][1] / det;
            g_inv[0][1] = -g[0][1] / det;
            g_inv[1][0] = -g[1][0] / det;
            g_inv[1][1] = g[0][0] / det;
        } else if constexpr (N == 3) {
            double det = g[0][0] * (g[1][1] * g[2][2] - g[1][2] * g[2][1])
                       - g[0][1] * (g[1][0] * g[2][2] - g[1][2] * g[2][0])
                       + g[0][2] * (g[1][0] * g[2][1] - g[1][1] * g[2][0]);
            if (std::abs(det) < 1e-15) throw std::runtime_error("Singular metric");

            g_inv[0][0] = (g[1][1] * g[2][2] - g[1][2] * g[2][1]) / det;
            g_inv[0][1] = (g[0][2] * g[2][1] - g[0][1] * g[2][2]) / det;
            g_inv[0][2] = (g[0][1] * g[1][2] - g[0][2] * g[1][1]) / det;
            g_inv[1][0] = (g[1][2] * g[2][0] - g[1][0] * g[2][2]) / det;
            g_inv[1][1] = (g[0][0] * g[2][2] - g[0][2] * g[2][0]) / det;
            g_inv[1][2] = (g[0][2] * g[1][0] - g[0][0] * g[1][2]) / det;
            g_inv[2][0] = (g[1][0] * g[2][1] - g[1][1] * g[2][0]) / det;
            g_inv[2][1] = (g[0][1] * g[2][0] - g[0][0] * g[2][1]) / det;
            g_inv[2][2] = (g[0][0] * g[1][1] - g[0][1] * g[1][0]) / det;
        }
        return g_inv;
    }

    double innerProduct(const Point<N>& p, const Vector<N>& v, const Vector<N>& w) const {
        auto g = at(p);
        double result = 0.0;
        for (int i = 0; i < N; ++i) {
            for (int j = 0; j < N; ++j) {
                result += g[i][j] * v[i] * w[j];
            }
        }
        return result;
    }

    double norm(const Point<N>& p, const Vector<N>& v) const {
        return std::sqrt(innerProduct(p, v, v));
    }

    double angle(const Point<N>& p, const Vector<N>& v, const Vector<N>& w) const {
        double ip = innerProduct(p, v, w);
        double norm_v = norm(p, v);
        double norm_w = norm(p, w);
        if (norm_v < 1e-15 || norm_w < 1e-15) return 0.0;
        return std::acos(std::clamp(ip / (norm_v * norm_w), -1.0, 1.0));
    }

    Covector<N> lowerIndex(const Point<N>& p, const Vector<N>& v) const {
        auto g = at(p);
        Covector<N> omega;
        for (int i = 0; i < N; ++i) {
            omega[i] = 0.0;
            for (int j = 0; j < N; ++j) {
                omega[i] += g[i][j] * v[j];
            }
        }
        return omega;
    }

    Vector<N> raiseIndex(const Point<N>& p, const Covector<N>& omega) const {
        auto g_inv = inverse(p);
        Vector<N> v;
        for (int i = 0; i < N; ++i) {
            v[i] = 0.0;
            for (int j = 0; j < N; ++j) {
                v[i] += g_inv[i][j] * omega[j];
            }
        }
        return v;
    }
};

template<int N>
class EuclideanMetric : public RiemannianMetric<N> {
public:
    Tensor2<N> at(const Point<N>& p) const override {
        Tensor2<N> g;
        for (int i = 0; i < N; ++i) {
            for (int j = 0; j < N; ++j) {
                g[i][j] = (i == j) ? 1.0 : 0.0;
            }
        }
        return g;
    }
};

template<int N>
class HyperbolicMetric : public RiemannianMetric<N> {
public:
    Tensor2<N> at(const Point<N>& p) const override {
        double y = p[N-1];
        if (y <= 0.0) throw std::runtime_error("Invalid point for hyperbolic metric");

        Tensor2<N> g;
        for (int i = 0; i < N; ++i) {
            for (int j = 0; j < N; ++j) {
                g[i][j] = (i == j) ? 1.0 / (y * y) : 0.0;
            }
        }
        return g;
    }
};

// ============================================================================
// LEVI-CIVITA CONNECTION
// ============================================================================

template<int N>
class LeviCivitaConnection {
private:
    const RiemannianMetric<N>& metric_;
    double epsilon_;

public:
    LeviCivitaConnection(const RiemannianMetric<N>& g, double eps = 1e-6)
        : metric_(g), epsilon_(eps) {}

    double christoffel(int k, int i, int j, const Point<N>& p) const {
        auto g_inv = metric_.inverse(p);
        double Gamma = 0.0;

        for (int l = 0; l < N; ++l) {
            Gamma += 0.5 * g_inv[k][l] * (
                partialDerivative(l, i, j, p) +
                partialDerivative(l, j, i, p) -
                partialDerivative(i, j, l, p)
            );
        }
        return Gamma;
    }

    Tensor3<N> christoffelSymbols(const Point<N>& p) const {
        Tensor3<N> Gamma;
        for (int k = 0; k < N; ++k) {
            for (int i = 0; i < N; ++i) {
                for (int j = 0; j < N; ++j) {
                    Gamma[k][i][j] = christoffel(k, i, j, p);
                }
            }
        }
        return Gamma;
    }

    Vector<N> covariantDerivative(const Point<N>& p,
                                  const std::function<Vector<N>(const Point<N>&)>& X,
                                  const Vector<N>& Y) const {
        Vector<N> nabla_Y_X;
        auto Gamma = christoffelSymbols(p);

        for (int k = 0; k < N; ++k) {
            nabla_Y_X[k] = 0.0;

            for (int i = 0; i < N; ++i) {
                Point<N> p_shift = p;
                p_shift[i] += epsilon_;
                nabla_Y_X[k] += Y[i] * (X(p_shift)[k] - X(p)[k]) / epsilon_;
            }

            for (int i = 0; i < N; ++i) {
                for (int j = 0; j < N; ++j) {
                    nabla_Y_X[k] += Gamma[k][i][j] * Y[i] * X(p)[j];
                }
            }
        }
        return nabla_Y_X;
    }

    Vector<N> parallelTransport(const Vector<N>& v0,
                               const std::vector<Point<N>>& curve) const {
        if (curve.size() < 2) return v0;

        Vector<N> v = v0;
        for (size_t i = 0; i < curve.size() - 1; ++i) {
            const auto& p = curve[i];
            const auto& p_next = curve[i + 1];

            Vector<N> tangent;
            for (int k = 0; k < N; ++k) {
                tangent[k] = p_next[k] - p[k];
            }

            auto Gamma = christoffelSymbols(p);
            Vector<N> correction;
            for (int k = 0; k < N; ++k) {
                correction[k] = 0.0;
                for (int i = 0; i < N; ++i) {
                    for (int j = 0; j < N; ++j) {
                        correction[k] += Gamma[k][i][j] * tangent[i] * v[j];
                    }
                }
            }

            for (int k = 0; k < N; ++k) {
                v[k] -= correction[k];
            }
        }
        return v;
    }

private:
    double partialDerivative(int alpha, int mu, int nu, const Point<N>& p) const {
        Point<N> p_plus = p, p_minus = p;
        p_plus[alpha] += epsilon_;
        p_minus[alpha] -= epsilon_;

        auto g_plus = metric_.at(p_plus);
        auto g_minus = metric_.at(p_minus);

        return (g_plus[mu][nu] - g_minus[mu][nu]) / (2.0 * epsilon_);
    }
};

// ============================================================================
// CURVATURE TENSORS
// ============================================================================

template<int N>
class RiemannCurvatureTensor {
private:
    LeviCivitaConnection<N> connection_;
    double epsilon_;

public:
    RiemannCurvatureTensor(const LeviCivitaConnection<N>& conn, double eps = 1e-6)
        : connection_(conn), epsilon_(eps) {}

    double component(int rho, int sigma, int mu, int nu, const Point<N>& p) const {
        Point<N> p_mu_plus = p, p_mu_minus = p;
        p_mu_plus[mu] += epsilon_;
        p_mu_minus[mu] -= epsilon_;
        double dGamma_mu = (connection_.christoffel(rho, nu, sigma, p_mu_plus)
                          - connection_.christoffel(rho, nu, sigma, p_mu_minus))
                          / (2.0 * epsilon_);

        Point<N> p_nu_plus = p, p_nu_minus = p;
        p_nu_plus[nu] += epsilon_;
        p_nu_minus[nu] -= epsilon_;
        double dGamma_nu = (connection_.christoffel(rho, mu, sigma, p_nu_plus)
                          - connection_.christoffel(rho, mu, sigma, p_nu_minus))
                          / (2.0 * epsilon_);

        double quadratic = 0.0;
        for (int lambda = 0; lambda < N; ++lambda) {
            quadratic += connection_.christoffel(rho, mu, lambda, p)
                       * connection_.christoffel(lambda, nu, sigma, p);
            quadratic -= connection_.christoffel(rho, nu, lambda, p)
                       * connection_.christoffel(lambda, mu, sigma, p);
        }

        return dGamma_mu - dGamma_nu + quadratic;
    }

    Tensor4<N> allComponents(const Point<N>& p) const {
        Tensor4<N> R;
        for (int rho = 0; rho < N; ++rho) {
            for (int sigma = 0; sigma < N; ++sigma) {
                for (int mu = 0; mu < N; ++mu) {
                    for (int nu = 0; nu < N; ++nu) {
                        R[rho][sigma][mu][nu] = component(rho, sigma, mu, nu, p);
                    }
                }
            }
        }
        return R;
    }

    double sectionalCurvature(const Point<N>& p, const Vector<N>& X, const Vector<N>& Y,
                             const RiemannianMetric<N>& metric) const {
        double num = 0.0;
        for (int rho = 0; rho < N; ++rho) {
            for (int sigma = 0; sigma < N; ++sigma) {
                for (int mu = 0; mu < N; ++mu) {
                    for (int nu = 0; nu < N; ++nu) {
                        num += component(rho, sigma, mu, nu, p) * X[rho] * Y[sigma] * X[mu] * Y[nu];
                    }
                }
            }
        }

        double XX = metric.innerProduct(p, X, X);
        double YY = metric.innerProduct(p, Y, Y);
        double XY = metric.innerProduct(p, X, Y);
        double denom = XX * YY - XY * XY;

        if (std::abs(denom) < 1e-15) return 0.0;
        return num / denom;
    }
};

template<int N>
class RicciTensor {
private:
    RiemannCurvatureTensor<N> riemann_;

public:
    RicciTensor(const RiemannCurvatureTensor<N>& R) : riemann_(R) {}

    Tensor2<N> at(const Point<N>& p) const {
        Tensor2<N> Ric;
        for (int mu = 0; mu < N; ++mu) {
            for (int nu = 0; nu < N; ++nu) {
                Ric[mu][nu] = 0.0;
                for (int lambda = 0; lambda < N; ++lambda) {
                    Ric[mu][nu] += riemann_.component(lambda, mu, lambda, nu, p);
                }
            }
        }
        return Ric;
    }

    double scalar(const Point<N>& p, const RiemannianMetric<N>& metric) const {
        auto Ric = at(p);
        auto g_inv = metric.inverse(p);
        double R = 0.0;
        for (int mu = 0; mu < N; ++mu) {
            for (int nu = 0; nu < N; ++nu) {
                R += g_inv[mu][nu] * Ric[mu][nu];
            }
        }
        return R;
    }

    bool isEinstein(const Point<N>& p, const RiemannianMetric<N>& metric,
                   double tol = 1e-6) const {
        auto Ric = at(p);
        double R = scalar(p, metric);
        auto g = metric.at(p);

        for (int mu = 0; mu < N; ++mu) {
            for (int nu = 0; nu < N; ++nu) {
                double expected = R / N * g[mu][nu];
                if (std::abs(Ric[mu][nu] - expected) > tol) {
                    return false;
                }
            }
        }
        return true;
    }
};

// ============================================================================
// GEODESICS
// ============================================================================

template<int N>
class GeodesicEquation {
private:
    LeviCivitaConnection<N> connection_;

public:
    struct State {
        Point<N> position;
        Vector<N> velocity;
    };

    GeodesicEquation(const LeviCivitaConnection<N>& conn) : connection_(conn) {}

    State derivative(const State& s) const {
        State ds;
        ds.position = s.velocity;

        auto Gamma = connection_.christoffelSymbols(s.position);
        for (int k = 0; k < N; ++k) {
            ds.velocity[k] = 0.0;
            for (int i = 0; i < N; ++i) {
                for (int j = 0; j < N; ++j) {
                    ds.velocity[k] -= Gamma[k][i][j] * s.velocity[i] * s.velocity[j];
                }
            }
        }
        return ds;
    }

    std::vector<State> integrate(const State& initial, double dt, int steps) const {
        std::vector<State> trajectory;
        trajectory.push_back(initial);

        State current = initial;
        for (int step = 0; step < steps; ++step) {
            State k1 = derivative(current);

            State temp;
            for (int i = 0; i < N; ++i) {
                temp.position[i] = current.position[i] + 0.5 * dt * k1.position[i];
                temp.velocity[i] = current.velocity[i] + 0.5 * dt * k1.velocity[i];
            }
            State k2 = derivative(temp);

            for (int i = 0; i < N; ++i) {
                temp.position[i] = current.position[i] + 0.5 * dt * k2.position[i];
                temp.velocity[i] = current.velocity[i] + 0.5 * dt * k2.velocity[i];
            }
            State k3 = derivative(temp);

            for (int i = 0; i < N; ++i) {
                temp.position[i] = current.position[i] + dt * k3.position[i];
                temp.velocity[i] = current.velocity[i] + dt * k3.velocity[i];
            }
            State k4 = derivative(temp);

            for (int i = 0; i < N; ++i) {
                current.position[i] += (dt / 6.0) * (k1.position[i] + 2*k2.position[i]
                                                    + 2*k3.position[i] + k4.position[i]);
                current.velocity[i] += (dt / 6.0) * (k1.velocity[i] + 2*k2.velocity[i]
                                                     + 2*k3.velocity[i] + k4.velocity[i]);
            }

            trajectory.push_back(current);
        }
        return trajectory;
    }

    double arclength(const std::vector<Point<N>>& curve,
                    const RiemannianMetric<N>& metric) const {
        double length = 0.0;
        for (size_t i = 0; i < curve.size() - 1; ++i) {
            Vector<N> tangent;
            for (int k = 0; k < N; ++k) {
                tangent[k] = curve[i+1][k] - curve[i][k];
            }
            length += metric.norm(curve[i], tangent);
        }
        return length;
    }
};

// ============================================================================
// EXPONENTIAL MAP
// ============================================================================

template<int N>
class ExponentialMap {
private:
    GeodesicEquation<N> geodesic_;
    const RiemannianMetric<N>& metric_;

public:
    ExponentialMap(const GeodesicEquation<N>& geo, const RiemannianMetric<N>& g)
        : geodesic_(geo), metric_(g) {}

    Point<N> exp(const Point<N>& p, const Vector<N>& v, double t = 1.0) const {
        typename GeodesicEquation<N>::State initial;
        initial.position = p;
        initial.velocity = v;

        double dt = t / 100.0;
        auto trajectory = geodesic_.integrate(initial, dt, 100);
        return trajectory.back().position;
    }

    Vector<N> log(const Point<N>& p, const Point<N>& q, int max_iter = 100,
                 double tol = 1e-8) const {
        Vector<N> v;
        for (int k = 0; k < N; ++k) {
            v[k] = q[k] - p[k];
        }

        for (int iter = 0; iter < max_iter; ++iter) {
            auto exp_v = exp(p, v);
            Vector<N> diff;
            double err = 0.0;
            for (int k = 0; k < N; ++k) {
                diff[k] = q[k] - exp_v[k];
                err += diff[k] * diff[k];
            }

            if (std::sqrt(err) < tol) break;

            for (int k = 0; k < N; ++k) {
                v[k] += 0.5 * diff[k];
            }
        }
        return v;
    }

    double distance(const Point<N>& p, const Point<N>& q) const {
        auto v = log(p, q);
        return metric_.norm(p, v);
    }
};

// ============================================================================
// VOLUME FORMS AND INTEGRATION
// ============================================================================

template<int N>
class VolumeForm {
private:
    const RiemannianMetric<N>& metric_;

public:
    VolumeForm(const RiemannianMetric<N>& g) : metric_(g) {}

    double volumeDensity(const Point<N>& p) const {
        auto g = metric_.at(p);
        double det = determinant(g);
        return std::sqrt(std::abs(det));
    }

    double integrate(std::function<double(const Point<N>&)> f,
                    const std::vector<std::pair<double, double>>& bounds,
                    int n_per_dim = 10) const {
        if (bounds.size() != N) throw std::invalid_argument("Bounds dimension mismatch");

        return integrateRecursive(f, bounds, 0, Point<N>(), n_per_dim);
    }

private:
    double determinant(const Tensor2<N>& g) const {
        if constexpr (N == 2) {
            return g[0][0] * g[1][1] - g[0][1] * g[1][0];
        } else if constexpr (N == 3) {
            return g[0][0] * (g[1][1] * g[2][2] - g[1][2] * g[2][1])
                 - g[0][1] * (g[1][0] * g[2][2] - g[1][2] * g[2][0])
                 + g[0][2] * (g[1][0] * g[2][1] - g[1][1] * g[2][0]);
        }
        return 1.0;
    }

    double integrateRecursive(std::function<double(const Point<N>&)> f,
                             const std::vector<std::pair<double, double>>& bounds,
                             int dim, Point<N> p, int n) const {
        if (dim == N) {
            return f(p) * volumeDensity(p);
        }

        double sum = 0.0;
        double dx = (bounds[dim].second - bounds[dim].first) / n;
        for (int i = 0; i < n; ++i) {
            p[dim] = bounds[dim].first + (i + 0.5) * dx;
            sum += integrateRecursive(f, bounds, dim + 1, p, n) * dx;
        }
        return sum;
    }
};

// ============================================================================
// LIE DERIVATIVES
// ============================================================================

template<int N>
class LieDerivative {
public:
    static Vector<N> ofVectorField(
        const Point<N>& p,
        const std::function<Vector<N>(const Point<N>&)>& X,
        const std::function<Vector<N>(const Point<N>&)>& Y,
        double epsilon = 1e-6) {

        Vector<N> L_X_Y;
        for (int k = 0; k < N; ++k) {
            L_X_Y[k] = 0.0;
            for (int i = 0; i < N; ++i) {
                Point<N> p_shift = p;
                p_shift[i] += epsilon;
                L_X_Y[k] += X(p)[i] * (Y(p_shift)[k] - Y(p)[k]) / epsilon;
                L_X_Y[k] -= Y(p)[i] * (X(p_shift)[k] - X(p)[k]) / epsilon;
            }
        }
        return L_X_Y;
    }

    static Tensor2<N> ofMetric(
        const Point<N>& p,
        const std::function<Vector<N>(const Point<N>&)>& X,
        const RiemannianMetric<N>& g,
        double epsilon = 1e-6) {

        Tensor2<N> L_X_g;
        auto g_p = g.at(p);

        for (int i = 0; i < N; ++i) {
            for (int j = 0; j < N; ++j) {
                L_X_g[i][j] = 0.0;
                for (int k = 0; k < N; ++k) {
                    Point<N> p_shift = p;
                    p_shift[k] += epsilon;
                    L_X_g[i][j] += X(p)[k] * (g.at(p_shift)[i][j] - g_p[i][j]) / epsilon;
                }

                Point<N> p_i = p;
                p_i[i] += epsilon;
                L_X_g[i][j] += g_p[i][j] * (X(p_i)[j] - X(p)[j]) / epsilon;

                Point<N> p_j = p;
                p_j[j] += epsilon;
                L_X_g[i][j] += g_p[i][j] * (X(p_j)[i] - X(p)[i]) / epsilon;
            }
        }
        return L_X_g;
    }

    static bool isKillingField(
        const Point<N>& p,
        const std::function<Vector<N>(const Point<N>&)>& X,
        const RiemannianMetric<N>& g,
        double tol = 1e-6) {

        auto L_X_g = ofMetric(p, X, g);
        for (int i = 0; i < N; ++i) {
            for (int j = 0; j < N; ++j) {
                if (std::abs(L_X_g[i][j]) > tol) {
                    return false;
                }
            }
        }
        return true;
    }
};

} // namespace differential_geometry
} // namespace maths

#endif // MATHS_DIFFERENTIAL_GEOMETRY_HPP
