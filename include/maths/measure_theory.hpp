#ifndef MATHS_MEASURE_THEORY_HPP
#define MATHS_MEASURE_THEORY_HPP

#include <vector>
#include <set>
#include <map>
#include <functional>
#include <cmath>
#include <algorithm>
#include <numeric>
#include <limits>
#include <stdexcept>

namespace maths {
namespace measure_theory {

template<typename T>
using Set = std::set<T>;

template<typename T>
using SetFamily = std::set<Set<T>>;

// ============================================================================
// SIGMA ALGEBRAS
// ============================================================================

template<typename T>
class SigmaAlgebra {
private:
    Set<T> universe_;
    SetFamily<T> algebra_;

public:
    SigmaAlgebra(const Set<T>& universe) : universe_(universe) {
        algebra_.insert(Set<T>{});
        algebra_.insert(universe_);
    }

    void addSet(const Set<T>& s) {
        algebra_.insert(s);
        Set<T> complement;
        std::set_difference(universe_.begin(), universe_.end(),
                          s.begin(), s.end(),
                          std::inserter(complement, complement.begin()));
        algebra_.insert(complement);
    }

    bool contains(const Set<T>& s) const {
        return algebra_.find(s) != algebra_.end();
    }

    Set<T> complement(const Set<T>& s) const {
        Set<T> comp;
        std::set_difference(universe_.begin(), universe_.end(),
                          s.begin(), s.end(),
                          std::inserter(comp, comp.begin()));
        return comp;
    }

    Set<T> setUnion(const Set<T>& a, const Set<T>& b) const {
        Set<T> result;
        std::set_union(a.begin(), a.end(), b.begin(), b.end(),
                      std::inserter(result, result.begin()));
        return result;
    }

    Set<T> intersection(const Set<T>& a, const Set<T>& b) const {
        Set<T> result;
        std::set_intersection(a.begin(), a.end(), b.begin(), b.end(),
                            std::inserter(result, result.begin()));
        return result;
    }

    const SetFamily<T>& getSets() const { return algebra_; }
};

// ============================================================================
// INTERVAL REPRESENTATION (for R)
// ============================================================================

struct Interval {
    double a, b;
    bool left_closed, right_closed;

    Interval(double a_, double b_, bool lc = true, bool rc = false)
        : a(a_), b(b_), left_closed(lc), right_closed(rc) {}

    double length() const { return b - a; }

    bool contains(double x) const {
        bool left_ok = left_closed ? (x >= a) : (x > a);
        bool right_ok = right_closed ? (x <= b) : (x < b);
        return left_ok && right_ok;
    }

    Interval intersect(const Interval& other) const {
        double new_a = std::max(a, other.a);
        double new_b = std::min(b, other.b);
        if (new_a >= new_b) {
            return Interval(0, 0, true, false);
        }
        bool lc = (new_a == a) ? left_closed : other.left_closed;
        bool rc = (new_b == b) ? right_closed : other.right_closed;
        return Interval(new_a, new_b, lc, rc);
    }
};

// ============================================================================
// MEASURE SPACES
// ============================================================================

class Measure {
public:
    virtual double operator()(const Interval& interval) const = 0;
    virtual ~Measure() = default;

    double operator()(const std::vector<Interval>& intervals) const {
        double sum = 0.0;
        for (const auto& I : intervals) {
            sum += (*this)(I);
        }
        return sum;
    }
};

class LebesgueMeasure : public Measure {
public:
    double operator()(const Interval& interval) const override {
        return interval.length();
    }
};

class CountingMeasure : public Measure {
public:
    double operator()(const Interval& interval) const override {
        if (interval.length() == 0.0) return 0.0;
        return std::numeric_limits<double>::infinity();
    }
};

class DiracMeasure : public Measure {
private:
    double point_;

public:
    DiracMeasure(double x0) : point_(x0) {}

    double operator()(const Interval& interval) const override {
        return interval.contains(point_) ? 1.0 : 0.0;
    }
};

// ============================================================================
// MEASURABLE FUNCTIONS
// ============================================================================

class MeasurableFunction {
protected:
    std::function<double(double)> f_;

public:
    MeasurableFunction(std::function<double(double)> f) : f_(f) {}

    double operator()(double x) const { return f_(x); }

    double supremum(double a, double b, int n = 1000) const {
        double sup = -std::numeric_limits<double>::infinity();
        double dx = (b - a) / n;
        for (int i = 0; i <= n; ++i) {
            double x = a + i * dx;
            sup = std::max(sup, f_(x));
        }
        return sup;
    }

    double infimum(double a, double b, int n = 1000) const {
        double inf = std::numeric_limits<double>::infinity();
        double dx = (b - a) / n;
        for (int i = 0; i <= n; ++i) {
            double x = a + i * dx;
            inf = std::min(inf, f_(x));
        }
        return inf;
    }

    MeasurableFunction compose(const MeasurableFunction& g) const {
        return MeasurableFunction([this, &g](double x) {
            return f_(g(x));
        });
    }
};

// ============================================================================
// SIMPLE FUNCTIONS
// ============================================================================

class SimpleFunction : public MeasurableFunction {
private:
    std::vector<double> values_;
    std::vector<Interval> sets_;

public:
    SimpleFunction() : MeasurableFunction([](double) { return 0.0; }) {}

    void addLevel(double value, const Interval& set) {
        values_.push_back(value);
        sets_.push_back(set);

        f_ = [this](double x) {
            for (size_t i = 0; i < sets_.size(); ++i) {
                if (sets_[i].contains(x)) {
                    return values_[i];
                }
            }
            return 0.0;
        };
    }

    const std::vector<double>& getValues() const { return values_; }
    const std::vector<Interval>& getSets() const { return sets_; }
};

// ============================================================================
// LEBESGUE INTEGRATION
// ============================================================================

class LebesgueIntegral {
private:
    const Measure& measure_;

public:
    LebesgueIntegral(const Measure& mu) : measure_(mu) {}

    double integrate(const SimpleFunction& f, double a, double b) const {
        double sum = 0.0;
        const auto& values = f.getValues();
        const auto& sets = f.getSets();

        for (size_t i = 0; i < values.size(); ++i) {
            Interval domain(a, b);
            Interval E_i = sets[i].intersect(domain);
            sum += values[i] * measure_(E_i);
        }
        return sum;
    }

    double integrate(const MeasurableFunction& f, double a, double b,
                    int n = 10000) const {
        double sum = 0.0;
        double dx = (b - a) / n;

        for (int i = 0; i < n; ++i) {
            double x = a + i * dx;
            double y = f(x);
            if (std::isfinite(y)) {
                sum += y * dx;
            }
        }
        return sum;
    }

    double integrateLowerBound(const MeasurableFunction& f, double a, double b,
                               int n = 1000) const {
        std::vector<SimpleFunction> approximations;
        double max_val = 0.0;

        for (int k = 1; k <= 10; ++k) {
            SimpleFunction s;
            double dx = (b - a) / (n * k);
            for (int i = 0; i < n * k; ++i) {
                double x_i = a + i * dx;
                double x_next = a + (i + 1) * dx;
                double val = f.infimum(x_i, x_next, 100);
                if (std::isfinite(val)) {
                    s.addLevel(val, Interval(x_i, x_next));
                    max_val = std::max(max_val, std::abs(val));
                }
            }

            double integral = integrate(s, a, b);
            if (k > 1 && std::abs(integral - max_val) < 1e-6) {
                return integral;
            }
            max_val = integral;
        }
        return max_val;
    }
};

// ============================================================================
// L^p SPACES
// ============================================================================

class LpSpace {
private:
    double p_;
    const Measure& measure_;

public:
    LpSpace(double p, const Measure& mu) : p_(p), measure_(mu) {
        if (p < 1.0) {
            throw std::invalid_argument("p must be >= 1");
        }
    }

    double norm(const MeasurableFunction& f, double a, double b,
                int n = 10000) const {
        if (std::isinf(p_)) {
            return f.supremum(a, b, n);
        }

        LebesgueIntegral integrator(measure_);
        auto f_p = MeasurableFunction([&f, this](double x) {
            return std::pow(std::abs(f(x)), p_);
        });

        double integral = integrator.integrate(f_p, a, b, n);
        return std::pow(integral, 1.0 / p_);
    }

    double distance(const MeasurableFunction& f, const MeasurableFunction& g,
                   double a, double b, int n = 10000) const {
        auto diff = MeasurableFunction([&f, &g](double x) {
            return f(x) - g(x);
        });
        return norm(diff, a, b, n);
    }

    double innerProduct(const MeasurableFunction& f,
                       const MeasurableFunction& g,
                       double a, double b, int n = 10000) const {
        if (p_ != 2.0) {
            throw std::logic_error("Inner product only defined for L^2");
        }

        LebesgueIntegral integrator(measure_);
        auto product = MeasurableFunction([&f, &g](double x) {
            return f(x) * g(x);
        });

        return integrator.integrate(product, a, b, n);
    }
};

// ============================================================================
// CONVERGENCE THEOREMS (Numerical Verification)
// ============================================================================

class ConvergenceTheorems {
public:
    static bool verifyMonotoneConvergence(
        const std::vector<MeasurableFunction>& sequence,
        const MeasurableFunction& limit,
        double a, double b,
        const Measure& mu,
        double tol = 1e-6) {

        LebesgueIntegral integrator(mu);
        int n = sequence.size();

        std::vector<double> integrals(n);
        for (int i = 0; i < n; ++i) {
            integrals[i] = integrator.integrate(sequence[i], a, b);
        }

        for (int i = 1; i < n; ++i) {
            if (integrals[i] < integrals[i-1] - tol) {
                return false;
            }
        }

        double limit_integral = integrator.integrate(limit, a, b);
        return std::abs(integrals.back() - limit_integral) < tol;
    }

    static bool verifyDominatedConvergence(
        const std::vector<MeasurableFunction>& sequence,
        const MeasurableFunction& limit,
        const MeasurableFunction& dominating,
        double a, double b,
        const Measure& mu,
        double tol = 1e-6) {

        LebesgueIntegral integrator(mu);

        double dx = (b - a) / 1000;
        for (double x = a; x < b; x += dx) {
            for (const auto& f_n : sequence) {
                if (std::abs(f_n(x)) > std::abs(dominating(x)) + tol) {
                    return false;
                }
            }
        }

        std::vector<double> integrals(sequence.size());
        for (size_t i = 0; i < sequence.size(); ++i) {
            integrals[i] = integrator.integrate(sequence[i], a, b);
        }

        double limit_integral = integrator.integrate(limit, a, b);
        return std::abs(integrals.back() - limit_integral) < tol;
    }

    static double computeLimitIntegral(
        std::function<MeasurableFunction(int)> sequence_generator,
        double a, double b,
        const Measure& mu,
        int max_n = 100,
        double tol = 1e-8) {

        LebesgueIntegral integrator(mu);

        double prev_integral = integrator.integrate(sequence_generator(1), a, b);
        for (int n = 2; n <= max_n; ++n) {
            double curr_integral = integrator.integrate(sequence_generator(n), a, b);
            if (std::abs(curr_integral - prev_integral) < tol) {
                return curr_integral;
            }
            prev_integral = curr_integral;
        }
        return prev_integral;
    }
};

// ============================================================================
// PRODUCT MEASURES
// ============================================================================

class ProductMeasure {
private:
    const Measure& mu1_;
    const Measure& mu2_;

public:
    ProductMeasure(const Measure& mu1, const Measure& mu2)
        : mu1_(mu1), mu2_(mu2) {}

    double operator()(const Interval& I1, const Interval& I2) const {
        return mu1_(I1) * mu2_(I2);
    }

    double fubiniIntegrate(
        std::function<double(double, double)> f,
        double a1, double b1,
        double a2, double b2,
        int n = 1000) const {

        double sum = 0.0;
        double dx1 = (b1 - a1) / n;
        double dx2 = (b2 - a2) / n;

        for (int i = 0; i < n; ++i) {
            for (int j = 0; j < n; ++j) {
                double x1 = a1 + i * dx1;
                double x2 = a2 + j * dx2;
                double val = f(x1, x2);
                if (std::isfinite(val)) {
                    sum += val * dx1 * dx2;
                }
            }
        }
        return sum;
    }

    bool verifyFubini(
        std::function<double(double, double)> f,
        double a1, double b1,
        double a2, double b2,
        double tol = 1e-6) const {

        auto integrate_x_first = [&]() {
            double sum = 0.0;
            int n = 500;
            double dx2 = (b2 - a2) / n;
            for (int j = 0; j < n; ++j) {
                double x2 = a2 + j * dx2;
                double inner_sum = 0.0;
                double dx1 = (b1 - a1) / n;
                for (int i = 0; i < n; ++i) {
                    double x1 = a1 + i * dx1;
                    inner_sum += f(x1, x2) * dx1;
                }
                sum += inner_sum * dx2;
            }
            return sum;
        };

        auto integrate_y_first = [&]() {
            double sum = 0.0;
            int n = 500;
            double dx1 = (b1 - a1) / n;
            for (int i = 0; i < n; ++i) {
                double x1 = a1 + i * dx1;
                double inner_sum = 0.0;
                double dx2 = (b2 - a2) / n;
                for (int j = 0; j < n; ++j) {
                    double x2 = a2 + j * dx2;
                    inner_sum += f(x1, x2) * dx2;
                }
                sum += inner_sum * dx1;
            }
            return sum;
        };

        double I1 = integrate_x_first();
        double I2 = integrate_y_first();
        return std::abs(I1 - I2) < tol;
    }
};

// ============================================================================
// SIGNED MEASURES AND HAHN DECOMPOSITION
// ============================================================================

class SignedMeasure {
private:
    std::function<double(const Interval&)> measure_;

public:
    SignedMeasure(std::function<double(const Interval&)> mu) : measure_(mu) {}

    double operator()(const Interval& I) const {
        return measure_(I);
    }

    std::pair<Interval, Interval> hahnDecomposition(double a, double b,
                                                     int n = 1000) const {
        std::vector<Interval> positive_sets;
        std::vector<Interval> negative_sets;

        double dx = (b - a) / n;
        for (int i = 0; i < n; ++i) {
            Interval I(a + i * dx, a + (i + 1) * dx);
            double mu_I = measure_(I);
            if (mu_I >= 0) {
                positive_sets.push_back(I);
            } else {
                negative_sets.push_back(I);
            }
        }

        if (positive_sets.empty()) {
            return {Interval(a, a), Interval(a, b)};
        }
        if (negative_sets.empty()) {
            return {Interval(a, b), Interval(a, a)};
        }

        Interval P(positive_sets.front().a, positive_sets.back().b);
        Interval N(negative_sets.front().a, negative_sets.back().b);
        return {P, N};
    }

    double totalVariation(double a, double b, int n = 1000) const {
        double pos = 0.0, neg = 0.0;
        double dx = (b - a) / n;

        for (int i = 0; i < n; ++i) {
            Interval I(a + i * dx, a + (i + 1) * dx);
            double mu_I = measure_(I);
            if (mu_I >= 0) {
                pos += mu_I;
            } else {
                neg += -mu_I;
            }
        }
        return pos + neg;
    }
};

// ============================================================================
// RADON-NIKODYM DERIVATIVE (Numerical)
// ============================================================================

class RadonNikodym {
public:
    static MeasurableFunction computeDerivative(
        const Measure& nu,
        const Measure& mu,
        double a, double b,
        int n = 1000) {

        std::vector<double> ratios;
        std::vector<double> points;

        double dx = (b - a) / n;
        for (int i = 0; i < n; ++i) {
            Interval I(a + i * dx, a + (i + 1) * dx);
            double mu_I = mu(I);
            double nu_I = nu(I);

            if (mu_I > 1e-10) {
                ratios.push_back(nu_I / mu_I);
                points.push_back(a + (i + 0.5) * dx);
            }
        }

        return MeasurableFunction([ratios, points, dx](double x) {
            for (size_t i = 0; i < points.size(); ++i) {
                if (std::abs(x - points[i]) < dx) {
                    return ratios[i];
                }
            }
            return 0.0;
        });
    }

    static bool verifyAbsoluteContinuity(
        const Measure& nu,
        const Measure& mu,
        double a, double b,
        int n = 1000,
        double epsilon = 1e-6) {

        double dx = (b - a) / n;
        for (int i = 0; i < n; ++i) {
            Interval I(a + i * dx, a + (i + 1) * dx);
            double mu_I = mu(I);
            double nu_I = nu(I);

            if (mu_I < epsilon && nu_I > epsilon) {
                return false;
            }
        }
        return true;
    }
};

} // namespace measure_theory
} // namespace maths

#endif // MATHS_MEASURE_THEORY_HPP
