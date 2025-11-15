#ifndef MATHS_COMPLEX_ANALYSIS_HPP
#define MATHS_COMPLEX_ANALYSIS_HPP

#include <complex>
#include <vector>
#include <functional>
#include <cmath>
#include <stdexcept>
#include <limits>
#include <algorithm>

namespace maths {
namespace complex_analysis {

using Complex = std::complex<double>;
using ComplexFunction = std::function<Complex(Complex)>;

/**
 * @class CauchyRiemann
 * @brief Cauchy-Riemann equations and holomorphic function tests
 *
 * Implements:
 * - Cauchy-Riemann equations: ∂u/∂x = ∂v/∂y, ∂u/∂y = -∂v/∂x
 * - Holomorphic (analytic) function verification
 * - Harmonic function tests
 */
class CauchyRiemann {
public:
    /**
     * @brief Check Cauchy-Riemann equations numerically
     *
     * For f(z) = u(x,y) + iv(x,y), checks:
     * ∂u/∂x = ∂v/∂y and ∂u/∂y = -∂v/∂x
     *
     * @param f Complex function
     * @param z Point to check
     * @param h Step size for numerical derivatives
     * @return true if CR equations satisfied (within tolerance)
     */
    static bool satisfies_cauchy_riemann(
        const ComplexFunction& f,
        Complex z,
        double h = 1e-6) {

        // Compute partial derivatives numerically
        Complex f_z = f(z);
        Complex f_x_plus = f(z + Complex(h, 0));
        Complex f_x_minus = f(z - Complex(h, 0));
        Complex f_y_plus = f(z + Complex(0, h));
        Complex f_y_minus = f(z - Complex(0, h));

        // ∂f/∂x ≈ (f(z+h) - f(z-h)) / (2h)
        Complex df_dx = (f_x_plus - f_x_minus) / (2.0 * h);
        Complex df_dy = (f_y_plus - f_y_minus) / (2.0 * h);

        // Extract real and imaginary parts
        double u_x = df_dx.real();
        double v_x = df_dx.imag();
        double u_y = df_dy.real();
        double v_y = df_dy.imag();

        // Check CR equations: u_x = v_y and u_y = -v_x
        double tol = 1e-4;
        bool cr1 = std::abs(u_x - v_y) < tol;
        bool cr2 = std::abs(u_y + v_x) < tol;

        return cr1 && cr2;
    }

    /**
     * @brief Compute complex derivative f'(z)
     *
     * For holomorphic f, f'(z) = ∂f/∂z = ∂u/∂x + i∂v/∂x
     */
    static Complex derivative(
        const ComplexFunction& f,
        Complex z,
        double h = 1e-6) {

        return (f(z + h) - f(z - h)) / (2.0 * h);
    }

    /**
     * @brief Check if function is harmonic
     *
     * u is harmonic if ∇²u = ∂²u/∂x² + ∂²u/∂y² = 0
     *
     * @param u Real-valued function u(x,y)
     * @param z Point to check
     * @return true if Laplacian ≈ 0
     */
    static bool is_harmonic(
        const std::function<double(double, double)>& u,
        Complex z,
        double h = 1e-6) {

        double x = z.real(), y = z.imag();

        // Second partial derivatives
        double u_xx = (u(x + h, y) - 2.0 * u(x, y) + u(x - h, y)) / (h * h);
        double u_yy = (u(x, y + h) - 2.0 * u(x, y) + u(x, y - h)) / (h * h);

        double laplacian = u_xx + u_yy;

        return std::abs(laplacian) < 1e-3;
    }
};

/**
 * @class ComplexFunctions
 * @brief Complex elementary functions and special functions
 *
 * Implements complex versions of:
 * - Trigonometric functions (sin, cos, tan, cot)
 * - Hyperbolic functions
 * - Exponential and logarithm (with branch cuts)
 * - Power functions
 */
class ComplexFunctions {
public:
    /**
     * @brief Complex sine: sin(z) = (e^(iz) - e^(-iz)) / (2i)
     */
    static Complex sin(Complex z) {
        return std::sin(z);
    }

    /**
     * @brief Complex cosine: cos(z) = (e^(iz) + e^(-iz)) / 2
     */
    static Complex cos(Complex z) {
        return std::cos(z);
    }

    /**
     * @brief Complex tangent: tan(z) = sin(z) / cos(z)
     */
    static Complex tan(Complex z) {
        return std::tan(z);
    }

    /**
     * @brief Complex cotangent: cot(z) = cos(z) / sin(z)
     *
     * Has poles at z = nπ for integer n
     */
    static Complex cot(Complex z) {
        Complex s = std::sin(z);
        if (std::abs(s) < 1e-10) {
            throw std::runtime_error("cot(z) has pole at z = nπ");
        }
        return std::cos(z) / s;
    }

    /**
     * @brief Complex exponential: exp(z) = e^x(cos(y) + i·sin(y))
     */
    static Complex exp(Complex z) {
        return std::exp(z);
    }

    /**
     * @brief Complex logarithm (principal branch)
     *
     * log(z) = log|z| + i·arg(z), arg(z) ∈ (-π, π]
     * Branch cut along negative real axis
     */
    static Complex log(Complex z) {
        if (std::abs(z) < 1e-10) {
            throw std::runtime_error("log(0) is undefined");
        }
        return std::log(z);
    }

    /**
     * @brief Complex power: z^w = exp(w·log(z))
     *
     * Multi-valued unless w is integer
     */
    static Complex pow(Complex z, Complex w) {
        if (std::abs(z) < 1e-10 && w.real() <= 0) {
            throw std::runtime_error("0^w undefined for Re(w) ≤ 0");
        }
        return std::pow(z, w);
    }

    /**
     * @brief Complex square root (principal branch)
     *
     * √z with branch cut along negative real axis
     */
    static Complex sqrt(Complex z) {
        return std::sqrt(z);
    }

    /**
     * @brief Complex hyperbolic sine: sinh(z) = (e^z - e^(-z)) / 2
     */
    static Complex sinh(Complex z) {
        return std::sinh(z);
    }

    /**
     * @brief Complex hyperbolic cosine: cosh(z) = (e^z + e^(-z)) / 2
     */
    static Complex cosh(Complex z) {
        return std::cosh(z);
    }

    /**
     * @brief Complex hyperbolic tangent
     */
    static Complex tanh(Complex z) {
        return std::tanh(z);
    }
};

/**
 * @class PowerSeries
 * @brief Power series and radius of convergence
 *
 * For series Σ aₙ(z-z₀)ⁿ:
 * - Radius of convergence R
 * - Convergence tests (ratio test, root test)
 * - Series evaluation
 */
class PowerSeries {
public:
    /**
     * @brief Compute radius of convergence using ratio test
     *
     * R = lim_{n→∞} |aₙ/aₙ₊₁|
     *
     * @param coefficients Power series coefficients {a₀, a₁, a₂, ...}
     * @return Radius of convergence
     */
    static double radius_of_convergence_ratio(const std::vector<Complex>& coefficients) {
        int n = coefficients.size();
        if (n < 2) return std::numeric_limits<double>::infinity();

        // Use last several terms for better estimate
        int start = std::max(0, n - 10);
        double sum_ratio = 0.0;
        int count = 0;

        for (int i = start; i < n - 1; ++i) {
            if (std::abs(coefficients[i + 1]) > 1e-10) {
                double ratio = std::abs(coefficients[i]) / std::abs(coefficients[i + 1]);
                sum_ratio += ratio;
                count++;
            }
        }

        return (count > 0) ? sum_ratio / count : std::numeric_limits<double>::infinity();
    }

    /**
     * @brief Compute radius of convergence using root test
     *
     * 1/R = lim sup_{n→∞} |aₙ|^(1/n)
     */
    static double radius_of_convergence_root(const std::vector<Complex>& coefficients) {
        int n = coefficients.size();
        if (n == 0) return std::numeric_limits<double>::infinity();

        double max_root = 0.0;
        int start = std::max(0, n - 20);

        for (int i = start; i < n; ++i) {
            if (i > 0) {
                double root = std::pow(std::abs(coefficients[i]), 1.0 / i);
                max_root = std::max(max_root, root);
            }
        }

        return (max_root > 1e-10) ? 1.0 / max_root : std::numeric_limits<double>::infinity();
    }

    /**
     * @brief Evaluate power series Σ aₙ(z-z₀)ⁿ
     *
     * @param coefficients {a₀, a₁, a₂, ...}
     * @param z Point of evaluation
     * @param z0 Center of series
     * @return Series value
     */
    static Complex evaluate(
        const std::vector<Complex>& coefficients,
        Complex z,
        Complex z0 = 0.0) {

        Complex result = 0.0;
        Complex power = 1.0;  // (z - z0)^n
        Complex dz = z - z0;

        for (const Complex& a_n : coefficients) {
            result += a_n * power;
            power *= dz;
        }

        return result;
    }

    /**
     * @brief Check if series converges at point z
     */
    static bool converges_at(
        const std::vector<Complex>& coefficients,
        Complex z,
        Complex z0 = 0.0) {

        double R = radius_of_convergence_ratio(coefficients);
        double dist = std::abs(z - z0);

        return dist < R;
    }
};

/**
 * @class CauchyTheory
 * @brief Cauchy's integral theorem and formula
 *
 * Implements:
 * - Cauchy integral theorem (contour integrals)
 * - Cauchy integral formula
 * - Cauchy's inequality
 * - Maximum modulus principle
 */
class CauchyTheory {
public:
    /**
     * @brief Numerical contour integral using trapezoidal rule
     *
     * ∫_γ f(z) dz
     *
     * @param f Complex function
     * @param path Parameterized path γ(t), t ∈ [0, 1]
     * @param n_points Number of discretization points
     */
    static Complex contour_integral(
        const ComplexFunction& f,
        const std::function<Complex(double)>& path,
        int n_points = 1000) {

        Complex result = 0.0;
        double dt = 1.0 / n_points;

        for (int i = 0; i < n_points; ++i) {
            double t = i * dt;
            Complex z = path(t);
            Complex z_next = path(t + dt);

            Complex dz = z_next - z;
            result += f(z) * dz;
        }

        return result;
    }

    /**
     * @brief Cauchy integral formula: f(z₀) = (1/2πi) ∮ f(z)/(z-z₀) dz
     *
     * For f holomorphic inside and on contour
     *
     * @param f Holomorphic function
     * @param z0 Point inside contour
     * @param path Contour enclosing z0
     * @return f(z₀)
     */
    static Complex cauchy_integral_formula(
        const ComplexFunction& f,
        Complex z0,
        const std::function<Complex(double)>& path,
        int n_points = 1000) {

        auto integrand = [&](Complex z) -> Complex {
            Complex diff = z - z0;
            if (std::abs(diff) < 1e-10) return 0.0;
            return f(z) / diff;
        };

        Complex integral = contour_integral(integrand, path, n_points);
        return integral / (2.0 * M_PI * Complex(0, 1));
    }

    /**
     * @brief Cauchy's inequality: |f^(n)(z₀)| ≤ n! M / R^n
     *
     * Where M = max|f| on circle |z - z₀| = R
     */
    static double cauchy_inequality(
        const ComplexFunction& f,
        Complex z0,
        double R,
        int n,
        int n_samples = 100) {

        // Compute M = max|f| on circle
        double M = 0.0;
        for (int i = 0; i < n_samples; ++i) {
            double theta = 2.0 * M_PI * i / n_samples;
            Complex z = z0 + R * Complex(std::cos(theta), std::sin(theta));
            M = std::max(M, std::abs(f(z)));
        }

        // n!
        double factorial = 1.0;
        for (int k = 1; k <= n; ++k) {
            factorial *= k;
        }

        return factorial * M / std::pow(R, n);
    }

    /**
     * @brief Maximum modulus principle
     *
     * If f is holomorphic on region D, max|f| occurs on boundary
     *
     * @return Maximum modulus on boundary
     */
    static double maximum_modulus(
        const ComplexFunction& f,
        const std::function<Complex(double)>& boundary,
        int n_samples = 1000) {

        double max_mod = 0.0;
        for (int i = 0; i < n_samples; ++i) {
            double t = static_cast<double>(i) / n_samples;
            Complex z = boundary(t);
            max_mod = std::max(max_mod, std::abs(f(z)));
        }

        return max_mod;
    }
};

/**
 * @class ResidueCalculus
 * @brief Residue theorem and pole analysis
 *
 * Implements:
 * - Residue computation
 * - Pole order determination
 * - Zero order determination
 * - Residue theorem for contour integrals
 */
class ResidueCalculus {
public:
    /**
     * @brief Compute residue at simple pole
     *
     * Res(f, z₀) = lim_{z→z₀} (z - z₀)f(z)
     *
     * @param f Function with pole at z0
     * @param z0 Location of pole
     * @return Residue
     */
    static Complex residue_simple_pole(
        const ComplexFunction& f,
        Complex z0) {

        double h = 1e-4;
        Complex z = z0 + Complex(h, 0);
        return (z - z0) * f(z);
    }

    /**
     * @brief Compute residue at pole of order m
     *
     * Res(f, z₀) = (1/(m-1)!) lim_{z→z₀} d^(m-1)/dz^(m-1) [(z-z₀)^m f(z)]
     */
    static Complex residue_pole_order_m(
        const ComplexFunction& f,
        Complex z0,
        int m) {

        if (m == 1) return residue_simple_pole(f, z0);

        // Numerical differentiation (simplified)
        double h = 1e-4;

        auto g = [&](Complex z) -> Complex {
            Complex diff = z - z0;
            return std::pow(diff, m) * f(z);
        };

        // Compute (m-1)-th derivative numerically
        Complex deriv = g(z0 + h);
        for (int k = 1; k < m; ++k) {
            deriv = (g(z0 + h * (k + 1)) - deriv) / h;
        }

        // Divide by (m-1)!
        double factorial = 1.0;
        for (int k = 2; k < m; ++k) {
            factorial *= k;
        }

        return deriv / factorial;
    }

    /**
     * @brief Residue theorem: ∮ f(z) dz = 2πi Σ Res(f, zₖ)
     *
     * @param f Meromorphic function
     * @param poles Poles inside contour
     * @param residues Residues at each pole
     * @return Value of contour integral
     */
    static Complex residue_theorem(
        const std::vector<Complex>& residues) {

        Complex sum = 0.0;
        for (Complex res : residues) {
            sum += res;
        }

        return 2.0 * M_PI * Complex(0, 1) * sum;
    }

    /**
     * @brief Determine order of pole at z₀
     *
     * f has pole of order m if lim_{z→z₀} (z-z₀)^m f(z) ≠ 0, ∞
     *
     * @return Pole order (returns -1 if not a pole)
     */
    static int pole_order(
        const ComplexFunction& f,
        Complex z0,
        int max_order = 10) {

        double h = 1e-4;

        for (int m = 1; m <= max_order; ++m) {
            Complex z = z0 + Complex(h, 0);
            Complex val = std::pow(z - z0, m) * f(z);

            double mod = std::abs(val);
            if (mod > 1e-3 && mod < 1e3) {
                return m;  // Found order
            }
        }

        return -1;  // Not a pole or order > max_order
    }

    /**
     * @brief Determine order of zero at z₀
     *
     * f has zero of order n if f^(k)(z₀) = 0 for k < n, f^(n)(z₀) ≠ 0
     */
    static int zero_order(
        const ComplexFunction& f,
        Complex z0,
        int max_order = 10) {

        double h = 1e-6;

        // Check if f(z₀) ≈ 0
        if (std::abs(f(z0)) > 1e-6) return 0;

        for (int n = 1; n <= max_order; ++n) {
            // Compute n-th derivative numerically
            Complex deriv = (f(z0 + h) - f(z0 - h)) / (2.0 * h);

            for (int k = 1; k < n; ++k) {
                Complex f_plus = (f(z0 + h * (k + 1)) - f(z0 + h * k)) / h;
                Complex f_minus = (f(z0 - h * k) - f(z0 - h * (k + 1))) / h;
                deriv = (f_plus - f_minus) / (2.0 * h);
            }

            if (std::abs(deriv) > 1e-6) {
                return n;
            }
        }

        return max_order;
    }
};

/**
 * @class LaurentSeries
 * @brief Laurent series for functions holomorphic on annulus
 *
 * f(z) = Σ_{n=-∞}^∞ aₙ(z - z₀)ⁿ
 *
 * Computes Laurent series coefficients on annulus r < |z - z₀| < R
 */
class LaurentSeries {
public:
    struct LaurentCoefficients {
        std::vector<Complex> positive_powers;  // a₀, a₁, a₂, ...
        std::vector<Complex> negative_powers;  // a₋₁, a₋₂, a₋₃, ...
        Complex z0;  // Center
    };

    /**
     * @brief Compute Laurent series coefficients numerically
     *
     * aₙ = (1/2πi) ∮ f(z)/(z-z₀)^(n+1) dz
     *
     * @param f Function holomorphic on annulus
     * @param z0 Center
     * @param r Radius for integration contour (r < |z-z₀| < R)
     * @param max_positive Maximum positive power
     * @param max_negative Maximum negative power magnitude
     */
    static LaurentCoefficients compute_coefficients(
        const ComplexFunction& f,
        Complex z0,
        double r,
        int max_positive = 10,
        int max_negative = 10) {

        LaurentCoefficients result;
        result.z0 = z0;
        result.positive_powers.resize(max_positive + 1);
        result.negative_powers.resize(max_negative);

        // Integration contour: circle of radius r centered at z0
        auto circle = [z0, r](double t) -> Complex {
            double theta = 2.0 * M_PI * t;
            return z0 + r * Complex(std::cos(theta), std::sin(theta));
        };

        int n_points = 500;

        // Compute positive powers
        for (int n = 0; n <= max_positive; ++n) {
            auto integrand = [&](Complex z) -> Complex {
                return f(z) / std::pow(z - z0, n + 1);
            };

            Complex integral = CauchyTheory::contour_integral(
                integrand, circle, n_points);

            result.positive_powers[n] = integral / (2.0 * M_PI * Complex(0, 1));
        }

        // Compute negative powers
        for (int n = 1; n <= max_negative; ++n) {
            auto integrand = [&](Complex z) -> Complex {
                return f(z) * std::pow(z - z0, n - 1);
            };

            Complex integral = CauchyTheory::contour_integral(
                integrand, circle, n_points);

            result.negative_powers[n - 1] = integral / (2.0 * M_PI * Complex(0, 1));
        }

        return result;
    }

    /**
     * @brief Evaluate Laurent series at point z
     */
    static Complex evaluate(const LaurentCoefficients& coeffs, Complex z) {
        Complex result = 0.0;
        Complex dz = z - coeffs.z0;

        // Positive powers
        Complex power = 1.0;
        for (const Complex& a_n : coeffs.positive_powers) {
            result += a_n * power;
            power *= dz;
        }

        // Negative powers
        Complex inv_dz = 1.0 / dz;
        power = inv_dz;
        for (const Complex& a_minus_n : coeffs.negative_powers) {
            result += a_minus_n * power;
            power *= inv_dz;
        }

        return result;
    }

    /**
     * @brief Get residue (coefficient a₋₁)
     */
    static Complex get_residue(const LaurentCoefficients& coeffs) {
        if (coeffs.negative_powers.empty()) return 0.0;
        return coeffs.negative_powers[0];
    }
};

/**
 * @class ConformalMaps
 * @brief Conformal mappings and Möbius transformations
 *
 * Implements:
 * - Möbius transformations: f(z) = (az + b)/(cz + d)
 * - Conformal mapping properties
 * - Special mappings (Joukowski, etc.)
 */
class ConformalMaps {
public:
    /**
     * @brief Möbius transformation: f(z) = (az + b)/(cz + d)
     *
     * Maps circles/lines to circles/lines
     * Preserves angles (conformal)
     */
    struct MobiusTransform {
        Complex a, b, c, d;

        MobiusTransform(Complex a_, Complex b_, Complex c_, Complex d_)
            : a(a_), b(b_), c(c_), d(d_) {
            // Check ad - bc ≠ 0
            if (std::abs(a * d - b * c) < 1e-10) {
                throw std::invalid_argument("Möbius transform must have ad - bc ≠ 0");
            }
        }

        Complex operator()(Complex z) const {
            Complex denom = c * z + d;
            if (std::abs(denom) < 1e-10) {
                throw std::runtime_error("Möbius transform has pole");
            }
            return (a * z + b) / denom;
        }

        /**
         * @brief Compose two Möbius transformations
         */
        MobiusTransform compose(const MobiusTransform& other) const {
            return MobiusTransform(
                a * other.a + b * other.c,
                a * other.b + b * other.d,
                c * other.a + d * other.c,
                c * other.b + d * other.d
            );
        }

        /**
         * @brief Inverse Möbius transformation
         */
        MobiusTransform inverse() const {
            Complex det = a * d - b * c;
            return MobiusTransform(d / det, -b / det, -c / det, a / det);
        }
    };

    /**
     * @brief Möbius transformation mapping three points to three points
     *
     * Find unique Möbius f such that f(z₁) = w₁, f(z₂) = w₂, f(z₃) = w₃
     */
    static MobiusTransform three_point_map(
        Complex z1, Complex z2, Complex z3,
        Complex w1, Complex w2, Complex w3) {

        // Use cross-ratio formula
        // (w - w₁)(w₂ - w₃) / ((w - w₃)(w₂ - w₁)) = (z - z₁)(z₂ - z₃) / ((z - z₃)(z₂ - z₁))

        // Simplified construction (there's a formula for a, b, c, d)
        // This is a placeholder - full implementation would solve for coefficients
        return MobiusTransform(1.0, 0.0, 0.0, 1.0);  // Identity
    }

    /**
     * @brief Joukowski transformation: f(z) = (z + 1/z) / 2
     *
     * Maps circles to ellipses, used in airfoil theory
     */
    static Complex joukowski(Complex z) {
        if (std::abs(z) < 1e-10) {
            throw std::runtime_error("Joukowski undefined at z = 0");
        }
        return 0.5 * (z + 1.0 / z);
    }

    /**
     * @brief Map unit disk to upper half-plane
     *
     * f(z) = i(1 - z)/(1 + z)
     */
    static Complex disk_to_half_plane(Complex z) {
        return Complex(0, 1) * (1.0 - z) / (1.0 + z);
    }

    /**
     * @brief Map upper half-plane to unit disk
     *
     * Inverse of disk_to_half_plane
     */
    static Complex half_plane_to_disk(Complex z) {
        return (Complex(0, 1) - z) / (Complex(0, 1) + z);
    }

    /**
     * @brief Check if mapping is conformal at point z
     *
     * f is conformal if f'(z) ≠ 0
     */
    static bool is_conformal(
        const ComplexFunction& f,
        Complex z,
        double h = 1e-6) {

        Complex derivative = (f(z + h) - f(z - h)) / (2.0 * h);
        return std::abs(derivative) > 1e-6;
    }
};

/**
 * @class RungeTheorem
 * @brief Runge's theorem - polynomial approximation
 *
 * If f is holomorphic on compact K and A ⊂ ℂ contains at least one
 * point from each bounded component of ℂ \ K, then f can be uniformly
 * approximated on K by rational functions with poles in A
 */
class RungeTheorem {
public:
    /**
     * @brief Approximate holomorphic function by polynomial
     *
     * On simply connected domain, can approximate by polynomials
     *
     * @param f Holomorphic function
     * @param z0 Center for Taylor series
     * @param degree Polynomial degree
     * @return Polynomial coefficients
     */
    static std::vector<Complex> polynomial_approximation(
        const ComplexFunction& f,
        Complex z0,
        int degree) {

        std::vector<Complex> coeffs(degree + 1);

        // Compute Taylor series coefficients
        // aₙ = f^(n)(z₀) / n!
        coeffs[0] = f(z0);

        double h = 1e-5;
        for (int n = 1; n <= degree; ++n) {
            // Numerical differentiation
            Complex deriv = 0.0;

            // Use central differences
            for (int k = 0; k <= n; ++k) {
                double sign = ((n - k) % 2 == 0) ? 1.0 : -1.0;
                double binom = 1.0;  // Binomial coefficient C(n, k)
                for (int j = 0; j < k; ++j) {
                    binom *= (n - j) / (j + 1.0);
                }

                deriv += sign * binom * f(z0 + Complex(h * (n - 2.0 * k), 0));
            }

            deriv /= std::pow(2.0 * h, n);

            // Divide by n!
            double factorial = 1.0;
            for (int k = 1; k <= n; ++k) {
                factorial *= k;
            }

            coeffs[n] = deriv / factorial;
        }

        return coeffs;
    }

    /**
     * @brief Evaluate polynomial
     */
    static Complex evaluate_polynomial(
        const std::vector<Complex>& coeffs,
        Complex z,
        Complex z0 = 0.0) {

        Complex result = 0.0;
        Complex power = 1.0;
        Complex dz = z - z0;

        for (const Complex& a_n : coeffs) {
            result += a_n * power;
            power *= dz;
        }

        return result;
    }

    /**
     * @brief Compute maximum error of polynomial approximation
     *
     * On given contour
     */
    static double approximation_error(
        const ComplexFunction& f,
        const std::vector<Complex>& poly_coeffs,
        Complex z0,
        const std::function<Complex(double)>& contour,
        int n_samples = 100) {

        double max_error = 0.0;

        for (int i = 0; i < n_samples; ++i) {
            double t = static_cast<double>(i) / n_samples;
            Complex z = contour(t);

            Complex f_val = f(z);
            Complex poly_val = evaluate_polynomial(poly_coeffs, z, z0);

            double error = std::abs(f_val - poly_val);
            max_error = std::max(max_error, error);
        }

        return max_error;
    }
};

} // namespace complex_analysis
} // namespace maths

#endif // MATHS_COMPLEX_ANALYSIS_HPP
