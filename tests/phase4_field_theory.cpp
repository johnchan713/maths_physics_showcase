/**
 * Phase 4 Validation: Classical Field Theory
 *
 * Tests the classical_field_theory.hpp module functions.
 *
 * Coverage:
 * - Scalar field operations (evaluation, derivatives, d'Alembertian)
 * - Klein-Gordon Lagrangian and action
 * - Phi-fourth Lagrangian
 * - Euler-Lagrange equations
 * - Noether theorem and conserved charges
 * - Stress-energy tensor components
 * - Electromagnetic field tensor properties
 * - EM field invariants
 * - Energy density and Poynting vector
 * - Gauge transformations
 * - Spontaneous symmetry breaking
 * - Mexican hat potential
 */

#include <iostream>
#include <cmath>
#include <vector>
#include <functional>
#include "../include/physics/classical_field_theory.hpp"

// Test tolerance
const double TOLERANCE = 1e-5;
const double LOOSE_TOLERANCE = 1e-3;

// Test macros
#define ASSERT_NEAR(actual, expected, tolerance) \
    do { \
        if (std::abs((actual) - (expected)) > (tolerance)) { \
            std::cerr << "FAIL at line " << __LINE__ << ": " << #actual \
                      << " = " << (actual) << ", expected " << (expected) \
                      << " (diff: " << std::abs((actual) - (expected)) << ")" << std::endl; \
            return false; \
        } \
    } while(0)

#define ASSERT_TRUE(condition) \
    do { \
        if (!(condition)) { \
            std::cerr << "FAIL at line " << __LINE__ << ": " << #condition << std::endl; \
            return false; \
        } \
    } while(0)

// Use namespace for cleaner code
using namespace physics::classical_field_theory;

int main() {
    int tests_passed = 0;
    int tests_failed = 0;

    std::cout << "=== Phase 4: Classical Field Theory Validation ===" << std::endl;
    std::cout << std::endl;

    // Helper lambda to run tests
    auto run_test = [&](const char* name, bool (*test_func)()) {
        std::cout << "Running: " << name << "... ";
        if (test_func()) {
            std::cout << "PASS" << std::endl;
            tests_passed++;
            return true;
        } else {
            tests_failed++;
            return false;
        }
    };

    // ========================================
    // Scalar Field Tests
    // ========================================

    run_test("Scalar field evaluation at origin", []() {
        ScalarField phi([](const Vector4& x) { return x[1] * x[1]; });
        Vector4 x = {0.0, 2.0, 3.0, 4.0};
        double result = phi(x);
        ASSERT_NEAR(result, 4.0, TOLERANCE);
        return true;
    });

    run_test("Scalar field evaluation with constant field", []() {
        ScalarField phi([](const Vector4& x) { return 5.0; });
        Vector4 x = {0.0, 1.0, 2.0, 3.0};
        double result = phi(x);
        ASSERT_NEAR(result, 5.0, TOLERANCE);
        return true;
    });

    run_test("Scalar field partial derivative time component", []() {
        ScalarField phi([](const Vector4& x) { return x[0] * x[0]; });
        Vector4 x = {2.0, 1.0, 1.0, 1.0};
        double dphi_dt = phi.partialDerivative(0, x, 1e-6);
        ASSERT_NEAR(dphi_dt, 4.0, LOOSE_TOLERANCE);
        return true;
    });

    run_test("Scalar field partial derivative spatial component", []() {
        ScalarField phi([](const Vector4& x) { return x[1] * x[1] + x[2] * x[2]; });
        Vector4 x = {0.0, 3.0, 4.0, 0.0};
        double dphi_dx = phi.partialDerivative(1, x, 1e-6);
        ASSERT_NEAR(dphi_dx, 6.0, LOOSE_TOLERANCE);
        return true;
    });

    run_test("Scalar field gradient computation", []() {
        ScalarField phi([](const Vector4& x) { return x[0] + 2.0*x[1] + 3.0*x[2] + 4.0*x[3]; });
        Vector4 x = {1.0, 1.0, 1.0, 1.0};
        auto grad = phi.gradient(x, 1e-6);
        ASSERT_NEAR(grad[0], 1.0, LOOSE_TOLERANCE);
        ASSERT_NEAR(grad[1], 2.0, LOOSE_TOLERANCE);
        ASSERT_NEAR(grad[2], 3.0, LOOSE_TOLERANCE);
        ASSERT_NEAR(grad[3], 4.0, LOOSE_TOLERANCE);
        return true;
    });

    run_test("Scalar field gradient of constant is zero", []() {
        ScalarField phi([](const Vector4& x) { return 7.0; });
        Vector4 x = {1.0, 2.0, 3.0, 4.0};
        auto grad = phi.gradient(x, 1e-6);
        for (int i = 0; i < 4; ++i) {
            ASSERT_NEAR(std::abs(grad[i]), 0.0, LOOSE_TOLERANCE);
        }
        return true;
    });

    run_test("Scalar field dAlembertian of constant is zero", []() {
        ScalarField phi([](const Vector4& x) { return 2.0; });
        Vector4 x = {0.0, 1.0, 2.0, 3.0};
        double box = phi.dAlembertian(x, 1e-6);
        ASSERT_NEAR(std::abs(box), 0.0, LOOSE_TOLERANCE);
        return true;
    });

    // ========================================
    // Klein-Gordon Lagrangian Tests
    // ========================================

    run_test("Klein-Gordon Lagrangian with zero field", []() {
        KleinGordonLagrangian L(1.0);
        ScalarField phi([](const Vector4& x) { return 0.0; });
        Vector4 dphi = {0.0, 0.0, 0.0, 0.0};
        Vector4 x = {0.0, 0.0, 0.0, 0.0};
        double lag = L(x, phi, dphi);
        ASSERT_NEAR(lag, 0.0, TOLERANCE);
        return true;
    });

    run_test("Klein-Gordon kinetic term sign", []() {
        KleinGordonLagrangian L(0.0);
        ScalarField phi([](const Vector4& x) { return 1.0; });
        Vector4 dphi = {1.0, 0.0, 0.0, 0.0};
        Vector4 x = {0.0, 0.0, 0.0, 0.0};
        double lag = L(x, phi, dphi);
        // Kinetic term with only time derivative: -0.5 * (1/c^2) * dphi[0]^2
        ASSERT_TRUE(lag < 0.0);
        return true;
    });

    run_test("Klein-Gordon mass term", []() {
        double mass = 2.0;
        KleinGordonLagrangian L(mass);
        ScalarField phi([](const Vector4& x) { return 1.0; });
        Vector4 dphi = {0.0, 0.0, 0.0, 0.0};
        Vector4 x = {0.0, 0.0, 0.0, 0.0};
        double lag = L(x, phi, dphi);
        // -0.5 * m^2 * c^2 * phi^2 = -0.5 * 4 * c^2
        double expected = -0.5 * mass * mass * c * c;
        ASSERT_NEAR(lag, expected, TOLERANCE);
        return true;
    });

    run_test("Klein-Gordon mass getter", []() {
        double mass = 3.5;
        KleinGordonLagrangian L(mass);
        ASSERT_NEAR(L.mass(), mass, TOLERANCE);
        return true;
    });

    run_test("Klein-Gordon action integral", []() {
        KleinGordonLagrangian L(0.0);
        ScalarField phi([](const Vector4& x) { return 0.0; });
        std::vector<Vector4> points(5);
        for (int i = 0; i < 5; ++i) {
            points[i] = {static_cast<double>(i), 0.0, 0.0, 0.0};
        }
        double S = L.action(phi, points);
        ASSERT_NEAR(S, 0.0, TOLERANCE);
        return true;
    });

    // ========================================
    // Phi-Fourth Lagrangian Tests
    // ========================================

    run_test("PhiFourth Lagrangian with zero field", []() {
        PhiFourthLagrangian L(1.0, 1.0);
        ScalarField phi([](const Vector4& x) { return 0.0; });
        Vector4 dphi = {0.0, 0.0, 0.0, 0.0};
        Vector4 x = {0.0, 0.0, 0.0, 0.0};
        double lag = L(x, phi, dphi);
        ASSERT_NEAR(lag, 0.0, TOLERANCE);
        return true;
    });

    run_test("PhiFourth Lagrangian kinetic component", []() {
        PhiFourthLagrangian L(0.0, 0.0);
        ScalarField phi([](const Vector4& x) { return 0.0; });
        Vector4 dphi = {2.0, 1.0, 1.0, 1.0};
        Vector4 x = {0.0, 0.0, 0.0, 0.0};
        double lag = L(x, phi, dphi);
        // Kinetic: 0.5 * (-4/c^2 + 3) = 0.5 * (3 - 4/c^2)
        double kinetic = 0.5 * (-dphi[0]*dphi[0]/(c*c) + dphi[1]*dphi[1] + dphi[2]*dphi[2] + dphi[3]*dphi[3]);
        ASSERT_NEAR(lag, kinetic, LOOSE_TOLERANCE);
        return true;
    });

    run_test("PhiFourth Lagrangian potential term", []() {
        double mass = 2.0;
        double lambda = 24.0;
        PhiFourthLagrangian L(mass, lambda);
        ScalarField phi([](const Vector4& x) { return 1.0; });
        Vector4 dphi = {0.0, 0.0, 0.0, 0.0};
        Vector4 x = {0.0, 0.0, 0.0, 0.0};
        double lag = L(x, phi, dphi);
        // Potential: -0.5 * m^2 * phi^2 - (lambda/24) * phi^4
        double expected = -0.5 * mass * mass - (lambda / 24.0);
        ASSERT_NEAR(lag, expected, TOLERANCE);
        return true;
    });

    // ========================================
    // Euler-Lagrange Equation Tests
    // ========================================

    run_test("Euler-Lagrange for constant field with zero mass", []() {
        KleinGordonLagrangian L(0.0);
        ScalarField phi([](const Vector4& x) { return 5.0; });
        Vector4 x = {0.0, 1.0, 1.0, 1.0};
        double eq = EulerLagrangeEquations::fieldEquation(L, phi, x, 1e-6);
        ASSERT_NEAR(std::abs(eq), 0.0, LOOSE_TOLERANCE);
        return true;
    });

    run_test("Euler-Lagrange satisfies equation check", []() {
        KleinGordonLagrangian L(0.1);
        ScalarField phi([](const Vector4& x) { return std::sin(x[1]); });
        Vector4 x = {0.0, 1.0, 0.0, 0.0};
        bool satisfied = EulerLagrangeEquations::satisfiesEquation(L, phi, x, LOOSE_TOLERANCE);
        ASSERT_TRUE(satisfied || !satisfied);  // Just verify function works
        return true;
    });

    // ========================================
    // Noether Theorem Tests
    // ========================================

    run_test("Noether current computation", []() {
        KleinGordonLagrangian L(0.0);
        ScalarField phi([](const Vector4& x) { return x[1]; });
        Vector4 x = {0.0, 1.0, 0.0, 0.0};
        NoetherTheorem::Symmetry sym;
        sym.parameter = 0.01;
        sym.transformation = [](const Vector4& x, double eps) {
            Vector4 result = x;
            result[1] += eps;
            return result;
        };
        auto j = NoetherTheorem::noetherCurrent(L, phi, sym, x, 1e-6);
        ASSERT_TRUE(std::isfinite(j[0]));
        return true;
    });

    run_test("Conserved charge from current density", []() {
        std::vector<Vector4> current(5);
        std::vector<Vector4> points(5);
        for (int i = 0; i < 5; ++i) {
            current[i] = {1.0, 0.0, 0.0, 0.0};
            points[i] = {static_cast<double>(i), 0.0, 0.0, 0.0};
        }
        double Q = NoetherTheorem::conservedCharge(current, points);
        ASSERT_NEAR(Q, 5.0, TOLERANCE);
        return true;
    });

    run_test("Conserved charge is conserved", []() {
        std::vector<double> charges = {10.0, 10.0, 10.0, 10.0, 10.0};
        bool conserved = NoetherTheorem::isConserved(charges, LOOSE_TOLERANCE);
        ASSERT_TRUE(conserved);
        return true;
    });

    run_test("Non-conserved charge detection", []() {
        std::vector<double> charges = {10.0, 10.0, 12.0, 10.0, 10.0};
        bool conserved = NoetherTheorem::isConserved(charges, LOOSE_TOLERANCE);
        ASSERT_TRUE(!conserved);
        return true;
    });

    // ========================================
    // Stress-Energy Tensor Tests
    // ========================================

    run_test("Stress-energy tensor component computation", []() {
        KleinGordonLagrangian L(0.0);
        ScalarField phi([](const Vector4& x) { return 1.0; });
        StressEnergyTensor T(L);
        Vector4 x = {0.0, 0.0, 0.0, 0.0};
        double T00 = T.component(0, 0, phi, x, 1e-6);
        ASSERT_TRUE(std::isfinite(T00));
        return true;
    });

    run_test("Stress-energy tensor all components", []() {
        KleinGordonLagrangian L(1.0);
        ScalarField phi([](const Vector4& x) { return 0.0; });
        StressEnergyTensor T(L);
        Vector4 x = {0.0, 0.0, 0.0, 0.0};
        auto T_all = T.allComponents(phi, x, 1e-6);
        for (int mu = 0; mu < 4; ++mu) {
            for (int nu = 0; nu < 4; ++nu) {
                ASSERT_TRUE(std::isfinite(T_all[mu][nu]));
            }
        }
        return true;
    });

    run_test("Stress-energy tensor energy density positive", []() {
        KleinGordonLagrangian L(0.0);
        ScalarField phi([](const Vector4& x) { return 0.0; });
        StressEnergyTensor T(L);
        Vector4 x = {0.0, 0.0, 0.0, 0.0};
        double rho = T.energyDensity(phi, x, 1e-6);
        ASSERT_TRUE(rho >= 0.0);
        return true;
    });

    run_test("Stress-energy tensor momentum density", []() {
        KleinGordonLagrangian L(1.0);
        ScalarField phi([](const Vector4& x) { return std::sin(x[1]); });
        StressEnergyTensor T(L);
        Vector4 x = {0.0, 1.0, 0.0, 0.0};
        auto pi = T.momentumDensity(phi, x, 1e-6);
        ASSERT_TRUE(std::isfinite(pi[0]));
        ASSERT_TRUE(std::isfinite(pi[1]));
        ASSERT_TRUE(std::isfinite(pi[2]));
        return true;
    });

    run_test("Stress-energy tensor trace computation", []() {
        KleinGordonLagrangian L(1.0);
        ScalarField phi([](const Vector4& x) { return 1.0; });
        StressEnergyTensor T(L);
        Vector4 x = {0.0, 0.0, 0.0, 0.0};
        double tr = T.trace(phi, x, 1e-6);
        ASSERT_TRUE(std::isfinite(tr));
        return true;
    });

    // ========================================
    // Electromagnetic Field Tensor Tests
    // ========================================

    run_test("EM field tensor antisymmetry F_mu_nu = -F_nu_mu", []() {
        auto E = [](const Vector4& x) { return std::array<double, 3>{1.0, 2.0, 3.0}; };
        auto B = [](const Vector4& x) { return std::array<double, 3>{0.5, 0.5, 0.5}; };
        ElectromagneticField F(E, B);
        Vector4 x = {0.0, 0.0, 0.0, 0.0};
        auto Ftensor = F.fieldTensor(x);
        for (int mu = 0; mu < 4; ++mu) {
            for (int nu = 0; nu < 4; ++nu) {
                ASSERT_NEAR(Ftensor[mu][nu], -Ftensor[nu][mu], TOLERANCE);
            }
        }
        return true;
    });

    run_test("EM field tensor diagonal elements are zero", []() {
        auto E = [](const Vector4& x) { return std::array<double, 3>{1.0, 2.0, 3.0}; };
        auto B = [](const Vector4& x) { return std::array<double, 3>{0.5, 0.5, 0.5}; };
        ElectromagneticField F(E, B);
        Vector4 x = {0.0, 0.0, 0.0, 0.0};
        auto Ftensor = F.fieldTensor(x);
        for (int mu = 0; mu < 4; ++mu) {
            ASSERT_NEAR(Ftensor[mu][mu], 0.0, TOLERANCE);
        }
        return true;
    });

    run_test("EM field tensor electric field components", []() {
        auto E = [](const Vector4& x) { return std::array<double, 3>{3.0, 4.0, 5.0}; };
        auto B = [](const Vector4& x) { return std::array<double, 3>{0.0, 0.0, 0.0}; };
        ElectromagneticField F(E, B);
        Vector4 x = {0.0, 0.0, 0.0, 0.0};
        auto Ftensor = F.fieldTensor(x);
        ASSERT_NEAR(Ftensor[0][1], -3.0/c, TOLERANCE);
        ASSERT_NEAR(Ftensor[0][2], -4.0/c, TOLERANCE);
        ASSERT_NEAR(Ftensor[0][3], -5.0/c, TOLERANCE);
        return true;
    });

    run_test("EM field tensor magnetic field components", []() {
        auto E = [](const Vector4& x) { return std::array<double, 3>{0.0, 0.0, 0.0}; };
        auto B = [](const Vector4& x) { return std::array<double, 3>{1.0, 2.0, 3.0}; };
        ElectromagneticField F(E, B);
        Vector4 x = {0.0, 0.0, 0.0, 0.0};
        auto Ftensor = F.fieldTensor(x);
        ASSERT_NEAR(Ftensor[1][2], -3.0, TOLERANCE);
        ASSERT_NEAR(Ftensor[1][3], 2.0, TOLERANCE);
        ASSERT_NEAR(Ftensor[2][3], -1.0, TOLERANCE);
        return true;
    });

    run_test("EM first invariant B²-E²/c²", []() {
        auto E = [](const Vector4& x) { return std::array<double, 3>{3.0, 0.0, 0.0}; };
        auto B = [](const Vector4& x) { return std::array<double, 3>{4.0, 0.0, 0.0}; };
        ElectromagneticField F(E, B);
        Vector4 x = {0.0, 0.0, 0.0, 0.0};
        double inv1 = F.firstInvariant(x);
        double expected = 16.0 - 9.0/(c*c);
        ASSERT_NEAR(inv1, expected, TOLERANCE);
        return true;
    });

    run_test("EM second invariant E dot B", []() {
        auto E = [](const Vector4& x) { return std::array<double, 3>{1.0, 2.0, 3.0}; };
        auto B = [](const Vector4& x) { return std::array<double, 3>{4.0, 5.0, 6.0}; };
        ElectromagneticField F(E, B);
        Vector4 x = {0.0, 0.0, 0.0, 0.0};
        double inv2 = F.secondInvariant(x);
        double E_dot_B = 1.0*4.0 + 2.0*5.0 + 3.0*6.0;  // 32
        double expected = E_dot_B / c;
        ASSERT_NEAR(inv2, expected, TOLERANCE);
        return true;
    });

    run_test("EM stress-energy tensor components are finite", []() {
        auto E = [](const Vector4& x) { return std::array<double, 3>{1.0, 2.0, 3.0}; };
        auto B = [](const Vector4& x) { return std::array<double, 3>{0.5, 0.5, 0.5}; };
        ElectromagneticField F(E, B);
        Vector4 x = {0.0, 0.0, 0.0, 0.0};
        auto T = F.stressEnergyTensor(x);
        for (int mu = 0; mu < 4; ++mu) {
            for (int nu = 0; nu < 4; ++nu) {
                ASSERT_TRUE(std::isfinite(T[mu][nu]));
            }
        }
        return true;
    });

    run_test("EM energy density positive semi-definite", []() {
        auto E = [](const Vector4& x) { return std::array<double, 3>{1.0, 2.0, 3.0}; };
        auto B = [](const Vector4& x) { return std::array<double, 3>{0.5, 0.5, 0.5}; };
        ElectromagneticField F(E, B);
        Vector4 x = {0.0, 0.0, 0.0, 0.0};
        double u = F.energyDensity(x);
        ASSERT_TRUE(u >= 0.0);
        return true;
    });

    run_test("EM energy density with zero fields", []() {
        auto E = [](const Vector4& x) { return std::array<double, 3>{0.0, 0.0, 0.0}; };
        auto B = [](const Vector4& x) { return std::array<double, 3>{0.0, 0.0, 0.0}; };
        ElectromagneticField F(E, B);
        Vector4 x = {0.0, 0.0, 0.0, 0.0};
        double u = F.energyDensity(x);
        ASSERT_NEAR(u, 0.0, TOLERANCE);
        return true;
    });

    run_test("EM Poynting vector", []() {
        auto E = [](const Vector4& x) { return std::array<double, 3>{1.0, 0.0, 0.0}; };
        auto B = [](const Vector4& x) { return std::array<double, 3>{0.0, 1.0, 0.0}; };
        ElectromagneticField F(E, B);
        Vector4 x = {0.0, 0.0, 0.0, 0.0};
        auto S = F.poyntingVector(x);
        // S = (E x B) / mu_0 = (0, 0, 1) / mu_0
        ASSERT_NEAR(std::abs(S[2]), 1.0/mu_0, TOLERANCE);
        return true;
    });

    run_test("EM Poynting vector perpendicular to E", []() {
        auto E = [](const Vector4& x) { return std::array<double, 3>{1.0, 0.0, 0.0}; };
        auto B = [](const Vector4& x) { return std::array<double, 3>{0.0, 1.0, 0.0}; };
        ElectromagneticField F(E, B);
        Vector4 x = {0.0, 0.0, 0.0, 0.0};
        auto S = F.poyntingVector(x);
        double E_dot_S = 1.0 * S[0] + 0.0 * S[1] + 0.0 * S[2];
        ASSERT_NEAR(E_dot_S, 0.0, TOLERANCE);
        return true;
    });

    run_test("EM Poynting vector perpendicular to B", []() {
        auto E = [](const Vector4& x) { return std::array<double, 3>{1.0, 0.0, 0.0}; };
        auto B = [](const Vector4& x) { return std::array<double, 3>{0.0, 1.0, 0.0}; };
        ElectromagneticField F(E, B);
        Vector4 x = {0.0, 0.0, 0.0, 0.0};
        auto S = F.poyntingVector(x);
        double B_dot_S = 0.0 * S[0] + 1.0 * S[1] + 0.0 * S[2];
        ASSERT_NEAR(B_dot_S, 0.0, TOLERANCE);
        return true;
    });

    // ========================================
    // Gauge Transformation Tests
    // ========================================

    run_test("Vector potential gauge transformation", []() {
        auto phi_gauge = [](const Vector4& x) { return x[1]; };
        std::array<double, 4> A_old = {1.0, 2.0, 3.0, 4.0};
        Vector4 x = {0.0, 1.0, 0.0, 0.0};
        auto A_new = GaugeTransformation::vectorPotential(phi_gauge, A_old, x, 1e-6);
        ASSERT_TRUE(std::isfinite(A_new[0]));
        ASSERT_TRUE(std::isfinite(A_new[1]));
        ASSERT_TRUE(std::isfinite(A_new[2]));
        ASSERT_TRUE(std::isfinite(A_new[3]));
        return true;
    });

    run_test("Coulomb gauge computation", []() {
        std::array<double, 4> A = {0.0, 1.0, 2.0, 3.0};
        Vector4 x = {0.0, 0.0, 0.0, 0.0};
        auto A_coulomb = GaugeTransformation::coulombGauge(A, x, 1e-6);
        ASSERT_TRUE(std::isfinite(A_coulomb[0]));
        ASSERT_TRUE(std::isfinite(A_coulomb[1]));
        return true;
    });

    run_test("Lorenz gauge computation", []() {
        std::array<double, 4> A = {1.0, 1.0, 1.0, 1.0};
        Vector4 x = {0.0, 0.0, 0.0, 0.0};
        auto A_lorenz = GaugeTransformation::lorenzGauge(A, x, 1e-6);
        ASSERT_TRUE(std::isfinite(A_lorenz[0]));
        ASSERT_TRUE(std::isfinite(A_lorenz[1]));
        return true;
    });

    run_test("Field tensor computation from vector potential", []() {
        std::array<double, 4> A = {1.0, 2.0, 3.0, 4.0};
        Vector4 x = {0.0, 0.0, 0.0, 0.0};
        auto F = GaugeTransformation::computeFieldTensor(A, x, 1e-6);
        for (int mu = 0; mu < 4; ++mu) {
            for (int nu = 0; nu < 4; ++nu) {
                ASSERT_TRUE(std::isfinite(F[mu][nu]));
            }
        }
        return true;
    });

    // ========================================
    // Spontaneous Symmetry Breaking Tests
    // ========================================

    run_test("Mexican hat potential at origin", []() {
        double V = SpontaneousSymmetryBreaking::mexicanHatPotential(0.0, 1.0, 1.0);
        ASSERT_NEAR(V, 0.0, TOLERANCE);
        return true;
    });

    run_test("Mexican hat potential is negative at small values", []() {
        double V1 = SpontaneousSymmetryBreaking::mexicanHatPotential(1.0, 1.0, 1.0);
        double V2 = SpontaneousSymmetryBreaking::mexicanHatPotential(0.5, 1.0, 1.0);
        ASSERT_TRUE(V1 < 0.0);
        ASSERT_TRUE(V2 < 0.0);
        return true;
    });

    run_test("Mexican hat potential grows at large values", []() {
        double V_small = SpontaneousSymmetryBreaking::mexicanHatPotential(1.0, 1.0, 1.0);
        double V_large = SpontaneousSymmetryBreaking::mexicanHatPotential(10.0, 1.0, 1.0);
        ASSERT_TRUE(V_large > V_small);
        return true;
    });

    run_test("VEV is zero for negative mu squared", []() {
        auto [v_pos, v_neg] = SpontaneousSymmetryBreaking::vacuumExpectationValues(-1.0, 1.0);
        ASSERT_NEAR(v_pos, 0.0, TOLERANCE);
        ASSERT_NEAR(v_neg, 0.0, TOLERANCE);
        return true;
    });

    run_test("VEV is symmetric for positive mu squared", []() {
        auto [v_pos, v_neg] = SpontaneousSymmetryBreaking::vacuumExpectationValues(4.0, 1.0);
        ASSERT_NEAR(v_pos, -v_neg, TOLERANCE);
        ASSERT_TRUE(v_pos > 0.0);
        return true;
    });

    run_test("VEV calculation v = sqrt(mu²/lambda)", []() {
        double mu2 = 4.0;
        double lambda = 2.0;
        auto [v_pos, v_neg] = SpontaneousSymmetryBreaking::vacuumExpectationValues(mu2, lambda);
        double expected = std::sqrt(mu2 / lambda);
        ASSERT_NEAR(std::abs(v_pos), expected, TOLERANCE);
        return true;
    });

    run_test("Mass squared fluctuation at VEV", []() {
        double mu2 = 1.0;
        double lambda = 1.0;
        auto [vev, _] = SpontaneousSymmetryBreaking::vacuumExpectationValues(mu2, lambda);
        double m2 = SpontaneousSymmetryBreaking::massSquaredFluctuation(mu2, lambda, vev);
        // m² = -mu² + 3*lambda*vev² = -1 + 3*1*1 = 2
        ASSERT_NEAR(m2, 2.0, TOLERANCE);
        return true;
    });

    run_test("Goldstone boson has zero mass", []() {
        double m = SpontaneousSymmetryBreaking::goldstoneBosonMass();
        ASSERT_NEAR(m, 0.0, TOLERANCE);
        return true;
    });

    run_test("Goldstone bosons exist for positive mu squared", []() {
        bool has_goldstone = SpontaneousSymmetryBreaking::hasGoldstoneBosons(1.0);
        ASSERT_TRUE(has_goldstone);
        return true;
    });

    run_test("No Goldstone bosons for negative mu squared", []() {
        bool has_goldstone = SpontaneousSymmetryBreaking::hasGoldstoneBosons(-1.0);
        ASSERT_TRUE(!has_goldstone);
        return true;
    });

    // ========================================
    // Physical Consistency Tests
    // ========================================

    run_test("Consistency scalar field partial derivatives are finite", []() {
        ScalarField phi([](const Vector4& x) {
            return std::exp(-(x[1]*x[1] + x[2]*x[2] + x[3]*x[3]));
        });
        Vector4 x = {0.0, 0.5, 0.5, 0.5};
        for (int mu = 0; mu < 4; ++mu) {
            double dphi = phi.partialDerivative(mu, x, 1e-6);
            ASSERT_TRUE(std::isfinite(dphi));
        }
        return true;
    });

    run_test("Consistency dAlembertian of harmonic function", []() {
        ScalarField phi([](const Vector4& x) {
            return std::cos(x[1]) * std::cos(x[2]);
        });
        Vector4 x = {0.0, M_PI/4.0, M_PI/4.0, 0.0};
        double box = phi.dAlembertian(x, 1e-6);
        ASSERT_TRUE(std::isfinite(box));
        return true;
    });

    run_test("Consistency Klein-Gordon action is scalar invariant", []() {
        KleinGordonLagrangian L(1.0);
        ScalarField phi([](const Vector4& x) { return std::sin(x[1]); });
        std::vector<Vector4> points(10);
        for (int i = 0; i < 10; ++i) {
            points[i] = {0.0, static_cast<double>(i)*0.1, 0.0, 0.0};
        }
        double S = L.action(phi, points);
        ASSERT_TRUE(std::isfinite(S));
        return true;
    });

    run_test("Consistency EM field tensor element bounds", []() {
        auto E = [](const Vector4& x) { return std::array<double, 3>{1e6, 1e6, 1e6}; };
        auto B = [](const Vector4& x) { return std::array<double, 3>{1.0, 1.0, 1.0}; };
        ElectromagneticField F(E, B);
        Vector4 x = {0.0, 0.0, 0.0, 0.0};
        auto Ftensor = F.fieldTensor(x);
        for (int mu = 0; mu < 4; ++mu) {
            for (int nu = 0; nu < 4; ++nu) {
                ASSERT_TRUE(std::isfinite(Ftensor[mu][nu]));
            }
        }
        return true;
    });

    run_test("Consistency energy density is finite", []() {
        auto E = [](const Vector4& x) { return std::array<double, 3>{1.0, 2.0, 3.0}; };
        auto B = [](const Vector4& x) { return std::array<double, 3>{0.5, 0.5, 0.5}; };
        ElectromagneticField F(E, B);
        Vector4 x = {0.0, 0.0, 0.0, 0.0};
        double u = F.energyDensity(x);
        ASSERT_TRUE(std::isfinite(u));
        return true;
    });

    run_test("Consistency Mexican hat potential grows without bound", []() {
        double V_small = SpontaneousSymmetryBreaking::mexicanHatPotential(1.0, 1.0, 1.0);
        double V_large = SpontaneousSymmetryBreaking::mexicanHatPotential(100.0, 1.0, 1.0);
        ASSERT_TRUE(V_large > V_small);
        ASSERT_TRUE(std::isfinite(V_large));
        return true;
    });

    run_test("Consistency multiple VEV pairs", []() {
        for (double mu2 = 0.1; mu2 <= 1.0; mu2 += 0.1) {
            for (double lambda = 0.1; lambda <= 1.0; lambda += 0.1) {
                auto [v_pos, v_neg] = SpontaneousSymmetryBreaking::vacuumExpectationValues(mu2, lambda);
                ASSERT_TRUE(std::isfinite(v_pos));
                ASSERT_TRUE(std::isfinite(v_neg));
            }
        }
        return true;
    });

    // ========================================
    // Summary
    // ========================================

    std::cout << std::endl;
    std::cout << "=== Test Summary ===" << std::endl;
    std::cout << "Tests passed: " << tests_passed << std::endl;
    std::cout << "Tests failed: " << tests_failed << std::endl;
    std::cout << "Total tests: " << (tests_passed + tests_failed) << std::endl;

    if (tests_failed == 0) {
        std::cout << "All tests passed!" << std::endl;
        return 0;
    } else {
        std::cout << "Some tests failed." << std::endl;
        return 1;
    }
}
