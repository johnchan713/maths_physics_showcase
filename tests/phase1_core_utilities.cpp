/**
 * @file phase1_core_utilities.cpp
 * @brief Phase 1 Validation: Core Utilities (Vectors, Matrices, Units)
 *
 * This file validates the correctness of:
 * - Vector operations (dot, cross, norm, normalize, projections)
 * - Matrix operations (multiply, transpose, determinant, inverse, rank, solve)
 * - Unit conversions (force, mass, length, velocity, energy, acceleration, angles)
 *
 * Validation approach:
 * - Test with known mathematical values
 * - Verify mathematical properties and identities
 * - Check edge cases and error handling
 * - Ensure numerical stability within tolerance (1e-6)
 */

#include <iostream>
#include <iomanip>
#include <cmath>
#include <cassert>
#include <string>
#include <vector>

#include "../include/maths/vectors.hpp"
#include "../include/maths/matrices.hpp"
#include "../include/physics/units.hpp"

using namespace maths::linear_algebra;
using namespace physics::units;

// Test tolerance
constexpr double TOLERANCE = 1e-6;

// Test counter
int tests_passed = 0;
int tests_failed = 0;

// Helper macros for testing
#define ASSERT_NEAR(actual, expected, tolerance) \
    do { \
        if (std::abs((actual) - (expected)) <= (tolerance)) { \
            tests_passed++; \
        } else { \
            tests_failed++; \
            std::cerr << "FAILED: " << __FILE__ << ":" << __LINE__ << std::endl; \
            std::cerr << "  Expected: " << (expected) << std::endl; \
            std::cerr << "  Actual:   " << (actual) << std::endl; \
            std::cerr << "  Diff:     " << std::abs((actual) - (expected)) << std::endl; \
        } \
    } while(0)

#define ASSERT_TRUE(condition) \
    do { \
        if (condition) { \
            tests_passed++; \
        } else { \
            tests_failed++; \
            std::cerr << "FAILED: " << __FILE__ << ":" << __LINE__ << std::endl; \
            std::cerr << "  Condition failed: " << #condition << std::endl; \
        } \
    } while(0)

#define ASSERT_FALSE(condition) ASSERT_TRUE(!(condition))

// ============================================================================
// VECTOR TESTS
// ============================================================================

void test_vector_basic_operations() {
    std::cout << "\n=== Testing Vector Basic Operations ===" << std::endl;

    // Test construction and access
    Vector v({1.0, 2.0, 3.0});
    ASSERT_NEAR(v[0], 1.0, TOLERANCE);
    ASSERT_NEAR(v[1], 2.0, TOLERANCE);
    ASSERT_NEAR(v[2], 3.0, TOLERANCE);
    ASSERT_TRUE(v.dimension() == 3);

    // Test addition
    Vector w({4.0, 5.0, 6.0});
    Vector sum = v + w;
    ASSERT_NEAR(sum[0], 5.0, TOLERANCE);
    ASSERT_NEAR(sum[1], 7.0, TOLERANCE);
    ASSERT_NEAR(sum[2], 9.0, TOLERANCE);

    // Test subtraction
    Vector diff = w - v;
    ASSERT_NEAR(diff[0], 3.0, TOLERANCE);
    ASSERT_NEAR(diff[1], 3.0, TOLERANCE);
    ASSERT_NEAR(diff[2], 3.0, TOLERANCE);

    // Test scalar multiplication
    Vector scaled = v * 2.0;
    ASSERT_NEAR(scaled[0], 2.0, TOLERANCE);
    ASSERT_NEAR(scaled[1], 4.0, TOLERANCE);
    ASSERT_NEAR(scaled[2], 6.0, TOLERANCE);

    // Test negation
    Vector neg = -v;
    ASSERT_NEAR(neg[0], -1.0, TOLERANCE);
    ASSERT_NEAR(neg[1], -2.0, TOLERANCE);
    ASSERT_NEAR(neg[2], -3.0, TOLERANCE);

    std::cout << "Basic vector operations: OK" << std::endl;
}

void test_vector_dot_product() {
    std::cout << "\n=== Testing Vector Dot Product ===" << std::endl;

    // Test dot product
    Vector v({1.0, 2.0, 3.0});
    Vector w({4.0, 5.0, 6.0});
    double dot = v.dot(w);
    ASSERT_NEAR(dot, 1*4 + 2*5 + 3*6, TOLERANCE); // = 4 + 10 + 18 = 32

    // Test orthogonal vectors (dot product = 0)
    Vector i({1.0, 0.0, 0.0});
    Vector j({0.0, 1.0, 0.0});
    ASSERT_NEAR(i.dot(j), 0.0, TOLERANCE);

    // Test dot product with self = norm^2
    Vector u({3.0, 4.0});
    double norm_squared = u.dot(u);
    ASSERT_NEAR(norm_squared, 25.0, TOLERANCE); // 3^2 + 4^2 = 25

    std::cout << "Dot product: OK" << std::endl;
}

void test_vector_norm_and_normalize() {
    std::cout << "\n=== Testing Vector Norm and Normalize ===" << std::endl;

    // Test norm (3-4-5 triangle)
    Vector v({3.0, 4.0});
    ASSERT_NEAR(v.norm(), 5.0, TOLERANCE);

    // Test norm (unit vectors)
    Vector i({1.0, 0.0, 0.0});
    ASSERT_NEAR(i.norm(), 1.0, TOLERANCE);

    // Test normalize
    Vector normalized = v.normalize();
    ASSERT_NEAR(normalized.norm(), 1.0, TOLERANCE);
    ASSERT_NEAR(normalized[0], 3.0/5.0, TOLERANCE);
    ASSERT_NEAR(normalized[1], 4.0/5.0, TOLERANCE);

    // Test distance
    Vector a({1.0, 2.0, 3.0});
    Vector b({4.0, 6.0, 8.0});
    double dist = a.distance(b);
    // (4-1)^2 + (6-2)^2 + (8-3)^2 = 9 + 16 + 25 = 50
    ASSERT_NEAR(dist, std::sqrt(50.0), TOLERANCE);

    std::cout << "Norm and normalize: OK" << std::endl;
}

void test_vector_cross_product() {
    std::cout << "\n=== Testing Vector Cross Product ===" << std::endl;

    // Test standard basis cross products
    Vector i({1.0, 0.0, 0.0});
    Vector j({0.0, 1.0, 0.0});
    Vector k({0.0, 0.0, 1.0});

    // i × j = k
    Vector cross_ij = CrossProduct::compute(i, j);
    ASSERT_NEAR(cross_ij[0], 0.0, TOLERANCE);
    ASSERT_NEAR(cross_ij[1], 0.0, TOLERANCE);
    ASSERT_NEAR(cross_ij[2], 1.0, TOLERANCE);

    // j × k = i
    Vector cross_jk = CrossProduct::compute(j, k);
    ASSERT_NEAR(cross_jk[0], 1.0, TOLERANCE);
    ASSERT_NEAR(cross_jk[1], 0.0, TOLERANCE);
    ASSERT_NEAR(cross_jk[2], 0.0, TOLERANCE);

    // k × i = j
    Vector cross_ki = CrossProduct::compute(k, i);
    ASSERT_NEAR(cross_ki[0], 0.0, TOLERANCE);
    ASSERT_NEAR(cross_ki[1], 1.0, TOLERANCE);
    ASSERT_NEAR(cross_ki[2], 0.0, TOLERANCE);

    // Test anticommutativity: v × w = -(w × v)
    Vector v({1.0, 2.0, 3.0});
    Vector w({4.0, 5.0, 6.0});
    Vector vxw = CrossProduct::compute(v, w);
    Vector wxv = CrossProduct::compute(w, v);
    ASSERT_NEAR(vxw[0], -wxv[0], TOLERANCE);
    ASSERT_NEAR(vxw[1], -wxv[1], TOLERANCE);
    ASSERT_NEAR(vxw[2], -wxv[2], TOLERANCE);

    // Test perpendicularity: (v × w) · v = 0
    ASSERT_NEAR(vxw.dot(v), 0.0, TOLERANCE);
    ASSERT_NEAR(vxw.dot(w), 0.0, TOLERANCE);

    // Test parallel vectors: v × v = 0
    Vector self_cross = CrossProduct::compute(v, v);
    ASSERT_TRUE(self_cross.isZero(TOLERANCE));

    std::cout << "Cross product: OK" << std::endl;
}

void test_vector_projections() {
    std::cout << "\n=== Testing Vector Projections ===" << std::endl;

    Vector v({3.0, 4.0});
    Vector w({1.0, 0.0});

    // Scalar projection of v onto w
    double scalarProj = VectorProjection::scalarProjection(v, w);
    ASSERT_NEAR(scalarProj, 3.0, TOLERANCE); // v·w / |w| = 3/1 = 3

    // Vector projection of v onto w
    Vector vectorProj = VectorProjection::vectorProjection(v, w);
    ASSERT_NEAR(vectorProj[0], 3.0, TOLERANCE);
    ASSERT_NEAR(vectorProj[1], 0.0, TOLERANCE);

    // Orthogonal component
    Vector perp = VectorProjection::orthogonalComponent(v, w);
    ASSERT_NEAR(perp[0], 0.0, TOLERANCE);
    ASSERT_NEAR(perp[1], 4.0, TOLERANCE);

    // Verify decomposition: v = proj + perp
    ASSERT_TRUE(VectorProjection::verifyDecomposition(v, w, TOLERANCE));

    // Verify orthogonality: proj · perp = 0
    ASSERT_NEAR(vectorProj.dot(perp), 0.0, TOLERANCE);

    // Verify Pythagorean theorem: |v|^2 = |proj|^2 + |perp|^2
    double v_norm_sq = v.norm() * v.norm();
    double proj_norm_sq = vectorProj.norm() * vectorProj.norm();
    double perp_norm_sq = perp.norm() * perp.norm();
    ASSERT_NEAR(v_norm_sq, proj_norm_sq + perp_norm_sq, TOLERANCE);

    std::cout << "Projections: OK" << std::endl;
}

void test_vector_orthogonality() {
    std::cout << "\n=== Testing Vector Orthogonality ===" << std::endl;

    // Standard basis is orthonormal
    Vector i({1.0, 0.0, 0.0});
    Vector j({0.0, 1.0, 0.0});
    Vector k({0.0, 0.0, 1.0});

    ASSERT_TRUE(Orthogonality::areOrthogonal(i, j, TOLERANCE));
    ASSERT_TRUE(Orthogonality::areOrthogonal(j, k, TOLERANCE));
    ASSERT_TRUE(Orthogonality::areOrthogonal(k, i, TOLERANCE));

    ASSERT_TRUE(Orthogonality::isUnitVector(i, TOLERANCE));
    ASSERT_TRUE(Orthogonality::isUnitVector(j, TOLERANCE));
    ASSERT_TRUE(Orthogonality::isUnitVector(k, TOLERANCE));

    std::vector<Vector> basis = {i, j, k};
    ASSERT_TRUE(Orthogonality::areOrthonormal(basis, TOLERANCE));

    // Test Gram-Schmidt
    Vector v1({1.0, 1.0, 0.0});
    Vector v2({1.0, 0.0, 1.0});
    Vector v3({0.0, 1.0, 1.0});

    std::vector<Vector> input = {v1, v2, v3};
    std::vector<Vector> orthogonal = Orthogonality::gramSchmidt(input);

    // Check orthogonality
    for (size_t i = 0; i < orthogonal.size(); ++i) {
        for (size_t j = i + 1; j < orthogonal.size(); ++j) {
            ASSERT_TRUE(Orthogonality::areOrthogonal(orthogonal[i], orthogonal[j], TOLERANCE));
        }
    }

    // Test orthonormal Gram-Schmidt
    std::vector<Vector> orthonormal = Orthogonality::gramSchmidtOrthonormal(input);
    ASSERT_TRUE(Orthogonality::areOrthonormal(orthonormal, TOLERANCE));

    std::cout << "Orthogonality: OK" << std::endl;
}

void test_linear_independence() {
    std::cout << "\n=== Testing Linear Independence ===" << std::endl;

    // Standard basis is independent
    std::vector<Vector> basis = Orthogonality::standardBasis(3);
    ASSERT_TRUE(LinearIndependence::areIndependent(basis, TOLERANCE));
    ASSERT_TRUE(LinearIndependence::spanDimension(basis) == 3);

    // Dependent vectors
    Vector v1({1.0, 0.0});
    Vector v2({2.0, 0.0}); // v2 = 2*v1, so dependent
    std::vector<Vector> dependent = {v1, v2};
    ASSERT_FALSE(LinearIndependence::areIndependent(dependent, TOLERANCE));

    // Independent vectors
    Vector u1({1.0, 0.0});
    Vector u2({0.0, 1.0});
    std::vector<Vector> independent = {u1, u2};
    ASSERT_TRUE(LinearIndependence::areIndependent(independent, TOLERANCE));

    std::cout << "Linear independence: OK" << std::endl;
}

// ============================================================================
// MATRIX TESTS
// ============================================================================

void test_matrix_basic_operations() {
    std::cout << "\n=== Testing Matrix Basic Operations ===" << std::endl;

    // Test construction
    Matrix A({{1.0, 2.0}, {3.0, 4.0}});
    ASSERT_TRUE(A.rows() == 2);
    ASSERT_TRUE(A.cols() == 2);
    ASSERT_NEAR(A(0, 0), 1.0, TOLERANCE);
    ASSERT_NEAR(A(0, 1), 2.0, TOLERANCE);
    ASSERT_NEAR(A(1, 0), 3.0, TOLERANCE);
    ASSERT_NEAR(A(1, 1), 4.0, TOLERANCE);

    // Test addition
    Matrix B({{5.0, 6.0}, {7.0, 8.0}});
    Matrix sum = A + B;
    ASSERT_NEAR(sum(0, 0), 6.0, TOLERANCE);
    ASSERT_NEAR(sum(0, 1), 8.0, TOLERANCE);
    ASSERT_NEAR(sum(1, 0), 10.0, TOLERANCE);
    ASSERT_NEAR(sum(1, 1), 12.0, TOLERANCE);

    // Test subtraction
    Matrix diff = B - A;
    ASSERT_NEAR(diff(0, 0), 4.0, TOLERANCE);
    ASSERT_NEAR(diff(0, 1), 4.0, TOLERANCE);
    ASSERT_NEAR(diff(1, 0), 4.0, TOLERANCE);
    ASSERT_NEAR(diff(1, 1), 4.0, TOLERANCE);

    // Test scalar multiplication
    Matrix scaled = A * 2.0;
    ASSERT_NEAR(scaled(0, 0), 2.0, TOLERANCE);
    ASSERT_NEAR(scaled(0, 1), 4.0, TOLERANCE);
    ASSERT_NEAR(scaled(1, 0), 6.0, TOLERANCE);
    ASSERT_NEAR(scaled(1, 1), 8.0, TOLERANCE);

    std::cout << "Basic matrix operations: OK" << std::endl;
}

void test_matrix_multiplication() {
    std::cout << "\n=== Testing Matrix Multiplication ===" << std::endl;

    // Test matrix multiplication
    Matrix A({{1.0, 2.0}, {3.0, 4.0}});
    Matrix B({{5.0, 6.0}, {7.0, 8.0}});
    Matrix product = A * B;

    // [1 2] * [5 6] = [1*5+2*7  1*6+2*8] = [19 22]
    // [3 4]   [7 8]   [3*5+4*7  3*6+4*8]   [43 50]
    ASSERT_NEAR(product(0, 0), 19.0, TOLERANCE);
    ASSERT_NEAR(product(0, 1), 22.0, TOLERANCE);
    ASSERT_NEAR(product(1, 0), 43.0, TOLERANCE);
    ASSERT_NEAR(product(1, 1), 50.0, TOLERANCE);

    // Test identity property: A * I = A
    Matrix I = Matrix::identity(2);
    Matrix AI = A * I;
    ASSERT_NEAR(AI(0, 0), A(0, 0), TOLERANCE);
    ASSERT_NEAR(AI(0, 1), A(0, 1), TOLERANCE);
    ASSERT_NEAR(AI(1, 0), A(1, 0), TOLERANCE);
    ASSERT_NEAR(AI(1, 1), A(1, 1), TOLERANCE);

    // Test matrix-vector multiplication
    Vector v({1.0, 2.0});
    Vector Av = A * v;
    // [1 2] * [1] = [5]
    // [3 4]   [2]   [11]
    ASSERT_NEAR(Av[0], 5.0, TOLERANCE);
    ASSERT_NEAR(Av[1], 11.0, TOLERANCE);

    std::cout << "Matrix multiplication: OK" << std::endl;
}

void test_matrix_transpose() {
    std::cout << "\n=== Testing Matrix Transpose ===" << std::endl;

    Matrix A({{1.0, 2.0, 3.0}, {4.0, 5.0, 6.0}});
    Matrix AT = A.transpose();

    ASSERT_TRUE(AT.rows() == 3);
    ASSERT_TRUE(AT.cols() == 2);
    ASSERT_NEAR(AT(0, 0), 1.0, TOLERANCE);
    ASSERT_NEAR(AT(0, 1), 4.0, TOLERANCE);
    ASSERT_NEAR(AT(1, 0), 2.0, TOLERANCE);
    ASSERT_NEAR(AT(1, 1), 5.0, TOLERANCE);
    ASSERT_NEAR(AT(2, 0), 3.0, TOLERANCE);
    ASSERT_NEAR(AT(2, 1), 6.0, TOLERANCE);

    // Test (A^T)^T = A
    Matrix ATT = AT.transpose();
    ASSERT_TRUE(ATT.rows() == A.rows());
    ASSERT_TRUE(ATT.cols() == A.cols());
    for (size_t i = 0; i < A.rows(); ++i) {
        for (size_t j = 0; j < A.cols(); ++j) {
            ASSERT_NEAR(ATT(i, j), A(i, j), TOLERANCE);
        }
    }

    std::cout << "Transpose: OK" << std::endl;
}

void test_matrix_determinant() {
    std::cout << "\n=== Testing Matrix Determinant ===" << std::endl;

    // Test 2x2 determinant
    Matrix A({{1.0, 2.0}, {3.0, 4.0}});
    double det = A.determinant();
    // det = 1*4 - 2*3 = -2
    ASSERT_NEAR(det, -2.0, TOLERANCE);

    // Test 3x3 determinant
    Matrix B({{1.0, 2.0, 3.0},
              {4.0, 5.0, 6.0},
              {7.0, 8.0, 9.0}});
    double det_B = B.determinant();
    // This matrix is singular (rows are linearly dependent)
    ASSERT_NEAR(det_B, 0.0, TOLERANCE);

    // Test 3x3 non-singular
    Matrix C({{1.0, 0.0, 0.0},
              {0.0, 2.0, 0.0},
              {0.0, 0.0, 3.0}});
    double det_C = C.determinant();
    // Diagonal matrix: det = product of diagonal = 6
    ASSERT_NEAR(det_C, 6.0, TOLERANCE);

    // Test identity determinant = 1
    Matrix I = Matrix::identity(3);
    ASSERT_NEAR(I.determinant(), 1.0, TOLERANCE);

    std::cout << "Determinant: OK" << std::endl;
}

void test_matrix_trace() {
    std::cout << "\n=== Testing Matrix Trace ===" << std::endl;

    Matrix A({{1.0, 2.0, 3.0},
              {4.0, 5.0, 6.0},
              {7.0, 8.0, 9.0}});
    double trace = A.trace();
    // trace = 1 + 5 + 9 = 15
    ASSERT_NEAR(trace, 15.0, TOLERANCE);

    // Test trace of identity
    Matrix I = Matrix::identity(4);
    ASSERT_NEAR(I.trace(), 4.0, TOLERANCE);

    std::cout << "Trace: OK" << std::endl;
}

void test_matrix_inverse() {
    std::cout << "\n=== Testing Matrix Inverse ===" << std::endl;

    // Test 2x2 inverse
    Matrix A({{4.0, 7.0}, {2.0, 6.0}});
    Matrix A_inv = A.inverse();

    // Verify A * A^(-1) = I
    Matrix product = A * A_inv;
    ASSERT_TRUE(product.isIdentity(TOLERANCE));

    // Test specific inverse values
    // A = [4 7], det(A) = 24-14 = 10
    //     [2 6]
    // A^(-1) = (1/10) * [6 -7] = [0.6 -0.7]
    //                    [-2 4]   [-0.2 0.4]
    ASSERT_NEAR(A_inv(0, 0), 0.6, TOLERANCE);
    ASSERT_NEAR(A_inv(0, 1), -0.7, TOLERANCE);
    ASSERT_NEAR(A_inv(1, 0), -0.2, TOLERANCE);
    ASSERT_NEAR(A_inv(1, 1), 0.4, TOLERANCE);

    // Test 3x3 inverse
    Matrix B({{1.0, 2.0, 3.0},
              {0.0, 1.0, 4.0},
              {5.0, 6.0, 0.0}});
    Matrix B_inv = B.inverse();
    Matrix B_product = B * B_inv;
    ASSERT_TRUE(B_product.isIdentity(TOLERANCE));

    std::cout << "Inverse: OK" << std::endl;
}

void test_matrix_rank() {
    std::cout << "\n=== Testing Matrix Rank ===" << std::endl;

    // Full rank matrix
    Matrix A({{1.0, 0.0}, {0.0, 1.0}});
    ASSERT_TRUE(A.rank(TOLERANCE) == 2);

    // Rank-deficient matrix
    Matrix B({{1.0, 2.0}, {2.0, 4.0}}); // Second row = 2 * first row
    ASSERT_TRUE(B.rank(TOLERANCE) == 1);

    // Zero matrix has rank 0
    Matrix Z(2, 2);
    ASSERT_TRUE(Z.rank(TOLERANCE) == 0);

    // 3x3 rank-2 matrix
    Matrix C({{1.0, 2.0, 3.0},
              {2.0, 4.0, 6.0},  // Row 2 = 2 * Row 1
              {0.0, 1.0, 1.0}});
    ASSERT_TRUE(C.rank(TOLERANCE) == 2);

    std::cout << "Rank: OK" << std::endl;
}

void test_matrix_solve() {
    std::cout << "\n=== Testing Matrix Solve ===" << std::endl;

    // Solve: [2 1] * [x] = [5]
    //        [1 3]   [y]   [7]
    Matrix A({{2.0, 1.0}, {1.0, 3.0}});
    Vector b({5.0, 7.0});
    Vector x = A.solve(b);

    // Solution: x = 1.6, y = 1.8 (verify: 2*1.6 + 1*1.8 = 5, 1*1.6 + 3*1.8 = 7)
    ASSERT_NEAR(x[0], 1.6, TOLERANCE);
    ASSERT_NEAR(x[1], 1.8, TOLERANCE);

    // Verify A*x = b
    Vector Ax = A * x;
    ASSERT_NEAR(Ax[0], b[0], TOLERANCE);
    ASSERT_NEAR(Ax[1], b[1], TOLERANCE);

    std::cout << "Solve: OK" << std::endl;
}

void test_matrix_properties() {
    std::cout << "\n=== Testing Matrix Properties ===" << std::endl;

    // Test symmetric matrix
    Matrix S({{1.0, 2.0, 3.0},
              {2.0, 4.0, 5.0},
              {3.0, 5.0, 6.0}});
    ASSERT_TRUE(S.isSymmetric(TOLERANCE));

    // Test non-symmetric
    Matrix N({{1.0, 2.0}, {3.0, 4.0}});
    ASSERT_FALSE(N.isSymmetric(TOLERANCE));

    // Test diagonal matrix
    Matrix D({{1.0, 0.0, 0.0},
              {0.0, 2.0, 0.0},
              {0.0, 0.0, 3.0}});
    ASSERT_TRUE(D.isDiagonal(TOLERANCE));

    // Test identity
    Matrix I = Matrix::identity(3);
    ASSERT_TRUE(I.isIdentity(TOLERANCE));
    ASSERT_TRUE(I.isDiagonal(TOLERANCE));
    ASSERT_TRUE(I.isSymmetric(TOLERANCE));

    std::cout << "Matrix properties: OK" << std::endl;
}

// ============================================================================
// UNITS TESTS
// ============================================================================

void test_force_conversions() {
    std::cout << "\n=== Testing Force Conversions ===" << std::endl;

    // Test Newtons to Dynes
    double force_N = 1.0;
    double force_dynes = newtonsToДynes(force_N);
    ASSERT_NEAR(force_dynes, 1.0e5, TOLERANCE);

    // Test round-trip conversion
    double back_to_N = dynesToNewtons(force_dynes);
    ASSERT_NEAR(back_to_N, force_N, TOLERANCE);

    // Test force calculation in Dynes
    double mass_g = 100.0;
    double accel_cm_s2 = 980.0; // ~1g
    double calculated_dynes = calculateForceDynes(mass_g, accel_cm_s2);
    ASSERT_NEAR(calculated_dynes, 98000.0, TOLERANCE);

    // Test SI to Dynes conversion
    double mass_kg = 1.0;
    double accel_m_s2 = 9.8;
    double force_dynes_from_SI = forceInDynesFromSI(mass_kg, accel_m_s2);
    ASSERT_NEAR(force_dynes_from_SI, 9.8e5, TOLERANCE);

    std::cout << "Force conversions: OK" << std::endl;
}

void test_mass_conversions() {
    std::cout << "\n=== Testing Mass Conversions ===" << std::endl;

    // Test kg to grams
    double mass_kg = 2.5;
    double mass_g = kilogramsToGrams(mass_kg);
    ASSERT_NEAR(mass_g, 2500.0, TOLERANCE);

    // Test round-trip
    double back_to_kg = gramsToKilograms(mass_g);
    ASSERT_NEAR(back_to_kg, mass_kg, TOLERANCE);

    std::cout << "Mass conversions: OK" << std::endl;
}

void test_length_conversions() {
    std::cout << "\n=== Testing Length Conversions ===" << std::endl;

    // Test meters to cm
    double length_m = 1.5;
    double length_cm = metersToCentimeters(length_m);
    ASSERT_NEAR(length_cm, 150.0, TOLERANCE);

    // Test round-trip
    double back_to_m = centimetersToMeters(length_cm);
    ASSERT_NEAR(back_to_m, length_m, TOLERANCE);

    // Test meters to feet
    double length_ft = metersToFeet(1.0);
    ASSERT_NEAR(length_ft, 3.28084, TOLERANCE);

    // Test feet to meters
    double back_to_m2 = feetToMeters(length_ft);
    ASSERT_NEAR(back_to_m2, 1.0, TOLERANCE);

    std::cout << "Length conversions: OK" << std::endl;
}

void test_velocity_conversions() {
    std::cout << "\n=== Testing Velocity Conversions ===" << std::endl;

    // Test m/s to km/h
    double vel_mps = 10.0;
    double vel_kmph = mpsToKmph(vel_mps);
    ASSERT_NEAR(vel_kmph, 36.0, TOLERANCE);

    // Test round-trip
    double back_to_mps = kmphToMps(vel_kmph);
    ASSERT_NEAR(back_to_mps, vel_mps, TOLERANCE);

    // Test m/s to mph
    double vel_mph = mpsToMph(10.0);
    ASSERT_NEAR(vel_mph, 22.3694, 1e-3);

    // Test mph to m/s
    double back_to_mps2 = mphToMps(vel_mph);
    ASSERT_NEAR(back_to_mps2, 10.0, 1e-3);

    std::cout << "Velocity conversions: OK" << std::endl;
}

void test_energy_conversions() {
    std::cout << "\n=== Testing Energy Conversions ===" << std::endl;

    // Test Joules to ergs
    double energy_J = 1.0;
    double energy_ergs = joulesToErgs(energy_J);
    ASSERT_NEAR(energy_ergs, 1.0e7, TOLERANCE);

    // Test round-trip
    double back_to_J = ergsToJoules(energy_ergs);
    ASSERT_NEAR(back_to_J, energy_J, TOLERANCE);

    // Test Joules to calories
    double energy_cal = joulesToCalories(4.184);
    ASSERT_NEAR(energy_cal, 1.0, 1e-3);

    // Test calories to Joules
    double back_to_J2 = caloriesToJoules(energy_cal);
    ASSERT_NEAR(back_to_J2, 4.184, 1e-3);

    std::cout << "Energy conversions: OK" << std::endl;
}

void test_angle_conversions() {
    std::cout << "\n=== Testing Angle Conversions ===" << std::endl;

    // Test degrees to radians
    double deg = 180.0;
    double rad = degreesToRadians(deg);
    ASSERT_NEAR(rad, M_PI, TOLERANCE);

    // Test round-trip
    double back_to_deg = radiansToDegrees(rad);
    ASSERT_NEAR(back_to_deg, deg, TOLERANCE);

    // Test common angles
    ASSERT_NEAR(degreesToRadians(90.0), M_PI/2.0, TOLERANCE);
    ASSERT_NEAR(degreesToRadians(360.0), 2.0*M_PI, TOLERANCE);
    ASSERT_NEAR(radiansToDegrees(M_PI/4.0), 45.0, TOLERANCE);

    std::cout << "Angle conversions: OK" << std::endl;
}

void test_acceleration_conversions() {
    std::cout << "\n=== Testing Acceleration Conversions ===" << std::endl;

    // Test acceleration in g's
    double accel_mps2 = 9.81;
    double accel_gs = accelerationInGs(accel_mps2);
    ASSERT_NEAR(accel_gs, 1.0, TOLERANCE);

    double accel_2g = accelerationInGs(19.62);
    ASSERT_NEAR(accel_2g, 2.0, 1e-2);

    std::cout << "Acceleration conversions: OK" << std::endl;
}

// ============================================================================
// MAIN TEST RUNNER
// ============================================================================

int main() {
    std::cout << "=====================================================" << std::endl;
    std::cout << "PHASE 1 VALIDATION: CORE UTILITIES" << std::endl;
    std::cout << "=====================================================" << std::endl;

    // Vector tests
    test_vector_basic_operations();
    test_vector_dot_product();
    test_vector_norm_and_normalize();
    test_vector_cross_product();
    test_vector_projections();
    test_vector_orthogonality();
    test_linear_independence();

    // Matrix tests
    test_matrix_basic_operations();
    test_matrix_multiplication();
    test_matrix_transpose();
    test_matrix_determinant();
    test_matrix_trace();
    test_matrix_inverse();
    test_matrix_rank();
    test_matrix_solve();
    test_matrix_properties();

    // Units tests
    test_force_conversions();
    test_mass_conversions();
    test_length_conversions();
    test_velocity_conversions();
    test_energy_conversions();
    test_angle_conversions();
    test_acceleration_conversions();

    // Summary
    std::cout << "\n=====================================================" << std::endl;
    std::cout << "VALIDATION SUMMARY" << std::endl;
    std::cout << "=====================================================" << std::endl;
    std::cout << "Tests Passed: " << tests_passed << std::endl;
    std::cout << "Tests Failed: " << tests_failed << std::endl;
    std::cout << "Total Tests:  " << (tests_passed + tests_failed) << std::endl;

    if (tests_failed == 0) {
        std::cout << "\n✓ ALL TESTS PASSED!" << std::endl;
        std::cout << "Phase 1 validation: COMPLETE" << std::endl;
        return 0;
    } else {
        std::cout << "\n✗ SOME TESTS FAILED" << std::endl;
        std::cout << "Phase 1 validation: INCOMPLETE" << std::endl;
        return 1;
    }
}
