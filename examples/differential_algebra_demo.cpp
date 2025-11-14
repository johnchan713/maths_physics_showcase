/**
 * @file differential_algebra_demo.cpp
 * @brief Comprehensive demonstration of differential algebra computational implementations
 *
 * Demonstrates practical usage of:
 * - Differential polynomials and fields
 * - Characteristic sets and Ritt's reduction algorithm
 * - Differential ideals and membership testing
 * - Algebraic differential manifolds
 * - Resolvents and dimension theory
 * - Elimination methods
 * - Low power theorem for singular solutions
 */

#include "maths/algebra/differential_algebra.hpp"
#include <iostream>
#include <iomanip>

using namespace maths::algebra;

void printSeparator(const std::string& title) {
    std::cout << "\n=== " << title << " ===" << std::endl;
}

// Demo 1: Basic Differential Polynomials
void demoDifferentialPolynomials() {
    printSeparator("Demo 1: Differential Polynomials");

    // Create y_0^(0) (i.e., y)
    auto y = DifferentialPolynomial::derivative(0, 0);
    std::cout << "Created differential polynomial y (variable 0, order 0)" << std::endl;

    // Create y_0^(1) (i.e., y')
    auto y_prime = DifferentialPolynomial::derivative(0, 1);
    std::cout << "Created differential polynomial y' (variable 0, order 1)" << std::endl;

    // Create y_0^(2) (i.e., y'')
    auto y_double_prime = DifferentialPolynomial::derivative(0, 2);
    std::cout << "Created differential polynomial y'' (variable 0, order 2)" << std::endl;

    // Create polynomial: y'' + y = 0 (harmonic oscillator)
    auto constant_one = DifferentialPolynomial::constant(1.0);
    auto harmonic = y_double_prime + y;
    std::cout << "\nHarmonic oscillator equation: y'' + y" << std::endl;
    std::cout << "  Order: " << harmonic.order() << std::endl;
    std::cout << "  Number of terms: " << harmonic.terms.size() << std::endl;

    // Create polynomial: (y')^2 - 2y (related to y = x^2)
    auto y_prime_squared = y_prime * y_prime;
    auto two_y = y * DifferentialPolynomial::constant(2.0);
    auto clairaut_like = y_prime_squared + two_y * DifferentialPolynomial::constant(-1.0);
    std::cout << "\nDifferential polynomial: (y')^2 - 2y" << std::endl;
    std::cout << "  Order: " << clairaut_like.order() << std::endl;
    std::cout << "  Degree: " << clairaut_like.degree() << std::endl;
}

// Demo 2: Differentiation in Differential Fields
void demoDifferentialFields() {
    printSeparator("Demo 2: Differential Field Operations");

    // Create polynomial y^2
    auto y = DifferentialPolynomial::derivative(0, 0);
    auto y_squared = y * y;

    std::cout << "Starting polynomial: y^2" << std::endl;
    std::cout << "  Order: " << y_squared.order() << std::endl;

    // Differentiate using Leibniz rule: D(y^2) = 2y·y'
    auto d_y_squared = DifferentialField::differentiate(y_squared);
    std::cout << "\nDifferentiated D(y^2) = 2y·y'" << std::endl;
    std::cout << "  Order: " << d_y_squared.order() << std::endl;
    std::cout << "  Number of terms: " << d_y_squared.terms.size() << std::endl;

    // Differentiate again: D^2(y^2) = 2(y')^2 + 2y·y''
    auto d2_y_squared = DifferentialField::differentiate(d_y_squared);
    std::cout << "\nSecond derivative D²(y^2) = 2(y')^2 + 2y·y''" << std::endl;
    std::cout << "  Order: " << d2_y_squared.order() << std::endl;

    // Check if constant
    auto constant = DifferentialPolynomial::constant(5.0);
    bool is_const = DifferentialField::isConstant(constant);
    std::cout << "\nIs constant(5.0) actually constant? " << (is_const ? "Yes" : "No") << std::endl;

    bool y_is_const = DifferentialField::isConstant(y);
    std::cout << "Is y constant? " << (y_is_const ? "Yes" : "No") << std::endl;
}

// Demo 3: Rankings and Leaders
void demoRankingsAndLeaders() {
    printSeparator("Demo 3: Rankings and Leaders");

    // Create differential polynomial y·y' + y''
    auto y = DifferentialPolynomial::derivative(0, 0);
    auto y_prime = DifferentialPolynomial::derivative(0, 1);
    auto y_double_prime = DifferentialPolynomial::derivative(0, 2);

    auto poly = y * y_prime + y_double_prime;

    std::cout << "Polynomial: y·y' + y''" << std::endl;

    // Test orderly ranking
    Ranking orderly(Ranking::ORDERLY);
    auto leader_orderly = CharacteristicSet::leader(poly, orderly);
    std::cout << "\nOrderly ranking - Leader: y_" << leader_orderly.variable_index
              << "^(" << leader_orderly.order << ")" << std::endl;
    std::cout << "(Orderly: compare by order first, then variable)" << std::endl;

    // Test eliminative ranking
    Ranking eliminative(Ranking::ELIMINATIVE);
    auto leader_elim = CharacteristicSet::leader(poly, eliminative);
    std::cout << "\nEliminative ranking - Leader: y_" << leader_elim.variable_index
              << "^(" << leader_elim.order << ")" << std::endl;
    std::cout << "(Eliminative: compare by variable first, then order)" << std::endl;

    // Test initial (leading coefficient)
    auto init = CharacteristicSet::initial(poly, orderly);
    std::cout << "\nInitial (leading coefficient w.r.t. orderly ranking):" << std::endl;
    std::cout << "  Is zero: " << init.isZero() << std::endl;
    std::cout << "  Number of terms: " << init.terms.size() << std::endl;
}

// Demo 4: Characteristic Sets and Reduction
void demoCharacteristicSets() {
    printSeparator("Demo 4: Characteristic Sets and Reduction");

    // Create system: {y'' + y = 0}
    auto y = DifferentialPolynomial::derivative(0, 0);
    auto y_double_prime = DifferentialPolynomial::derivative(0, 2);
    auto harmonic = y_double_prime + y;

    std::cout << "System: y'' + y = 0 (harmonic oscillator)" << std::endl;

    // Create characteristic set
    CharacteristicSet char_set(Ranking::ORDERLY);
    char_set.addPolynomial(harmonic);

    std::cout << "\nCharacteristic set size: " << char_set.polynomials.size() << std::endl;

    // Test reduction
    auto y_prime = DifferentialPolynomial::derivative(0, 1);
    auto test_poly = y_double_prime + y * DifferentialPolynomial::constant(2.0);

    std::cout << "\nTest polynomial: y'' + 2y" << std::endl;
    auto reduced = char_set.reduce(test_poly);
    std::cout << "After reduction: " << (reduced.isZero() ? "Zero" : "Non-zero") << std::endl;
    std::cout << "Reduced polynomial terms: " << reduced.terms.size() << std::endl;

    // Test ideal membership
    bool in_ideal = char_set.inIdeal(harmonic);
    std::cout << "\nIs y'' + y in the ideal? " << (in_ideal ? "Yes" : "No") << std::endl;
}

// Demo 5: Differential Ideals
void demoDifferentialIdeals() {
    printSeparator("Demo 5: Differential Ideals");

    // Create ideal generated by y' - y = 0 (exponential function)
    auto y = DifferentialPolynomial::derivative(0, 0);
    auto y_prime = DifferentialPolynomial::derivative(0, 1);
    auto exponential_eq = y_prime + y * DifferentialPolynomial::constant(-1.0);

    std::cout << "Generator: y' - y = 0 (exponential growth)" << std::endl;

    std::vector<DifferentialPolynomial> generators = {exponential_eq};
    auto ideal = DifferentialIdeal::generate(generators);

    std::cout << "Created differential ideal" << std::endl;
    std::cout << "Number of generators: " << ideal.generators.size() << std::endl;

    // Compute characteristic set
    auto char_set = ideal.characteristicSet(Ranking::ORDERLY, 3);
    std::cout << "\nCharacteristic set size: " << char_set.polynomials.size() << std::endl;

    // Test membership
    auto y_double_prime = DifferentialPolynomial::derivative(0, 2);
    auto test = y_double_prime + y * DifferentialPolynomial::constant(-1.0);
    bool contains = ideal.contains(test, 3);
    std::cout << "\nDoes ideal contain y'' - y? " << (contains ? "Yes" : "No") << std::endl;
    std::cout << "(Should be Yes since y'' = (y')' = y' = y for exponential)" << std::endl;
}

// Demo 6: Algebraic Differential Manifolds
void demoManifolds() {
    printSeparator("Demo 6: Algebraic Differential Manifolds");

    // System: y'' + y = 0
    auto y = DifferentialPolynomial::derivative(0, 0);
    auto y_double_prime = DifferentialPolynomial::derivative(0, 2);
    auto harmonic = y_double_prime + y;

    std::vector<DifferentialPolynomial> system = {harmonic};
    auto manifold = AlgebraicDifferentialManifold::fromSystem(system);

    std::cout << "System: y'' + y = 0" << std::endl;
    std::cout << "Defining system size: " << manifold.defining_system.size() << std::endl;

    // Compute dimension
    int dim = manifold.dimension(3);
    std::cout << "\nDimension of solution manifold: " << dim << std::endl;
    std::cout << "(Number of arbitrary constants in general solution)" << std::endl;
    std::cout << "General solution: y = C₁·cos(x) + C₂·sin(x) has 2 constants" << std::endl;

    // Check irreducibility
    bool irreducible = manifold.isIrreducible();
    std::cout << "\nIs irreducible? " << (irreducible ? "Yes" : "No") << std::endl;

    // Decompose
    auto components = manifold.decompose();
    std::cout << "Number of irreducible components: " << components.size() << std::endl;
}

// Demo 7: Resolvents
void demoResolvents() {
    printSeparator("Demo 7: Resolvents");

    // Create characteristic set for y'' + y = 0
    auto y = DifferentialPolynomial::derivative(0, 0);
    auto y_double_prime = DifferentialPolynomial::derivative(0, 2);
    auto harmonic = y_double_prime + y;

    CharacteristicSet char_set(Ranking::ORDERLY);
    char_set.addPolynomial(harmonic);

    std::cout << "System: y'' + y = 0" << std::endl;

    // Construct resolvent
    Derivative target(0, 0);  // Eliminate to get equation in y only
    auto resolvent = Resolvent::construct(char_set, target);

    std::cout << "\nResolvent constructed" << std::endl;
    std::cout << "Order: " << resolvent.order() << std::endl;

    // Evaluate resolvent
    double test_val = 1.0;
    double result = resolvent.evaluate(test_val);
    std::cout << "\nResolvent evaluated at x=1.0: " << result << std::endl;

    // Find roots
    auto roots = resolvent.roots();
    std::cout << "\nRoots of resolvent: ";
    for (const auto& root : roots) {
        std::cout << std::fixed << std::setprecision(4) << root << " ";
    }
    std::cout << std::endl;
}

// Demo 8: Elimination Theory
void demoEliminationTheory() {
    printSeparator("Demo 8: Elimination Theory");

    // System with two variables: y₀' - y₁ = 0, y₁' - y₀ = 0
    auto y0 = DifferentialPolynomial::derivative(0, 0);
    auto y0_prime = DifferentialPolynomial::derivative(0, 1);
    auto y1 = DifferentialPolynomial::derivative(1, 0);
    auto y1_prime = DifferentialPolynomial::derivative(1, 1);

    auto eq1 = y0_prime + y1 * DifferentialPolynomial::constant(-1.0);
    auto eq2 = y1_prime + y0 * DifferentialPolynomial::constant(-1.0);

    std::vector<DifferentialPolynomial> system = {eq1, eq2};

    std::cout << "System:" << std::endl;
    std::cout << "  y₀' - y₁ = 0" << std::endl;
    std::cout << "  y₁' - y₀ = 0" << std::endl;

    // Test consistency
    bool consistent = EliminationTheory::isConsistent(system);
    std::cout << "\nSystem is consistent? " << (consistent ? "Yes" : "No") << std::endl;

    // Compute solution dimension
    int sol_dim = EliminationTheory::solutionDimension(system);
    std::cout << "Solution space dimension: " << sol_dim << std::endl;

    // Eliminate variable y₁
    auto eliminated = EliminationTheory::eliminate(system, 1);
    std::cout << "\nAfter eliminating y₁: " << eliminated.size()
              << " equations remain" << std::endl;
    std::cout << "(Should give y₀'' - y₀ = 0)" << std::endl;
}

// Demo 9: Low Power Theorem and Singular Solutions
void demoLowPowerTheorem() {
    printSeparator("Demo 9: Low Power Theorem - Singular Solutions");

    // Clairaut equation: y = xy' + f(y')
    // Simplified form: y - xy' - (y')² = 0
    auto y = DifferentialPolynomial::derivative(0, 0);
    auto y_prime = DifferentialPolynomial::derivative(0, 1);

    // Create y - (y')² (simplified version)
    auto y_prime_squared = y_prime * y_prime;
    auto clairaut = y + y_prime_squared * DifferentialPolynomial::constant(-1.0);

    std::cout << "Differential polynomial: y - (y')²" << std::endl;
    std::cout << "Order: " << clairaut.order() << std::endl;
    std::cout << "Degree: " << clairaut.degree() << std::endl;

    // Extract low power terms
    Ranking rank(Ranking::ORDERLY);
    auto low_terms = LowPowerTheorem::lowPowerTerms(clairaut, rank);
    std::cout << "\nLow power terms extracted" << std::endl;
    std::cout << "  Number of terms: " << low_terms.terms.size() << std::endl;

    // Check for singular solutions
    bool has_singular = LowPowerTheorem::hasSingularSolutions(clairaut);
    std::cout << "\nHas singular solutions? " << (has_singular ? "Yes" : "No") << std::endl;
    std::cout << "(Clairaut equations typically have singular solutions)" << std::endl;
}

// Demo 10: Practical Example - Pendulum Equation
void demoPendulumEquation() {
    printSeparator("Demo 10: Application - Pendulum Equation");

    // Linearized pendulum: θ'' + (g/L)θ = 0
    // Using y for θ, and coefficient k = g/L = 1
    auto y = DifferentialPolynomial::derivative(0, 0);
    auto y_double_prime = DifferentialPolynomial::derivative(0, 2);
    auto pendulum = y_double_prime + y;

    std::cout << "Linearized pendulum equation: θ'' + (g/L)θ = 0" << std::endl;
    std::cout << "(Using θ = y and g/L = 1)" << std::endl;

    // Create differential ideal
    std::vector<DifferentialPolynomial> system = {pendulum};
    auto ideal = DifferentialIdeal::generate(system);

    std::cout << "\nAnalysis:" << std::endl;

    // Characteristic set
    auto char_set = ideal.characteristicSet(Ranking::ORDERLY, 4);
    std::cout << "  Characteristic set size: " << char_set.polynomials.size() << std::endl;

    // Manifold properties
    auto manifold = AlgebraicDifferentialManifold::fromSystem(system);
    int dim = manifold.dimension(4);
    std::cout << "  Solution dimension: " << dim << " (arbitrary constants)" << std::endl;

    // Check if y'' - y is in the ideal
    auto test = y_double_prime + y * DifferentialPolynomial::constant(-1.0);
    bool related = ideal.contains(test, 4);
    std::cout << "  Is θ'' - θ related to original? " << (related ? "Yes" : "No") << std::endl;

    std::cout << "\nGeneral solution: θ(t) = A·cos(√(g/L)·t) + B·sin(√(g/L)·t)" << std::endl;
    std::cout << "Physical meaning: Simple harmonic motion with period T = 2π√(L/g)" << std::endl;
}

int main() {
    std::cout << "========================================" << std::endl;
    std::cout << "Differential Algebra Computational Demo" << std::endl;
    std::cout << "Ritt-Kolchin Theory Implementation" << std::endl;
    std::cout << "========================================" << std::endl;

    try {
        demoDifferentialPolynomials();
        demoDifferentialFields();
        demoRankingsAndLeaders();
        demoCharacteristicSets();
        demoDifferentialIdeals();
        demoManifolds();
        demoResolvents();
        demoEliminationTheory();
        demoLowPowerTheorem();
        demoPendulumEquation();

        std::cout << "\n========================================" << std::endl;
        std::cout << "All demos completed successfully!" << std::endl;
        std::cout << "========================================" << std::endl;

    } catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << std::endl;
        return 1;
    }

    return 0;
}
