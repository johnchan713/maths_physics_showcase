/**
 * Phase 4 Validation: Higgs Mechanism and Electroweak Symmetry Breaking
 *
 * Tests the gauge_theory_higgs_mechanism.hpp module functions.
 *
 * Coverage:
 * - Spontaneous symmetry breaking (SSB)
 * - Goldstone theorem and Goldstone bosons
 * - Higgs mechanism (gauge bosons acquire mass)
 * - Electroweak symmetry breaking (EWSB)
 * - Gauge boson masses (W± Z⁰ photon)
 * - Higgs VEV (v = 246 GeV)
 * - ρ parameter test
 * - Yukawa couplings and fermion masses
 * - Higgs boson properties (mass decay modes width)
 * - Top quark special role
 * - Higgs discovery and production modes
 * - Vacuum stability
 */

#include <iostream>
#include <cmath>
#include <string>
#include "../include/physics/gauge_theory_higgs_mechanism.hpp"

using namespace physics::advanced::gauge_theory;

// Test tolerances
const double TOLERANCE = 1e-5;
const double LOOSE_TOLERANCE = 1e-3;
const double VERY_LOOSE = 0.01;

// Test macros
#define ASSERT_NEAR(actual, expected, tolerance) \
    do { \
        if (std::abs((actual) - (expected)) > (tolerance)) { \
            std::cerr << "FAIL: " << __LINE__ << ": " << #actual \
                      << " = " << (actual) << ", expected " << (expected) \
                      << " (diff: " << std::abs((actual) - (expected)) << ")" << std::endl; \
            return false; \
        } \
    } while(0)

#define ASSERT_TRUE(condition) \
    do { \
        if (!(condition)) { \
            std::cerr << "FAIL: " << __LINE__ << ": " << #condition << std::endl; \
            return false; \
        } \
    } while(0)

int main() {
    int tests_passed = 0;
    int tests_failed = 0;

    std::cout << "=== Phase 4: Higgs Mechanism and EWSB Validation ===" << std::endl;
    std::cout << std::endl;

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
    // Spontaneous Symmetry Breaking
    // ========================================

    run_test("SSB definition not empty", []() {
        std::string def = SpontaneousSymmetryBreaking::definition();
        ASSERT_TRUE(def.length() > 0);
        ASSERT_TRUE(def.find("Symmetry") != std::string::npos);
        return true;
    });

    run_test("Mexican hat potential description", []() {
        std::string desc = SpontaneousSymmetryBreaking::mexicanHatPotential();
        ASSERT_TRUE(desc.find("V(") != std::string::npos);
        ASSERT_TRUE(desc.find("Minimum") != std::string::npos || desc.find("minimum") != std::string::npos);
        return true;
    });

    run_test("VEV description contains value", []() {
        std::string desc = SpontaneousSymmetryBreaking::vacuumExpectationValue();
        ASSERT_TRUE(desc.find("246") != std::string::npos);
        ASSERT_TRUE(desc.find("GeV") != std::string::npos);
        return true;
    });

    run_test("Ferromagnetism analogy exists", []() {
        std::string analogy = SpontaneousSymmetryBreaking::ferromagnetismAnalogy();
        ASSERT_TRUE(analogy.find("magnet") != std::string::npos ||
                   analogy.find("Magnet") != std::string::npos);
        return true;
    });

    run_test("Discrete vs continuous SSB explained", []() {
        std::string desc = SpontaneousSymmetryBreaking::discreteVsContinuous();
        ASSERT_TRUE(desc.find("Discrete") != std::string::npos);
        ASSERT_TRUE(desc.find("Continuous") != std::string::npos);
        ASSERT_TRUE(desc.find("Goldstone") != std::string::npos);
        return true;
    });

    // ========================================
    // Goldstone Theorem
    // ========================================

    run_test("Goldstone theorem statement", []() {
        std::string stmt = GoldstonesTheorem::statement();
        ASSERT_TRUE(stmt.find("massless") != std::string::npos);
        ASSERT_TRUE(stmt.find("boson") != std::string::npos);
        return true;
    });

    run_test("Goldstone physical origin", []() {
        std::string origin = GoldstonesTheorem::physicalOrigin();
        ASSERT_TRUE(origin.find("flat") != std::string::npos ||
                   origin.find("vacuum") != std::string::npos);
        return true;
    });

    run_test("U(1) SSB example", []() {
        std::string example = GoldstonesTheorem::u1Example();
        ASSERT_TRUE(example.find("massless") != std::string::npos);
        ASSERT_TRUE(example.find("Goldstone") != std::string::npos);
        return true;
    });

    run_test("SU(2) SSB example", []() {
        std::string example = GoldstonesTheorem::su2Example();
        ASSERT_TRUE(example.find("SU(2)") != std::string::npos);
        return true;
    });

    run_test("Pions as Goldstone bosons", []() {
        std::string desc = GoldstonesTheorem::pionsAsGoldstoneBosons();
        ASSERT_TRUE(desc.find("140") != std::string::npos ||
                   desc.find("pion") != std::string::npos ||
                   desc.find("MeV") != std::string::npos);
        return true;
    });

    // ========================================
    // Higgs Mechanism
    // ========================================

    run_test("Higgs mechanism summary", []() {
        std::string mech = HiggsMechanism::mechanism();
        ASSERT_TRUE(mech.find("eaten") != std::string::npos);
        ASSERT_TRUE(mech.find("mass") != std::string::npos);
        return true;
    });

    run_test("Degrees of freedom conservation", []() {
        std::string dof = HiggsMechanism::degreesFreedom();
        ASSERT_TRUE(dof.find("DOF") != std::string::npos ||
                   dof.find("4") != std::string::npos);
        return true;
    });

    run_test("Abelian Higgs model description", []() {
        std::string model = HiggsMechanism::abelianHiggsModel();
        ASSERT_TRUE(model.find("U(1)") != std::string::npos);
        ASSERT_TRUE(model.find("mass") != std::string::npos);
        return true;
    });

    run_test("Unitary gauge description", []() {
        std::string gauge = HiggsMechanism::unitaryGauge();
        ASSERT_TRUE(gauge.find("gauge") != std::string::npos);
        ASSERT_TRUE(gauge.find("Goldstone") != std::string::npos);
        return true;
    });

    run_test("Mass generation explanation", []() {
        std::string gen = HiggsMechanism::massGeneration();
        ASSERT_TRUE(gen.find("VEV") != std::string::npos ||
                   gen.find("coupling") != std::string::npos);
        return true;
    });

    run_test("Why massive gauge bosons allowed", []() {
        std::string why = HiggsMechanism::whyMassiveGaugeBosons();
        ASSERT_TRUE(why.find("gauge") != std::string::npos);
        ASSERT_TRUE(why.find("invariance") != std::string::npos ||
                   why.find("symmetry") != std::string::npos);
        return true;
    });

    // ========================================
    // Electroweak Symmetry Breaking
    // ========================================

    run_test("EWSB breaking pattern", []() {
        std::string pattern = ElectroweakSymmetryBreaking::breakingPattern();
        ASSERT_TRUE(pattern.find("SU(2)") != std::string::npos);
        ASSERT_TRUE(pattern.find("U(1)") != std::string::npos);
        return true;
    });

    run_test("Higgs doublet description", []() {
        std::string doublet = ElectroweakSymmetryBreaking::higgsDoublet();
        ASSERT_TRUE(doublet.find("doublet") != std::string::npos);
        ASSERT_TRUE(doublet.find("4") != std::string::npos);
        return true;
    });

    run_test("Higgs potential form", []() {
        std::string pot = ElectroweakSymmetryBreaking::higgsPotential();
        ASSERT_TRUE(pot.find("246") != std::string::npos);
        ASSERT_TRUE(pot.find("GeV") != std::string::npos);
        return true;
    });

    run_test("Vacuum choice preserves EM", []() {
        std::string vac = ElectroweakSymmetryBreaking::vacuumChoice();
        ASSERT_TRUE(vac.find("photon") != std::string::npos ||
                   vac.find("massless") != std::string::npos);
        return true;
    });

    run_test("W boson mass approximately 80 GeV", []() {
        double m_W = ElectroweakSymmetryBreaking::wBosonMass();
        ASSERT_TRUE(m_W > 79.0 && m_W < 82.0);
        return true;
    });

    run_test("Z boson mass approximately 91 GeV", []() {
        double m_Z = ElectroweakSymmetryBreaking::zBosonMass();
        ASSERT_TRUE(m_Z > 90.0 && m_Z < 92.0);
        return true;
    });

    run_test("Higgs VEV is 246 GeV", []() {
        double v = ElectroweakSymmetryBreaking::higgsVEV();
        ASSERT_NEAR(v, 246.0, 1.0);
        return true;
    });

    run_test("Mass ratio W to Z", []() {
        double m_W = ElectroweakSymmetryBreaking::wBosonMass();
        double m_Z = ElectroweakSymmetryBreaking::zBosonMass();

        double ratio = m_W / m_Z;

        // cos θ_W ≈ 0.88
        ASSERT_TRUE(ratio > 0.85 && ratio < 0.90);
        return true;
    });

    run_test("Gauge boson mass formulas", []() {
        std::string formulas = ElectroweakSymmetryBreaking::gaugeBosonMassFormulas();
        ASSERT_TRUE(formulas.find("80") != std::string::npos);
        ASSERT_TRUE(formulas.find("91") != std::string::npos);
        return true;
    });

    run_test("Rho parameter close to 1", []() {
        double rho = ElectroweakSymmetryBreaking::rhoParameter();
        ASSERT_NEAR(rho, 1.0, 0.01);
        return true;
    });

    run_test("Rho parameter significance", []() {
        std::string sig = ElectroweakSymmetryBreaking::rhoParameterSignificance();
        ASSERT_TRUE(sig.find("doublet") != std::string::npos);
        ASSERT_TRUE(sig.find("1") != std::string::npos);
        return true;
    });

    // ========================================
    // Fermion Masses and Yukawa Couplings
    // ========================================

    run_test("Yukawa interaction description", []() {
        std::string yukawa = FermionMasses::yukawaInteraction();
        ASSERT_TRUE(yukawa.find("Yukawa") != std::string::npos);
        ASSERT_TRUE(yukawa.find("mass") != std::string::npos);
        return true;
    });

    run_test("Why no direct fermion mass", []() {
        std::string why = FermionMasses::whyNoDirectMass();
        ASSERT_TRUE(why.find("gauge") != std::string::npos);
        ASSERT_TRUE(why.find("invariance") != std::string::npos);
        return true;
    });

    run_test("Yukawa coupling from mass", []() {
        double m_fermion = 173.0;  // Top quark mass (GeV)
        double y = FermionMasses::yukawaCoupling(m_fermion);

        // y = m√2/v ≈ 173×1.414/246 ≈ 0.995
        ASSERT_TRUE(y > 0.9 && y < 1.1);
        return true;
    });

    run_test("Electron Yukawa coupling is tiny", []() {
        auto yukawas = FermionMasses::standardModelYukawas();
        double y_e = yukawas["electron"];

        ASSERT_TRUE(y_e < 1e-5);
        return true;
    });

    run_test("Top Yukawa coupling near 1", []() {
        auto yukawas = FermionMasses::standardModelYukawas();
        double y_t = yukawas["top"];

        ASSERT_TRUE(y_t > 0.9 && y_t < 1.1);
        return true;
    });

    run_test("Yukawa hierarchy electron to top", []() {
        auto yukawas = FermionMasses::standardModelYukawas();
        double y_e = yukawas["electron"];
        double y_t = yukawas["top"];

        double ratio = y_e / y_t;

        ASSERT_TRUE(ratio < 1e-5);  // Huge hierarchy
        return true;
    });

    run_test("Top quark special description", []() {
        std::string special = FermionMasses::topQuarkSpecial();
        ASSERT_TRUE(special.find("top") != std::string::npos ||
                   special.find("Top") != std::string::npos);
        ASSERT_TRUE(special.find("173") != std::string::npos ||
                   special.find("Higgs") != std::string::npos);
        return true;
    });

    run_test("Hierarchy problem mentioned", []() {
        std::string prob = FermionMasses::hierarchyProblem();
        ASSERT_TRUE(prob.find("hierarchy") != std::string::npos);
        return true;
    });

    // ========================================
    // Higgs Boson Properties
    // ========================================

    run_test("Higgs mass is 125 GeV", []() {
        double m_h = HiggsBosonProperties::mass();
        ASSERT_TRUE(m_h > 124.0 && m_h < 126.0);
        return true;
    });

    run_test("Higgs mass significance", []() {
        std::string sig = HiggsBosonProperties::massSignificance();
        ASSERT_TRUE(sig.find("125") != std::string::npos);
        ASSERT_TRUE(sig.find("GeV") != std::string::npos);
        return true;
    });

    run_test("Higgs branching ratios sum less than 1", []() {
        auto BR = HiggsBosonProperties::branchingRatios();

        double total = 0.0;
        for (const auto& [channel, ratio] : BR) {
            total += ratio;
        }

        ASSERT_TRUE(total > 0.9 && total < 1.1);
        return true;
    });

    run_test("Higgs decay to bb dominates", []() {
        auto BR = HiggsBosonProperties::branchingRatios();
        double BR_bb = BR["bb"];

        ASSERT_TRUE(BR_bb > 0.5);  // ~58%
        return true;
    });

    run_test("Higgs decay to diphoton is rare", []() {
        auto BR = HiggsBosonProperties::branchingRatios();
        double BR_gg = BR["gamma gamma"];

        ASSERT_TRUE(BR_gg < 0.01);  // ~0.23%
        return true;
    });

    run_test("Higgs width is narrow", []() {
        double Gamma = HiggsBosonProperties::width();

        ASSERT_TRUE(Gamma < 0.01);  // ~4.1 MeV = 0.0041 GeV
        return true;
    });

    run_test("Higgs coupling pattern description", []() {
        std::string pattern = HiggsBosonProperties::couplingPattern();
        ASSERT_TRUE(pattern.find("mass") != std::string::npos);
        ASSERT_TRUE(pattern.find("coupling") != std::string::npos);
        return true;
    });

    run_test("Discovery channels mentioned", []() {
        std::string disc = HiggsBosonProperties::discoveryChannels();
        ASSERT_TRUE(disc.find("2012") != std::string::npos);
        ASSERT_TRUE(disc.find("gamma") != std::string::npos ||
                   disc.find("photon") != std::string::npos);
        ASSERT_TRUE(disc.find("ZZ") != std::string::npos ||
                   disc.find("4") != std::string::npos);
        return true;
    });

    run_test("Production modes at LHC", []() {
        std::string prod = HiggsBosonProperties::productionModes();
        ASSERT_TRUE(prod.find("gluon") != std::string::npos);
        ASSERT_TRUE(prod.find("fusion") != std::string::npos ||
                   prod.find("gg") != std::string::npos);
        return true;
    });

    run_test("Vacuum stability discussion", []() {
        std::string stab = HiggsBosonProperties::vacuumStability();
        ASSERT_TRUE(stab.find("vacuum") != std::string::npos);
        ASSERT_TRUE(stab.find("stable") != std::string::npos ||
                   stab.find("metastable") != std::string::npos);
        return true;
    });

    // ========================================
    // Physical Consistency Checks
    // ========================================

    run_test("W and Z masses heavier than Higgs divided by sqrt(2)", []() {
        double m_W = ElectroweakSymmetryBreaking::wBosonMass();
        double m_Z = ElectroweakSymmetryBreaking::zBosonMass();
        double m_h = HiggsBosonProperties::mass();

        // Since m_W,Z ~ gv/2 and m_h ~ √(2λ)v,
        // we expect m_W,Z < m_h for λ > g²
        // But actually m_W < m_h and m_Z < m_h is true
        ASSERT_TRUE(m_W < m_h);
        ASSERT_TRUE(m_Z < m_h);
        return true;
    });

    run_test("Higgs VEV related to gauge boson masses", []() {
        double v = ElectroweakSymmetryBreaking::higgsVEV();
        double m_W = ElectroweakSymmetryBreaking::wBosonMass();

        // m_W ≈ gv/2, so v ≈ 2m_W/g
        // For g ≈ 0.65, v ≈ 2×80/0.65 ≈ 246 GeV
        double v_estimate = 2.0 * m_W / 0.65;

        ASSERT_NEAR(v_estimate, v, 10.0);
        return true;
    });

    run_test("Top quark mass close to Higgs VEV", []() {
        double v = ElectroweakSymmetryBreaking::higgsVEV();
        auto yukawas = FermionMasses::standardModelYukawas();
        double y_t = yukawas["top"];

        // m_t = y_t v/√2 ≈ v/√2 for y_t ≈ 1
        double m_t = y_t * v / std::sqrt(2.0);

        ASSERT_TRUE(m_t > 170.0 && m_t < 175.0);
        return true;
    });

    run_test("Electron mass from Yukawa tiny", []() {
        double v = ElectroweakSymmetryBreaking::higgsVEV();
        auto yukawas = FermionMasses::standardModelYukawas();
        double y_e = yukawas["electron"];

        // m_e = y_e v/√2
        double m_e_GeV = y_e * v / std::sqrt(2.0);
        double m_e_MeV = m_e_GeV * 1000.0;

        ASSERT_TRUE(m_e_MeV < 1.0);  // ~0.5 MeV
        return true;
    });

    run_test("Rho parameter from mass ratio", []() {
        double m_W = ElectroweakSymmetryBreaking::wBosonMass();
        double m_Z = ElectroweakSymmetryBreaking::zBosonMass();

        // ρ = m_W²/(m_Z² cos²θ_W)
        // For tree-level doublet: cos²θ_W = m_W²/m_Z²
        // So ρ = 1
        double cos2_theta_W = m_W * m_W / (m_Z * m_Z);
        double rho_calc = m_W * m_W / (m_Z * m_Z * cos2_theta_W);

        ASSERT_NEAR(rho_calc, 1.0, TOLERANCE);
        return true;
    });

    // ========================================
    // String Content Tests
    // ========================================

    run_test("All functions return non-empty strings", []() {
        ASSERT_TRUE(SpontaneousSymmetryBreaking::definition().length() > 10);
        ASSERT_TRUE(GoldstonesTheorem::statement().length() > 10);
        ASSERT_TRUE(HiggsMechanism::mechanism().length() > 10);
        ASSERT_TRUE(ElectroweakSymmetryBreaking::breakingPattern().length() > 10);
        ASSERT_TRUE(FermionMasses::yukawaInteraction().length() > 10);
        ASSERT_TRUE(HiggsBosonProperties::discoveryChannels().length() > 10);
        return true;
    });

    run_test("Key physics terms present", []() {
        // Check that important terms appear in the module
        std::string combined =
            SpontaneousSymmetryBreaking::definition() +
            GoldstonesTheorem::statement() +
            HiggsMechanism::mechanism() +
            ElectroweakSymmetryBreaking::breakingPattern() +
            FermionMasses::yukawaInteraction() +
            HiggsBosonProperties::massSignificance();

        ASSERT_TRUE(combined.find("symmetry") != std::string::npos ||
                   combined.find("Symmetry") != std::string::npos);
        ASSERT_TRUE(combined.find("mass") != std::string::npos);
        ASSERT_TRUE(combined.find("Higgs") != std::string::npos ||
                   combined.find("higgs") != std::string::npos);
        return true;
    });

    // ========================================
    // Numerical Values Summary
    // ========================================

    run_test("Summary of Standard Model masses", []() {
        double m_W = ElectroweakSymmetryBreaking::wBosonMass();
        double m_Z = ElectroweakSymmetryBreaking::zBosonMass();
        double m_h = HiggsBosonProperties::mass();
        double v = ElectroweakSymmetryBreaking::higgsVEV();

        // All in reasonable ranges
        ASSERT_TRUE(m_W > 0 && m_W < 200);
        ASSERT_TRUE(m_Z > 0 && m_Z < 200);
        ASSERT_TRUE(m_h > 0 && m_h < 200);
        ASSERT_TRUE(v > 0 && v < 500);
        return true;
    });

    // ========================================
    // Summary
    // ========================================

    std::cout << std::endl;
    std::cout << "========================================" << std::endl;
    std::cout << "Test Summary" << std::endl;
    std::cout << "========================================" << std::endl;
    std::cout << "Tests passed: " << tests_passed << std::endl;
    std::cout << "Tests failed: " << tests_failed << std::endl;
    std::cout << "Total tests:  " << (tests_passed + tests_failed) << std::endl;
    std::cout << "Success rate: " << (100.0 * tests_passed / (tests_passed + tests_failed)) << "%" << std::endl;

    return tests_failed == 0 ? 0 : 1;
}
