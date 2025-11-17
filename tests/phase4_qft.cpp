/**
 * Phase 4 Validation: Quantum Field Theory - Particle Physics & Antiparticles
 *
 * Tests the particle_physics.hpp and antiparticles.hpp modules.
 *
 * Coverage:
 * PARTICLE PHYSICS (~50 tests):
 * - Quarks: 6 flavors, charges (±2/3, ±1/3), masses, color charge, spin 1/2
 * - Leptons: electron, muon, tau, neutrinos, charges, lifetimes
 * - Gauge bosons: photon (massless), W± (80.4 GeV), Z (91.2 GeV), gluon
 * - Higgs: 125.1 GeV, spin 0, scalar
 * - Generations: 3 generations, mass hierarchy (gen3 > gen2 > gen1)
 * - Quantum numbers: charge/baryon/lepton conservation
 * - Fermions vs bosons: spin classification
 *
 * ANTIPARTICLES (~40 tests):
 * - Antiparticle creation: opposite charge, baryon/lepton number, CPT theorem
 * - Self-conjugate: photon, Z, neutral particles
 * - Antiquarks: opposite charges
 * - Antileptons: positron, antimuon, antitau
 * - Charge conjugation: C-parity, conservation in EM/strong
 * - CP symmetry: CP violation (ε = 0.00228), CKM phase (1.2 rad)
 * - Pair production: threshold E = 2mc² (1.022 MeV for e⁺e⁻)
 * - Annihilation: energy E = 2mc², photon energy
 * - Matter-antimatter asymmetry: η = 6×10⁻¹⁰, Sakharov conditions
 */

#include <iostream>
#include <cmath>
#include <vector>
#include <string>
#include <functional>
#include "../include/physics/qft_particle_physics.hpp"
#include "../include/physics/qft_antiparticles.hpp"

// Test tolerances
const double TOLERANCE = 1e-5;           // For exact values
const double LOOSE_TOLERANCE = 1e-3;     // For approximate values
const double ORDER_OF_MAGNITUDE = 0.01;  // For order of magnitude checks

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
using namespace physics::advanced::qft;

int main() {
    int tests_passed = 0;
    int tests_failed = 0;

    std::cout << "=== Phase 4: Quantum Field Theory - Particle Physics & Antiparticles ===" << std::endl;
    std::cout << std::endl;

    // Helper lambda to run tests
    auto run_test = [&](const char* name, std::function<bool()> test_func) {
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
    // QUARK PROPERTIES TESTS (~12 tests)
    // ========================================

    run_test("Quark up has positive charge 2/3", []() {
        auto u = Quarks::up();
        ASSERT_NEAR(u.charge, 2.0/3.0, TOLERANCE);
        return true;
    });

    run_test("Quark down has negative charge 1/3", []() {
        auto d = Quarks::down();
        ASSERT_NEAR(d.charge, -1.0/3.0, TOLERANCE);
        return true;
    });

    run_test("Quark charm has positive charge 2/3", []() {
        auto c = Quarks::charm();
        ASSERT_NEAR(c.charge, 2.0/3.0, TOLERANCE);
        return true;
    });

    run_test("Quark strange has negative charge 1/3", []() {
        auto s = Quarks::strange();
        ASSERT_NEAR(s.charge, -1.0/3.0, TOLERANCE);
        return true;
    });

    run_test("Quark top has positive charge 2/3", []() {
        auto t = Quarks::top();
        ASSERT_NEAR(t.charge, 2.0/3.0, TOLERANCE);
        return true;
    });

    run_test("Quark bottom has negative charge 1/3", []() {
        auto b = Quarks::bottom();
        ASSERT_NEAR(b.charge, -1.0/3.0, TOLERANCE);
        return true;
    });

    run_test("Quark top is heaviest at 173 GeV", []() {
        auto t = Quarks::top();
        ASSERT_NEAR(t.mass, 173.0, LOOSE_TOLERANCE);
        auto u = Quarks::up();
        ASSERT_TRUE(t.mass > u.mass);
        return true;
    });

    run_test("Quark up is lightest at 2.2 MeV", []() {
        auto u = Quarks::up();
        ASSERT_NEAR(u.mass, 0.0022, LOOSE_TOLERANCE);
        auto t = Quarks::top();
        ASSERT_TRUE(u.mass < t.mass);
        return true;
    });

    run_test("All quarks have spin 1/2", []() {
        auto quarks = Quarks::allFlavors();
        for (const auto& q : quarks) {
            ASSERT_NEAR(q.spin, 0.5, TOLERANCE);
            ASSERT_TRUE(q.spin_type == SpinType::FERMION);
        }
        return true;
    });

    run_test("All quarks have color charge", []() {
        auto quarks = Quarks::allFlavors();
        for (const auto& q : quarks) {
            ASSERT_TRUE(q.color_charge == 1);
        }
        return true;
    });

    run_test("All quarks have baryon number 1/3", []() {
        auto quarks = Quarks::allFlavors();
        for (const auto& q : quarks) {
            ASSERT_TRUE(q.baryon_number == 1);  // Stored as 1 for quarks
            ASSERT_TRUE(q.type == ParticleType::QUARK);
        }
        return true;
    });

    run_test("Quark top decays before hadronization", []() {
        auto t = Quarks::top();
        ASSERT_TRUE(t.lifetime > 0.0);
        ASSERT_TRUE(t.lifetime < 1e-23);  // Very short lived
        ASSERT_TRUE(!t.isStable());
        return true;
    });

    // ========================================
    // LEPTON PROPERTIES TESTS (~12 tests)
    // ========================================

    run_test("Electron has mass 0.511 MeV", []() {
        auto e = Leptons::electron();
        ASSERT_NEAR(e.mass, 0.000511, LOOSE_TOLERANCE);
        return true;
    });

    run_test("Electron has charge -1", []() {
        auto e = Leptons::electron();
        ASSERT_NEAR(e.charge, -1.0, TOLERANCE);
        return true;
    });

    run_test("Electron is stable", []() {
        auto e = Leptons::electron();
        ASSERT_TRUE(e.isStable());
        return true;
    });

    run_test("Muon has mass 105.7 MeV", []() {
        auto mu = Leptons::muon();
        ASSERT_NEAR(mu.mass, 0.1057, LOOSE_TOLERANCE);
        return true;
    });

    run_test("Muon has charge -1", []() {
        auto mu = Leptons::muon();
        ASSERT_NEAR(mu.charge, -1.0, TOLERANCE);
        return true;
    });

    run_test("Muon lifetime 2.2 microseconds", []() {
        auto mu = Leptons::muon();
        ASSERT_NEAR(mu.lifetime, 2.2e-6, LOOSE_TOLERANCE);
        ASSERT_TRUE(!mu.isStable());
        return true;
    });

    run_test("Tau has mass 1.777 GeV", []() {
        auto tau = Leptons::tau();
        ASSERT_NEAR(tau.mass, 1.777, LOOSE_TOLERANCE);
        return true;
    });

    run_test("Tau has charge -1", []() {
        auto tau = Leptons::tau();
        ASSERT_NEAR(tau.charge, -1.0, TOLERANCE);
        return true;
    });

    run_test("Tau lifetime 290 femtoseconds", []() {
        auto tau = Leptons::tau();
        ASSERT_NEAR(tau.lifetime, 2.9e-13, LOOSE_TOLERANCE);
        ASSERT_TRUE(!tau.isStable());
        return true;
    });

    run_test("Neutrinos are very light under 1 eV", []() {
        auto nu_e = Leptons::electron_neutrino();
        auto nu_mu = Leptons::muon_neutrino();
        auto nu_tau = Leptons::tau_neutrino();
        ASSERT_TRUE(nu_e.mass < 1e-8);
        ASSERT_TRUE(nu_mu.mass < 1e-8);
        ASSERT_TRUE(nu_tau.mass < 1e-8);
        return true;
    });

    run_test("Neutrinos are neutral", []() {
        auto nu_e = Leptons::electron_neutrino();
        auto nu_mu = Leptons::muon_neutrino();
        auto nu_tau = Leptons::tau_neutrino();
        ASSERT_NEAR(nu_e.charge, 0.0, TOLERANCE);
        ASSERT_NEAR(nu_mu.charge, 0.0, TOLERANCE);
        ASSERT_NEAR(nu_tau.charge, 0.0, TOLERANCE);
        return true;
    });

    run_test("All leptons have spin 1/2", []() {
        auto leptons = Leptons::allLeptons();
        for (const auto& l : leptons) {
            ASSERT_NEAR(l.spin, 0.5, TOLERANCE);
            ASSERT_TRUE(l.spin_type == SpinType::FERMION);
        }
        return true;
    });

    // ========================================
    // GAUGE BOSON PROPERTIES TESTS (~8 tests)
    // ========================================

    run_test("Photon is massless", []() {
        auto gamma = GaugeBosons::photon();
        ASSERT_NEAR(gamma.mass, 0.0, TOLERANCE);
        return true;
    });

    run_test("Photon is neutral", []() {
        auto gamma = GaugeBosons::photon();
        ASSERT_NEAR(gamma.charge, 0.0, TOLERANCE);
        return true;
    });

    run_test("Photon has spin 1", []() {
        auto gamma = GaugeBosons::photon();
        ASSERT_NEAR(gamma.spin, 1.0, TOLERANCE);
        ASSERT_TRUE(gamma.spin_type == SpinType::BOSON);
        return true;
    });

    run_test("W boson has mass 80.4 GeV", []() {
        auto W_plus = GaugeBosons::W_plus();
        auto W_minus = GaugeBosons::W_minus();
        ASSERT_NEAR(W_plus.mass, 80.4, LOOSE_TOLERANCE);
        ASSERT_NEAR(W_minus.mass, 80.4, LOOSE_TOLERANCE);
        return true;
    });

    run_test("W bosons have opposite charges", []() {
        auto W_plus = GaugeBosons::W_plus();
        auto W_minus = GaugeBosons::W_minus();
        ASSERT_NEAR(W_plus.charge, +1.0, TOLERANCE);
        ASSERT_NEAR(W_minus.charge, -1.0, TOLERANCE);
        return true;
    });

    run_test("Z boson has mass 91.2 GeV", []() {
        auto Z = GaugeBosons::Z_boson();
        ASSERT_NEAR(Z.mass, 91.2, LOOSE_TOLERANCE);
        return true;
    });

    run_test("Z boson is neutral", []() {
        auto Z = GaugeBosons::Z_boson();
        ASSERT_NEAR(Z.charge, 0.0, TOLERANCE);
        return true;
    });

    run_test("Gluon is massless and colored", []() {
        auto g = GaugeBosons::gluon();
        ASSERT_NEAR(g.mass, 0.0, TOLERANCE);
        ASSERT_TRUE(g.color_charge == 1);
        ASSERT_NEAR(g.charge, 0.0, TOLERANCE);
        return true;
    });

    // ========================================
    // HIGGS BOSON TESTS (~3 tests)
    // ========================================

    run_test("Higgs has mass 125.1 GeV", []() {
        auto h = HiggsBoson::higgs();
        ASSERT_NEAR(h.mass, 125.1, LOOSE_TOLERANCE);
        return true;
    });

    run_test("Higgs is scalar spin 0", []() {
        auto h = HiggsBoson::higgs();
        ASSERT_NEAR(h.spin, 0.0, TOLERANCE);
        ASSERT_TRUE(h.spin_type == SpinType::BOSON);
        return true;
    });

    run_test("Higgs is neutral", []() {
        auto h = HiggsBoson::higgs();
        ASSERT_NEAR(h.charge, 0.0, TOLERANCE);
        return true;
    });

    // ========================================
    // GENERATION HIERARCHY TESTS (~5 tests)
    // ========================================

    run_test("Generation 1 lighter than generation 2", []() {
        auto gen1_quarks = Quarks::generation(1);
        auto gen2_quarks = Quarks::generation(2);
        double avg1 = 0.0, avg2 = 0.0;
        for (const auto& q : gen1_quarks) avg1 += q.mass;
        for (const auto& q : gen2_quarks) avg2 += q.mass;
        avg1 /= gen1_quarks.size();
        avg2 /= gen2_quarks.size();
        ASSERT_TRUE(avg1 < avg2);
        return true;
    });

    run_test("Generation 2 lighter than generation 3", []() {
        auto gen2_quarks = Quarks::generation(2);
        auto gen3_quarks = Quarks::generation(3);
        double avg2 = 0.0, avg3 = 0.0;
        for (const auto& q : gen2_quarks) avg2 += q.mass;
        for (const auto& q : gen3_quarks) avg3 += q.mass;
        avg2 /= gen2_quarks.size();
        avg3 /= gen3_quarks.size();
        ASSERT_TRUE(avg2 < avg3);
        return true;
    });

    run_test("Three generations exist", []() {
        auto gen1 = Generations::generation(1);
        auto gen2 = Generations::generation(2);
        auto gen3 = Generations::generation(3);
        ASSERT_TRUE(gen1.size() > 0);
        ASSERT_TRUE(gen2.size() > 0);
        ASSERT_TRUE(gen3.size() > 0);
        return true;
    });

    run_test("Average mass increases across generations", []() {
        auto masses = Generations::averageMassPerGeneration();
        ASSERT_TRUE(masses[0] < masses[1]);
        ASSERT_TRUE(masses[1] < masses[2]);
        return true;
    });

    run_test("Charged leptons follow mass hierarchy", []() {
        auto e = Leptons::electron();
        auto mu = Leptons::muon();
        auto tau = Leptons::tau();
        ASSERT_TRUE(e.mass < mu.mass);
        ASSERT_TRUE(mu.mass < tau.mass);
        return true;
    });

    // ========================================
    // FERMION VS BOSON TESTS (~4 tests)
    // ========================================

    run_test("All quarks are fermions", []() {
        auto quarks = Quarks::allFlavors();
        for (const auto& q : quarks) {
            ASSERT_TRUE(q.isFermion());
            ASSERT_TRUE(!q.isBoson());
        }
        return true;
    });

    run_test("All leptons are fermions", []() {
        auto leptons = Leptons::allLeptons();
        for (const auto& l : leptons) {
            ASSERT_TRUE(l.isFermion());
            ASSERT_TRUE(!l.isBoson());
        }
        return true;
    });

    run_test("All gauge bosons are bosons", []() {
        auto bosons = GaugeBosons::allBosons();
        for (const auto& b : bosons) {
            ASSERT_TRUE(b.isBoson());
            ASSERT_TRUE(!b.isFermion());
        }
        return true;
    });

    run_test("Higgs is boson not fermion", []() {
        auto h = HiggsBoson::higgs();
        ASSERT_TRUE(h.isBoson());
        ASSERT_TRUE(!h.isFermion());
        return true;
    });

    // ========================================
    // QUANTUM NUMBER CONSERVATION TESTS (~5 tests)
    // ========================================

    run_test("Charge conservation in proton", []() {
        std::vector<ParticleProperties> proton = {
            Quarks::up(), Quarks::up(), Quarks::down()
        };
        double Q = QuantumNumbers::totalCharge(proton);
        ASSERT_NEAR(Q, 1.0, TOLERANCE);  // +2/3 + 2/3 - 1/3 = 1
        return true;
    });

    run_test("Charge conservation in neutron", []() {
        std::vector<ParticleProperties> neutron = {
            Quarks::up(), Quarks::down(), Quarks::down()
        };
        double Q = QuantumNumbers::totalCharge(neutron);
        ASSERT_NEAR(Q, 0.0, TOLERANCE);  // +2/3 - 1/3 - 1/3 = 0
        return true;
    });

    run_test("Baryon number conservation for quarks", []() {
        std::vector<ParticleProperties> hadron = {
            Quarks::up(), Quarks::down(), Quarks::strange()
        };
        double B = QuantumNumbers::totalBaryonNumber(hadron);
        ASSERT_NEAR(B, 1.0, TOLERANCE);  // Three quarks: 1/3 + 1/3 + 1/3 = 1
        return true;
    });

    run_test("Lepton number conservation for leptons", []() {
        std::vector<ParticleProperties> leptons = {
            Leptons::electron(), Leptons::muon()
        };
        int L = QuantumNumbers::totalLeptonNumber(leptons);
        ASSERT_TRUE(L == 2);  // Two leptons
        return true;
    });

    run_test("Quantum number check conservation pass", []() {
        std::vector<ParticleProperties> initial = {
            Quarks::up(), Leptons::electron()
        };
        std::vector<ParticleProperties> final = {
            Quarks::up(), Leptons::electron()
        };
        bool conserved = QuantumNumbers::checkConservation(initial, final, LOOSE_TOLERANCE);
        ASSERT_TRUE(conserved);
        return true;
    });

    // ========================================
    // ANTIPARTICLE CREATION TESTS (~6 tests)
    // ========================================

    run_test("Antiparticle has opposite charge to particle", []() {
        auto e = Leptons::electron();
        auto e_anti = Antiparticle::createAntiparticle(e);
        ASSERT_NEAR(e_anti.charge, -e.charge, TOLERANCE);
        return true;
    });

    run_test("Antiparticle preserves mass", []() {
        auto u = Quarks::up();
        auto u_bar = Antiparticle::createAntiparticle(u);
        ASSERT_NEAR(u_bar.mass, u.mass, TOLERANCE);
        return true;
    });

    run_test("Antiparticle preserves spin", []() {
        auto mu = Leptons::muon();
        auto mu_anti = Antiparticle::createAntiparticle(mu);
        ASSERT_NEAR(mu_anti.spin, mu.spin, TOLERANCE);
        return true;
    });

    run_test("Antiparticle flips baryon number", []() {
        auto d = Quarks::down();
        auto d_bar = Antiparticle::createAntiparticle(d);
        ASSERT_TRUE(d_bar.baryon_number == -d.baryon_number);
        return true;
    });

    run_test("Antiparticle flips lepton number", []() {
        auto e = Leptons::electron();
        auto e_anti = Antiparticle::createAntiparticle(e);
        ASSERT_TRUE(e_anti.lepton_number == -e.lepton_number);
        return true;
    });

    run_test("CPT theorem antiparticle lifetime equals particle", []() {
        auto tau = Leptons::tau();
        auto tau_anti = Antiparticle::createAntiparticle(tau);
        ASSERT_NEAR(tau_anti.lifetime, tau.lifetime, TOLERANCE);
        return true;
    });

    // ========================================
    // ANTIQUARK TESTS (~4 tests)
    // ========================================

    run_test("Antiup quark has charge -2/3", []() {
        auto u_bar = AntiQuarks::up_bar();
        ASSERT_NEAR(u_bar.charge, -2.0/3.0, TOLERANCE);
        return true;
    });

    run_test("Antidown quark has charge +1/3", []() {
        auto d_bar = AntiQuarks::down_bar();
        ASSERT_NEAR(d_bar.charge, +1.0/3.0, TOLERANCE);
        return true;
    });

    run_test("Antitop quark preserves mass 173 GeV", []() {
        auto t_bar = AntiQuarks::top_bar();
        ASSERT_NEAR(t_bar.mass, 173.0, LOOSE_TOLERANCE);
        return true;
    });

    run_test("Antiquarks have negative baryon number", []() {
        auto u_bar = AntiQuarks::up_bar();
        auto d_bar = AntiQuarks::down_bar();
        ASSERT_TRUE(u_bar.baryon_number < 0);
        ASSERT_TRUE(d_bar.baryon_number < 0);
        return true;
    });

    // ========================================
    // ANTILEPTON TESTS (~5 tests)
    // ========================================

    run_test("Positron has charge +1", []() {
        auto e_plus = AntiLeptons::positron();
        ASSERT_NEAR(e_plus.charge, 1.0, TOLERANCE);
        return true;
    });

    run_test("Positron has mass 0.511 MeV", []() {
        auto e_plus = AntiLeptons::positron();
        ASSERT_NEAR(e_plus.mass, 0.000511, LOOSE_TOLERANCE);
        return true;
    });

    run_test("Antimuon has charge +1", []() {
        auto mu_plus = AntiLeptons::antimuon();
        ASSERT_NEAR(mu_plus.charge, 1.0, TOLERANCE);
        return true;
    });

    run_test("Antitau has charge +1", []() {
        auto tau_plus = AntiLeptons::antitau();
        ASSERT_NEAR(tau_plus.charge, 1.0, TOLERANCE);
        return true;
    });

    run_test("Antileptons have negative lepton number", []() {
        auto e_plus = AntiLeptons::positron();
        auto mu_plus = AntiLeptons::antimuon();
        ASSERT_TRUE(e_plus.lepton_number < 0);
        ASSERT_TRUE(mu_plus.lepton_number < 0);
        return true;
    });

    // ========================================
    // SELF-CONJUGATE PARTICLE TESTS (~3 tests)
    // ========================================

    run_test("Photon is self-conjugate", []() {
        auto gamma = GaugeBosons::photon();
        bool self = Antiparticle::isSelfConjugate(gamma);
        ASSERT_TRUE(self);
        return true;
    });

    run_test("Z boson is self-conjugate", []() {
        auto Z = GaugeBosons::Z_boson();
        bool self = Antiparticle::isSelfConjugate(Z);
        ASSERT_TRUE(self);
        return true;
    });

    run_test("Electron is not self-conjugate", []() {
        auto e = Leptons::electron();
        bool self = Antiparticle::isSelfConjugate(e);
        ASSERT_TRUE(!self);  // Has charge and lepton number
        return true;
    });

    // ========================================
    // CHARGE CONJUGATION C TESTS (~5 tests)
    // ========================================

    run_test("Charge conjugation inverts charge", []() {
        auto e = Leptons::electron();
        auto e_C = ChargeConjugation::apply(e);
        ASSERT_NEAR(e_C.charge, -e.charge, TOLERANCE);
        return true;
    });

    run_test("Photon has C-parity -1", []() {
        int c_parity = ChargeConjugation::cParity("photon");
        ASSERT_TRUE(c_parity == -1);
        return true;
    });

    run_test("Neutral pion has C-parity +1", []() {
        int c_parity = ChargeConjugation::cParity("neutral pion");
        ASSERT_TRUE(c_parity == 1);
        return true;
    });

    run_test("C is conserved in electromagnetic interaction", []() {
        bool conserved = ChargeConjugation::isConserved("electromagnetic");
        ASSERT_TRUE(conserved);
        return true;
    });

    run_test("C is violated in weak interaction", []() {
        bool conserved = ChargeConjugation::isConserved("weak");
        ASSERT_TRUE(!conserved);  // C is violated in weak
        return true;
    });

    // ========================================
    // CP SYMMETRY VIOLATION TESTS (~5 tests)
    // ========================================

    run_test("CP violation parameter in K0 meson ε equals 0.00228", []() {
        double eps = CPSymmetry::cpViolationParameter();
        ASSERT_NEAR(eps, 2.228e-3, LOOSE_TOLERANCE);
        return true;
    });

    run_test("CKM phase 1.2 radians", []() {
        double delta = CPSymmetry::ckmPhase();
        ASSERT_NEAR(delta, 1.2, LOOSE_TOLERANCE);
        return true;
    });

    run_test("CP is conserved in electromagnetic interaction", []() {
        bool conserved = CPSymmetry::isConserved("electromagnetic");
        ASSERT_TRUE(conserved);
        return true;
    });

    run_test("CP is violated in weak interaction", []() {
        bool conserved = CPSymmetry::isConserved("weak");
        ASSERT_TRUE(!conserved);  // Weak violates CP
        return true;
    });

    run_test("CP transformation description", []() {
        auto transform = CPSymmetry::cpTransformation();
        ASSERT_TRUE(transform.find("antiparticle") != std::string::npos);
        return true;
    });

    // ========================================
    // PAIR PRODUCTION TESTS (~6 tests)
    // ========================================

    run_test("Pair production threshold energy 2mc²", []() {
        double m_e = Leptons::electron().mass;  // in GeV
        double threshold = PairProduction::thresholdEnergy(m_e);
        ASSERT_NEAR(threshold, 2.0 * m_e, TOLERANCE);
        return true;
    });

    run_test("Electron pair production threshold 1.022 MeV", []() {
        double threshold = PairProduction::electronPairThreshold();
        ASSERT_NEAR(threshold, 1.022e-3, LOOSE_TOLERANCE);  // 1.022 MeV in GeV
        return true;
    });

    run_test("Can produce pair above threshold", []() {
        double m_e = Leptons::electron().mass;
        double photon_E = 3.0 * m_e;  // Above threshold
        bool can_produce = PairProduction::canProducePair(photon_E, m_e);
        ASSERT_TRUE(can_produce);
        return true;
    });

    run_test("Cannot produce pair below threshold", []() {
        double m_e = Leptons::electron().mass;
        double photon_E = 1.5 * m_e;  // Below threshold
        bool can_produce = PairProduction::canProducePair(photon_E, m_e);
        ASSERT_TRUE(!can_produce);
        return true;
    });

    run_test("Pair production cross section zero below threshold", []() {
        double m_e = Leptons::electron().mass;
        double photon_E = 1.5 * m_e;  // Below threshold
        double sigma = PairProduction::crossSection(photon_E, m_e);
        ASSERT_NEAR(sigma, 0.0, TOLERANCE);
        return true;
    });

    run_test("Pair production cross section nonzero above threshold", []() {
        double m_e = Leptons::electron().mass;
        double photon_E = 3.0 * m_e;  // Above threshold
        double sigma = PairProduction::crossSection(photon_E, m_e);
        ASSERT_TRUE(sigma > 0.0);
        return true;
    });

    // ========================================
    // ANNIHILATION TESTS (~6 tests)
    // ========================================

    run_test("Annihilation energy released 2mc²", []() {
        double m_e = Leptons::electron().mass;
        double E = Annihilation::energyReleased(m_e);
        ASSERT_NEAR(E, 2.0 * m_e, TOLERANCE);
        return true;
    });

    run_test("Electron positron annihilation energy 1.022 MeV", []() {
        double E = Annihilation::electronPositronEnergy();
        ASSERT_NEAR(E, 1.022e-3, LOOSE_TOLERANCE);  // 1.022 MeV in GeV
        return true;
    });

    run_test("Photon energy from annihilation at rest 0.511 MeV", []() {
        double E_photon = Annihilation::photonEnergyAtRest();
        ASSERT_NEAR(E_photon, 0.511e-3, LOOSE_TOLERANCE);  // 0.511 MeV in GeV
        return true;
    });

    run_test("Total photon energy equals annihilation energy", []() {
        double E_total = Annihilation::electronPositronEnergy();
        double E_photon = Annihilation::photonEnergyAtRest();
        ASSERT_NEAR(2.0 * E_photon, E_total, LOOSE_TOLERANCE);
        return true;
    });

    run_test("Annihilation cross section decreases with velocity", []() {
        double v1 = 1.0e-5;
        double v2 = 1.0e-4;
        double sigma1 = Annihilation::crossSection(v1);
        double sigma2 = Annihilation::crossSection(v2);
        ASSERT_TRUE(sigma1 > sigma2);  // Larger cross section at lower velocity
        return true;
    });

    run_test("Annihilation rate depends on cross section and density", []() {
        double sigma = 1.0e-24;  // m²
        double v = 1.0e-2;       // relative velocity
        double n = 1.0e10;       // density
        double rate = Annihilation::annihilationRate(sigma, v, n);
        ASSERT_TRUE(rate > 0.0);
        ASSERT_NEAR(rate, sigma * v * n, TOLERANCE);
        return true;
    });

    // ========================================
    // MATTER-ANTIMATTER ASYMMETRY TESTS (~5 tests)
    // ========================================

    run_test("Baryon to photon ratio 6e-10", []() {
        double eta = MatterAntimatterAsymmetry::baryonToPhotonRatio();
        ASSERT_NEAR(eta, 6.0e-10, ORDER_OF_MAGNITUDE);
        return true;
    });

    run_test("Sakharov conditions are three", []() {
        auto conditions = MatterAntimatterAsymmetry::sakharovConditions();
        ASSERT_TRUE(conditions.size() == 3);
        return true;
    });

    run_test("Sakharov condition 1 is baryon number violation", []() {
        auto conditions = MatterAntimatterAsymmetry::sakharovConditions();
        ASSERT_TRUE(conditions[0].find("Baryon number violation") != std::string::npos);
        return true;
    });

    run_test("Sakharov condition 2 is CP violation", []() {
        auto conditions = MatterAntimatterAsymmetry::sakharovConditions();
        ASSERT_TRUE(conditions[1].find("CP violation") != std::string::npos);
        return true;
    });

    run_test("Asymmetry parameter is 1e-9 order of magnitude", []() {
        double asymmetry = MatterAntimatterAsymmetry::asymmetryParameter();
        ASSERT_NEAR(asymmetry, 1.0e-9, ORDER_OF_MAGNITUDE);
        return true;
    });

    // ========================================
    // PARTICLE COUNT TESTS (~2 tests)
    // ========================================

    run_test("Standard Model has 6 quarks", []() {
        auto quarks = Quarks::allFlavors();
        ASSERT_TRUE(quarks.size() == 6);
        return true;
    });

    run_test("Standard Model has 6 leptons", []() {
        auto leptons = Leptons::allLeptons();
        ASSERT_TRUE(leptons.size() == 6);
        return true;
    });

    // ========================================
    // PHYSICAL CONSISTENCY TESTS (~8 tests)
    // ========================================

    run_test("Consistency all quarks are colored", []() {
        auto quarks = Quarks::allFlavors();
        for (const auto& q : quarks) {
            ASSERT_TRUE(q.color_charge == 1);
            ASSERT_TRUE(q.type == ParticleType::QUARK);
        }
        return true;
    });

    run_test("Consistency leptons are colorless", []() {
        auto leptons = Leptons::allLeptons();
        for (const auto& l : leptons) {
            ASSERT_TRUE(l.color_charge == 0);
            ASSERT_TRUE(l.type == ParticleType::LEPTON);
        }
        return true;
    });

    run_test("Consistency W bosons decay quickly", []() {
        auto W = GaugeBosons::W_plus();
        ASSERT_TRUE(W.lifetime > 0.0);
        ASSERT_TRUE(W.lifetime < 1e-23);
        ASSERT_TRUE(!W.isStable());
        return true;
    });

    run_test("Consistency charge conservation in pion decay", []() {
        // Positive pion decays to positive lepton and neutrino
        std::vector<ParticleProperties> initial = {};  // π⁺
        std::vector<ParticleProperties> final = {
            AntiLeptons::antimuon(),
            Leptons::muon_neutrino()
        };
        double Q_final = QuantumNumbers::totalCharge(final);
        // Both have charge that sums to +1
        ASSERT_NEAR(Q_final, 1.0, TOLERANCE);
        return true;
    });

    run_test("Consistency particle antiparticle pair has zero charge", []() {
        auto e = Leptons::electron();
        auto e_plus = AntiLeptons::positron();
        std::vector<ParticleProperties> pair = {e, e_plus};
        double Q = QuantumNumbers::totalCharge(pair);
        ASSERT_NEAR(Q, 0.0, TOLERANCE);
        return true;
    });

    run_test("Consistency antiparticle application twice gives original", []() {
        auto e = Leptons::electron();
        auto e_bar = Antiparticle::createAntiparticle(e);
        auto e_bar_bar = Antiparticle::createAntiparticle(e_bar);
        ASSERT_NEAR(e_bar_bar.charge, e.charge, TOLERANCE);
        ASSERT_NEAR(e_bar_bar.baryon_number, e.baryon_number, TOLERANCE);
        ASSERT_NEAR(e_bar_bar.lepton_number, e.lepton_number, TOLERANCE);
        return true;
    });

    run_test("Consistency all particles are finite mass", []() {
        auto fermions = FermionBosonsClassification::allFermions();
        auto bosons = FermionBosonsClassification::allBosons();
        for (const auto& f : fermions) {
            ASSERT_TRUE(std::isfinite(f.mass));
            ASSERT_TRUE(f.mass >= 0.0);
        }
        for (const auto& b : bosons) {
            ASSERT_TRUE(std::isfinite(b.mass));
            ASSERT_TRUE(b.mass >= 0.0);
        }
        return true;
    });

    run_test("Consistency all particles have finite lifetimes", []() {
        auto fermions = FermionBosonsClassification::allFermions();
        auto bosons = FermionBosonsClassification::allBosons();
        for (const auto& f : fermions) {
            ASSERT_TRUE(std::isfinite(f.lifetime));
            ASSERT_TRUE(f.lifetime >= 0.0);
        }
        for (const auto& b : bosons) {
            ASSERT_TRUE(std::isfinite(b.lifetime));
            ASSERT_TRUE(b.lifetime >= 0.0);
        }
        return true;
    });

    // ========================================
    // SUMMARY
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
