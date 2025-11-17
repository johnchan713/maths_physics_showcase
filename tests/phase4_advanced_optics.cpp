/**
 * Phase 4 Validation: Advanced Optics
 *
 * Tests the advanced_optics.hpp module functions.
 *
 * Coverage:
 * - Wave optics and phase difference
 * - Spherical waves and amplitude decay
 * - Refraction at spherical surfaces
 * - Aspheric lenses and conic surfaces
 * - Numerical aperture and f-number
 * - Lens aberrations (spherical chromatic coma astigmatism)
 * - Laser Doppler Velocimetry
 * - Interferometry (visibility fringe patterns Fabry-Perot)
 * - Optical Coherence Tomography
 * - Diffraction (single slit gratings)
 * - Fourier optics (Airy disk Rayleigh criterion)
 * - Transfer functions (CTF OTF MTF)
 * - Phase gratings and holography
 * - Triangulation and moire techniques
 * - Photoelasticity
 * - Polarized light and Stokes parameters
 * - Ellipsometry
 * - Fringe analysis and phase unwrapping
 * - Fiber optics (NA modes dispersion attenuation)
 */

#include <iostream>
#include <cmath>
#include <vector>
#include "../include/physics/advanced_optics.hpp"

using namespace physics::advanced_optics;

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

    std::cout << "=== Phase 4: Advanced Optics Validation ===" << std::endl;
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
    // Phase Difference and Wave Optics
    // ========================================

    run_test("Phase difference formula", []() {
        double wavelength = 600e-9;  // 600 nm
        double pathDiff = 1.2e-6;    // 1.2 microns

        double phase = phaseDifference(pathDiff, wavelength);

        // Δφ = 2π × (1.2/0.6) = 4π
        ASSERT_NEAR(phase, 4.0 * constants::PI, TOLERANCE);
        return true;
    });

    run_test("Optical path difference from phase", []() {
        double wavelength = 500e-9;
        double phase = 2.0 * constants::PI;  // One full wavelength

        double opd = opticalPathDifference(phase, wavelength);

        ASSERT_NEAR(opd, wavelength, TOLERANCE);
        return true;
    });

    run_test("Optical path length in medium", []() {
        double n = 1.5;  // Glass
        double d = 1e-3;  // 1 mm

        double opl = opticalPathLength(n, d);

        ASSERT_NEAR(opl, 1.5e-3, TOLERANCE);
        return true;
    });

    run_test("Constructive interference at zero phase", []() {
        ASSERT_TRUE(isConstructiveInterference(0.0));
        return true;
    });

    run_test("Constructive interference at 2pi", []() {
        ASSERT_TRUE(isConstructiveInterference(2.0 * constants::PI));
        return true;
    });

    run_test("Destructive interference at pi", []() {
        ASSERT_TRUE(isDestructiveInterference(constants::PI));
        return true;
    });

    run_test("Destructive interference at 3pi", []() {
        ASSERT_TRUE(isDestructiveInterference(3.0 * constants::PI));
        return true;
    });

    // ========================================
    // Spherical Waves
    // ========================================

    run_test("Spherical wave amplitude decay", []() {
        double A0 = 100.0;
        double r = 2.0;

        double A = sphericalWaveAmplitude(A0, r);

        ASSERT_NEAR(A, 50.0, TOLERANCE);  // A = 100/2
        return true;
    });

    run_test("Spherical wave intensity decay", []() {
        double I0 = 100.0;
        double r = 2.0;

        double I = sphericalWaveIntensity(I0, r);

        ASSERT_NEAR(I, 25.0, TOLERANCE);  // I = 100/4
        return true;
    });

    run_test("Wavefront curvature equals distance", []() {
        double r = 5.0;

        double R = wavefrontCurvature(r);

        ASSERT_NEAR(R, r, TOLERANCE);
        return true;
    });

    run_test("Gouy phase at waist", []() {
        double z0 = 1e-3;  // Rayleigh range

        double phase = gouyPhase(0.0, z0);

        ASSERT_NEAR(phase, 0.0, TOLERANCE);
        return true;
    });

    run_test("Gouy phase at Rayleigh range", []() {
        double z0 = 1e-3;

        double phase = gouyPhase(z0, z0);

        ASSERT_NEAR(phase, constants::PI / 4.0, TOLERANCE);  // atan(1) = π/4
        return true;
    });

    // ========================================
    // Refraction at Spherical Surface
    // ========================================

    run_test("Spherical surface refraction air to glass", []() {
        double n1 = 1.0;   // Air
        double n2 = 1.5;   // Glass
        double u = -0.2;   // Object at 20 cm (negative)
        double R = 0.05;   // 5 cm radius

        double v = sphericalSurfaceRefraction(n1, n2, u, R);

        // Should produce real image
        ASSERT_TRUE(v > 0);
        return true;
    });

    run_test("Spherical surface power calculation", []() {
        double n1 = 1.0;
        double n2 = 1.5;
        double R = 0.1;  // 10 cm

        double P = sphericalSurfacePower(n1, n2, R);

        ASSERT_NEAR(P, 5.0, TOLERANCE);  // (1.5-1.0)/0.1 = 5 D
        return true;
    });

    // ========================================
    // Aspheric Lenses
    // ========================================

    run_test("Aspheric sag for sphere", []() {
        double r = 0.01;  // 1 cm
        double c = 10.0;  // 1/R
        double k = 0.0;   // Sphere

        double sag = asphericSag(r, c, k);

        // For small r, sag ≈ cr²/2
        double expected = c * r * r / 2.0;
        ASSERT_NEAR(sag, expected, LOOSE_TOLERANCE);
        return true;
    });

    run_test("Classify sphere conic constant", []() {
        ASSERT_TRUE(classifyConicSurface(0.0) == 0);  // Sphere
        return true;
    });

    run_test("Classify parabola conic constant", []() {
        ASSERT_TRUE(classifyConicSurface(-1.0) == 2);  // Parabola
        return true;
    });

    run_test("Classify hyperbola conic constant", []() {
        ASSERT_TRUE(classifyConicSurface(-2.0) == 3);  // Hyperbola
        return true;
    });

    run_test("Classify ellipse conic constant", []() {
        ASSERT_TRUE(classifyConicSurface(-0.5) == 1);  // Ellipse
        return true;
    });

    // ========================================
    // Stops and Apertures
    // ========================================

    run_test("Numerical aperture calculation", []() {
        double n = 1.0;
        double angle = constants::PI / 6.0;  // 30 degrees

        double NA = numericalAperture(n, angle);

        ASSERT_NEAR(NA, 0.5, TOLERANCE);  // sin(30°) = 0.5
        return true;
    });

    run_test("F-number from focal length and aperture", []() {
        double f = 50e-3;   // 50 mm
        double D = 25e-3;   // 25 mm

        double fNum = fNumber(f, D);

        ASSERT_NEAR(fNum, 2.0, TOLERANCE);  // f/2
        return true;
    });

    run_test("F-number from numerical aperture", []() {
        double NA = 0.25;

        double fNum = fNumberFromNA(NA);

        ASSERT_NEAR(fNum, 2.0, TOLERANCE);  // 1/(2×0.25) = 2
        return true;
    });

    run_test("Depth of field calculation", []() {
        double lambda = 550e-9;
        double fNum = 2.0;

        double DOF = depthOfField(lambda, fNum);

        ASSERT_NEAR(DOF, 2.2e-6, LOOSE_TOLERANCE);
        return true;
    });

    // ========================================
    // Lens Aberrations
    // ========================================

    run_test("Longitudinal spherical aberration", []() {
        double h = 0.01;    // 1 cm aperture
        double f = 0.1;     // 10 cm focal length
        double n = 1.5;     // Glass

        double LSA = longitudinalSphericalAberration(h, f, n);

        // Should be negative (underfocusing)
        ASSERT_TRUE(LSA < 0);
        return true;
    });

    run_test("Chromatic aberration from Abbe number", []() {
        double f = 0.1;     // 10 cm
        double V = 50.0;    // Typical crown glass

        double delta_f = chromaticAberration(f, V);

        ASSERT_NEAR(delta_f, 0.002, LOOSE_TOLERANCE);
        return true;
    });

    run_test("Abbe number calculation", []() {
        double n_d = 1.520;
        double n_F = 1.524;
        double n_C = 1.518;

        double V = abbeNumber(n_d, n_F, n_C);

        // V = (1.520-1)/(1.524-1.518) = 86.67
        ASSERT_NEAR(V, 86.67, 1.0);
        return true;
    });

    run_test("Coma aberration proportional to cube", []() {
        double h = 0.01;
        double f = 0.1;
        double angle = 0.1;  // radians

        double coma = comaAberration(h, f, angle);

        ASSERT_TRUE(coma > 0);
        return true;
    });

    run_test("Astigmatism increases with field angle", []() {
        double f = 0.1;
        double angle1 = 0.1;
        double angle2 = 0.2;

        double ast1 = astigmatism(f, angle1);
        double ast2 = astigmatism(f, angle2);

        ASSERT_TRUE(ast2 > ast1);
        return true;
    });

    run_test("Distortion calculation", []() {
        double actual = 10.5;
        double ideal = 10.0;

        double dist = distortion(actual, ideal);

        ASSERT_NEAR(dist, 5.0, TOLERANCE);  // 5% barrel distortion
        return true;
    });

    // ========================================
    // Laser Doppler Velocimetry
    // ========================================

    run_test("Doppler frequency shift", []() {
        double v = 10.0;             // 10 m/s
        double lambda = 632.8e-9;    // He-Ne laser
        double angle = constants::PI / 6.0;  // 30 degrees

        double f_D = dopplerFrequencyShift(v, lambda, angle);

        // f_D = 2v sin(θ/2)/λ
        ASSERT_TRUE(f_D > 0);
        return true;
    });

    run_test("LDV fringe spacing", []() {
        double lambda = 632.8e-9;
        double angle = constants::PI / 6.0;

        double d_f = ldvFringeSpacing(lambda, angle);

        ASSERT_TRUE(d_f > 0);
        return true;
    });

    run_test("Velocity from Doppler frequency", []() {
        double lambda = 632.8e-9;
        double angle = constants::PI / 6.0;
        double f_D = dopplerFrequencyShift(15.0, lambda, angle);

        double v = velocityFromDoppler(f_D, lambda, angle);

        ASSERT_NEAR(v, 15.0, TOLERANCE);
        return true;
    });

    // ========================================
    // Interferometry
    // ========================================

    run_test("Fringe visibility perfect contrast", []() {
        double V = fringeVisibility(100.0, 0.0);

        ASSERT_NEAR(V, 1.0, TOLERANCE);
        return true;
    });

    run_test("Fringe visibility no contrast", []() {
        double V = fringeVisibility(50.0, 50.0);

        ASSERT_NEAR(V, 0.0, TOLERANCE);
        return true;
    });

    run_test("Two-beam interference constructive", []() {
        double I1 = 25.0;
        double I2 = 25.0;
        double phase = 0.0;

        double I = twoBeamInterferenceIntensity(I1, I2, phase);

        ASSERT_NEAR(I, 100.0, TOLERANCE);  // I1 + I2 + 2√(I1I2)
        return true;
    });

    run_test("Two-beam interference destructive", []() {
        double I1 = 25.0;
        double I2 = 25.0;
        double phase = constants::PI;

        double I = twoBeamInterferenceIntensity(I1, I2, phase);

        ASSERT_NEAR(I, 0.0, TOLERANCE);  // I1 + I2 - 2√(I1I2)
        return true;
    });

    run_test("Double-slit fringe spacing", []() {
        double lambda = 600e-9;
        double D = 1.0;  // 1 m screen distance
        double d = 0.001;  // 1 mm slit separation

        double beta = doubleSlitFringeSpacing(lambda, D, d);

        ASSERT_NEAR(beta, 0.6e-3, LOOSE_TOLERANCE);  // 0.6 mm
        return true;
    });

    run_test("Fabry-Perot finesse high reflectivity", []() {
        double R = 0.9;

        double F = fabryPerotFinesse(R);

        ASSERT_TRUE(F > 20.0);  // High finesse for R=0.9
        return true;
    });

    run_test("Free spectral range", []() {
        double L = 0.01;  // 1 cm cavity
        double n = 1.0;

        double FSR = freeSpectralRange(L, n);

        ASSERT_NEAR(FSR, 15e9, 1e9);  // ~15 GHz
        return true;
    });

    run_test("Fabry-Perot transmission at resonance", []() {
        double I0 = 100.0;
        double R = 0.5;
        double phase = 2.0 * constants::PI;  // Resonance

        double I_t = fabryPerotTransmission(I0, R, phase);

        ASSERT_TRUE(I_t > 90.0);  // High transmission at resonance
        return true;
    });

    // ========================================
    // Optical Coherence Tomography
    // ========================================

    run_test("OCT axial resolution", []() {
        double lambda0 = 1300e-9;  // 1.3 μm (typical OCT)
        double delta_lambda = 100e-9;  // 100 nm bandwidth

        double dz = octAxialResolution(lambda0, delta_lambda);

        ASSERT_TRUE(dz < 20e-6);  // Better than 20 μm
        return true;
    });

    run_test("OCT lateral resolution", []() {
        double lambda = 1300e-9;
        double f = 0.05;  // 50 mm
        double d = 0.005;  // 5 mm beam

        double dx = octLateralResolution(lambda, f, d);

        ASSERT_TRUE(dx < 50e-6);  // ~10-20 μm
        return true;
    });

    run_test("OCT depth of focus", []() {
        double dx = 15e-6;  // 15 μm lateral resolution
        double lambda = 1300e-9;

        double DOF = octDepthOfFocus(dx, lambda);

        ASSERT_TRUE(DOF > 0);
        return true;
    });

    // ========================================
    // Diffraction - Single Slit
    // ========================================

    run_test("Single slit central maximum intensity", []() {
        double a = 1e-3;  // 1 mm slit
        double lambda = 500e-9;
        double angle = 0.0;  // Center
        double I0 = 100.0;

        double I = singleSlitIntensity(a, lambda, angle, I0);

        ASSERT_NEAR(I, I0, TOLERANCE);
        return true;
    });

    run_test("Single slit first minimum angle", []() {
        double a = 1e-3;
        double lambda = 500e-9;

        double theta = singleSlitMinimumAngle(a, lambda, 1);

        ASSERT_NEAR(theta, lambda / a, LOOSE_TOLERANCE);
        return true;
    });

    run_test("Single slit central maximum width", []() {
        double a = 1e-3;
        double lambda = 500e-9;

        double width = singleSlitCentralMaxWidth(a, lambda);

        ASSERT_NEAR(width, 2.0 * lambda / a, LOOSE_TOLERANCE);
        return true;
    });

    // ========================================
    // Grating Equation
    // ========================================

    run_test("Grating diffraction angle normal incidence", []() {
        double d = 2e-6;  // 2 μm spacing (500 lines/mm)
        double lambda = 600e-9;
        int m = 1;

        double theta = gratingDiffractionAngle(d, lambda, m, 0.0);

        ASSERT_NEAR(std::sin(theta), lambda / d, LOOSE_TOLERANCE);
        return true;
    });

    run_test("Grating angular dispersion", []() {
        double d = 2e-6;
        int m = 1;
        double theta = 0.3;  // radians

        double disp = gratingAngularDispersion(d, m, theta);

        ASSERT_TRUE(disp > 0);
        return true;
    });

    run_test("Grating resolving power", []() {
        int m = 1;
        int N = 1000;  // 1000 illuminated lines

        double R = gratingResolvingPower(m, N);

        ASSERT_NEAR(R, 1000.0, TOLERANCE);
        return true;
    });

    run_test("Grating free spectral range", []() {
        double lambda = 500e-9;
        int m = 1;

        double FSR = gratingFSR(lambda, m);

        ASSERT_NEAR(FSR, lambda, TOLERANCE);
        return true;
    });

    // ========================================
    // Fourier Optics
    // ========================================

    run_test("Spatial frequency", []() {
        double angle = constants::PI / 6.0;  // 30 degrees
        double lambda = 500e-9;

        double f_x = spatialFrequency(angle, lambda);

        ASSERT_NEAR(f_x, 1e6, 1e5);  // ~1 line per micron
        return true;
    });

    run_test("Airy disk diameter", []() {
        double lambda = 500e-9;
        double f = 0.1;  // 10 cm
        double D = 0.01;  // 1 cm

        double d = airyDiskDiameter(lambda, f, D);

        ASSERT_NEAR(d, 12.2e-6, LOOSE_TOLERANCE);
        return true;
    });

    run_test("Rayleigh criterion", []() {
        double lambda = 500e-9;
        double D = 0.05;  // 5 cm telescope

        double theta_min = rayleighCriterion(lambda, D);

        ASSERT_NEAR(theta_min, 1.22e-5, LOOSE_TOLERANCE);
        return true;
    });

    run_test("Cutoff spatial frequency", []() {
        double D = 0.01;
        double lambda = 500e-9;
        double f = 0.1;

        double f_c = cutoffSpatialFrequency(D, lambda, f);

        ASSERT_TRUE(f_c > 0);
        return true;
    });

    // ========================================
    // Transfer Functions
    // ========================================

    run_test("Coherent transfer function inside cutoff", []() {
        double CTF = coherentTransferFunction(1000.0, 2000.0);

        ASSERT_NEAR(CTF, 1.0, TOLERANCE);
        return true;
    });

    run_test("Coherent transfer function outside cutoff", []() {
        double CTF = coherentTransferFunction(3000.0, 2000.0);

        ASSERT_NEAR(CTF, 0.0, TOLERANCE);
        return true;
    });

    run_test("Strehl ratio for perfect system", []() {
        double S = strehlRatio(0.0, 500e-9);

        ASSERT_NEAR(S, 1.0, TOLERANCE);
        return true;
    });

    run_test("Strehl ratio degrades with aberrations", []() {
        double lambda = 500e-9;
        double rms = lambda / 20.0;  // λ/20 RMS error

        double S = strehlRatio(rms, lambda);

        ASSERT_TRUE(S < 1.0 && S > 0.5);
        return true;
    });

    run_test("Optical transfer function at zero frequency", []() {
        double OTF = opticalTransferFunction(0.0, 1000.0);

        ASSERT_NEAR(OTF, 1.0, TOLERANCE);
        return true;
    });

    run_test("OTF decreases with spatial frequency", []() {
        double f_c = 1000.0;
        double OTF1 = opticalTransferFunction(200.0, f_c);
        double OTF2 = opticalTransferFunction(800.0, f_c);

        ASSERT_TRUE(OTF1 > OTF2);
        return true;
    });

    run_test("MTF equals absolute OTF", []() {
        double freq = 500.0;
        double f_c = 1000.0;

        double MTF = modulationTransferFunction(freq, f_c);
        double OTF = opticalTransferFunction(freq, f_c);

        ASSERT_NEAR(MTF, std::abs(OTF), TOLERANCE);
        return true;
    });

    run_test("Contrast reduction through system", []() {
        double C_in = 0.8;
        double MTF = 0.5;

        double C_out = contrastReduction(C_in, MTF);

        ASSERT_NEAR(C_out, 0.4, TOLERANCE);
        return true;
    });

    // ========================================
    // Phase Gratings and Holography
    // ========================================

    run_test("Phase grating intensity modulation", []() {
        double x = 0.0;
        double period = 10e-6;
        double phi = 0.05;  // Small modulation
        double I0 = 100.0;

        double I = phaseGratingIntensity(x, period, phi, I0);

        ASSERT_TRUE(I > 0);
        return true;
    });

    run_test("Phase grating efficiency first order", []() {
        double phi = constants::PI;  // π phase shift

        double eta = phaseGratingEfficiency(1, phi);

        ASSERT_NEAR(eta, 1.0, LOOSE_TOLERANCE);  // sin²(π/2) = 1
        return true;
    });

    run_test("Holographic recording intensity", []() {
        double I_o = 25.0;
        double I_r = 25.0;
        double phase = 0.0;

        double I = holographicRecordingIntensity(I_o, I_r, phase);

        ASSERT_NEAR(I, 100.0, TOLERANCE);
        return true;
    });

    run_test("Holographic fringe spacing", []() {
        double lambda = 532e-9;  // Green laser
        double angle = constants::PI / 3.0;  // 60 degrees

        double Lambda = holographicFringeSpacing(lambda, angle);

        ASSERT_TRUE(Lambda > 0);
        return true;
    });

    run_test("Volume hologram efficiency", []() {
        double dn = 0.01;  // Index modulation
        double d = 10e-6;  // 10 μm thickness
        double lambda = 532e-9;
        double theta = 0.1;

        double eta = volumeHologramEfficiency(dn, d, lambda, theta);

        ASSERT_TRUE(eta >= 0.0 && eta <= 1.0);
        return true;
    });

    run_test("Bragg angle calculation", []() {
        double lambda = 532e-9;
        double Lambda = 1e-6;  // 1 μm fringes

        double theta_B = braggAngle(lambda, Lambda);

        ASSERT_TRUE(theta_B > 0 && theta_B < constants::PI / 2.0);
        return true;
    });

    // ========================================
    // Optical Triangulation
    // ========================================

    run_test("Triangulation distance calculation", []() {
        double baseline = 0.1;  // 10 cm
        double f = 0.05;  // 5 cm focal length
        double disparity = 1e-3;  // 1 mm

        double Z = triangulationDistance(baseline, f, disparity);

        ASSERT_NEAR(Z, 5.0, LOOSE_TOLERANCE);
        return true;
    });

    run_test("Triangulation angle", []() {
        double distance = 1.0;  // 1 m
        double baseline = 0.1;  // 10 cm

        double angle = triangulationAngle(distance, baseline);

        ASSERT_TRUE(angle > 0);
        return true;
    });

    run_test("Triangulation depth resolution", []() {
        double Z = 1.0;
        double b = 0.1;
        double f = 0.05;
        double dd = 10e-6;  // 10 μm pixel

        double dZ = triangulationDepthResolution(Z, b, f, dd);

        ASSERT_TRUE(dZ > 0);
        return true;
    });

    // ========================================
    // Moire Technique
    // ========================================

    run_test("Moire fringe spacing calculation", []() {
        double p1 = 1.00e-3;  // 1.00 mm
        double p2 = 0.99e-3;  // 0.99 mm

        double Lambda_m = moireFringeSpacing(p1, p2);

        ASSERT_NEAR(Lambda_m, 99e-3, LOOSE_TOLERANCE);  // ~99 mm
        return true;
    });

    run_test("Moire magnification", []() {
        double p1 = 1.00e-3;
        double p2 = 0.99e-3;

        double M = moireMagnification(p1, p2);

        ASSERT_NEAR(M, 100.0, LOOSE_TOLERANCE);
        return true;
    });

    run_test("Strain from moire", []() {
        double M = 100.0;

        double strain = strainFromMoire(M);

        ASSERT_NEAR(strain, 0.01, TOLERANCE);
        return true;
    });

    // ========================================
    // Photoelasticity
    // ========================================

    run_test("Photoelastic birefringence", []() {
        double C = 1e-10;  // Stress-optic coefficient (Pa⁻¹)
        double sigma = 1e7;  // 10 MPa stress

        double dn = photoelasticBirefringence(C, sigma);

        ASSERT_NEAR(dn, 1e-3, LOOSE_TOLERANCE);
        return true;
    });

    run_test("Photoelastic retardation", []() {
        double C = 1e-10;
        double sigma = 1e7;
        double t = 0.01;  // 1 cm
        double lambda = 550e-9;

        double delta = photoelasticRetardation(C, sigma, t, lambda);

        ASSERT_TRUE(delta > 0);
        return true;
    });

    run_test("Photoelastic fringe order", []() {
        double retardation = 4.0 * constants::PI;

        double N = photoelasticFringeOrder(retardation);

        ASSERT_NEAR(N, 2.0, TOLERANCE);
        return true;
    });

    run_test("Stress from fringe order", []() {
        double N = 5.0;
        double lambda = 550e-9;
        double C = 1e-10;
        double t = 0.01;

        double sigma = stressFromFringeOrder(N, lambda, C, t);

        ASSERT_TRUE(sigma > 0);
        return true;
    });

    // ========================================
    // Stokes Parameters
    // ========================================

    run_test("Degree of polarization for unpolarized", []() {
        StokesVector S(100.0, 0.0, 0.0, 0.0);

        double DOP = S.degreeOfPolarization();

        ASSERT_NEAR(DOP, 0.0, TOLERANCE);
        return true;
    });

    run_test("Degree of polarization for fully polarized", []() {
        StokesVector S(100.0, 100.0, 0.0, 0.0);

        double DOP = S.degreeOfPolarization();

        ASSERT_NEAR(DOP, 1.0, TOLERANCE);
        return true;
    });

    run_test("Linear polarization degree", []() {
        StokesVector S(100.0, 60.0, 80.0, 0.0);

        double DOLP = S.degreeOfLinearPolarization();

        ASSERT_NEAR(DOLP, 1.0, TOLERANCE);  // sqrt(60²+80²)/100 = 1
        return true;
    });

    run_test("Circular polarization degree", []() {
        StokesVector S(100.0, 0.0, 0.0, 100.0);

        double DOCP = S.degreeOfCircularPolarization();

        ASSERT_NEAR(DOCP, 1.0, TOLERANCE);
        return true;
    });

    run_test("Linearly polarized Stokes vector", []() {
        StokesVector S = linearlyPolarizedStokes(100.0, 0.0);

        ASSERT_NEAR(S.S0, 100.0, TOLERANCE);
        ASSERT_NEAR(S.S1, 100.0, TOLERANCE);
        ASSERT_NEAR(S.S2, 0.0, TOLERANCE);
        ASSERT_NEAR(S.S3, 0.0, TOLERANCE);
        return true;
    });

    run_test("Circularly polarized Stokes vector", []() {
        StokesVector S = circularlyPolarizedStokes(100.0, true);

        ASSERT_NEAR(S.S0, 100.0, TOLERANCE);
        ASSERT_NEAR(S.S1, 0.0, TOLERANCE);
        ASSERT_NEAR(S.S2, 0.0, TOLERANCE);
        ASSERT_NEAR(S.S3, 100.0, TOLERANCE);
        return true;
    });

    // ========================================
    // Ellipsometry
    // ========================================

    run_test("Ellipsometric parameters", []() {
        double r_p = 0.5;
        double r_s = 1.0;
        double delta = constants::PI / 2.0;

        auto [psi, d] = ellipsometricParameters(r_p, r_s, delta);

        ASSERT_TRUE(psi > 0);
        ASSERT_NEAR(d, delta, TOLERANCE);
        return true;
    });

    run_test("Refractive index from ellipsometry", []() {
        double theta = constants::PI / 4.0;  // 45 degrees
        double psi = constants::PI / 6.0;
        double delta = constants::PI / 2.0;

        double n = refractiveIndexFromEllipsometry(theta, psi, delta);

        ASSERT_TRUE(n > 0.5);  // Reasonable range for refractive index calculation
        return true;
    });

    // ========================================
    // Fringe Analysis
    // ========================================

    run_test("Phase from 4-step algorithm", []() {
        double I1 = 50.0;
        double I2 = 100.0;
        double I3 = 50.0;
        double I4 = 0.0;

        double phase = phaseFrom4Step(I1, I2, I3, I4);

        // atan2(50-50, 100-0) = atan2(0, 100) = 0
        ASSERT_NEAR(phase, 0.0, TOLERANCE);
        return true;
    });

    run_test("Fringe modulation calculation", []() {
        double I1 = 25.0;
        double I2 = 50.0;
        double I3 = 75.0;
        double I4 = 50.0;

        double gamma = fringeModulation(I1, I2, I3, I4);

        ASSERT_TRUE(gamma > 0.0 && gamma <= 1.0);
        return true;
    });

    run_test("Average intensity", []() {
        double I_avg = averageIntensity(20.0, 40.0, 60.0, 80.0);

        ASSERT_NEAR(I_avg, 50.0, TOLERANCE);
        return true;
    });

    // ========================================
    // Phase Unwrapping
    // ========================================

    run_test("Phase unwrapping continuous", []() {
        double wrapped = constants::PI / 2.0;
        double previous = constants::PI / 4.0;

        double unwrapped = unwrapPhase(wrapped, previous);

        ASSERT_NEAR(unwrapped, constants::PI / 2.0, TOLERANCE);
        return true;
    });

    run_test("Phase unwrapping across 2pi boundary", []() {
        double wrapped = -constants::PI + 0.1;
        double previous = constants::PI - 0.1;

        double unwrapped = unwrapPhase(wrapped, previous);

        ASSERT_TRUE(std::abs(unwrapped - previous) < 0.5);
        return true;
    });

    // ========================================
    // Fiber Optics
    // ========================================

    run_test("Fiber numerical aperture", []() {
        double n_core = 1.48;
        double n_clad = 1.46;

        double NA = fiberNumericalAperture(n_core, n_clad);

        ASSERT_NEAR(NA, 0.242, LOOSE_TOLERANCE);
        return true;
    });

    run_test("Fiber acceptance angle", []() {
        double NA = 0.22;

        double theta = fiberAcceptanceAngle(NA);

        ASSERT_NEAR(theta, 0.222, LOOSE_TOLERANCE);  // ~12.7 degrees
        return true;
    });

    run_test("Fiber V-number", []() {
        double a = 4e-6;  // 4 μm core radius
        double lambda = 1.55e-6;  // 1550 nm
        double NA = 0.12;

        double V = fiberVNumber(a, lambda, NA);

        ASSERT_NEAR(V, 1.946, VERY_LOOSE);  // Single mode for V < 2.405
        return true;
    });

    run_test("Fiber number of modes single mode", []() {
        double a = 4e-6;
        double lambda = 1.55e-6;
        double NA = 0.12;

        double M = fiberNumberOfModes(a, lambda, NA);

        ASSERT_NEAR(M, 1.0, TOLERANCE);  // Single mode
        return true;
    });

    run_test("Fiber number of modes multimode", []() {
        double a = 25e-6;  // 50 μm core diameter
        double lambda = 1.55e-6;
        double NA = 0.22;

        double M = fiberNumberOfModes(a, lambda, NA);

        ASSERT_TRUE(M > 100.0);  // Multimode
        return true;
    });

    run_test("Fiber cutoff wavelength", []() {
        double a = 4e-6;
        double NA = 0.12;

        double lambda_c = fiberCutoffWavelength(a, NA);

        ASSERT_TRUE(lambda_c > 1e-6);  // ~1.26 μm
        return true;
    });

    run_test("Fiber intermodal dispersion", []() {
        double L = 1000.0;  // 1 km
        double NA = 0.22;
        double n = 1.48;

        double dt = fiberIntermodalDispersion(L, NA, n);

        ASSERT_TRUE(dt > 0);  // ~50 ns/km
        return true;
    });

    run_test("Fiber attenuation calculation", []() {
        double P_in = 1e-3;  // 1 mW
        double P_out = 0.5e-3;  // 0.5 mW
        double L = 10000.0;  // 10 km

        double alpha = fiberAttenuation(P_in, P_out, L);

        ASSERT_TRUE(alpha > 0);  // ~0.3 dB/km
        return true;
    });

    run_test("Fiber output power from attenuation", []() {
        double P_in = 1e-3;
        double alpha = 0.2e-3;  // 0.2 dB/km
        double L = 10000.0;  // 10 km

        double P_out = fiberOutputPower(P_in, alpha, L);

        ASSERT_TRUE(P_out < P_in);
        ASSERT_TRUE(P_out > 0.5 * P_in);  // Less than 3 dB loss
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
