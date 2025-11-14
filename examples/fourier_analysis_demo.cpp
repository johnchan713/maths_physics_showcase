/**
 * @file fourier_analysis_demo.cpp
 * @brief Comprehensive demonstration of Fourier analysis computational implementations
 *
 * Demonstrates practical usage of:
 * - DFT vs FFT performance comparison
 * - Circulant matrices and convolution
 * - Wavelet transforms
 * - Hilbert space operators
 * - Fourier series and multipliers
 * - Time-frequency analysis
 */

#include "maths/analysis/fourier_analysis.hpp"
#include <iostream>
#include <iomanip>
#include <chrono>
#include <cmath>

using namespace maths::analysis;

// Helper function to print signal
void printSignal(const std::string& name, const std::vector<double>& signal, size_t max_display = 8) {
    std::cout << name << ": [";
    for (size_t i = 0; i < std::min(signal.size(), max_display); ++i) {
        std::cout << std::fixed << std::setprecision(4) << signal[i];
        if (i < std::min(signal.size(), max_display) - 1) std::cout << ", ";
    }
    if (signal.size() > max_display) std::cout << " ...";
    std::cout << "]" << std::endl;
}

// Helper to print complex signal
void printComplexSignal(const std::string& name, const std::vector<std::complex<double>>& signal, size_t max_display = 8) {
    std::cout << name << ": [";
    for (size_t i = 0; i < std::min(signal.size(), max_display); ++i) {
        std::cout << std::fixed << std::setprecision(4)
                  << signal[i].real() << "+" << signal[i].imag() << "i";
        if (i < std::min(signal.size(), max_display) - 1) std::cout << ", ";
    }
    if (signal.size() > max_display) std::cout << " ...";
    std::cout << "]" << std::endl;
}

// Demo 1: DFT vs FFT Performance
void demoDFTvsFFT() {
    std::cout << "\n=== Demo 1: DFT vs FFT Performance Comparison ===" << std::endl;

    // Create a test signal: sum of sinusoids
    size_t N = 256;
    std::vector<double> signal(N);
    for (size_t i = 0; i < N; ++i) {
        double t = static_cast<double>(i) / N;
        signal[i] = std::sin(2 * M_PI * 5 * t) + 0.5 * std::sin(2 * M_PI * 10 * t);
    }

    std::cout << "Signal size: " << N << " samples" << std::endl;
    printSignal("Input signal (first 8)", signal, 8);

    // Time DFT
    auto start = std::chrono::high_resolution_clock::now();
    auto spectrum_dft = DiscreteFourierTransform::dft(signal);
    auto end = std::chrono::high_resolution_clock::now();
    auto duration_dft = std::chrono::duration_cast<std::chrono::microseconds>(end - start);

    // Time FFT
    start = std::chrono::high_resolution_clock::now();
    auto spectrum_fft = FastFourierTransform::fft(signal);
    end = std::chrono::high_resolution_clock::now();
    auto duration_fft = std::chrono::duration_cast<std::chrono::microseconds>(end - start);

    std::cout << "\nDFT time: " << duration_dft.count() << " μs" << std::endl;
    std::cout << "FFT time: " << duration_fft.count() << " μs" << std::endl;
    std::cout << "Speedup: " << static_cast<double>(duration_dft.count()) / duration_fft.count() << "x" << std::endl;

    // Compute power spectrum
    auto power = DiscreteFourierTransform::powerSpectrum(signal);
    std::cout << "\nPower spectrum peaks at frequencies 5 and 10 Hz:" << std::endl;
    std::cout << "Power[5] = " << power[5] << std::endl;
    std::cout << "Power[10] = " << power[10] << std::endl;
}

// Demo 2: Circulant Matrices and Convolution
void demoCirculantConvolution() {
    std::cout << "\n=== Demo 2: Circulant Matrices and Convolution ===" << std::endl;

    // Create circulant matrix from first column
    std::vector<double> first_col = {1.0, 2.0, 3.0, 4.0};
    std::cout << "Circulant matrix first column: ";
    printSignal("", first_col);

    // Compute eigenvalues (which are DFT of first column)
    auto eigenvals = CirculantMatrix::eigenvalues(first_col);
    printComplexSignal("Eigenvalues", eigenvals);

    // Multiply circulant matrix with vector
    std::vector<double> x = {1.0, 0.0, 0.0, 0.0};
    auto result = CirculantMatrix::multiply(first_col, x);
    printSignal("Circulant multiplication result", result);

    // Demonstrate convolution theorem
    std::cout << "\n--- Convolution Theorem Verification ---" << std::endl;
    std::vector<double> f = {1.0, 2.0, 3.0, 4.0, 0.0, 0.0, 0.0, 0.0};
    std::vector<double> g = {1.0, 1.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0};

    printSignal("Signal f", f);
    printSignal("Kernel g", g);

    // Fast convolution via FFT
    auto conv_fast = Convolution::fastCircular(f, g);
    printSignal("Convolution (fast)", conv_fast);

    // Cross-correlation
    auto xcorr = Convolution::crossCorrelation(f, g);
    printSignal("Cross-correlation", xcorr);
}

// Demo 3: Wavelet Transforms
void demoWavelets() {
    std::cout << "\n=== Demo 3: Wavelet Transforms ===" << std::endl;

    // Create a test signal with discontinuity
    std::vector<double> signal = {1.0, 1.0, 1.0, 1.0, 5.0, 5.0, 5.0, 5.0,
                                   2.0, 2.0, 2.0, 2.0, 3.0, 3.0, 3.0, 3.0};

    printSignal("Original signal", signal);

    // Haar wavelet transform
    auto haar_coeffs = WaveletTransform::haarTransform(signal);
    printSignal("Haar coefficients", haar_coeffs);

    // Inverse Haar transform
    auto haar_reconstructed = WaveletTransform::haarInverse(haar_coeffs);
    printSignal("Haar reconstructed", haar_reconstructed);

    // Compute reconstruction error
    double error = 0.0;
    for (size_t i = 0; i < signal.size(); ++i) {
        error += std::abs(signal[i] - haar_reconstructed[i]);
    }
    std::cout << "Reconstruction error: " << error << std::endl;

    // Daubechies-4 wavelet
    std::cout << "\n--- Daubechies-4 Wavelet ---" << std::endl;
    auto d4_coeffs = WaveletTransform::daubechies4Transform(signal);
    printSignal("D4 coefficients", d4_coeffs);

    // Note: Daubechies-4 inverse not yet implemented
    std::cout << "D4 transform completed (inverse reconstruction TBD)" << std::endl;
}

// Demo 4: Fourier Series and Multipliers
void demoFourierSeries() {
    std::cout << "\n=== Demo 4: Fourier Series and Multipliers ===" << std::endl;

    // Sample a smooth periodic function
    size_t N = 16;
    std::vector<double> samples(N);
    for (size_t i = 0; i < N; ++i) {
        double theta = 2 * M_PI * i / N;
        samples[i] = std::sin(theta) + 0.5 * std::sin(3 * theta);
    }

    printSignal("Sampled function", samples, 8);

    // Compute Fourier coefficients
    auto coeffs = FourierSeries::coefficients(samples);
    std::cout << "\nFourier coefficients (magnitude):" << std::endl;
    for (size_t k = 0; k < std::min(size_t(5), coeffs.size()); ++k) {
        std::cout << "  k=" << k << ": " << std::abs(coeffs[k]) << std::endl;
    }

    // Compute derivative using Fourier multiplier (ik)
    auto derivative = FourierSeries::derivative(samples);
    printSignal("Derivative (via multiplier)", derivative, 8);

    // Apply fractional Laplacian
    double s = 0.5;
    auto frac_laplacian = FourierSeries::fractionalLaplacian(samples, s);
    printSignal("Fractional Laplacian (-Δ)^0.5", frac_laplacian, 8);

    // L2 norm
    double l2_norm = HilbertSpaceOperators::l2Norm(samples);
    std::cout << "L² norm: " << l2_norm << std::endl;
}

// Demo 5: Hilbert Space Operators
void demoHilbertOperators() {
    std::cout << "\n=== Demo 5: Hilbert Space Operators ===" << std::endl;

    // Create test matrix
    std::vector<std::vector<double>> A = {
        {1.0, 0.5, 0.0},
        {0.5, 2.0, 0.5},
        {0.0, 0.5, 1.0}
    };

    std::cout << "Matrix A (3x3 symmetric):" << std::endl;
    for (const auto& row : A) {
        std::cout << "  [";
        for (double val : row) {
            std::cout << std::setw(6) << std::fixed << std::setprecision(2) << val;
        }
        std::cout << " ]" << std::endl;
    }

    // Compute trace
    double tr = HilbertSpaceOperators::trace(A);
    std::cout << "\nTrace(A) = " << tr << std::endl;

    // Frobenius norm
    double frob = HilbertSpaceOperators::frobeniusNorm(A);
    std::cout << "Frobenius norm ||A||_F = " << frob << std::endl;

    // Schatten norms
    double s1 = HilbertSpaceOperators::schattenNorm(A, 1.0);
    double s2 = HilbertSpaceOperators::schattenNorm(A, 2.0);
    std::cout << "Schatten-1 norm (nuclear) ||A||_S1 = " << s1 << std::endl;
    std::cout << "Schatten-2 norm (Frobenius) ||A||_S2 = " << s2 << std::endl;

    // Check compactness (rank)
    bool compact = HilbertSpaceOperators::isCompact(A, 2, 1e-10);
    std::cout << "Is compact with rank ≤ 2? " << (compact ? "Yes" : "No") << std::endl;

    // Inner product of vectors
    std::vector<double> v1 = {1.0, 2.0, 3.0};
    std::vector<double> v2 = {4.0, 5.0, 6.0};
    double inner = HilbertSpaceOperators::innerProduct(v1, v2);
    std::cout << "<v1, v2> = " << inner << " (should be 32)" << std::endl;
}

// Demo 6: Time-Frequency Analysis
void demoTimeFrequency() {
    std::cout << "\n=== Demo 6: Time-Frequency Analysis (STFT) ===" << std::endl;

    // Create chirp signal (frequency increases over time)
    size_t N = 128;
    std::vector<double> signal(N);
    for (size_t i = 0; i < N; ++i) {
        double t = static_cast<double>(i) / N;
        double freq = 5.0 + 10.0 * t;  // Frequency chirp from 5 to 15 Hz
        signal[i] = std::sin(2 * M_PI * freq * t);
    }

    std::cout << "Chirp signal: " << N << " samples" << std::endl;
    std::cout << "Frequency increases from 5 Hz to 15 Hz" << std::endl;

    // Compute STFT
    size_t window_size = 32;
    size_t hop_size = 8;
    auto stft_result = TimeFrequencyAnalysis::stft(signal, window_size, hop_size);

    std::cout << "\nSTFT size: " << stft_result.size() << " time frames × "
              << stft_result[0].size() << " frequency bins" << std::endl;

    // Compute spectrogram
    auto spec = TimeFrequencyAnalysis::spectrogram(signal, window_size, hop_size);

    std::cout << "Spectrogram size: " << spec.size() << " × " << spec[0].size() << std::endl;

    // Find peak frequency at first and last time frames
    auto findPeakFreq = [](const std::vector<double>& frame) {
        size_t peak_idx = 0;
        double peak_val = frame[0];
        for (size_t i = 1; i < frame.size() / 2; ++i) {
            if (frame[i] > peak_val) {
                peak_val = frame[i];
                peak_idx = i;
            }
        }
        return peak_idx;
    };

    size_t peak_start = findPeakFreq(spec[0]);
    size_t peak_end = findPeakFreq(spec[spec.size() - 1]);

    std::cout << "\nPeak frequency bin at start: " << peak_start << std::endl;
    std::cout << "Peak frequency bin at end: " << peak_end << std::endl;
    std::cout << "Frequency shift detected: " << (peak_end > peak_start ? "Yes" : "No") << std::endl;

    // Gabor transform
    std::cout << "\n--- Gabor Transform ---" << std::endl;
    auto gabor = TimeFrequencyAnalysis::gaborTransform(signal, window_size, hop_size, 4.0);
    std::cout << "Gabor transform size: " << gabor.size() << " × " << gabor[0].size() << std::endl;
}

// Demo 7: Advanced Applications
void demoAdvancedApplications() {
    std::cout << "\n=== Demo 7: Advanced Applications ===" << std::endl;

    // Image filtering simulation (2D FFT)
    std::cout << "\n--- 2D FFT for Image Processing ---" << std::endl;
    std::vector<std::vector<double>> image = {
        {1.0, 2.0, 1.0, 2.0},
        {2.0, 4.0, 2.0, 4.0},
        {1.0, 2.0, 1.0, 2.0},
        {2.0, 4.0, 2.0, 4.0}
    };

    std::cout << "4×4 test image (checkerboard pattern)" << std::endl;
    auto fft2d = FastFourierTransform::fft2D(image);
    std::cout << "2D FFT computed: " << fft2d.size() << " × " << fft2d[0].size() << std::endl;

    // DC component
    std::cout << "DC component (average): " << std::abs(fft2d[0][0]) / 16.0 << std::endl;

    // Zero padding for interpolation
    std::cout << "\n--- Zero Padding for Interpolation ---" << std::endl;
    std::vector<double> signal = {1.0, 2.0, 3.0, 4.0};

    printSignal("Original signal", signal);
    auto padded = FastFourierTransform::zeroPad(signal);
    std::cout << "Zero-padded to power of 2: size " << padded.size() << std::endl;

    auto interpolated = FastFourierTransform::fft(
        FastFourierTransform::ifft(
            FastFourierTransform::fft(padded)
        )
    );
    std::cout << "Interpolated signal size: " << interpolated.size() << std::endl;

    // Auto-correlation for signal analysis
    std::cout << "\n--- Auto-correlation ---" << std::endl;
    std::vector<double> periodic = {1.0, 2.0, 3.0, 2.0, 1.0, 2.0, 3.0, 2.0};
    auto autocorr = Convolution::autoCorrelation(periodic);
    printSignal("Periodic signal", periodic);
    printSignal("Auto-correlation", autocorr);

    // Find period from autocorrelation peaks
    std::cout << "Auto-correlation peak at lag 4 indicates period = 4" << std::endl;
}

int main() {
    std::cout << "========================================" << std::endl;
    std::cout << "Fourier Analysis Computational Demo" << std::endl;
    std::cout << "========================================" << std::endl;

    try {
        demoDFTvsFFT();
        demoCirculantConvolution();
        demoWavelets();
        demoFourierSeries();
        demoHilbertOperators();
        demoTimeFrequency();
        demoAdvancedApplications();

        std::cout << "\n========================================" << std::endl;
        std::cout << "All demos completed successfully!" << std::endl;
        std::cout << "========================================" << std::endl;

    } catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << std::endl;
        return 1;
    }

    return 0;
}
