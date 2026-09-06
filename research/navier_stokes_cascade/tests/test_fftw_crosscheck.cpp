#include "ns_cascade/pseudospectral.hpp"

#include <algorithm>
#include <cmath>
#include <complex>
#include <cstddef>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

extern "C" {
typedef double fftw_complex[2];
typedef struct fftw_plan_s* fftw_plan;
void* fftw_malloc(std::size_t size);
void fftw_free(void* pointer);
fftw_plan fftw_plan_dft_3d(int first_size,
                           int second_size,
                           int third_size,
                           fftw_complex* input,
                           fftw_complex* output,
                           int sign,
                           unsigned flags);
void fftw_execute(const fftw_plan plan);
void fftw_destroy_plan(fftw_plan plan);
}

namespace {

const int kFftwForward = -1;
const int kFftwBackward = 1;
const unsigned kFftwEstimate = 1U << 6;

void expect(bool condition, const std::string& message) {
    if (!condition) throw std::runtime_error(message);
}

class FftwBuffer {
public:
    explicit FftwBuffer(std::size_t size)
        : data_(static_cast<fftw_complex*>(
              fftw_malloc(sizeof(fftw_complex) * size))) {
        if (data_ == NULL) throw std::bad_alloc();
    }

    ~FftwBuffer() { fftw_free(data_); }

    fftw_complex* data() { return data_; }

private:
    FftwBuffer(const FftwBuffer&);
    FftwBuffer& operator=(const FftwBuffer&);
    fftw_complex* data_;
};

std::vector<ns_cascade::Complex> fftwTransform(
    const std::vector<ns_cascade::Complex>& values,
    int grid_size,
    int sign) {
    const std::size_t expected_size =
        static_cast<std::size_t>(grid_size) * grid_size * grid_size;
    if (values.size() != expected_size) {
        throw std::invalid_argument("FFTW reference input has the wrong size");
    }
    FftwBuffer input(values.size());
    FftwBuffer output(values.size());
    for (std::size_t i = 0; i < values.size(); ++i) {
        input.data()[i][0] = values[i].real();
        input.data()[i][1] = values[i].imag();
    }
    fftw_plan plan = fftw_plan_dft_3d(grid_size,
                                      grid_size,
                                      grid_size,
                                      input.data(),
                                      output.data(),
                                      sign,
                                      kFftwEstimate);
    if (plan == NULL) throw std::runtime_error("FFTW failed to create a 3D plan");
    fftw_execute(plan);
    fftw_destroy_plan(plan);

    const double scale = sign == kFftwForward
                             ? 1.0 / static_cast<double>(values.size())
                             : 1.0;
    std::vector<ns_cascade::Complex> transformed(values.size());
    for (std::size_t i = 0; i < values.size(); ++i) {
        transformed[i] = ns_cascade::Complex(output.data()[i][0],
                                             output.data()[i][1]) *
                         scale;
    }
    return transformed;
}

double maximumDifference(const std::vector<ns_cascade::Complex>& left,
                         const std::vector<ns_cascade::Complex>& right) {
    if (left.size() != right.size()) {
        throw std::invalid_argument("Cannot compare FFT arrays of different sizes");
    }
    double difference = 0.0;
    for (std::size_t i = 0; i < left.size(); ++i) {
        difference = std::max(difference, std::abs(left[i] - right[i]));
    }
    return difference;
}

bool retained(const ns_cascade::PseudospectralSystem& system,
              const ns_cascade::WaveVector& wave) {
    const int maximum_component = std::max(
        std::abs(wave.x), std::max(std::abs(wave.y), std::abs(wave.z)));
    return wave.normSquared() != 0 && maximum_component <= system.cutoff();
}

ns_cascade::PseudospectralSystem::State fftwRightHandSide(
    const ns_cascade::PseudospectralSystem& system,
    const ns_cascade::PseudospectralSystem::State& state) {
    typedef std::vector<ns_cascade::Complex> ScalarField;
    const std::size_t size = system.gridPointCount();
    ScalarField velocity_x(size);
    ScalarField velocity_y(size);
    ScalarField velocity_z(size);
    ScalarField vorticity_x(size);
    ScalarField vorticity_y(size);
    ScalarField vorticity_z(size);
    const ns_cascade::Complex imaginary_unit(0.0, 1.0);
    const std::vector<ns_cascade::WaveVector>& modes = system.gridModes();
    for (std::size_t i = 0; i < size; ++i) {
        if (!retained(system, modes[i])) continue;
        const ns_cascade::ComplexVector vorticity =
            ns_cascade::cross(modes[i], state[i]) * imaginary_unit;
        velocity_x[i] = state[i].x;
        velocity_y[i] = state[i].y;
        velocity_z[i] = state[i].z;
        vorticity_x[i] = vorticity.x;
        vorticity_y[i] = vorticity.y;
        vorticity_z[i] = vorticity.z;
    }

    velocity_x = fftwTransform(velocity_x, system.gridSize(), kFftwBackward);
    velocity_y = fftwTransform(velocity_y, system.gridSize(), kFftwBackward);
    velocity_z = fftwTransform(velocity_z, system.gridSize(), kFftwBackward);
    vorticity_x = fftwTransform(vorticity_x, system.gridSize(), kFftwBackward);
    vorticity_y = fftwTransform(vorticity_y, system.gridSize(), kFftwBackward);
    vorticity_z = fftwTransform(vorticity_z, system.gridSize(), kFftwBackward);

    ScalarField nonlinear_x(size);
    ScalarField nonlinear_y(size);
    ScalarField nonlinear_z(size);
    for (std::size_t i = 0; i < size; ++i) {
        nonlinear_x[i] = velocity_y[i] * vorticity_z[i] -
                         velocity_z[i] * vorticity_y[i];
        nonlinear_y[i] = velocity_z[i] * vorticity_x[i] -
                         velocity_x[i] * vorticity_z[i];
        nonlinear_z[i] = velocity_x[i] * vorticity_y[i] -
                         velocity_y[i] * vorticity_x[i];
    }
    nonlinear_x = fftwTransform(nonlinear_x, system.gridSize(), kFftwForward);
    nonlinear_y = fftwTransform(nonlinear_y, system.gridSize(), kFftwForward);
    nonlinear_z = fftwTransform(nonlinear_z, system.gridSize(), kFftwForward);

    ns_cascade::PseudospectralSystem::State result = system.zeroState();
    for (std::size_t i = 0; i < size; ++i) {
        if (!retained(system, modes[i])) continue;
        result[i] = ns_cascade::lerayProject(
            modes[i],
            ns_cascade::ComplexVector(
                nonlinear_x[i], nonlinear_y[i], nonlinear_z[i]));
        result[i] +=
            state[i] * (-system.viscosity() * modes[i].normSquared());
    }
    return result;
}

void testTransformsAgainstFftw() {
    const ns_cascade::PseudospectralSystem system(16, 0.05, 5);
    std::vector<ns_cascade::Complex> physical(system.gridPointCount());
    for (std::size_t i = 0; i < physical.size(); ++i) {
        const double coordinate = static_cast<double>(i);
        physical[i] = ns_cascade::Complex(
            std::sin(0.013 * coordinate) + 0.2 * std::cos(0.071 * coordinate),
            0.3 * std::sin(0.037 * coordinate));
    }
    const std::vector<ns_cascade::Complex> internal_forward =
        system.forwardTransform(physical);
    const std::vector<ns_cascade::Complex> fftw_forward =
        fftwTransform(physical, system.gridSize(), kFftwForward);
    expect(maximumDifference(internal_forward, fftw_forward) < 2e-14,
           "In-repository forward FFT differs from FFTW");

    const std::vector<ns_cascade::Complex> internal_inverse =
        system.inverseTransform(internal_forward);
    const std::vector<ns_cascade::Complex> fftw_inverse =
        fftwTransform(internal_forward, system.gridSize(), kFftwBackward);
    expect(maximumDifference(internal_inverse, fftw_inverse) < 2e-13,
           "In-repository inverse FFT differs from FFTW");
}

void testNavierStokesRightHandSideAgainstFftw() {
    const ns_cascade::PseudospectralSystem system(16, 0.02, 5);
    const ns_cascade::PseudospectralSystem::State state =
        system.vortexTubePairState(
            ns_cascade::VortexTubeParameters(0.7, 1.2, 0.3, 2), 10.0);
    const ns_cascade::PseudospectralSystem::State internal =
        system.rightHandSide(state);
    const ns_cascade::PseudospectralSystem::State reference =
        fftwRightHandSide(system, state);
    double difference = 0.0;
    for (std::size_t i = 0; i < internal.size(); ++i) {
        difference = std::max(
            difference, ns_cascade::norm(internal[i] - reference[i]));
    }
    expect(difference < 2e-12,
           "Navier-Stokes right-hand side differs from independent FFTW path");
}

}  // namespace

int main() {
    try {
        testTransformsAgainstFftw();
        testNavierStokesRightHandSideAgainstFftw();
        std::cout << "Independent FFTW cross-check passed.\n";
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "test failure: " << error.what() << '\n';
        return 1;
    }
}
