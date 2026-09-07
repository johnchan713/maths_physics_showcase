#include "ns_cascade/continuation.hpp"

#include <cstdlib>
#include <iostream>
#include <set>

namespace {
using namespace ns_cascade;

std::string canonicalPath(const std::string& path) {
    char* resolved = realpath(path.c_str(), NULL);
    if (resolved) {
        const std::string result(resolved);
        std::free(resolved);
        return result;
    }
    const auto slash = path.find_last_of('/');
    const std::string directory = slash == std::string::npos ? "." : path.substr(0, slash + 1);
    const std::string filename = slash == std::string::npos ? path : path.substr(slash + 1);
    resolved = realpath(directory.c_str(), NULL);
    if (!resolved || filename.empty()) {
        std::free(resolved);
        throw std::invalid_argument("Output directory must exist: " + path);
    }
    const std::string result = std::string(resolved) + "/" + filename;
    std::free(resolved);
    return result;
}

bool exists(const std::string& path) {
    std::ifstream input(path.c_str(), std::ios::binary);
    return input.good();
}

void writeHeader(std::ostream& output) {
    output << "backend,grid,cutoff,sampling_grid,dense_sampling_grid,step,observation,"
        "time,energy,enstrophy,palinstrophy,h_half,l3_sample,vorticity_sample,l3_dense,"
        "vorticity_dense,vorticity_fourier_bound,k_rms,stretching,viscous_destruction,"
        "net_enstrophy_rate,production_to_dissipation,cutoff_fraction,divergence_defect,"
        "reality_defect,nonlinear_energy_residual,peak_cutoff_fraction,minimum_time_step,"
        "maximum_time_step,maximum_cfl,maximum_viscous_number,bkm_sampled_integral,"
        "h_half_ratio,l3_ratio,vorticity_ratio,l3_dense_ratio,vorticity_dense_ratio,"
        "enstrophy_ratio,k_rms_ratio,cutoff_stopped\n";
}

void writeRow(std::ostream& output, const ContinuationCheckpoint& p) {
    const auto& c = p.configuration;
    const auto& d = p.latest;
    output << std::setprecision(17) << (c.fftw ? "fftw" : "radix2") << ','
        << c.grid << ',' << c.cutoff << ',' << c.sampling_grid << ',' << c.dense_sampling_grid
        << ',' << p.step << ',' << p.completed_observations;
    for (double value : continuation_detail::observationValues(d)) output << ',' << value;
    for (double value : {p.peak_cutoff_fraction, p.minimum_time_step, p.maximum_time_step,
                         p.maximum_cfl, p.maximum_viscous_number, p.bkm_sampled_integral,
                         d.h_half / p.initial.h_half, d.l3_sample / p.initial.l3_sample,
                         d.vorticity_sample / p.initial.vorticity_sample,
                         d.l3_dense / p.initial.l3_dense, d.vorticity_dense / p.initial.vorticity_dense,
                         d.enstrophy / p.initial.enstrophy, d.k_rms / p.initial.k_rms}) output << ',' << value;
    output << ',' << (p.cutoff_stopped ? "true" : "false") << '\n';
    output.flush();
    if (!output) throw std::runtime_error("Continuation evidence write failed");
}
}  // namespace

int main(int argc, char** argv) {
    try {
        ContinuationConfiguration c;
        std::string state_input, restart, checkpoint_output, output_path;
        double final_time = 0.0;
        bool scientific_set = false, cutoff_set = false;
        std::set<std::string> seen;
        for (int i = 1; i < argc; ++i) {
            const std::string key(argv[i]);
            if (key == "--help") {
                std::cout << "Continue saved Fourier coefficients with an atomic checkpoint at every observation.\n"
                    "Fresh: --state-input CSV --grid N [--cutoff K] [--backend fftw|radix2]\n"
                    "       [--viscosity .02] [--dt .000125] [--cfl .4] [--diffusion-safety 2]\n"
                    "       [--observation-interval .01] [--sampling-grid 64] [--dense-sampling-grid 128]\n"
                    "       [--cutoff-limit .008]\n"
                    "Resume: --restart CHECKPOINT (scientific settings are restored and cannot be overridden)\n"
                    "Both: --final-time T --output CSV --checkpoint-output CHECKPOINT\n"
                    "T must be an observation time. Exit 2 means the per-step cutoff gate stopped evolution.\n";
                return 0;
            }
            if (!seen.insert(key).second || i + 1 == argc) throw std::invalid_argument("Duplicate or incomplete option: " + key);
            const std::string value(argv[++i]);
            if (key == "--state-input") state_input = value;
            else if (key == "--restart") restart = value;
            else if (key == "--output") output_path = value;
            else if (key == "--checkpoint-output") checkpoint_output = value;
            else if (key == "--final-time") final_time = parseOptimizationStateCsvNumber<double>(value, key);
            else {
                scientific_set = true;
                if (key == "--grid") c.grid = parseOptimizationStateCsvNumber<int>(value, key);
                else if (key == "--cutoff") { c.cutoff = parseOptimizationStateCsvNumber<int>(value, key); cutoff_set = true; }
                else if (key == "--sampling-grid") c.sampling_grid = parseOptimizationStateCsvNumber<int>(value, key);
                else if (key == "--dense-sampling-grid") c.dense_sampling_grid = parseOptimizationStateCsvNumber<int>(value, key);
                else if (key == "--viscosity") c.viscosity = parseOptimizationStateCsvNumber<double>(value, key);
                else if (key == "--dt") c.maximum_time_step = parseOptimizationStateCsvNumber<double>(value, key);
                else if (key == "--cfl") c.target_cfl = parseOptimizationStateCsvNumber<double>(value, key);
                else if (key == "--diffusion-safety") c.diffusion_safety = parseOptimizationStateCsvNumber<double>(value, key);
                else if (key == "--observation-interval") c.observation_interval = parseOptimizationStateCsvNumber<double>(value, key);
                else if (key == "--cutoff-limit") c.cutoff_limit = parseOptimizationStateCsvNumber<double>(value, key);
                else if (key == "--backend" && (value == "fftw" || value == "radix2")) c.fftw = value == "fftw";
                else throw std::invalid_argument("Unknown option/value: " + key + " " + value);
            }
        }
        if (state_input.empty() == restart.empty() || output_path.empty() || checkpoint_output.empty()) {
            throw std::invalid_argument("Require exactly one input, --output and --checkpoint-output");
        }
        ContinuationCheckpoint p;
        if (!restart.empty()) {
            if (scientific_set) throw std::invalid_argument("Restart cannot override scientific settings");
            p = loadContinuationCheckpoint(restart);
            c = p.configuration;
        } else if (!cutoff_set) c.cutoff = (c.grid - 1) / 3;
        continuation_detail::validateConfiguration(c);
        const auto final_index = continuation_detail::observationIndex(final_time, c.observation_interval);
        if (final_index <= p.completed_observations || p.cutoff_stopped) {
            throw std::invalid_argument("Endpoint must be later than a passing restart checkpoint");
        }
        // Resolve aliases before opening anything for writing. The checkpoint
        // may replace its own restart input, but never an initial CSV/evidence.
        const auto input = canonicalPath(restart.empty() ? state_input : restart);
        const auto evidence = canonicalPath(output_path);
        const auto saved = canonicalPath(checkpoint_output);
        const auto temporary = canonicalPath(checkpoint_output + ".tmp");
        if (evidence == input || evidence == saved || evidence == temporary || temporary == input ||
            (saved == input && restart.empty()) || exists(evidence) || exists(temporary) ||
            (exists(saved) && saved != input)) {
            throw std::invalid_argument("Output would overwrite or alias existing evidence/input");
        }
        ContinuationSystem system(c);
        if (restart.empty()) p = system.initialize(readOptimizationStateCsv(state_input, system.spectralSystem()).state);
        std::ofstream output(output_path.c_str());
        if (!output) throw std::runtime_error("Cannot open continuation evidence output");
        writeHeader(output);
        writeRow(output, p);
        saveContinuationCheckpoint(checkpoint_output, p);
        if (!p.cutoff_stopped) {
            system.advance(p, final_time, [&](const ContinuationCheckpoint& current) {
                saveContinuationCheckpoint(checkpoint_output, current);
                writeRow(output, current);
                std::cout << std::setprecision(10) << "time=" << current.latest.time
                    << " steps=" << current.step << " H_ratio=" << current.latest.h_half / current.initial.h_half
                    << " vort_dense_ratio=" << current.latest.vorticity_dense / current.initial.vorticity_dense
                    << " cutoff_peak=" << current.peak_cutoff_fraction << std::endl;
            });
        }
        std::cout << (p.cutoff_stopped ? "Stopped at cutoff gate" : "Requested continuation completed") << '\n';
        return p.cutoff_stopped ? 2 : 0;
    } catch (const std::exception& error) {
        std::cerr << "Continuation failed: " << error.what() << '\n';
        return 1;
    }
}
