/**
 * @file greens_functions_demo.cpp
 * @brief Demonstration of Green's Functions for PDEs
 *
 * Topics covered:
 * - Green's functions for parabolic equations (heat/diffusion)
 * - Green's functions for elliptic equations (Laplace/Poisson)
 * - Green's functions for hyperbolic equations (wave)
 * - Method of images for boundary conditions
 * - Solution via convolution with source terms
 */

#include <iostream>
#include <iomanip>
#include <cmath>
#include <vector>
#include <functional>
#include "maths/advanced/pde/pde_classification_solutions.hpp"

using namespace maths::pde;

void print_section(const std::string& title) {
    std::cout << "\n" << std::string(80, '=') << "\n";
    std::cout << title << "\n";
    std::cout << std::string(80, '=') << "\n";
}

void print_subsection(const std::string& title) {
    std::cout << "\n--- " << title << " ---\n";
}

/**
 * Demonstrate heat equation Green's functions
 */
void demo_heat_green_functions() {
    print_section("PARABOLIC PDEs: HEAT EQUATION GREEN'S FUNCTIONS");

    print_subsection("1D Heat Kernel on Infinite Domain");

    double alpha = 0.1;  // Thermal diffusivity
    double x = 1.0, t = 1.0;
    double xi = 0.0, tau = 0.0;

    double G = GreensFunctions::HeatGreenFunction::infinite1D(x, t, xi, tau, alpha);

    std::cout << "Heat equation: u_t = α u_xx\n";
    std::cout << "Green's function: G(x,t;ξ,τ) = 1/√(4πα(t-τ)) exp(-(x-ξ)²/(4α(t-τ)))\n\n";
    std::cout << "Parameters: α = " << alpha << "\n";
    std::cout << "Source at (ξ,τ) = (" << xi << ", " << tau << ")\n";
    std::cout << "Observation at (x,t) = (" << x << ", " << t << ")\n";
    std::cout << "G(x,t;ξ,τ) = " << G << "\n";

    print_subsection("Heat Diffusion from Point Source");

    std::cout << "\nTemperature distribution at different times:\n";
    std::cout << std::setw(10) << "x"
              << std::setw(15) << "t=0.1"
              << std::setw(15) << "t=0.5"
              << std::setw(15) << "t=1.0"
              << std::setw(15) << "t=2.0\n";
    std::cout << std::string(70, '-') << "\n";

    for (double x_eval = -2.0; x_eval <= 2.0; x_eval += 0.5) {
        std::cout << std::setw(10) << x_eval;
        for (double t_eval : {0.1, 0.5, 1.0, 2.0}) {
            double G_val = GreensFunctions::HeatGreenFunction::infinite1D(
                x_eval, t_eval, 0.0, 0.0, alpha);
            std::cout << std::setw(15) << G_val;
        }
        std::cout << "\n";
    }

    print_subsection("Method of Images: Half-Space with Dirichlet BC");

    std::cout << "\nHalf-space x > 0 with G(0,t) = 0:\n";
    std::cout << "Using method of images: G = G_original - G_image\n\n";

    x = 1.0; t = 1.0;
    xi = 0.5; tau = 0.0;

    double G_dirichlet = GreensFunctions::HeatGreenFunction::halfSpaceDirichlet(
        x, t, xi, tau, alpha);

    std::cout << "Source at ξ = " << xi << ", image at ξ = " << -xi << "\n";
    std::cout << "G_Dirichlet(" << x << ", " << t << ") = " << G_dirichlet << "\n";

    print_subsection("2D Heat Equation");

    double x2d = 1.0, y2d = 1.0, t2d = 1.0;
    double xi2d = 0.0, eta2d = 0.0, tau2d = 0.0;

    double G_2d = GreensFunctions::HeatGreenFunction::infinite2D(
        x2d, y2d, t2d, xi2d, eta2d, tau2d, alpha);

    std::cout << "\n2D heat kernel:\n";
    std::cout << "G(x,y,t;ξ,η,τ) = 1/(4πα(t-τ)) exp(-r²/(4α(t-τ)))\n";
    std::cout << "r² = (x-ξ)² + (y-η)²\n\n";
    std::cout << "At (x,y,t) = (" << x2d << ", " << y2d << ", " << t2d << ")\n";
    std::cout << "Source at (ξ,η,τ) = (" << xi2d << ", " << eta2d << ", " << tau2d << ")\n";
    std::cout << "G = " << G_2d << "\n";

    print_subsection("Solving with Source Term");

    // Source term: f(x,t) = exp(-(x²)/2)
    auto source = [](double x, double t) {
        return std::exp(-x*x / 2.0);
    };

    double x_eval = 0.0, t_eval = 0.5;
    double u_solution = GreensFunctions::HeatGreenFunction::solveWithSource(
        source, x_eval, t_eval, alpha);

    std::cout << "\nSolving u_t = αu_xx + f(x,t) with f(x,t) = exp(-x²/2)\n";
    std::cout << "Solution u(" << x_eval << ", " << t_eval << ") = " << u_solution << "\n";
}

/**
 * Demonstrate Poisson equation Green's functions
 */
void demo_poisson_green_functions() {
    print_section("ELLIPTIC PDEs: LAPLACE/POISSON GREEN'S FUNCTIONS");

    print_subsection("2D Laplacian on Unbounded Domain");

    double x = 1.0, y = 1.0;
    double xi = 0.0, eta = 0.0;

    double G_2d = GreensFunctions::PoissonGreenFunction::infinite2D(x, y, xi, eta);

    std::cout << "Poisson equation: ∇²u = f\n";
    std::cout << "Green's function: ∇²G = δ(x-ξ)\n";
    std::cout << "G(x,y;ξ,η) = -(1/2π) ln(r) where r = √((x-ξ)² + (y-η)²)\n\n";
    std::cout << "Source at (ξ,η) = (" << xi << ", " << eta << ")\n";
    std::cout << "Observation at (x,y) = (" << x << ", " << y << ")\n";
    std::cout << "G(x,y;ξ,η) = " << G_2d << "\n";

    print_subsection("Radial Distribution");

    std::cout << "\nGreen's function vs distance:\n";
    std::cout << std::setw(15) << "Distance r"
              << std::setw(20) << "G (2D)"
              << std::setw(20) << "G (3D)\n";
    std::cout << std::string(55, '-') << "\n";

    for (double r = 0.5; r <= 3.0; r += 0.5) {
        double G2d = GreensFunctions::PoissonGreenFunction::infinite2D(r, 0, 0, 0);
        double G3d = GreensFunctions::PoissonGreenFunction::infinite3D(r, 0, 0, 0, 0, 0);
        std::cout << std::setw(15) << r
                  << std::setw(20) << G2d
                  << std::setw(20) << G3d << "\n";
    }

    print_subsection("3D Laplacian");

    double x3d = 1.0, y3d = 1.0, z3d = 1.0;
    double xi3d = 0.0, eta3d = 0.0, zeta3d = 0.0;

    double G_3d = GreensFunctions::PoissonGreenFunction::infinite3D(
        x3d, y3d, z3d, xi3d, eta3d, zeta3d);

    std::cout << "\n3D Green's function: G(r) = -1/(4πr)\n";
    std::cout << "At (x,y,z) = (" << x3d << ", " << y3d << ", " << z3d << ")\n";
    std::cout << "Source at origin\n";
    std::cout << "G = " << G_3d << "\n";

    print_subsection("Half-Space with Dirichlet BC");

    z3d = 1.0;
    zeta3d = 0.5;

    double G_halfspace = GreensFunctions::PoissonGreenFunction::halfSpace3DDirichlet(
        x3d, y3d, z3d, xi3d, eta3d, zeta3d);

    std::cout << "\nHalf-space z > 0 with G(x,y,0) = 0\n";
    std::cout << "Method of images: G = G_original - G_image\n";
    std::cout << "Source at ζ = " << zeta3d << ", image at ζ = " << -zeta3d << "\n";
    std::cout << "G(" << x3d << ", " << y3d << ", " << z3d << ") = " << G_halfspace << "\n";

    print_subsection("Rectangle Domain with Eigenfunction Expansion");

    double a = 1.0, b = 1.0;  // Rectangle [0,a] × [0,b]
    x = 0.5; y = 0.5;
    xi = 0.3; eta = 0.3;

    double G_rect = GreensFunctions::PoissonGreenFunction::rectangle2D(
        x, y, xi, eta, a, b, 10);

    std::cout << "\nRectangle [0," << a << "] × [0," << b << "] with Dirichlet BC\n";
    std::cout << "Eigenfunction expansion (10 terms)\n";
    std::cout << "G(" << x << ", " << y << "; " << xi << ", " << eta << ") = " << G_rect << "\n";

    print_subsection("Solving Poisson Equation on Rectangle");

    // Source: f(x,y) = sin(πx)sin(πy)
    auto source_2d = [](double x, double y) {
        return std::sin(M_PI * x) * std::sin(M_PI * y);
    };

    double u_poisson = GreensFunctions::PoissonGreenFunction::solvePoissonRectangle(
        source_2d, 0.5, 0.5, a, b, 10);

    std::cout << "\nSolving ∇²u = sin(πx)sin(πy) on [0,1]×[0,1]\n";
    std::cout << "u(0.5, 0.5) = " << u_poisson << "\n";
    std::cout << "Analytical: u = -1/(2π²) sin(πx)sin(πy) = "
              << -1.0/(2.0*M_PI*M_PI) * std::sin(M_PI*0.5) * std::sin(M_PI*0.5) << "\n";
}

/**
 * Demonstrate wave equation Green's functions
 */
void demo_wave_green_functions() {
    print_section("HYPERBOLIC PDEs: WAVE EQUATION GREEN'S FUNCTIONS");

    print_subsection("1D Wave Equation - Retarded Green's Function");

    double c = 1.0;  // Wave speed
    double x = 2.0, t = 3.0;
    double xi = 0.0, tau = 0.0;

    double G_1d = GreensFunctions::WaveGreenFunction::retarded1D(x, t, xi, tau, c);

    std::cout << "Wave equation: u_tt = c² u_xx\n";
    std::cout << "Green's function: G_tt - c²G_xx = δ(x-ξ)δ(t-τ)\n";
    std::cout << "G(x,t;ξ,τ) = 1/(2c) H(t-τ-|x-ξ|/c)\n\n";
    std::cout << "Wave speed c = " << c << "\n";
    std::cout << "Source at (ξ,τ) = (" << xi << ", " << tau << ")\n";
    std::cout << "Observation at (x,t) = (" << x << ", " << t << ")\n";
    std::cout << "G(x,t;ξ,τ) = " << G_1d << "\n";

    print_subsection("Causality and Light Cone");

    std::cout << "\nCausality check - Is (x,t) in past light cone of (ξ,τ)?\n";
    std::cout << std::setw(10) << "x"
              << std::setw(10) << "t"
              << std::setw(15) << "In Cone?"
              << std::setw(20) << "G(x,t;0,0)\n";
    std::cout << std::string(55, '-') << "\n";

    for (double x_test : {0.5, 1.0, 2.0, 3.0}) {
        double t_test = 2.0;
        bool in_cone = GreensFunctions::WaveGreenFunction::isInPastLightCone(
            x_test, t_test, 0.0, 0.0, c);
        double G_val = GreensFunctions::WaveGreenFunction::retarded1D(
            x_test, t_test, 0.0, 0.0, c);

        std::cout << std::setw(10) << x_test
                  << std::setw(10) << t_test
                  << std::setw(15) << (in_cone ? "YES" : "NO")
                  << std::setw(20) << G_val << "\n";
    }

    print_subsection("2D Wave Equation");

    double x2d = 1.0, y2d = 1.0, t2d = 2.0;
    double xi2d = 0.0, eta2d = 0.0, tau2d = 0.0;

    double G_2d = GreensFunctions::WaveGreenFunction::retarded2D(
        x2d, y2d, t2d, xi2d, eta2d, tau2d, c);

    std::cout << "\n2D wave Green's function (odd dimensions):\n";
    std::cout << "G(x,y,t;ξ,η,τ) = H(c(t-τ) - r) / (2π√(c²(t-τ)² - r²))\n";
    std::cout << "where r = √((x-ξ)² + (y-η)²)\n\n";
    std::cout << "At (x,y,t) = (" << x2d << ", " << y2d << ", " << t2d << ")\n";
    std::cout << "Source at (ξ,η,τ) = (" << xi2d << ", " << eta2d << ", " << tau2d << ")\n";
    std::cout << "G = " << G_2d << "\n";

    print_subsection("Wavefront Propagation");

    std::cout << "\nWave amplitude at different distances (t = 2.0):\n";
    std::cout << std::setw(15) << "Distance r"
              << std::setw(20) << "G (2D)\n";
    std::cout << std::string(35, '-') << "\n";

    t2d = 2.0;
    for (double r = 0.0; r <= 2.0; r += 0.25) {
        xi2d = r; eta2d = 0.0;
        double G_val = GreensFunctions::WaveGreenFunction::retarded2D(
            0.0, 0.0, t2d, xi2d, eta2d, 0.0, c);
        std::cout << std::setw(15) << r
                  << std::setw(20) << G_val << "\n";
    }

    print_subsection("Solving with Source Term - Duhamel's Principle");

    // Source: f(x,t) = exp(-x²) * sin(t)
    auto source_wave = [](double x, double t) {
        return std::exp(-x*x) * std::sin(t);
    };

    double x_eval = 0.0, t_eval = 2.0;
    double u_wave = GreensFunctions::WaveGreenFunction::solveWithSource1D(
        source_wave, x_eval, t_eval, c);

    std::cout << "\nSolving u_tt = c²u_xx + f(x,t) with f(x,t) = exp(-x²)sin(t)\n";
    std::cout << "Using Duhamel's principle: u = ∫∫ G(x,t;ξ,τ) f(ξ,τ) dξdτ\n";
    std::cout << "u(" << x_eval << ", " << t_eval << ") = " << u_wave << "\n";
}

/**
 * Demonstrate method of images
 */
void demo_method_of_images() {
    print_section("METHOD OF IMAGES FOR BOUNDARY CONDITIONS");

    print_subsection("Concept");

    std::cout << "Method of images satisfies boundary conditions by adding\n";
    std::cout << "image sources with appropriate signs:\n\n";
    std::cout << "• Dirichlet BC (u = 0): G = G_original - G_image (opposite signs)\n";
    std::cout << "• Neumann BC (∂u/∂n = 0): G = G_original + G_image (same signs)\n";

    print_subsection("1D Image Locations");

    double xi = 0.7;
    double boundary = 0.0;
    double image = GreensFunctions::MethodOfImages::imageLocation1D(xi, boundary);

    std::cout << "\nBoundary at x = " << boundary << "\n";
    std::cout << "Source at ξ = " << xi << "\n";
    std::cout << "Image at ξ' = " << image << "\n";

    print_subsection("2D Image Locations");

    xi = 0.5;
    double eta = 0.7;
    double boundary_x = 0.0;
    double boundary_y = 0.0;

    auto [image_x, image_y_x] = GreensFunctions::MethodOfImages::imageLocation2D_x(
        xi, eta, boundary_x);
    auto [image_x_y, image_y] = GreensFunctions::MethodOfImages::imageLocation2D_y(
        xi, eta, boundary_y);

    std::cout << "\nSource at (ξ,η) = (" << xi << ", " << eta << ")\n";
    std::cout << "Image across x = " << boundary_x << ": ("
              << image_x << ", " << image_y_x << ")\n";
    std::cout << "Image across y = " << boundary_y << ": ("
              << image_x_y << ", " << image_y << ")\n";
}

/**
 * Main demonstration
 */
int main() {
    std::cout << std::fixed << std::setprecision(8);

    std::cout << "\n";
    std::cout << "╔═══════════════════════════════════════════════════════════════════════════╗\n";
    std::cout << "║                   GREEN'S FUNCTIONS FOR PDEs                              ║\n";
    std::cout << "╚═══════════════════════════════════════════════════════════════════════════╝\n";

    try {
        demo_heat_green_functions();
        demo_poisson_green_functions();
        demo_wave_green_functions();
        demo_method_of_images();

        std::cout << "\n" << std::string(80, '=') << "\n";
        std::cout << "All demonstrations completed successfully!\n";
        std::cout << std::string(80, '=') << "\n\n";

    } catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << "\n";
        return 1;
    }

    return 0;
}
