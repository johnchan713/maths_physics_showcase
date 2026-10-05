"""Exact outward arithmetic for the conditional annular operator constants."""
from fractions import Fraction as F


def ledger():
    # Axial second derivative bound: Hessian, f_X, f_eta, f terms.
    blocks = [[9, 0, 0], [10, 8, 0], [3, 4, 0], [0, 6, 4]]
    summed = [sum(row[j] for row in blocks) for j in range(3)]
    proposed = [32, 64, 32]
    difference = [a-b for a, b in zip(proposed, summed)]
    # Separate terms in Omega0, Gtheta and Gz, after replacing K^j by K^4.
    omega = F(2)+F(3, 2)+F(3)+F(2)
    theta = F(7, 2)+F(2)+F(6)+F(6)
    axial = F(3)+F(1)+F(5)+F(7)+F(4)
    return dict(conditional_on_final_profile_K=True, actual_K_evaluated=False,
                required_inputs=['F,U,Pi,V ordinary mixed C2 jets', 'X and X^-1', 'L^-1', 'nu'],
                axial_second_polynomial=summed, proposed_polynomial=proposed,
                nonnegative_polynomial_difference=difference,
                omega_coefficient=str(omega), theta_coefficient=str(theta), axial_coefficient=str(axial),
                proposed_inertial_coefficients=[10, 20, 20],
                viscous_coefficients=[32, 200, 128],
                physical_component_coefficients=[42, 440, 148],
                checks=dict(second_operator_dominated=all(v >= 0 for v in difference),
                            omega_outward=omega <= 10, theta_outward=theta <= 20, axial_outward=axial <= 20,
                            radial_component=42 == 10+32,
                            theta_component=440 == 2*(20+200),
                            axial_component=148 == 20+128))
