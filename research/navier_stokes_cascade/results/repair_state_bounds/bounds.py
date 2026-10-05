"""Check the constants derived in README without materializing G, H or N."""
from fractions import Fraction
import importlib.util
from pathlib import Path
import sys
import mpmath as mp

HERE = Path(__file__).resolve().parent
PROJECT = HERE.parent.parent
AMIN = 2**16


def load(name,path):
    spec = importlib.util.spec_from_file_location(name,path)
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


compact = load('repair_compact_bounds',HERE.parent/'compact_jet_envelope'/'bounds.py')
previous_paths=list(sys.path)
try:
    sys.path.insert(0,str(HERE.parent/'stress_realization_audit'))
    stress = load('repair_previous_stress_bounds',HERE.parent/'stress_realization_audit'/'bounds.py')
    repair = load('repair_previous_map',HERE.parent/'stress_realization_audit'/'repair.py')
finally:
    sys.path[:]=previous_paths
point,lower,upper = compact.p,compact.lower,compact.upper


def radius(A):
    """Exact rational radius; binary floats must not erase smallness."""
    if type(A) is not int or A < AMIN:
        raise ValueError('The analytic envelope requires an integer A>=65536 here')
    return Fraction(1,A**128)


def repair_budget(A,s):
    """A concrete bound on the coefficient neighborhood, not a field solver."""
    r = radius(A)
    s = Fraction(s)
    if not 0 <= s <= r:
        raise ValueError('Coefficient norm outside the proved repair neighborhood')
    d,q = 1280*A*s,4*10**6*A*s
    return dict(coefficient_norm=s,radius=r,field_C1=d,radial_field_C1=q,
                relative_swirl_error=1280*A*A*s,
                moment_C1=8*A**3*d,
                shear_error=4*A*q+4*A**3*d,
                state_error=327680*A**9*s)


def scalar_checks(digits=80):
    mp.mp.dps = mp.iv.dps = digits
    a = point(AMIN)
    r = a**(-128)
    d = 1280*a*r
    shear_ratio = 16*10**6/a**4+5120/a**2
    step = compact.c2.matching.certificate(digits)
    previous = stress.repair_certificate(digits)
    checks = dict(
        exact_bump_value_bound=previous['checks']['step_derivative_below_64'],
        exact_bump_radial_bound=step['checks']['second_step_derivative_below_20000'],
        exact_five_row_inverse=previous['checks']['U_inverse_below_2'] and previous['checks']['E_inverse_below_1000'],
        exact_C1_quadratic_map=previous['checks']['quadratic_C1_below_100000'],
        inherited_loop_cap_and_delta_comparisons=all(compact.previous.finite_integer_comparison().values()),
        positive_swirl_neighborhood=upper(1280*a*a*r) < mp.mpf('.5'),
        field_perturbation_below_one=upper(d) < 1,
        shear_error_below_A6=upper(shear_ratio) < 1,
        pressure_state_prefactor_below_A11=327680 < AMIN**2,
        chosen_Ccorr_dominates_raw_bound=11 < 128,
        active_third_gap_above_half_inverse_A=Fraction(45,64) > Fraction(1,2),
        active_fourth_gap_above_Gminus4=225*AMIN**2 > 256,
        epsilon_below_first_coordinate_half=AMIN**61 > 2,
        epsilon_below_cone_Lipschitz_tolerance=AMIN**30 > 512*4**10,
        state_threshold_log_constant=upper(mp.iv.log(4002)) < AMIN,
        coefficient_threshold_log_constant=upper(mp.iv.log(4000)) < AMIN,
        contraction_threshold_log_constant=upper(mp.iv.log(16*10**11)) < AMIN,
        incoming_memory_log_constant=upper(mp.iv.log(512)) < AMIN,
        state_frequency_comparison=193 < 16*AMIN**3,
        contraction_frequency_comparison=1 < 16*AMIN**3,
        radius_frequency_comparison=129 < 16*AMIN**3,
        incoming_memory_frequency_comparison=7 < 31*AMIN**3,
    )
    return dict(interval_digits=digits,checks=checks,
                relative_swirl_error_at_endpoint=1280*a*a*r,
                shear_to_A6_ratio=shear_ratio,
                exact_repair_inverse_bounds={k:previous[k] for k in
                    ('U_inverse_upper','E_inverse_upper','quadratic_C1_upper')},
                huge_frequency_materialized=False)


def ledger():
    """The arithmetic following the independently derived inequalities."""
    return dict(coefficient_norm='max_j max(sup|c_j|,sup|c_j_eta|)',
                field_angular_norm='sup|f|+sup|f_eta|',
                field_error='1280*A*s',radial_field_error='4000000*A*s',
                moment_error='10240*A^4*s',
                p1_coefficient=10*8+32,p1_power=7,
                p2_coefficient=24*8+8,p2_power=8,
                state_coefficient=256*1280,state_power=9,
                Ccorr='A^128',r='A^(-128)',epsilon='G^(-64)',
                Cstate='H^16',D='H^16',N='1+floor(H^32)',
                pressure_state_derivative_order=0,moment_angular_order=1,
                inherited_input_review_completed=False,
                checks=dict(p1_below_common=10*8+32 <= 256,
                            p2_below_common=24*8+8 <= 256,
                            coefficient_norm_conversion=1280 == 2*640,
                            radial_width_loss=4*10**6 == 2*20000*100))
