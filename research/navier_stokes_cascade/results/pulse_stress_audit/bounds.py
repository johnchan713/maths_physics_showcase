"""Outward constants for the corrected pulse, not a global PDE certificate.

README.md derives every estimate over the entire radial/angular rectangle.
The universal smallness checks use T=64 and monotonicity, not a sampled
profile or a floating-point representation of the actual huge parameters.
"""
import importlib.util
from pathlib import Path
import mpmath as mp

HERE = Path(__file__).resolve().parent
PULSE = HERE.parent/'pulse_moment_audit'
AXIAL = HERE.parent/'axial_stress_audit'
spec = importlib.util.spec_from_file_location('stress_interval_helpers', AXIAL/'bounds.py')
helpers = importlib.util.module_from_spec(spec)
spec.loader.exec_module(helpers)
point, lo, hi = helpers.point, helpers.lo, helpers.hi


def pulse_bound(md=64, digits=60, radius=100):
    """Finite corrected-pulse cone; earlier axial stage still needs Md=64.

    Products containing tiny positive lambda stay separated. In particular,
    v-2 >= 2 lambda is never established by subtracting rounded v and two.
    """
    if not isinstance(md, int) or md < 4 or digits < 50 or radius < 100:
        raise ValueError('Require integer Md>=4, >=50 interval digits and XR>=100')
    mp.iv.dps = digits
    p, exp = point, mp.iv.exp
    T = p(64)
    lam = exp(-4*T)
    h = lam**2
    eb = exp(-p('119.5')*T+p('.3'))
    memory = 24*exp(-p('120.5')*T+p('.7'))
    beta = p('.5')-lam
    amp = p('1.0100504')
    # The C1 end coefficients are <=2 lambda^20 on disjoint width-.3 bumps.
    end_R = 128*lam**20
    end_Ry = p(40000)/p('.3')*lam**20
    Reta = 33000*T*lam+end_R
    q_error = 100*eb
    # Positive finite remainder in w=2R-Cd R_xi+error, R=A R0 >=0.
    w_parts = dict(
        radial_convolution=8000000*lam,
        incoming_memory=4*memory/lam,
        end_corrections=1600*lam**19,
        angular_derivative=51000*T*mp.iv.sqrt(lam),
        energy_pressure_Q=10000*exp(26)*eb/lam**2,
    )
    w_error = sum(w_parts.values(), p(0))
    bs_error = 1800*lam
    # Cd is decreasing as c_eta increases. At zero angle, h=lambda^2
    # cancels lambda*(1-lambda), leaving D/beta^2 exactly.
    Cd = (p('.5')-h)/beta**2
    Cd_claim = p('2.001')
    linear = Cd_claim*amp
    wmax = p(1800)
    cross_error = 14*w_error+bs_error*wmax
    bsw = linear**2/8+cross_error
    second_error = 2*cross_error+14*bs_error+bs_error**2/2+2*lam*wmax**2
    second = 2*linear**2/7+second_error
    # Coarse compact ranges used to turn ratios into the actual finite cone.
    cmin, vmax, cvmax, pmin = p('.74'), p(101), p(1300000), p(20000000)
    finite_gap = p('.82')-4*cvmax/pmin
    log_p_floor = mp.iv.log(radius)+237*T+2-mp.iv.log(4)
    checks = {
        'T_floor': lo(exp(md)+10) >= 64,
        'lambda_below_one_millionth': hi(lam)<p('0.000001'),
        'beta_above_0p49': lo(beta)>p('.49'),
        'step_bounds_inherited': 256 <= 4*64 and 256**2+3072 < 4*20000,
        'main_R_below_14': hi(p('1.2')*11+end_R)<14,
        'main_first_derivative_below_846': 6*(1+11*64)<=5*846,
        'main_second_derivative_below_270000': p('1.2')*(3200+128+11*20000)<270000,
        'R_eta_below_one': hi(Reta)<1,
        'm_below_30': hi(14/beta+memory)<30,
        'm_eta_below_one': hi(Reta/beta+memory)<1,
        'power_memory_growth_below_two': hi(lam*240*T)<lo(mp.iv.log(2)),
        'incoming_Q_memory_below_eb': hi(32*exp(-p('120.5')*T-p('.3')))<1,
        'angular_source_error_below_32E': hi(61*lam+(1+2*h)*14)<32,
        'Q_error_below_100eb': hi(1+32/(1-lam))<100,
        'q0_above_angle_floor': lo(1-lam)>p(1)/3 and lo((1-2*h)/2)>p(1)/3,
        'Q_above_quarter_angle_floor': hi(q_error/lam)<p(1)/12,
        'Tlambda_below_one': hi(T*lam)<1,
        'S_below_600e26_over_lambda': hi(500*T*lam+197/2)<600,
        'S_eta_below_15e26_over_lambda': hi(lam+14)<15,
        'N_remainder_below_2000e26E_over_lambda': hi((61*14+3)*lam/exp(26)+4*h*600+15+1200)<2000,
        'Q_division_error_budget': hi(8000+54000/exp(26))<10000,
        'main_convolution_remainder': hi(270000/beta**3)<2400000,
        'end_convolution_below_300lambda20': hi(128/beta)<300,
        'angular_division_constant': hi(p(3)/4*33000/beta)<51000,
        'combined_radial_constant': 3*2400000+42<8000000,
        'combined_end_constant': 900+384+225<1600,
        'all_smallness_envelopes_decrease': 1-2*64<0 and 1-4*64<0 and -p('111.5')<0,
        'w_remainder_below_claim': 0<lo(w_error) and hi(w_error)<p('1e-48'),
        'Cd_below_2p001': hi(Cd)<lo(Cd_claim),
        'w_compact_range': hi(28+Cd_claim*846+w_error)<wmax,
        'bs_remainder_below_1800lambda': hi(1720+(end_R*(1+2*lam)+2*end_Ry)/lam)<1800,
        'first_ratio_claim': hi(bsw)<p('.52'),
        'second_ratio_claim': hi(second)<p('1.18'),
        'v_compact_range': hi(2+2*lam+(14+bs_error)**2/2)<vmax,
        'cv_compact_range': hi((1+(14+bs_error)*wmax/2)*vmax)<cvmax,
        'finite_ps1_floor': lo(log_p_floor)>hi(mp.iv.log(pmin)),
        'finite_Pc_above_v': lo(pmin*cmin-vmax)>0,
        'finite_quadratic_cone_gap': lo(finite_gap)>0,
        'strict_v_excess_retained': lo(lam)>0,
    }
    actual_T = exp(md)+10
    return dict(status='corrected-reference-pulse-cone-bounded' if all(checks.values()) else 'bound-inconclusive',
        Md=md, interval_digits=digits, XR_sufficient_lower=radius,
        checks=checks, T_interval=[lo(actual_T),hi(actual_T)],
        log_lambda_interval=[lo(-4*actual_T),hi(-4*actual_T)],
        log_h_interval=[lo(-8*actual_T),hi(-8*actual_T)],
        log_v_minus_two_lower_interval=[lo(mp.iv.log(2)-4*actual_T),hi(mp.iv.log(2)-4*actual_T)],
        Q_lower='(lambda+eta^2)/4', Q_absolute_error_bound='100 exp(-119.5T+.3)',
        Q_relative_error_universal_upper=hi(3*q_error/lam),
        w_remainder_parts={k:hi(v) for k,v in w_parts.items()},
        w_remainder_universal_upper=hi(w_error), w_remainder_claim='1e-48',
        Cd_universal_upper=hi(Cd), Cd_claim='2.001',
        bs_remainder_universal_upper=hi(bs_error),
        maximum_bsw_upper=hi(bsw), maximum_second_expression_upper=hi(second),
        first_A24_margin_lower=lo(2-bsw), second_A24_margin_lower=lo(2-second),
        first_A24_claim='1.48', second_A24_claim='0.82',
        w_absolute_upper=1800, bs_absolute_upper='14+1800lambda',
        finite_ps1_sufficient_lower=20000000, finite_Pc_over_ps1_lower='0.74',
        finite_normalized_quadratic_gap_lower=lo(finite_gap),
        universal_log_ps1_lower=lo(log_p_floor),
        scope='Only the corrected reference pulse; full outer profile, heat, attachment, PDE and smooth forcing are not certified.')
