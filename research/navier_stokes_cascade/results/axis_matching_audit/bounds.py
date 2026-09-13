"""Outward constants for the conditional five-moment annulus and axis data.

These are bounds on exact integrals and on an exact small root. The numerical
quadratures in moments.py are diagnostics, not replacements for these bounds.
"""
import importlib.util
from pathlib import Path
import mpmath as mp

HERE = Path(__file__).resolve().parent
RESULTS = HERE.parent
spec = importlib.util.spec_from_file_location('matching_interval_helpers',
    RESULTS/'axial_stress_audit'/'bounds.py')
helpers = importlib.util.module_from_spec(spec)
spec.loader.exec_module(helpers)
point,lo,hi = helpers.point,helpers.lo,helpers.hi


def inverse_bounds(rates, first, spacing):
    """Lagrange coefficients bound the inverse of an exponential bump matrix.

    B[k,j]=mass[k]*node[k]**j. Positivity and exact bump mass .08
    bound mass[k] below without numerical integration of the flat step.
    """
    p,exp = point,mp.iv.exp
    rates = list(map(p,rates))
    nodes = [exp(p(spacing)*r) for r in rates]
    rows = [p(0) for _ in rates]
    for k,r in enumerate(rates):
        other = [nodes[i] for i in range(len(rates)) if i!=k]
        denominator = p('.08')*exp((p(first)-p('.04'))*r)
        for node in other:
            denominator *= abs(nodes[k]-node)
        numerators = ([other[0],p(1)] if len(other)==1 else
                      [other[0]*other[1],other[0]+other[1],p(1)])
        for j,numerator in enumerate(numerators):
            rows[j] += numerator/denominator
    return [hi(v) for v in rows]


def certificate(digits=90):
    if not isinstance(digits,int) or digits<50:
        raise ValueError('At least 50 interval digits required')
    mp.iv.dps = digits
    p,exp = point,mp.iv.exp
    eps, beta, quadratic = p('1e-16'),p(1000),p(100)
    target = 3*eps
    coefficient = 2*beta*target
    ledger = 4*coefficient
    U_inverse = inverse_bounds(['1','1.6'],'.3','.4')
    E_inverse = inverse_bounds(['1.6','1.2','.2'],'.15','.35')
    zeta_C1 = 12*exp(p('1.2'))/16**2
    quadratic_bound = 32*p('.24')*exp(p('1.6'))+64*p('.16')*exp(1)*zeta_C1
    linear_bound = p('.24')*exp(p('1.6'))
    # Elementary derivative bounds: central interval and flat endpoint tails.
    g = 2/p('.25')**3+2/p('.75')**3
    gp = 6/p('.25')**4+6/p('.75')**4
    step_second_middle = (g*g+gp)/4
    step_second_edge = p(22368)*exp(-12)
    # Unedited pressure is a positive mixture of (1+eta^2)^(-2 theta).
    # Every later pressure-preserving correction keeps this axis datum.
    pressure_mass = (5+exp(p('.2'))+4*exp(p('-.4')))/2
    h = p('.01')
    q_ideal = p(9)/8-p(9)*h/8
    e_ideal = exp(p('-.8'))/2
    n_ideal = (20+32*(1+h))/16**2+4*(p('.5')+h)*5+10 \
        +(5+p(25)*h/3+p(25)/3)*exp(-1)
    profile = 64*coefficient
    slope = p(250000)*coefficient
    q_error = p('1e5')*coefficient
    n_error = p('1e4')*coefficient
    a_error = p('1e6')*coefficient
    W_error = 2*exp(2)*ledger
    Q_error_derived = W_error+78*exp(p('3.2'))*ledger/(1-profile)+10*profile/(1-profile)
    N_error_derived = (3*profile+W_error*(4+profile)+exp(2)*ledger)/16**2 \
        +p('1.04')*(16*exp(2)/16**2+4*exp(p('.8')))*ledger \
        +p('3.04')*4*exp(p('-1.2'))*ledger
    # P*b_s and w/P are used together: the large pressure scale cancels.
    b_scaled = p('3e6')*coefficient
    w_scaled = p(160)
    G = 1-b_scaled*32/(p('.7')*p('.2'))
    vs = p('.9')+b_scaled*b_scaled/(16**2*p('.7'))
    Pc = 10000*exp(-8)*p('.99')
    # Quantitative axis separation (B.2), using the same pressure envelope.
    small_Z_eta_ratio = (p(1)/100+5/p(256))/(p('1.25')-40/p(256))
    H_ratio = 1-p('4.5')/30-4*p('.05')**2/30**3-p('.05')**2/30**2
    chi_floor = p('.8')**2/(p('.8')**2+p('.01')**2)
    j = p('1e-18')
    sigma = p('1e-20')
    complex_radius = sigma/(4*(p('.5')+4+48+4*p('.05')))
    zeta_complex = 4*(1+8*h)*(1+5*(8+p('.05')))/sigma**2
    checks = dict(
        U_inverse_below_100=max(U_inverse)<100,
        E_inverse_below_1000=max(E_inverse)<1000,
        quadratic_below_100=hi(quadratic_bound)<100,
        linear_below_two=hi(linear_bound)<2,
        zeta_C1_below_one=hi(zeta_C1)<1,
        second_step_derivative_below_20000=max(hi(step_second_middle),hi(step_second_edge))<20000,
        exact_root_contraction=hi(8*beta**2*quadratic*target)<1,
        cumulative_ledger=hi(target+2*coefficient+100*coefficient**2)<lo(ledger),
        pressure_value_below_five=hi(pressure_mass)<5,
        pressure_derivative_below_ten=hi(2*pressure_mass)<10,
        ideal_Q_above_1p1=lo(q_ideal)>p('1.1'),
        ideal_N_over_P2_below_30=hi(n_ideal)<30,
        moment_formula_Q_budget=hi(Q_error_derived)<lo(q_error),
        moment_formula_N_budget=hi(N_error_derived)<lo(n_error),
        logarithmic_slope_budget=hi(2*slope/(1-profile))<lo(a_error),
        positive_E=lo(e_ideal*(1-profile))>p('.2'),
        positive_Q=lo(q_ideal-q_error)>1,
        N_over_P2_below_32=hi(n_ideal+n_error)<32,
        shear_a_between_p7_p9=hi(a_error)<p('.1'),
        scaled_shear_bound=hi(2*slope/p('.2'))<lo(b_scaled),
        relaxed_vs_below_one=hi(vs)<1,
        relaxed_G_above_p99=lo(G)>p('.99'),
        finite_Pc_above_two=lo(Pc)>2,
        old_XR100_fails_Pc_test=hi(100*exp(-8)*p(9)/8)<2,
        axis_j_in_matching_budget=hi(j)<lo(eps/12),
        small_Z_localizes_eta=hi(small_Z_eta_ratio)<p(1)/30,
        small_Z_keeps_H_away_from_zero=lo(H_ratio)>p('.8'),
        axis_complement_chi_above_p99=lo(chi_floor)>p('.99'),
        axis_Hzero_source_positive=lo(p('.25')-p(1)/256)>p('.2'),
        axis_negative_W_above_2p8=lo(3-8*h-p('.05'))>p('2.8'),
        common_complex_rectangle=0<lo(complex_radius) and hi(complex_radius)<p('1e-20'),
        complex_L_away_from_zero=lo(1-8*h)>p('.9'),
    )
    return dict(status='conditional-five-moment-annulus-bounded' if all(checks.values()) else 'inconclusive',
        interval_digits=digits,checks=checks,entry_C1_bound=hi(eps),target_C1_bound=hi(target),
        inverse_bound=1000,U_inverse_row_upper=U_inverse,E_inverse_row_upper=E_inverse,
        quadratic_bound=100,coefficient_C1_bound=hi(coefficient),ledger_C1_bound=hi(ledger),
        contraction_upper=hi(8*beta**2*quadratic*target),pressure_mass_upper=hi(pressure_mass),
        Q_error_upper=hi(q_error),N_over_P2_error_upper=hi(n_error),a_error_upper=hi(a_error),
        E_over_P_lower='.2',Q_lower=1,N_over_P2_upper=32,P_times_bs_upper=hi(b_scaled),
        w_over_P_upper=160,G_lower='.99',vs_upper=hi(vs),XR_min=10000,Pc_lower=lo(Pc),
        old_XR100_Pc_upper=hi(100*exp(-8)*p(9)/8),
        axis_j='1e-18',axis_sigma='1e-20',axis_delta_star='j P^2/100',
        axis_small_Z_eta_over_j_upper=hi(small_Z_eta_ratio),axis_chi_lower=lo(chi_floor),
        axis_complex_radius_lower=lo(complex_radius),
        log_C_complex_floor_over_Lambda=hi(2*zeta_complex+1),
        native_axial_C1_budget=lo(eps/12),activation_axial_C1_budget=lo(eps/12),
        core_ledger_C1_budget=lo(eps/4),
        actual_axis_entry_verified=False,full_admissible_cone_realized=False,blowup_verified=False)
