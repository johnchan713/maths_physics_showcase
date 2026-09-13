"""Outward constants and conditional frequency acceptance, not sampled minima."""
import mpmath as mp
from loop import point, lower, upper, variance_iv


def fixture_cap_certificate(digits=80, cells=128):
    """Cover every p in [-1,1], not merely the seven numerical fixtures.

    Evenness reduces to [0,1]; dyadic interval endpoints are exact. This is
    manufactured a=.8, d0=1, mu_cap=64 data, NOT the actual joined profile.
    """
    if digits < 50 or cells != 128:
        raise ValueError('This certificate uses 128 dyadic cells and >=50 digits')
    mp.iv.dps=digits
    floors=[]
    for j in range(cells):
        p=mp.iv.mpf([mp.mpf(j)/cells,mp.mpf(j+1)/cells])
        floors.append(lower(variance_iv(64,p,1)))
    return dict(digits=digits,cells=cells,pressure_range=['-1','1'],
                mu_cap=64,minimum_variance_lower=min(floors),
                required_variance_upper='3.75',
                full_pressure_range_covered=min(floors)>mp.mpf('3.75'),
                manufactured_input=True)


def repair_certificate(digits=80):
    if not isinstance(digits,int) or digits < 50:
        raise ValueError('At least 50 interval digits required')
    mp.iv.dps = digits
    p, exp = point, mp.iv.exp
    # Unit-mass positive bumps, width .1, centers U=1,2 and E=.5,1.5,2.5.
    du = p('.9')*exp(p('.95')+p('1.95')-p('.0205'))
    bu = p('3.05')*exp(p('2.05'))/du
    de = exp(-p('.51')*p('4.65')) \
        *(exp(p('1.45'))-exp(p('.55'))) \
        *(exp(p('2.45'))-exp(p('.55'))) \
        *(exp(p('2.45'))-exp(p('1.55')))
    be = 6*exp(3*p('2.55'))/de
    quadratic_c1 = 2*640*(2*exp(p('2.05'))+p('1.5')*exp(p('2.55')))
    middle_step = (2/p('.25')**3+2/p('.75')**3)/4
    edge_step = 192*exp(p(16)/9-16)
    checks = dict(U_inverse_below_2=upper(bu)<2,
                  E_inverse_below_1000=upper(be)<1000,
                  quadratic_C1_below_100000=upper(quadratic_c1)<100000,
                  step_derivative_below_64=max(upper(middle_step),upper(edge_step))<64,
                  Bessel_ratio_lower_constant=lower(8/(exp(1)*(p(22)/7)**2))>p('.25'),
                  excursion_small_z_constant=upper(2*exp(2))<16,
                  excursion_large_z_constant=upper((p(22)/7)*exp(p('.5'))+1)<16)
    return dict(digits=digits, U_determinant_lower=du, E_determinant_lower=de,
                U_inverse_upper=bu, E_inverse_upper=be,
                quadratic_C1_upper=quadratic_c1, checks=checks,
                lambda_domain=['0 (excluded)','.01'], inverse_bound=1000,
                quadratic_C1_bound=100000,
                original_discrepancy_lambda_loss_retained=True)


def uniform_mu_log_bound(a_min, d0, p_max):
    """A deliberately coarse explicit cap valid for ALL |p|<=p_max.

    Me(2z)/Me(z)^2 >= sqrt(z)/4 for z>=1, while
    V>=d0^2*mu^2*exp(-2|z|)/2 for all z. The split at Z gives V>=8/a_min.
    """
    a_min,d0,p_max = map(mp.mpf,(a_min,d0,p_max))
    if a_min <= 0 or d0 <= 0 or p_max < 0:
        raise ValueError('Invalid compact input bounds')
    Z = 16*(1+8*p_max*p_max/(a_min*d0*d0))**2
    return dict(Z=Z,log_mu_max=Z+mp.log(4)-mp.log(d0)-mp.log(a_min)/2,
                guaranteed_variance_lower=8/a_min, required_variance_upper=3/a_min)


def active_gap_bounds(a_min, ts_max, p1_max, p2_max, pc_input_gap,
                      d0, mu_cap, delta):
    """All-phase inequalities on the active cutoff region; no theta grid.

    Preconditions: the supplied cap reaches variance 3/a_min on the full
    parameter set, and the input data satisfy Pc(ts)>=2+pc_input_gap.
    The cap hypothesis is checked separately, not silently accepted as a proof.
    """
    a_min,ts_max,p1_max,p2_max,pc_input_gap,d0,mu_cap,delta = map(point,
        (a_min,ts_max,p1_max,p2_max,pc_input_gap,d0,mu_cap,delta))
    if min(lower(a_min),lower(d0),lower(mu_cap),lower(delta)) <= 0 \
            or min(lower(ts_max),lower(p1_max),lower(p2_max)) < 0 \
            or lower(pc_input_gap) <= 2*upper(d0):
        raise ValueError('Invalid active-loop bounds')
    T = ts_max+16*d0*mu_cap
    J = p2_max+p1_max*T
    gamma = pc_input_gap-d0
    gap = gamma-delta/2
    qgap = 2*gap*gap-delta*J*J/2
    checks = dict(delta_below_half=upper(delta)<=mp.mpf('.5'),
                  delta_below_gamma_half=upper(delta)<=lower(gamma/2),
                  delta_below_quadratic_budget=upper(delta)<=lower(gamma*gamma/(4*(1+J*J))),
                  fourth_gap_positive=lower(qgap)>0)
    if not all(checks.values()):
        raise ValueError('Delta is not certified inside the all-phase cone')
    # The cutoff transition has v-2>=delta/8, not necessarily delta/2.
    return dict(T_upper=T,J_upper=J,
                gaps=[2/(1+T*T),delta/8,gap,qgap], checks=checks,
                cap_reachability_is_separate_obligation=True)


def cone_tolerance(a_min, coordinate_bound, gap_min):
    """A coordinate sup-norm error preserving >=half of each cone gap.

    On the containing segment box, ||grad Psi_j||_1 <= 256 R^10.
    The first-coordinate floor also prevents crossing the singularity a=0.
    """
    a_min,coordinate_bound,gap_min=map(mp.mpf,(a_min,coordinate_bound,gap_min))
    if a_min <= 0 or coordinate_bound < a_min or gap_min <= 0:
        raise ValueError('Need valid coordinate bounds and positive margins')
    R=max(2,coordinate_bound+1,2/a_min)
    return min(1,a_min/2,gap_min/(512*R**10))


def frequency_requirement(*, epsilon, state_constant, correction_constant,
                          discrepancy_constant, coefficient_tolerance,
                          beta=1000, quadratic=100000):
    """Sufficient real lower bound for N; choose any larger finite integer.

    All supplied constants must already be certified for the actual compact
    profile. This helper does not certify them or select an actual-profile N.
    discrepancy_constant includes dimensional, angular, and lambda losses.
    """
    names=('epsilon','state_constant','correction_constant','discrepancy_constant',
           'coefficient_tolerance','beta','quadratic')
    values=map(mp.mpf,(epsilon,state_constant,correction_constant,discrepancy_constant,
                       coefficient_tolerance,beta,quadratic))
    d=dict(zip(names,values))
    if any(not mp.isfinite(v) or v<=0 for v in d.values()):
        raise ValueError('Every conditional bound must be finite and positive')
    e,c,k,D,r,b,q=[d[n] for n in names]
    terms=dict(profile_and_cone=2*(c+2*k*b*D)/e,
               contraction=16*b*b*q*D, small_coefficients=4*b*D/r,
               minimum_frequency=mp.mpf(1))
    return dict(terms=terms, sufficient_N_lower=max(terms.values()),
                actual_profile_inputs_verified=False)
