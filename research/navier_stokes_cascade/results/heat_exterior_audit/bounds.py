"""Explicit constants for the three-moment heat edit of the reference flow.

The exact heat kernel and exact small root define the mathematical edit.
Numerical Taylor and graded-root surrogates are separate diagnostics.
"""
import importlib.util
from pathlib import Path
import mpmath as mp

HERE = Path(__file__).resolve().parent
RESULTS = HERE.parent
spec = importlib.util.spec_from_file_location('heat_interval_helpers',RESULTS/'axial_stress_audit'/'bounds.py')
helpers = importlib.util.module_from_spec(spec)
spec.loader.exec_module(helpers)
point,lo,hi = helpers.point,helpers.lo,helpers.hi


def heat_bound(md=64,digits=90,radius=100):
    if not isinstance(md,int) or md<4 or digits<60 or radius<100:
        raise ValueError('Require Md>=4, >=60 digits, and XR>=100')
    mp.iv.dps = digits
    p,exp,log = point,mp.iv.exp,mp.iv.log
    lam_interval = mp.iv.mpf(['0','.01'])
    rates = [-1-2*lam_interval,-2*lam_interval,1-lam_interval]
    nodes = [exp(2*a) for a in rates]
    inverse_rows = [p(0) for _ in range(3)]
    for k in range(3):
        i,j = [n for n in range(3) if n!=k]
        denominator = abs((nodes[k]-nodes[i])*(nodes[k]-nodes[j]))
        factor = p('.3')*exp(p('.5')*rates[k]-p('.15')*abs(rates[k]))
        for degree,numerator in enumerate((nodes[i]*nodes[j],nodes[i]+nodes[j],p(1))):
            inverse_rows[degree] += numerator/(denominator*factor)
    # T>=64 bounds the whole chosen family without forming exp(-1/lambda).
    T,lam,h = p(64),exp(-256),exp(-512)
    xstar_log = log(radius)+241*T-18
    hxstar_log = log(radius)+233*T-18
    mu = exp(-xstar_log)
    target = 64*h**4*mu
    coefficient = 6400*h**4*mu
    eps_q = p('1e12')*exp(-hxstar_log)
    eps_w = p('1e17')*exp(-hxstar_log)
    eps_a = h*eps_q
    e_star = exp(-p('119.5')*T+p('10.3'))
    late_w = 60000*T*e_star/lam
    gap1_loss = eps_a+15*eps_w
    gap2_loss = p('1e7')*(eps_a+eps_w)
    finite_gap = p('.81')-4*14000*120/p('1e12')
    outer_ratio = 4*h**4
    outer_loss = 32*h**9
    checks = dict(
        T_floor=lo(exp(md)+10)>=64,
        lambda_range=0<lo(lam) and hi(lam)<p('.01'),
        h_range=0<lo(h) and hi(h)<p('.01'),
        independent_moment_weights=hi(nodes[0])<lo(nodes[1]) and hi(nodes[1])<lo(nodes[2]),
        inverse_below_50=max(hi(v) for v in inverse_rows)<50,
        quadratic_norm_below_100=p(32)*p('.9')<100,
        contraction=hi(8*50**2*100*target)<1,
        coefficient_small=0<lo(coefficient) and hi(coefficient)<p('1e-4'),
        relative_patch_positive=hi(64*coefficient)<p('.01'),
        patch_slope_small=hi(p(20000)/p('.3')*coefficient/(1-64*coefficient))<lo(lam/2),
        e_star_below_h14=hi(e_star)<lo(h**14),
        late_power_w_below_one=hi(late_w)<1,
        heat_factor_positive=hi(3*h/100)<p('.01'),
        heat_slope_excess=hi(4/ p(100))<p(1)/4,
        terminal_cutoff_slope=hi(p('.001')*64/(2*(1-p('.001')*h)))<p(1)/4,
        relative_Q_error=0<lo(eps_q) and hi(eps_q)<p('1e-6000'),
        absolute_w_error=0<lo(eps_w) and hi(eps_w)<p('1e-6000'),
        positive_shear=hi(eps_q)<p('.5'),
        first_margin=hi(gap1_loss)<p('.01'),
        second_margin=hi(gap2_loss)<p('.01'),
        large_finite_ps1=lo(hxstar_log-log(4000))>hi(log(p('1e12'))),
        finite_Pc_above_vs=p('.29')*p('1e12')>120,
        finite_cone_gap=lo(finite_gap)>p('.80'),
        outer_direction=0<lo(outer_loss) and hi(outer_loss)<p('.01'),
    )
    actual_T = exp(md)+10
    return dict(status='heat-compensated-outer-reference-bounded' if all(checks.values()) else 'inconclusive',
        Md=md,interval_digits=digits,XR=radius,checks=checks,
        inverse_row_upper=[hi(v) for v in inverse_rows],inverse_bound=50,quadratic_bound=100,
        target_C1_bound='64 h^4/Xstar',coefficient_C1_bound='6400 h^4/Xstar',
        target_universal_upper=hi(target),coefficient_universal_upper=hi(coefficient),
        log_Xstar_interval=[lo(log(radius)+241*actual_T-18),hi(log(radius)+241*actual_T-18)],
        log_h_Xstar_interval=[lo(log(radius)+233*actual_T-18),hi(log(radius)+233*actual_T-18)],
        Q_floor='h/4000',relative_Q_error_upper=hi(eps_q),absolute_w_error_upper=hi(eps_w),
        a_minus_two_floor='h',shear_error_over_h_upper=hi(eps_q),bs_unchanged=True,
        first_ratio_margin_lower='1.47',second_ratio_margin_lower='.81',
        finite_ps1_lower='1e12',finite_normalized_cone_gap_lower=lo(finite_gap),
        outer_direction_ratio_bound='4 h^4',outer_direction_ratio_upper=hi(outer_ratio),
        outer_direction_gap_lower='1.99',outer_direction_loss_upper=hi(outer_loss),
        outer_theta_coefficient='16 co h Epow(Xb) H(2d/Xb) exp(1)/sqrt(2 Xb)',
        outer_ratio_over_delta6_limit='eta Xb Epow(Xb) H(2d/Xb)/(64 L)',
        exact_moments=['M','J','S(infinity)','Cp(infinity)','renormalized I'],
        scope='Reference field for X>0; regular axis, annular attachment, full PDE corrections and smooth forcing remain unverified.')
