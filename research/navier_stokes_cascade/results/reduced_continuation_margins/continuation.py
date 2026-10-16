"""Exact taper identities and labelled diagnostics of the original sources."""
from fractions import Fraction as Q
import importlib.util
from pathlib import Path
import mpmath as mp

HERE = Path(__file__).resolve().parent
spec = importlib.util.spec_from_file_location('reduced_frozen_core',HERE.parent/'core_contraction_audit/core.py')
core = importlib.util.module_from_spec(spec)
spec.loader.exec_module(core)


def taper(p, q, kappa, beta, ratio, dp1=0, dp2=0):
    """Exact rational four-state map before the final a interpolation."""
    p,q,kappa,beta,ratio,dp1,dp2 = map(Q,(p,q,kappa,beta,ratio,dp1,dp2))
    if p <= 0 or kappa <= 0 or ratio <= 0 or not 0 <= beta <= 1:
        raise ValueError('Positive denominators and beta in [0,1] required')
    a, b = kappa*p,-beta*kappa*q*ratio
    t = -b/a
    v, Pc = a*(1+t*t),p+dp1+(q+dp2)*t
    main_Pc = p+beta*q*q/p
    main_v = kappa*(p+beta*beta*q*q/p)
    return dict(a=a,b=b,t=t,v=v,Pc=Pc,
                Pc_error=dp1+beta*q*q/p*(ratio-1)+beta*q/p*ratio*dp2,
                v_error=kappa*beta*beta*q*q/p*(ratio*ratio-1),
                main_Pc=main_Pc,main_v=main_v,
                literal_Pc_error=Pc-main_Pc,literal_v_error=v-main_v)


def taper_review():
    rows = []
    for q in (0,7,10**8):
        for beta in (Q(0),Q(1,3),Q(1)):
            for kappa in (Q(1,10**40),Q(1,10**80)):
                r=taper(3,q,kappa,beta,Q(1000001,1000000),Q(1,10**12),Q(-1,10**12))
                rows.append(r)
    invariant=all(rows[i]['t']==rows[i+1]['t'] and rows[i]['Pc']==rows[i+1]['Pc']
                  for i in range(0,len(rows),2))
    wrong=taper(3,7,Q(1,10**40),Q(1,3),1)
    return dict(samples=len(rows),checks=dict(
        exact_Pc_error_identity=all(r['Pc_error']==r['literal_Pc_error'] for r in rows),
        exact_v_error_identity=all(r['v_error']==r['literal_v_error'] for r in rows),
        common_kappa_cancels_before_estimation=invariant,
        axial_taper_main_Pc_not_below_p1=all(r['main_Pc']>=3 for r in rows),
        beta_not_beta_squared_in_Pc=wrong['main_Pc'] != 3+Q(1,9)*49/3,
        small_growth_fixture=all(r['v']<Q(1,2) for r in rows),
        manufactured_fixtures_not_actual_profile=True))


def ode_review():
    rows=[]
    for chi in (Q(99,100),Q(1)):
        A=Q(47,25)*chi
        # Replacing p1(0) by 2.3 gives A*r+(2.3-A)/r.
        for r in (Q(1),Q(2),Q(100)):
            p=A*r+(Q(23,10)-A)/r
            derivative=A-(Q(23,10)-A)/(r*r)
            rows.append((p,derivative))
    # Positive Duhamel weights preserve a source/initial sign, even for tiny j.
    j=Q(1,10**18)
    sign_controls=[]
    for sign in (-1,1):
        for r in (Q(1),Q(3),Q(1000)):
            n=sign*(Q(128,100)*j/r+Q(255,100)*j*(1-1/r))
            sign_controls.append(sign*n>j)
    return dict(checks=dict(
        angular_lower_comparison_increases=all(p>=Q(23,10) and dp>0 for p,dp in rows),
        positive_axial_integrating_factor_preserves_sign=all(sign_controls),
        crossing_barrier_points_upward=100-Q(3,5)*Q(5,2)>0,
        source_scalar_absorption=Q(524)<Q(2,100)*10**44*Q(98,100),
        Young_constant_sufficient=Q(1)/(Q(3,100)*Q(98,100))<35,
        universal_source_remainder_prefactor=Q(2196,10**6)+Q(1000,10**12)+Q(35,10**24)<Q(3,1000)))


def power_state(y, eta, alpha):
    """Manufactured power fields and all five primitives on X>=1.

    Incoming moment functions are arbitrary smooth algebraic fixtures.
    They are not the exact core, and no actual cone bound is inferred.
    """
    x=mp.exp(y)
    e=mp.mpf('1.2')+eta/10+eta*eta/20
    u=4*eta+mp.mpf('.01')
    E=e*x**alpha
    def integral(rate):
        return mp.expm1(rate*y)/rate if rate else y
    M=1+eta/5+u*mp.expm1(y)
    I=2+eta**2/10+mp.sqrt(2)*e*integral(alpha+mp.mpf('1.5'))
    J=eta/5+u*mp.sqrt(2)*e*integral(alpha+mp.mpf('1.5'))
    S=mp.mpf('.3')+eta**2/7+u*u*mp.expm1(y)-e*e*integral(2*alpha+1)/2
    Pi=-5-eta**2+e*e*integral(2*alpha)/2
    return dict(X=x,E=E,U=u,H=mp.sqrt(2*x)*E,M=M,I=I,J=J,S=S,Pi=Pi)


def power_coordinates(y,eta,alpha):
    r=power_state(y,eta,alpha)
    h=mp.mpf('.005')
    A,D,d,L=mp.mpf('.5')+h,mp.mpf('.5')-h,1-eta**2,1-2*h*eta**2
    M,I,J,S,Pi=[r[k] for k in ('M','I','J','S','Pi')]
    Me,Ie,Je,Se,Pie=[mp.diff(lambda t,k=k:power_state(y,t,alpha)[k],eta)
                     for k in ('M','I','J','S','Pi')]
    W=1-(2*D*eta*M+d*Me)/r['X']
    Qs=-W+((1-h)*I-D*eta*Ie-d*Je+2*(h-D)*eta*J)/(r['X']*r['H'])
    Ns=-W*r['U']+(D*(M-eta*Me)+4*h*eta*S-d*Se)/r['X']+4*A*eta*Pi-d*Pie
    return dict(**r,W=W,Ns=Ns,p1=r['X']*Qs/L,p2=r['X']*Ns/(L*r['E']),Pie=Pie)


def source_review():
    with mp.workdps(75):
        errors,pressure_value,pressure_eta,wrong_ell,wrong_Hc=[],[],[],[],[]
        for alpha in map(mp.mpf,('.5','.1')):
            for eta in map(mp.mpf,('-1','.3','1')):
                y=mp.mpf('.37')
                r=power_coordinates(y,eta,alpha)
                h=mp.mpf('.005')
                A,D,d,L=mp.mpf('.5')+h,mp.mpf('.5')-h,1-eta**2,1-2*h*eta**2
                Ueta=mp.diff(lambda t:power_state(y,t,alpha)['U'],eta)
                logEeta=mp.diff(lambda t:mp.log(power_state(y,t,alpha)['E']),eta)
                ell=mp.mpf('.5')+alpha
                Hc=D*eta+d*r['U']
                Sq=-r['W']*ell-h*(1-2*eta*r['U'])-Hc*logEeta
                Sn=-A*(1-2*eta*r['U'])*r['U']-Hc*Ueta-d*r['Pie']+4*A*eta*r['Pi']+eta*r['E']**2
                py=mp.diff(lambda z:power_coordinates(z,eta,alpha)['p1'],y)
                ny=mp.diff(lambda z:power_coordinates(z,eta,alpha)['Ns'],y)
                errors += [abs(py+ell*r['p1']-r['X']*Sq/L),abs(ny+r['Ns']-Sn)]
                pressure_value.append(abs(4*A*eta*r['Pi']))
                pressure_eta.append(abs(d*r['Pie']))
                # Wrong convention ell=alpha loses a fixed half, not a width error.
                wrong_ell.append(abs(py+alpha*r['p1']-r['X']*Sq/L))
                Sq_bad=-r['W']*ell-h*(1-2*eta*r['U'])-(D*eta+d*r['M']/r['X'])*logEeta
                wrong_Hc.append(abs(py+ell*r['p1']-r['X']*Sq_bad/L))
        # The post-taper log ratio persists at fixed X/Xend independently of width.
        defect=-mp.mpf(2)/5*mp.log(mp.mpf('1.1'))
        ratios=[abs(defect)/mp.exp(-n) for n in (20,40)]
        return dict(maximum_original_ODE_error=float(max(errors)),
                    post_taper_fixed_ratio_log_error=float(defect),
                    error_divided_by_width=[float(x) for x in ratios],
                    checks=dict(
                        complete_original_forward_source_ODEs=max(errors)<mp.mpf('1e-65'),
                        pressure_value_term_cannot_be_dropped=max(pressure_value)>1,
                        pressure_eta_term_cannot_be_dropped=max(pressure_eta)>mp.mpf('.01'),
                        ell_is_half_plus_E_logarithmic_slope=max(wrong_ell)>mp.mpf('.01'),
                        Hc_uses_local_U_not_radial_average=max(wrong_Hc)>mp.mpf('.01'),
                        final_field_error_not_width_small=abs(defect)>mp.mpf('.03') and ratios[1]>ratios[0]*10**8,
                        fixture_not_actual_continuum_profile=True))


def budget(scale,digits):
    """Outward logarithmic checks, without materializing C or tiny widths."""
    if digits<100:
        raise ValueError('Nested parameter budgets need at least 100 digits')
    b=core.core_budget(scale,digits)
    p,exp,log=mp.iv.mpf,mp.iv.exp,mp.iv.log
    upper=lambda z:mp.make_mpf(z._mpi_[1])
    lower=lambda z:mp.make_mpf(z._mpi_[0])
    enclosures=dict(wide_intervals=0,negligible_positive_tails=0)
    def add(*values):
        # A log-C interval can have huge absolute width even with tiny
        # relative width. Subtracting its upper endpoint then exponentiating
        # the lower gap would materialize a nested tiny number. Use the
        # elementary max <= logsumexp <= max+log(n) enclosure instead.
        if any(upper(v)-lower(v)>1 for v in values):
            enclosures['wide_intervals']+=1
            low=max(lower(v) for v in values)
            high=upper(p(max(upper(v) for v in values))+log(p(len(values))))
            return p([low,high])
        shift=p(max(upper(v) for v in values))
        # Never exponentiate a gap of size log C: even its exponent integer
        # would be impractical. A positive outward tail bound suffices.
        cutoff=log(p('1e-200'))
        terms=[]
        for v in values:
            if upper(v-shift)<lower(cutoff):
                enclosures['negligible_positive_tails']+=1
                terms.append(p([0,'1e-200']))
            else:
                terms.append(exp(v-shift))
        return shift+log(sum(terms,p(0)))
    lk,ll=b['logs']['K'],b['logs']['Lambda']
    lu=b['logs']['real_evaluation_bound']+b['logs']['ball_norm_bound']-ll
    le=b['logs']['six_parameter_two_radial_error']
    lJ=add(p(0),lk,ll)
    B=exp(log(p('1e20'))+64*lJ)
    logC=10000*B-100*log(p('1e-16'))
    ld=-20*logC
    lsigma=log(p('1e-20')); lj=log(p('1e-18'))
    source=add(log(p(2192))+le,log(p(4))+lu,log(p(1000))+lu-lsigma,
               log(p(35))+2*lu+ll-2*lsigma,log(p(240))+ld+2*lk)
    axial=add(log(p(100))+lu,log(p(5))-2*logC+B)
    comparison=log(p(10))+4*logC+100*B+ld
    logs=dict(reference_source_remainder=source,
              axial_source_error_over_j_over_100=axial-lj+log(p(100)),
              weak_sector_growth_lower_log=log(p(16))+2*lj+2*logC-3*B-2*ll,
              upper_reference_growth_over_half_C3=log(p(2))+6*B-logC,
              combined_taper_comparison_error=comparison,
              reduced_growth_upper=add(-17*logC-log(p(2)),comparison),
              final_source_change=ld+10*B,
              cutoff_inside_native=ld-log(log(p('4.1')/4)/2),
              taper_fits_before_Xi=ld-log(log(p('1.1'))/2),
              amplitude_above_K=lk-logC)
    checks=dict(reference_source_remainder_below_p1=upper(source)<lower(log(p('.1'))),
                axial_error_below_j_over_100=upper(logs['axial_source_error_over_j_over_100'])<0,
                weak_sector_growth_above_2p3=lower(logs['weak_sector_growth_lower_log'])>upper(log(p('2.3'))),
                upper_reference_growth_below_half_C3=upper(logs['upper_reference_growth_over_half_C3'])<0,
                taper_comparison_error_below_p01=upper(comparison)<lower(log(p('.01'))),
                reduced_growth_below_half=upper(logs['reduced_growth_upper'])<lower(log(p('.5'))),
                final_source_change_below_p1=upper(logs['final_source_change'])<lower(log(p('.1'))),
                reference_cutoff_inside_native=upper(logs['cutoff_inside_native'])<0,
                both_tapers_fit=upper(logs['taper_fits_before_Xi'])<0,
                amplitude_above_K=upper(logs['amplitude_above_K'])<0)
    # The full width itself is not returned or sampled.
    return dict(logs=logs,checks=checks,enclosures=enclosures)


def budget_review():
    rows,errors=[],[]
    with mp.workdps(170):
        for scale in ('P16','Md4','Md64'):
            a,b=budget(scale,100),budget(scale,160)
            for name in a['logs']:
                for i in (0,1):
                    x,y=[mp.make_mpf(z['logs'][name]._mpi_[i]) for z in (a,b)]
                    errors.append(abs(x-y)/max(1,abs(y)))
            encode=lambda z:{key:mp.nstr(mp.make_mpf(z._mpi_[i]),130) for i,key in enumerate(('lower','upper'))}
            rows.append(dict(scale=scale,checks=b['checks'],enclosures=b['enclosures'],logs={k:encode(v) for k,v in b['logs'].items()}))
        return dict(rows=rows,maximum_relative_refinement_error=float(max(errors)),checks=dict(
            all_outward_frozen_parameter_checks=all(all(r['checks'].values()) for r in rows),
            precision_refinement_consistent=max(errors)<mp.mpf('1e-55'),
            wide_intervals_and_positive_tails_exercised=all(all(n>0 for n in r['enclosures'].values()) for r in rows),
            universal_arguments_are_separate_from_finite_scales=True,
            no_full_amplitude_or_actual_profile_materialized=True))
