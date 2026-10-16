"""Shape/join bounds and labelled pressure-ledger controls."""
from fractions import Fraction as Q
import importlib.util
from pathlib import Path
import sys
import mpmath as mp

HERE=Path(__file__).resolve().parent
def load(name,path):
    spec=importlib.util.spec_from_file_location(name,path)
    module=importlib.util.module_from_spec(spec)
    sys.modules[name]=module
    spec.loader.exec_module(module)
    return module

previous=load('shape_previous_continuation',HERE.parent/'reduced_continuation_margins/continuation.py')
jets=load('shape_frozen_Jet2',HERE.parent/'angular_c2_transfer/jets.py')
matching=load('shape_frozen_annulus_bounds',HERE.parent/'axis_matching_audit/bounds.py')
Jet2=jets.Jet2


def scalar_review():
    eps=Q(1,10**16);gb=3*eps/100;d2=Q(286,1000)*eps
    coefficient=2000*d2
    # A tiny C2 ball gives the same root as the older C1 uniqueness ball.
    entry=eps/4+Q(7,50)*gb
    return dict(coefficient_C2_bound=str(coefficient),target_C2_bound=str(d2),checks=dict(
        source_shape_margin=Q(55,100)*Q(28,10)-Q(1,10)-Q(2,100)>1,
        barrier_crossing_points_upward=110-Q(65,100)*Q(5,2)>0,
        shape_slope_budget_at_B10=Q(64)*(10+3)/(2000*10)<Q(5,100),
        entry_below_p255_epsilon=entry<Q(255,1000)*eps,
        exp_minus_two_below_p14_from_positive_series=Q(109,15)>Q(50,7),
        target_below_p286_epsilon=Q(255,1000)*eps+gb+gb*gb<d2,
        C2_root_contraction=8*1000**2*100*d2<1,
        same_root_inside_C1_ball=coefficient<Q(6,10**13),
        raw_second_derivative_is_twice_factorial_bound=2*coefficient==Q(1144,10**16)))


def jet_review():
    checks=[]
    for eta in (Q(-1),Q(0),Q(1),Q(2,5)):
        e=Jet2(eta,Q(1),Q(0))
        g=Q(1,10**16)*(1+e+e*e)
        U=4*e+g
        exact=U*U-16*e*e-8*e*(U-4*e)
        checks.append(exact==g*g)
        inverse_f_squared=(1+e*e)*(1+e*e)
        if eta==1:
            norm=abs(inverse_f_squared.v)+abs(inverse_f_squared.d)+abs(inverse_f_squared.dd)/2
            checks.append(norm==20 and norm>12)
    with mp.workdps(80):
        # Independent differentiation of all five normalized manufactured rows.
        P,R=mp.mpf(16),mp.mpf(10000)
        def defects(e):
            return [1+e+e*e,2-e+e**3,1+2*e-e**3,3-e+2*e*e,1+e*e]
        def rows(e):
            M,I,J,S,Cp=defects(e);xc=mp.exp(-6);f=1/(1+e*e)
            ang=mp.sqrt(2)*R**mp.mpf('1.5')*P*f*xc**mp.mpf('1.6')
            return [M/(R*xc),(J-4*e*I)/ang,I/ang,-(S-8*e*M)/(R*P*P*f*f*xc**mp.mpf('1.2')),Cp/(P*P*f*f*xc**mp.mpf('.2'))]
        errors=[]
        for eta in map(mp.mpf,('-1','0','.3','1')):
            inputs=[Jet2(defects(eta)[i],mp.diff(lambda e: defects(e)[i],eta),mp.diff(lambda e: defects(e)[i],eta,2)) for i in range(5)]
            result=jets.physical_to_rows(inputs,eta,P,R)
            for i,r in enumerate(result):
                oracle=[rows(eta)[i],mp.diff(lambda e:rows(e)[i],eta),mp.diff(lambda e:rows(e)[i],eta,2)]
                errors.extend(abs(a-b) for a,b in zip((r.v,r.d,r.dd),oracle))
        zero=[Jet2(mp.mpf(0)) for _ in range(4)]+[Jet2(mp.mpf(1))]
        actual=jets.physical_to_rows(zero,mp.mpf(0),P,R)[4]
        bad=jets.physical_to_rows(zero,mp.mpf(0),P,R,True)[4]
        return dict(maximum_normalized_jet_error=float(max(errors)),checks=dict(
            energy_cancellation_and_inverse_norm_exact=all(checks),
            all_five_rows_keep_normalization_derivatives=max(errors)<mp.mpf('1e-65'),
            frozen_pressure_normalization_detected=actual.dd>mp.mpf('.01') and bad.dd==0,
            manufactured_inputs_not_actual_entry_bounds=True))


def pressure_review():
    checks=[]
    for eta in (Q(-1),Q(0),Q(1),Q(1,3)):
        e=Jet2(eta,Q(1),Q(0));A=1+e*e;total=A*A/Q(2)
        for X in (Q(1,10),Q(1),Q(10)):
            tail=A*A/(2*(1+X)**2)
            Cp=total-tail;Pi0=-total
            checks.append(Pi0+Cp==-tail)
            c=Q(1,10**30)*e
            checks.append(Pi0+c+Cp+tail==c)
        # Signed memory cancels only after BOTH separated edits.
        heat=(1+e*e)/Q(100)
        prefix=-heat
        checks.append(prefix.v!=0 and prefix+heat==Jet2(Q(0),Q(0),Q(0)))
        # On a toy fixed log-radius patch, beta=t(1-t): int beta=1/6,
        # int beta^2=1/30. Omitting the quadratic Cp term leaves c^2/60.
        c=Q(1,10**8)*(1+e)
        linear=c/Q(6);complete=linear+c*c/Q(60)
        if eta==0:
            quadratic_defect=complete-linear
            point_defect=Q(1,10**30)*e
    # I for H=A X^-h is not finite. Here h=1/2 is just an exact toy.
    I=[Q(2)*k for k in (1,2,4,8)]
    renormalized=[v-2*k for v,k in zip(I,(1,2,4,8))]
    return dict(checks=dict(
        forward_equals_backward_when_exact_total_restored=all(checks),
        quadratic_pressure_remainder_nonzero=quadratic_defect.v>0 and quadratic_defect.dd>0,
        one_angle_zero_not_functional_zero=point_defect.v==0 and point_defect.d!=0,
        raw_I_diverges_instead_of_having_zero_total=I==[2,4,8,16],
        renormalized_I_can_vanish=all(v==0 for v in renormalized),
        manufactured_pressure_not_actual_core_or_heat_profile=True),
        toy_raw_I=[str(v) for v in I],quadratic_pressure_residual=str(quadratic_defect.v))


def budget(scale,digits):
    if digits<100:
        raise ValueError('Nested parameters require at least 100 digits')
    old=previous.budget(scale,digits)
    b=previous.core.core_budget(scale,digits)
    p,log,exp=mp.iv.mpf,mp.iv.log,mp.iv.exp
    upper=lambda z:mp.make_mpf(z._mpi_[1]);lower=lambda z:mp.make_mpf(z._mpi_[0])
    def add(*v):
        shift=p(max(upper(x) for x in v));limit=log(p('1e-200'))
        # First discard only rigorously negligible positive tails. A wide
        # but wholly negligible term must not inflate the dominant j term.
        if any(upper(x-shift)>=lower(limit) and upper(x)-lower(x)>1 for x in v):
            return p([max(lower(x) for x in v),upper(p(max(upper(x) for x in v))+log(p(len(v))))])
        return shift+log(sum((p([0,'1e-200']) if upper(x-shift)<lower(limit) else exp(x-shift) for x in v),p(0)))
    lk,ll,lp=[b['logs'][k] for k in ('K','Lambda','P')]
    lu=b['logs']['axial_entry_error'];lJ=add(p(0),lk,ll)
    B=exp(log(p('1e20'))+64*lJ)
    lc=10000*B-100*log(p('1e-16'));ld=-20*lc;lj=log(p('1e-18'))
    ell=add(log(p(10))+12*lJ,p(0),ld+3*B)
    old_error=add(old['logs']['reference_source_remainder'],ld+10*B)
    new_error=add(log(p(2))+2*lj,log(p(4))+lu,log(p(2))+ld+3*B)
    axial=add(lj,log(p(3))+lu,ld+3*B)
    lq=2000*B-10*(lc+lp)
    row=log(p('2e6'))+412*B-2*lc
    eps=p('1e-16');gb=3*eps/100;d2=p('.286')*eps
    Uinv=matching.inverse_bounds(['1','1.6'],'.3','.4')
    Einv=matching.inverse_bounds(['1.6','1.2','.2'],'.15','.35')
    zeta=20*exp(p('1.2'))/256
    quadratic=32*p('.24')*exp(p('1.6'))+64*p('.16')*exp(1)*zeta
    coefficients=[18*exp(6),p(60)/16*exp(p('9.6')),
                  (p(10)/48+p(5)/8)*exp(p('9.6')),
                  (p(6500)/256+p(5)/12)*exp(p('7.2')),p('2.5')*exp(p('1.2'))+p(20)/256*exp(p('1.2'))]
    entry_increments=[exp(-2)*gb,p(5)/8*exp(-p('3.2'))*gb,p(0),p(20)/256*exp(-p('.8'))*gb*gb,p(0)]
    T=p(64) if scale=='P16' else exp(p(int(scale[2:])))+10
    llam=-4*T;lh=2*llam;lXR=log(p(110))+10*(lc+lp)
    logs=dict(endpoint_log_CE_C2_over_B=ell-log(B),endpoint_axial_offset=axial,
              old_shape_gradient_error=old_error,new_shape_gradient_error=new_error,
              W_error=add(log(p(4))+lu,log(p(2))+ld+3*B),
              shape_E_C2=204*B-lc,shape_Cp_C2=412*B-2*lc,
              pressure_length_factor=log((ll+log(p('27.5'))+2000*B)/2+2)-4*B,
              entry_separation=lq+8,normalized_rows=row,
              radius_above_10000=lXR-log(p(10000)),
              angular_Cp_root_contraction=log(p(8*25*32*700))+30*llam,
              heat_Cp_root_contraction=log(p(8*50**2*100*64))+4*lh-(lXR+241*T-18))
    checks=dict(inherited_continuation_budgets=all(old['checks'].values()),
                endpoint_log_CE_C2_below_B=upper(logs['endpoint_log_CE_C2_over_B'])<0,
                axial_offset_below_3epsilon_over100=upper(axial)<lower(log(gb)),
                old_shape_gradient_error_below_p01=upper(old_error)<lower(log(p('.01'))),
                new_shape_gradient_error_below_p01=upper(new_error)<lower(log(p('.01'))),
                W_error_below_p1=upper(logs['W_error'])<lower(log(p('.1'))),
                E_C2_below_one=upper(logs['shape_E_C2'])<0,
                pressure_length_majorant=upper(logs['pressure_length_factor'])<0,
                endpoint_before_entry=upper(logs['entry_separation'])<0,
                actual_normalized_rows_below_epsilon_over4=upper(row)<lower(log(eps/4)),
                actual_radius_above_10000=lower(logs['radius_above_10000'])>0,
                row_coefficients_below_one_million=max(upper(x) for x in coefficients)<10**6,
                entry_rows_below_p255epsilon=max(upper(eps/4+x) for x in entry_increments)<lower(p('.255')*eps),
                C2_zeta_below_one=upper(zeta)<1,
                exact_U_and_E_inverse_below1000=max(Uinv+Einv)<1000,
                C2_quadratic_bound_below100=upper(quadratic)<100,
                C2_root_contraction=upper(8*p(1000)**2*100*d2)<1,
                exact_angular_Cp_root_contraction=upper(logs['angular_Cp_root_contraction'])<0,
                exact_heat_Cp_root_contraction=upper(logs['heat_Cp_root_contraction'])<0)
    return dict(logs=logs,checks=checks,inverse_upper=float(max(Uinv+Einv)),
                entry_upper=float(max(upper(eps/4+x) for x in entry_increments)))


def budget_review():
    rows=[];errors=[]
    with mp.workdps(170):
        for scale in ('P16','Md4','Md64'):
            a,b=budget(scale,100),budget(scale,160)
            for k in a['logs']:
                for i in (0,1):
                    x,y=[mp.make_mpf(z['logs'][k]._mpi_[i]) for z in (a,b)]
                    errors.append(abs(x-y)/max(1,abs(x),abs(y)))
            encode=lambda z:{n:mp.nstr(mp.make_mpf(z._mpi_[i]),130) for i,n in enumerate(('lower','upper'))}
            rows.append(dict(scale=scale,checks=b['checks'],inverse_upper=b['inverse_upper'],entry_upper=b['entry_upper'],logs={k:encode(v) for k,v in b['logs'].items()}))
        return dict(rows=rows,maximum_relative_refinement_error=float(max(errors)),checks=dict(
            all_outward_shape_join_and_Cp_acceptance_checks=all(all(r['checks'].values()) for r in rows),
            precision_refinement_consistent=max(errors)<mp.mpf('1e-55'),
            outer_Md4_cone_failure_not_relabelled=True,finite_scales_not_continuum_proofs=True))
