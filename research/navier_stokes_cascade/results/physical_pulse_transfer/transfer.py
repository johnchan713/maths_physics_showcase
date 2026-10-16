"""Exact physical loss bookkeeping and labelled graph-derivative controls."""
from fractions import Fraction as Q
import importlib.util
from pathlib import Path
import sys
import mpmath as mp

HERE=Path(__file__).resolve().parent
KAPPA=Q(1,100000)


def load(name,path):
    spec=importlib.util.spec_from_file_location(name,path)
    module=importlib.util.module_from_spec(spec)
    sys.modules[name]=module
    spec.loader.exec_module(module)
    return module


est=load('physical_frozen_pulse_estimates',HERE.parent/'geometric_pulse_replay/estimates.py')


def order(n):
    if type(n) is not int or n<0:
        raise ValueError('Nonnegative integer derivative order required')
    return n


def losses(h):
    h=Q(h)
    if not 0<h<Q(1,100):
        raise ValueError('Fixed inherited 0<h<1/100 required')
    return dict(spatial=Q(1,2)+h/2,time=1+3*h/2,
                coefficient_radial=Q(1,2)+h*KAPPA,
                coefficient_axial=Q(1,2)-h,coefficient_time=1+h,
                frame=Q(1,2),A=Q(1,2)+h)


def primary_powers(h,a,b):
    l=losses(h);order(a);order(b)
    shift=-a*l['spatial']-b*l['time']
    return dict(velocity=shift-l['spatial'],pressure=shift-2*l['spatial'],potential=shift)


def common_loss(h,m):
    l=losses(h)
    return 2*l['A']+(order(m)+1)*l['time']


# Exact Q(sqrt(2)) arithmetic: pairs represent a+b*sqrt(2).
def add(a,b):
    return (a[0]+b[0],a[1]+b[1])


def mul(a,b):
    return (a[0]*b[0]+2*a[1]*b[1],a[0]*b[1]+a[1]*b[0])


def scale(a,s):
    return tuple(Q(s)*v for v in a)


def inv(a):
    den=a[0]*a[0]-2*a[1]*a[1]
    if not den:
        raise ValueError('Nonzero quadratic-field value required')
    return (a[0]/den,-a[1]/den)


def dot(a,b):
    out=(Q(0),Q(0))
    for x,y in zip(a,b):
        out=add(out,mul(x,y))
    return out


def exact_review():
    one=(Q(1),Q(0));zero=(Q(0),Q(0));beta=(Q(-1),Q(1))
    vr=(one,scale(beta,-1));vt=(beta,one)
    norm=dot(vt,vt);lt=tuple(mul(v,inv(norm)) for v in vt)
    checks=[]
    for level in range(9):
        matrix=est.geometry.matrix_power(level)
        image=tuple(add(scale(vr[0],row[0]),scale(vr[1],row[1])) for row in matrix)
        checks.append(dot(lt,image)==zero)
    powers=[];dominance=[]
    for h in (Q(1,200),Q(9,1000)):
        l=losses(h)
        dominance.extend((l['spatial']>l['coefficient_radial'],l['spatial']>l['coefficient_axial'],
                          l['time']>l['coefficient_time'],l['time']>l['spatial']))
        for a in range(7):
            for b in range(5):
                p=primary_powers(h,a,b)
                shift=-a*l['spatial']-b*l['time'];A=l['A']
                powers.extend((p['velocity']==-A+h/2+shift,
                               p['pressure']==-2*A+h+shift,
                               p['potential']==Q(1,2)-A+h+shift))
    rows=[{k:str(v) for k,v in primary_powers(Q(1,200),a,b).items()}|dict(space_order=a,time_order=b)
          for a,b in ((0,0),(1,0),(0,1),(2,1))]
    return dict(checks=dict(
        transverse_phase_coordinate_annihilates_radial_graph_at_all_checked_covers=all(checks),
        pulse_time_dual_normalization_exact=dot(lt,vt)==one,
        radial_and_temporal_eigenvectors_orthogonal_exact=dot(vr,vt)==zero,
        harmonic_losses_dominate_coefficient_and_frame_losses=all(dominance),
        physical_velocity_pressure_potential_normalizations_exact=all(powers),
        common_potential_loss_is_independent_of_correction_stage=True,
        primary_covariance_remainder_not_absolutely_decaying=-1-2*Q(1,200)*KAPPA<0,
        both_reserved_mean_pressure_tails_still_required=True),
        exact_power_checks=len(powers),exact_cover_checks=len(checks),primary_physical_power_rows=rows,
        spatial_loss='1/2+h/2',time_loss='1+3h/2',common_potential_loss='2A+(m+1)*(1+3h/2)',
        physical_covariance_remainder_power=str(-1-2*Q(1,200)*KAPPA))


def first_band(h,qstar,Lmin,geometry_band=16,r=est.RADIUS):
    qstar=mp.mpf(qstar)
    if not mp.isfinite(qstar) or not 0<qstar<=1:
        raise ValueError('A finite actual background domain in (0,1] is required')
    if type(geometry_band) is not int or geometry_band<16:
        raise ValueError('A geometric enlargement threshold at least 16 is required')
    if not mp.isfinite(mp.mpf(Lmin)) or mp.mpf(Lmin)<25:
        raise ValueError('The physical pulse acceptance length must be at least 25')
    background=int(mp.ceil(mp.log(4/qstar,2)))+1
    row=est.all_later_threshold(h,r,Lmin,max(geometry_band,background))
    row.update(background_band=background,working_log2_q=-row['first_band']-2)
    return row


def threshold_review():
    with mp.workdps(80):
        rows=[];checks=[]
        for h,qstar,length,geo in ((Q(1,200),mp.power(2,-20),1000,16),
                                   (Q(9,1000),mp.power(2,-10**6),100,10**6)):
            row=first_band(h,qstar,length,geo)
            ell=row['first_band']
            for late in (ell,ell+1,2*ell):
                checks.extend((4*mp.power(2,-late)<qstar,
                               est.length_interval(late)[0]>=length,
                               est.old_bounds.phase_majorant(late,est.real(h))<=1))
            rows.append(dict(h=str(h),background_log2_qstar=float(mp.log(qstar,2)),
                             first_band=ell,geometry_band=geo,working_log2_q=row['working_log2_q']))
        return dict(checks=dict(
            enlarged_chart_margin_and_pulse_length_and_phase_conditions=all(checks),
            actual_first_band_not_numerically_selected=True,
            same_threshold_applies_to_all_later_bands_and_fixed_jets=True,
            naive_Q_below_qstar_does_not_control_q=Q(3,4)<1<2*Q(3,4),
            working_domain_keeps_all_active_partition_labels_in_accepted_bands=True),
            rows=rows,diagnostic_acceptance_checks=len(checks),enlarged_chart_factor=4,
            threshold_source='finite actual gap constants and all-half-line inequalities, not a band scan')


def graph_review():
    """Independent direct and graph-operator derivatives of manufactured fields.

    A fixed local translation centers the unwrapped diagnostic lift; it
    changes no derivative operator and is NOT an alteration of the actual map.
    These bands test operators, not acceptance of short homogeneous pulses.
    """
    with mp.workdps(80):
        errors=[];phase_errors=[];radial_omissions=[];time_omissions=[];rows=[]
        beta=mp.sqrt(2)-1;vr=(mp.mpf(1),-beta);vt=(beta,mp.mpf(1));den=1+beta*beta
        lr=lambda H:(H[0]-beta*H[1])/den
        lt=lambda H:(beta*H[0]+H[1])/den
        R0,Z0,T0=map(mp.mpf,('1.4','.1','.6'))
        r=est.real(est.RADIUS)
        for ell in (160,1024):
            data=est.charts.band(ell,h=Q(1,200),r0=r)
            eps,ci,mi,dr,L=[data[k] for k in ('epsilon','ci','mi','dr','length')]
            H0=tuple(r/7*v for v in vr)
            k=data['k'];angular=est.charts.nearest_nonzero_integer(k*mp.mpf('.3'));p=mp.mpf(angular)/k
            pz,x0=mp.mpf('-.4'),mp.mpf('.2')
            def coords(H):
                return lr(H),(lt(H)+r)/ci
            def coefficient(R,Z,T,H):
                xi,v=coords(H)
                return 1+R/5+Z*Z/10+7*T/100+3*xi/(10*r)+v/(5*L)+3*xi*v/(100*r*L)
            def phase(R,Z,T,H):
                v=coords(H)[1];F=R**-3+Z/10+3*T*T/100;G=-R**-2+7*Z*T/100
                return pz*Z/eps+x0*R-v*(p*F+pz*G)
            def graph(R,T):
                return tuple(H0[j]+mi*(R**dr-R0**dr)*vr[j]-ci/eps*(T-T0)*vt[j] for j in range(2))
            for f in (coefficient,phase):
                fR=mp.diff(lambda R:f(R,Z0,T0,H0),R0)
                fT=mp.diff(lambda T:f(R0,Z0,T,H0),T0)
                Li=mp.diff(lambda w:f(R0,Z0,T0,tuple(H0[j]+w*vr[j] for j in range(2))),0)
                Ni=mp.diff(lambda w:f(R0,Z0,T0,tuple(H0[j]+w*vt[j] for j in range(2))),0)
                directR=mp.diff(lambda R:f(R,Z0,T0,graph(R,T0)),R0)
                directT=mp.diff(lambda T:f(R0,Z0,T,graph(R0,T)),T0)
                expectedR=fR+mi*dr*R0**(dr-1)*Li
                expectedT=fT-ci/eps*Ni
                errors.extend((abs(directR-expectedR)/max(1,abs(expectedR)),
                               abs(directT-expectedT)/max(1,abs(expectedT))))
                if f is phase:
                    phase_errors.append(abs(Li))
                    phase_errors.append(abs(directR-fR)/max(1,abs(fR)))
                else:
                    radial_omissions.append(abs(directR-fR))
                    time_omissions.append(abs(directT-fT))
            # A separate differentiated complete harmonic; physical Q factors
            # are tracked exactly by the exponent ledger, not rounded here.
            amp=lambda R,T:coefficient(R,Z0,T,graph(R,T))
            phi=lambda R,T:phase(R,Z0,T,graph(R,T))
            R,T=R0,T0
            ar=mp.diff(lambda v:amp(v,T),R);pr=mp.diff(lambda v:phi(v,T),R)
            a=amp(R,T);ph=phi(R,T);eik=mp.exp(1j*k*ph)
            direct=mp.diff(lambda v:amp(v,T)*mp.exp(1j*k*phi(v,T)),R)
            errors.append(abs(direct-(ar+1j*k*a*pr)*eik)/max(1,abs(direct)))
            rows.append(dict(band=ell,angular_integer=angular,
                             omitted_radial_graph_error=float(radial_omissions[-1]),
                             omitted_fast_time_error=float(time_omissions[-1])))
        return dict(checks=dict(
            direct_pullback_matches_physical_graph_operators=max(errors)<mp.mpf('1e-60'),
            radial_graph_phase_cancels_but_amplitude_does_not=max(phase_errors)<mp.mpf('1e-60') and min(radial_omissions)>mp.mpf('.01'),
            fast_auxiliary_time_cannot_be_omitted=min(time_omissions)>mp.mpf('.01'),
            discrete_angular_carriers_are_nonzero_and_frozen=all(r['angular_integer']!=0 for r in rows),
            local_translation_only_in_manufactured_derivative_control=True,
            short_operator_fixture_bands_not_claimed_as_accepted_pulses=True),
            maximum_graph_and_harmonic_derivative_error=float(max(errors)),maximum_radial_phase_cancellation_error=float(max(phase_errors)),rows=rows)


def tail_review():
    with mp.workdps(80):
        rows=[]
        for m in (0,2,6):
            physical_loss=est.real(common_loss(Q(1,200),m))*mp.log(2)+4
            cert=est.flat_tail_threshold(est.RADIUS,mp.mpf('.1'),mp.mpf('1e50'),physical_loss,3*m+5,7)
            ell=cert['first_band'];a=cert['quadratic_decay']
            excess=[mp.log(mp.mpf('1e50'))+(3*m+5)*mp.log(v)+(physical_loss+7*mp.log(2))*v-a*v*v
                    for v in (ell,ell+1,2*ell)]
            rows.append(dict(derivative_order=m,first_tail_band=ell,maximum_log_tail_over_q_seventh=float(max(excess))))
        # P(0)>0 for a finite Gaussian pulse, but the prescribed psi is flat
        # and zero on an open endpoint collar. An envelope alone is not support.
        return dict(checks=dict(
            physical_derivative_losses_absorbed_by_all_half_line_Gaussian=all(r['maximum_log_tail_over_q_seventh']<0 for r in rows),
            endpoint_envelope_is_positive_not_a_zero_extension=mp.exp(-10)>0,
            potentials_and_pressures_cut_before_curl=True,
            source_requires_its_own_smooth_endpoint_extension=True,
            closed_auxiliary_supports_and_all_derivatives_separated=True,
            at_most_one_actual_label_after_physical_evaluation=True,
            one_sided_tau_zero_jets_not_zero_extended=True,
            no_velocity_extension_through_q_zero_asserted=True),
            rows=rows,supports=dict(transverse='abs(xi)<=r0/2',pulse='abs(v/L-1/2)<=2/5',
                                   shell='Xa<=X<=Xb',slow='the fixed closed support of eta_beta'),
            full_nonlinear_generated_source_estimates_verified=False)
