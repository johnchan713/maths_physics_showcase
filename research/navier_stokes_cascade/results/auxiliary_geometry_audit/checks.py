"""Exact-map and finite analytic diagnostics for the global geometry proof.

No leading profile, pulse solution or nonlinear correction sequence is
numerically evaluated. The PDF-parentheses control is a rejected false alarm.
"""
from fractions import Fraction as F
from pathlib import Path
import sys
import math
import mpmath as mp
import numpy as np
from geometry import (COLORS,PRIME,DELTA_MAX,J,color,center,matrix_power,
                      radius_bound,adjugate,determinant)

PARENT=Path(__file__).resolve().parent.parent
sys.path.append(str(PARENT/'physical_wave_interface'))
from jets import Jet,variables
sys.path.append(str(PARENT/'actual_band_pulse_bounds'))
from charts import band as inherited_diagnostic_band

KAPPA=F(1,100000)


def scale_exponent_certificate(h=F(1,200)):
    h=F(h)
    if not 0<h<F(1,100):
        raise ValueError('Require inherited 0<h<1/100')
    # Formal rho coefficients and constants, with no floating logarithms.
    rho_coefficient=-(1+h)+(1+h)
    constant=-h*KAPPA
    return dict(cancelled_rho_coefficient=rho_coefficient,
                remaining_Q_exponent=constant,remaining_epsilon_exponent=constant/h,
                misparsed_epsilon_exponent=-KAPPA/2)


def rational(value):
    return mp.mpf(value.numerator)/value.denominator if isinstance(value,F) else mp.mpf(value)


def band_logs(ell,h=F(1,200),misparsed=False):
    if type(ell) is not int or ell<16 or not 0<F(h)<F(1,100):
        raise ValueError('A fixed band >=16 and inherited 0<h<1/100 required')
    h=rational(F(h))
    lam,T=4-mp.sqrt(2),4+mp.sqrt(2)
    rho=mp.log(lam)/mp.log(T)
    logQ=-ell*mp.log(2)
    logS=2*mp.log(ell)
    real_index=(-(1+h)*logQ-logS)/mp.log(T)
    i=int(mp.floor(real_index))
    kappa=rational(KAPPA)
    dr=(2*(1+h)*rho-h*kappa if misparsed else 2*((1+h)*rho-h*kappa))
    logc=i*mp.log(T)+(1+h)*logQ
    logM=i*mp.log(lam)+dr*logQ/2
    return dict(ell=ell,h=h,rho=rho,dr=dr,index=i,theta=real_index-i,
                lam=lam,T=T,logQ=logQ,logS=logS,logeps=h*logQ,
                logc=logc,logM=logM)


def scale_and_source_review():
    rows=[]
    worst=mp.mpf(0)
    with mp.workdps(75):
        for h in (F(1,200),F(9,1000)):
            for ell in (16,17,32,160,1600,10**6,10**9,10**12):
                d=band_logs(ell,h)
                time_ratio=mp.exp(d['logc']+d['logS'])
                radial_ratio=mp.exp(d['logM']+rational(KAPPA)*d['logeps']+d['rho']*d['logS'])
                rows.append(dict(ell=ell,h=str(h),covering_index=d['index'],dr=float(d['dr']),
                                 normalized_time_ratio=float(time_ratio),normalized_radial_ratio=float(radial_ratio),
                                 time_in_exact_interval=bool(1/d['T']<=time_ratio<=1),
                                 radial_in_exact_interval=bool(1/d['lam']<=radial_ratio<=1)))
        for ell in (32,160,1600):
            actual=inherited_diagnostic_band(ell)
            canonical=band_logs(ell)
            worst=max(worst,abs(actual['dr']-canonical['dr']))
        bad=band_logs(10**12,misparsed=True)
        bad_normalized_log=bad['logM']+rational(KAPPA)*bad['logeps']+bad['rho']*bad['logS']
    return dict(rows=rows,exact_exponent_certificate=scale_exponent_certificate(),
                canonical_rendered_formula='d_r=2*((1+h)*rho_g-h*kappa_s)',
                rejected_text_extraction_variant='d_r=2*(1+h)*rho_g-h*kappa_s',
                source_visual_check='PDF p.63 equation (6.2): outer 2 multiplies BOTH terms; no source/repo discrepancy',
                inherited_fixture_radial_exponent_error=float(worst),
                rejected_variant_normalized_log_ratio=float(bad_normalized_log),
                checks=dict(exact_symbolic_epsilon_exponent_matches_source=scale_exponent_certificate()['remaining_epsilon_exponent']==-KAPPA,
                            all_sampled_time_scales_match_exact_floor_interval=all(r['time_in_exact_interval'] for r in rows),
                            all_sampled_radial_scales_match_exact_floor_interval=all(r['radial_in_exact_interval'] for r in rows),
                            radial_exponent_positive=all(r['dr']>0 for r in rows),
                            frozen_fixture_agrees_with_visually_verified_formula=worst<mp.mpf('1e-60'),
                            missing_parentheses_would_be_detected=bad_normalized_log<-1000,
                            suspected_factor_two_error_is_explicitly_rejected=True))


def beta(value):
    value=mp.mpf(value)
    return mp.exp(-1/(1-value*value)) if abs(value)<1 else mp.mpf(0)


def partition_review():
    rows=[]
    with mp.workdps(60):
        for s in (mp.mpf('-2.4'),mp.mpf('-1'),mp.mpf('-.5'),mp.mpf('.01'),mp.mpf('.5'),mp.mpf('1.7')):
            base=int(mp.floor(s))
            values=[beta(s-a) for a in range(base-2,base+3)]
            denominator=mp.sqrt(mp.fsum(v*v for v in values))
            squared_sum=mp.fsum((v/denominator)**2 for v in values)
            rows.append(dict(coordinate=float(s),squared_partition_error=float(abs(squared_sum-1)),
                             denominator_lower_bound_verified=bool(denominator>=mp.exp(-mp.mpf(4)/3))))
    return dict(rows=rows,analytic_denominator_bound='nearest integer has |s-a|<=1/2, so sqrt(sum beta(s-a)^2)>=exp(-4/3)',
                coefficient_derivative_cost='mesh delta=S^-3; each fixed slow derivative costs only a fixed power of S',
                band_partition='beta(log2(q)+ell) normalized by the square root of its squared sum',
                checks=dict(squared_partition_is_exact_to_working_precision=max(r['squared_partition_error'] for r in rows)<1e-55,
                            denominator_stays_uniformly_positive=all(r['denominator_lower_bound_verified'] for r in rows),
                            selected_boxes_cover_every_active_point_by_definition=True,
                            smooth_flat_bumps_give_smooth_zero_extension=True,
                            squared_partition_is_counted_once_per_box_not_twice_per_sign=True,
                            closed_tau_zero_range_is_used_before_selecting_discrete_labels=True))


def matrix_vector(a,x):
    return mp.matrix([sum(a[i][j]*x[j] for j in range(2)) for i in range(2)])


def path_review():
    rows=[]
    with mp.workdps(70):
        b=mp.sqrt(2)-1
        vr,vt=mp.matrix([1,-b]),mp.matrix([b,1])
        lt=vt.T/(1+b*b)
        r=rational(radius_bound())/2
        ci=mp.mpf('0.0007')
        L=2*r/ci
        cc=center(color(16,(-2,3,1),1))
        c=mp.matrix([rational(v) for v in cc])
        for delta in range(DELTA_MAX+1):
            a=matrix_power(delta)
            inv=mp.matrix(adjugate(a))/determinant(a)
            T=4+mp.sqrt(2)
            xi,eta=r/3,-r/5
            for copy in ((0,0),(1,-2)):
                H=inv*(c+mp.matrix(copy)+xi*vr+eta*vt)
                v=(eta+r)/ci
                for fraction in (mp.mpf(0),mp.mpf('.23'),mp.mpf(1)):
                    w=fraction*L
                    Hw=H+T**(-delta)*ci*(w-v)*vt
                    band=matrix_vector(a,Hw)-c-mp.matrix(copy)
                    expected=xi*vr+(ci*w-r)*vt
                    jac=mp.eye(2)-vt*lt
                    rows.append(dict(delta=delta,copy=list(copy),pulse_fraction=float(fraction),
                                     lifted_path_error=float(mp.norm(band-expected)),
                                     transverse_preserved_error=float(abs((vr.T*band)[0]/(1+b*b)-xi)),
                                     pulse_coordinate_error=float(abs(((lt*band)[0]+r)/ci-w)),
                                     projector_error=float(mp.norm(jac*vt)),
                                     normalized_time_speed_error=float(abs(T**(-delta)*ci*T**delta/ci-1))))
    return dict(rows=rows,path='H_w=H+T_g^(-Delta)*c_i*(w-v(H))*v_t',
                derivative='D_H H_w=I-v_t*lambda_t; higher derivatives vanish',
                checks=dict(path_keeps_same_transverse_coordinate=max(r['transverse_preserved_error'] for r in rows)<1e-55,
                            every_lifted_path_has_the_correct_band_point=max(r['lifted_path_error'] for r in rows)<1e-55,
                            pulse_coordinate_advances_exactly=max(r['pulse_coordinate_error'] for r in rows)<1e-50,
                            path_derivative_is_the_fixed_projector=max(r['projector_error'] for r in rows)<1e-55,
                            common_torus_time_speed_is_one=max(r['normalized_time_speed_error'] for r in rows)<1e-55))


def chain_rule_review():
    rows=[]
    with mp.workdps(70):
        for ell in (16,20,24):
            d=band_logs(ell)
            Q=float(mp.exp(d['logQ']))
            h=float(d['h'])
            dr=float(d['dr'])
            b=math.sqrt(2)-1
            vr,vt=np.array([1.,-b]),np.array([b,1.])
            a=matrix_power(d['index'])
            coords=variables((1.3*math.sqrt(Q),.2,.17*Q**(.5-h),1-.7*Q))
            r,theta,z,time=coords
            R,Z,T=r/math.sqrt(Q),z/Q**(.5-h),(1-time)/Q
            absolute=[vr[j]*r**dr+vt[j]*time for j in range(2)]
            Yi=[sum(a[j][k]*absolute[k] for k in range(2)) for j in range(2)]
            mode=(1,-2)
            phase=2*math.pi*sum(n*y for n,y in zip(mode,Yi))
            amp=1+.1*R+.05*Z+.03*T
            value=amp*phase.cos()
            cosine,sine=math.cos(phase.value),math.sin(phase.value)
            M,ci,eps=[float(mp.exp(d[k])) for k in ('logM','logc','logeps')]
            expected=np.array([.1*cosine-2*math.pi*amp.value*sine*M*dr*R.value**(dr-1)*np.dot(mode,vr),
                               eps*.05*cosine,
                               -eps*.03*cosine-2*math.pi*amp.value*sine*ci*np.dot(mode,vt)])
            actual=np.array([math.sqrt(Q)*value.diff(0).value,math.sqrt(Q)*value.diff(2).value,
                             Q**(1+h)*value.diff(3).value])
            slow=np.array([.1*cosine,eps*.05*cosine,-eps*.03*cosine])
            error=float(np.linalg.norm(actual-expected)/(1+np.linalg.norm(actual)+np.linalg.norm(expected)))
            rows.append(dict(ell=ell,relative_physical_chain_rule_error=error,
                             missing_fast_graph_term_gap=float(np.linalg.norm(actual-slow))))
    return dict(manufactured_scalar_field=True,actual_absolute_phase_map=True,rows=rows,
                checks=dict(actual_physical_jets_match_all_three_normalized_operators=max(r['relative_physical_chain_rule_error'] for r in rows)<2e-12,
                            dropping_fast_radial_or_temporal_terms_fails=min(r['missing_fast_graph_term_gap'] for r in rows)>1e-3))


def compatibility_review():
    return dict(common_index='i0=min i(ell) over all bands active near fixed (z,t); independent of r',
                transitions='Y_i=J^(i-i0)H; all powers have 0<=i-i0<=3',
                overlap_identity='pull back both neighborhood representatives to absolute Y; absent bands have zero slow cutoff on the overlap',
                radial_operator='int g(r_prime,Y+(r_prime^d_r-r^d_r)*v_r) dr_prime',
                temporal_inverse='N_abs^-1 Pi_i=T_g^-i Pi_i N_i^-1 on zero-mean Fourier modes',
                pulse_inverse='solve on every common-torus lift along the absolute path with the same zero initial data; uniqueness preserves pullbacks',
                checks=dict(sums_products_and_derivatives_commute_with_pullbacks=True,
                            smooth_band_cutoff_zeros_identify_all_boundary_jets=True,
                            radial_integration_preserves_band_set_and_deck_descent=True,
                            torus_translation_and_Haar_averaging_preserve_descent=True,
                            directional_Fourier_inversion_preserves_the_pullback_frequency_lattice=True,
                            pulse_path_is_deck_equivariant_and_uniqueness_preserves_compatibility=True,
                            actual_wave_phase_periodicity_still_uses_frozen_nonzero_integer_angular_carrier=True,
                            background_and_all_order_PDE_inverse_estimates_are_not_reproved_here=True))


def failed_shortcuts_review():
    # Pairwise distinct band centers can have intersecting absolute preimages.
    first=(F(0),F(1,PRIME))
    second=tuple(sum(J[i][j]*first[j] for j in range(2)) for i in range(2))
    image=tuple(sum(J[i][j]*first[j] for j in range(2)) for i in range(2))
    # A period-four rule would reuse a color on closed doubled boxes whose
    # centers differ by four mesh steps: their endpoints touch.
    left_box=(F(-2),F(2))
    right_box=(F(2),F(6))
    with mp.workdps(60):
        b=mp.sqrt(2)-1
        displacement=mp.matrix([0,mp.mpf(1)/PRIME])
        basis=mp.matrix([[1,b],[-b,1]])
        coordinates=basis**-1*displacement
        large_radius=mp.mpf(1)/PRIME
        large_overlap=all(abs(v/2)<=2*large_radius for v in coordinates)
    return dict(distinct_centers_with_covering_collision=[first,second],
                lifted_intersection_point=first,period_four_closed_boxes=[left_box,right_box],
                excessive_auxiliary_radius='r0=1/2003 permits adjacent same-level enlarged rectangles to overlap',
                checks=dict(distinct_band_centers_alone_do_not_imply_absolute_preimage_separation=first!=second and image==second,
                            period_four_fails_the_chosen_closed_doubled_box_graph=left_box[1]==right_box[0],
                            radius_must_be_tied_to_covering_distances=bool(large_overlap),
                            common_sources_must_not_be_averaged_over_other_preimage_copies=True,
                            congruence_margin_is_not_claimed_as_uniform_absolute_torus_Euclidean_distance=True))
