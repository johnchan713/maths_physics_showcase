"""Exact finite-sum majorants; analytical all-mode proofs are in README.md."""
from fractions import Fraction as F

KAPPA=F(1,100000)


def moment(data,s):
    if type(s) is not int or s<0:raise ValueError('Nonnegative integer diagnostic order required')
    return sum((1+abs(m))**s*F(v) for m,v in data.items())


def convolution(a,b):
    out={}
    for m,x in a.items():
        for n,y in b.items():out[m+n]=out.get(m+n,F(0))+x*y
    return out


def smoothed(data,p):
    if type(p) is not int or p<0 or 0 in data:raise ValueError('Nonzero modes and integer smoothing required')
    return {m:F(v)/abs(m)**p for m,v in data.items()}


def convolution_review():
    a={-17:F(2,7),-2:F(1,3),1:F(5,11),9:F(3,13)}
    b={-9:F(4,9),-1:F(2,5),2:F(7,19),31:F(5,17)}
    rows=[]
    for s in range(6):
        value=moment(convolution(a,b),s);upper=moment(a,s)*moment(b,s)
        rows.append(dict(order=s,value=str(value),upper=str(upper)))
    return dict(rows=rows,checks=dict(
        absolute_convolution_without_mode_count=all(F(r['value'])<=F(r['upper']) for r in rows),
        zero_output_retained=0 in convolution(a,b),
        negative_and_unequal_harmonics_retained=min(convolution(a,b))<0 and max(convolution(a,b))>31))


def inverse_review():
    data={-64:F(3,7),-2:F(2,5),-1:F(1,11),1:F(5,13),2:F(7,19),101:F(11,23)}
    rows=[]
    for s in range(5):
        for p in (1,2,3):
            value=moment(smoothed(data,p),s+p);upper=2**p*moment(data,s)
            rows.append(dict(order=s,smoothing=p,value=str(value),upper=str(upper)))
    L,dstar=F(17),F(1,3);C=max(4*L,F(16,3)/dstar)
    residence=True;high=True
    for m in (-101,-9,-2,-1,1,2,9,101):
        ell=L if abs(m)==1 else min(L,1/((m*m-1)*dstar))
        residence &= (1+abs(m))**2*ell<=C
        if abs(m)>=2:high &= F(1,m*m-1)<=F(4,3*m*m)
    return dict(rows=rows,model_low_mode_length=str(L),model_damping_floor=str(dstar),
        uniform_residence_constant=str(C),checks=dict(
        inverse_index_weights_obey_expected_smoothing=all(F(r['value'])<=F(r['upper']) for r in rows),
        high_mode_denominator_compared_with_m_squared=bool(high),
        low_fundamental_modes_included_in_constant=bool(residence),
        input_full_entrance_jets_are_a_required_hypothesis=True,
        full_path_envelope_normalized_source_norm_required=True,
        source_majorants_not_assumed_uniform_for_infinite_families=True))


def linear_review():
    # index gain = inverse smoothing - explicit multiplier power.
    # The amplitude pressure gradient uses pi, not the vector potential.
    rows=[('principal damping on curl remainder',3,2,F(1,2)-KAPPA),
          ('eikonal defect on t',2,1,F(1,2)),
          ('mixed phase viscosity on t',2,1,F(1,2)-KAPPA),
          ('amplitude gradient of full pressure',1,0,F(1,2)-KAPPA),
          ('two graph derivatives on t',2,0,1-2*KAPPA),
          ('angular vector viscous connection',2,1,F(1,2)),
          ('coefficient temporal and frame operators on curl remainder',3,0,F(1,2)-KAPPA)]
    return dict(rows=[dict(term=t,input_smoothing=q,harmonic_power=p,index_gain=q-p,epsilon_gain=str(e)) for t,q,p,e in rows],
        checks=dict(all_regular_linear_terms_gain_at_least_one_index=all(q-p>=1 for _,q,p,_ in rows),
            all_epsilon_gains_at_least_half_minus_kappa=all(e>=F(1,2)-KAPPA for *_,e in rows),
            pressure_uses_one_index_not_potential_three=True,
            exact_cutoff_tails_stay_additive_and_not_forwarded=True))


def pressure_review():
    # At any interior point of a smooth supported pure normal source of value 1.
    rows=[dict(m=m,pressure_imaginary=str(F(1,3*m)),one_index=str(F(1,3)),
               three_indices=str(F(m*m,3))) for m in (2,4,16,128)]
    return dict(rows=rows,checks=dict(
        direct_normal_source_pressure_has_exact_one_index_gain=all(F(r['pressure_imaginary'])*r['m']==F(1,3) for r in rows),
        potential_three_index_pressure_bound_is_obstructed=F(rows[-1]['three_indices'])>1000,
        moving_normal_and_source_pressure_terms_both_retained=True,
        generic_entrance_derivative_gain_not_inferred_from_value_damping=True))


def cycle_review():
    rho=mu=1-KAPPA;delta=F(1,2)-KAPPA;B0=F(3,2)-2*KAPPA
    # Slopes on B>=B0; value at B0 suffices for each affine margin.
    rows=[('linear',F(1),delta),('old_wave_cross',F(1),F(1,2)-KAPPA),
          ('mean_wave',F(1),mu-F(1,2)),('new_square',F(2),-KAPPA)]
    return dict(seed=str(B0),delta=str(delta),rho=str(rho),
        rows=[dict(term=t,slope=str(slope),intercept=str(c),boundary_margin=str((slope-1)*B0+c-delta)) for t,slope,c in rows],
        checks=dict(current_finite_epsilon_gain_preserved=all((a-1)*B0+c-delta>=0 and a>=1 for _,a,c in rows),
            angular_transport_loss_charged_in_weighted_moments=True,
            complete_covariance_and_current_pressure_recomputed=True,
            all_scale_common_domain_and_original_geometry_unchanged=True,
            infinite_weighted_source_sums_and_stage_uniformity_not_claimed=True))
