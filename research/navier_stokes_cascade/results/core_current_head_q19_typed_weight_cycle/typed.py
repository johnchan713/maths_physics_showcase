"""Radial weight invariants only; no epsilon-rate or inverse certification."""
from fractions import Fraction as Q
HALF=Q(1,2)

def cycle(wave=HALF,mean=Q(1),stress=Q(1),target_floor=Q(1)):
    wave,mean,stress,target_floor=map(Q,(wave,mean,stress,target_floor))
    if wave<HALF or mean<1 or stress<1 or target_floor<1:
        raise ValueError('half-wave/full-mean/full-target premises required')
    covariance=2*wave
    # Exact mean equations: linear mean, quadratic mean, covariance, assigned stress.
    rows=min(mean,2*mean,covariance,stress)
    pressure=min(mean,2*mean,covariance)
    target=min(rows,pressure,target_floor)
    signed=target-HALF
    new_wave=min(wave,signed)
    changed_covariance=min(wave+signed,2*signed)
    new_mean=min(mean,changed_covariance,rows)
    harmonic_source=min(new_wave,2*new_wave,new_mean+new_wave)
    return dict(covariance=covariance,mean_rows=rows,mean_pressure=pressure,
                signed_target=target,signed_wave=signed,
                covariance_change=changed_covariance,mean=new_mean,
                wave=new_wave,harmonic_source=harmonic_source,
                harmonic_pressure=harmonic_source)

def multiply(a,b):
    out={}
    for k,x in a.items():
        for l,y in b.items():out[k+l]=out.get(k+l,Q(0))+x*y
    return {k:v for k,v in out.items() if v}

def derivative(a):return {k:Q(k)*v for k,v in a.items() if k} # theta derivative / i

def record():
    return {'radial_invariant':{k:str(v) for k,v in cycle().items()},
            'angular_mean_distinct_from_auxiliary_Haar':True,
            'compact_maps_only_receive_full_weight_mean_inputs':True,
            'harmonic_inverse_requires_whole_path_f_over_P_jets':True,
            'full_pressure_can_have_half_weight':True,
            'conditional_radial_cycle_closed':True,
            'actual_inverse_root_and_rate_contracts_certified':False,
            'unconditional_residual_or_blowup':False}
