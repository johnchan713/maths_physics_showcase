"""First-real residual rate and stage ledger, not a continuum certificate."""
from fractions import Fraction as Q

def rates(h=Q(1,200),k=Q(1,100000),m=0):
    h,k=Q(h),Q(k)
    if not 0<h<Q(1,100) or not 0<k<Q(1,4) or type(m) is not int or m<0:raise ValueError('Retained parameter ranges required')
    B=1-k;t=3-5*k;nu=3-6*k;K=Q(5,2)+2*h+2*m
    return dict(nonzero_regular=B,mean_rows=t,target=nu,missing_source_gain=nu-B,
                normalized_total=min(B,t),physical_total=h*B-K)

SOURCES={'angular_mean':'final_first_real','axial_regular':'final_first_real',
         'axial_compact':'post_signed_pre_new_fast_first_real',
         'radial_compact':'final_radial_first_real','compensator':'final_first_real'}

def check_ledger(stage,sources,gaussian_kinds):
    if stage!='first_real' or sources!=SOURCES:raise ValueError('Do not transplant first-real tails or swap their source states')
    if set(gaussian_kinds)-{'original_temporal','signed_temporal'}:
        raise ValueError('Spatial, Borel and localization cutoffs need their separate estimates')
    return dict(stage=stage,angular_mean_rate='3-5*kappa',nonzero_regular_rate='1-kappa',
                continuum_inputs_certified=False)

def record():
    return dict(status='first-real-regular-and-additive-residual-ledger',
                fixture_rates={k:str(v) for k,v in rates().items()},source_states=SOURCES,
                ledger=check_ledger('first_real',dict(SOURCES),['original_temporal','signed_temporal']),
                mean_row_gain_is_complete_residual_gain=False,
                unconditional_pulse_residual_or_blowup=False)
