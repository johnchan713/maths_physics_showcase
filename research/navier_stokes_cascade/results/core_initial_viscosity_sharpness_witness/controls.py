from fractions import Fraction as Q
from witness import Mode,exponents,local_witness,logical_consequence


def review():
    checks={};k=Q(1,100000);w=local_witness(k);o=w['orders'];c=logical_consequence(k)
    checks['principal_cross_realizes_nonzero_covariance_mode_at_alpha']=w['covariance']['mode']==Mode(0,1) and w['covariance']['order']==o['alpha']
    checks['one_graph_derivative_attains_complete_source_order_D']=w['source']['mode']==Mode(1,0) and w['source']['order']==o['D']
    checks['connection_part_has_strict_kappa_source_margin']=w['source']['connection_order']-o['D']==k
    checks['temporal_inverse_preserves_D_and_cancels_source_exactly']=w['fast']['order']==o['D'] and w['fast']['mode']==Mode(0,-1) and w['exact_fast_time_cancellation']==Mode()
    checks['two_graph_viscosity_derivatives_attain_nu']=w['viscosity']['mode']==Mode(0,-1) and w['viscosity']['order']==o['nu']
    checks['one_and_zero_graph_laplacian_parts_are_strictly_higher']=w['viscosity']['one_graph_order']-o['nu']==k and w['viscosity']['zero_graph_order']-o['nu']==2*k
    checks['witness_respects_zero_Haar_requirement']=w['all_modes_zero_Haar']
    checks['both_missing_gains_can_be_zero_in_same_local_jet']=c['source_rate_attained'] and c['radial_viscosity_rate_attained'] and w['strict_margin']==0
    checks['positive_gain_is_not_a_consequence_of_inherited_size_hypotheses']=not c['inherited_upper_bounds_force_positive_gain']
    checks['actual_special_identity_is_left_open']=not c['actual_pulse_special_cancellation_excluded'] and bool(c['required_actual_identity'])
    try:Mode(1,0)+1
    except Exception:checks['mode_algebra_accepts_only_exact_scalars']=False
    else:checks['mode_algebra_accepts_only_exact_scalars']=True
    return dict(checks={name:bool(value) for name,value in checks.items()},
      parent_commit='a05a5a2a91220e0e080f54e558b8c2fd73e36cba',
      scope=('exact local active-interior Fourier-jet sharpness witness for the inherited source and radial-viscosity '
             'bounds; not an actual-pulse lower bound, high-jet enclosure, complete residual estimate or blow-up proof'))
