"""Exact controls for support localization; no reference inverse is synthesized."""
from fractions import Fraction as Q

def length_in_U(lo,hi,edge):
    if any(isinstance(x,(float,bool)) or not isinstance(x,(int,Q)) for x in (lo,hi,edge)):
        raise TypeError('Exact bounds required')
    lo,hi,edge=map(Q,(lo,hi,edge))
    if not 0<=lo<=hi or edge<0:raise ValueError('Ordered positive radial interval required')
    return max(Q(0),min(hi,edge)-lo)

def record():
    return {'status':'principal-comparison-error-supported-on-positive-order-velocity-hull',
      'interval_controls':[str(length_in_U(2,8,e)) for e in (1,2,5,8,10)],
      'common_target_required':True,'common_reference_domain_required':True,
      'high_order_leading_integral_localized':False,'failed_branch_improved':False,
      'pressure_output_support_inferred':False,'physical_wave_replaced':False,
      'unconditional_pulse_residual_or_blowup':False}
