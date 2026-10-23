"""Conditional finite-source smooth closure, not actual norm certification."""
from fractions import Fraction as Q
from pathlib import Path
import importlib.util

HERE=Path(__file__).resolve().parent
def load(name,path):
    spec=importlib.util.spec_from_file_location(name,path)
    out=importlib.util.module_from_spec(spec);spec.loader.exec_module(out);return out
trace=load('finite_support_trace_parent',HERE.parent/'core_compact_trace_gate/bounds.py')
source=load('finite_support_edge_source',HERE.parent/'source_envelope_edge_audit/source_bounds.py')

def order(n):
    if type(n) is not int or n<0:raise ValueError('Nonnegative integer order required')
    return n

def flat_radius_recipe(inverse_power,vanishing_order):
    """x^-p exp(-c/x^2) <= L! c^-L x^(2L-p); no giant L! built."""
    p,K=order(inverse_power),order(vanishing_order);L=(p+K)//2+1
    return dict(inverse_power=p,requested_vanishing_order=K,series_index=L,
        resulting_power=2*L-p,strict_margin=2*L-p-K,
        constant='L! * c^(-L), original positive c retained',materialized=False)

def signed_half_weight(n):
    order(n)
    return dict(output_jet=n,formula='delta_a=(Hc^-1 Xi)/(2 epsilon ac), ac=sqrt(y)',
        denominator_gate='y>=c*zeta/s0>0 on the interior, same ORIGINAL root',
        derivative_gate='each required D^I Xi has full zeta; D^I y has zeta/s0 times finite edge poles',
        derivative_weight=Q(1,2),lower_exponential_left='C_setup^-40/2',lower_exponential_right='2',
        reciprocal_root_alone_extended=False,product_has_flat_zero_extension_conditionally=True,
        new_root_total_wave_or_finer_target_callback_used=False)

def signed_support_threshold(theta):
    if isinstance(theta,bool) or not isinstance(theta,(int,Q)):raise TypeError('Exact positive target weight required')
    theta=Q(theta)
    if theta<=0:raise ValueError('Positive target exponential weight required')
    power=theta-Q(1,2)
    return dict(target_weight=theta,remaining_signed_weight=power,
        flat_extension_for_all_inputs_in_weighted_class=power>0,
        equality_requires_additional_vanishing=power==0,
        actual_canonical_target_failure_claimed=False)

def support_recipe(n,*,physical_power=None):
    order(n)
    req=trace.A1_trace_request(physical_power=physical_power)
    return dict(requested_mixed_order=n,scope='finite CURRENT prefix only, same physical labels/covers',
        derivation='weighted edge jets -> complete cut waves -> full covariance -> compact actual means/pressure -> Q Omega,Q K1,Q Phi',
        node_contracts=[
            'original smooth fixed background/profile and current normal/root gaps on unchanged accepted bands',
            'EVERY earlier signed target and source has its retained weighted mixed jet bounds on the whole closed lifted path',
            'same temporal/transverse cutoff on full potential AND pressure, complete cylindrical curl',
            'zero entrance forced solve and smooth whole-path coefficients; raw outgoing response is not assumed zero',
            'shifted T_e with SAME compact chi, full potential hull and edge-directed bounds; actual gd-A1 gd retained',
            'original finite smooth reserved response, same moving rho and full current pressure',
            'physical global periodic pullbacks/copy sums before charts; no independent patch formulas or inverse branch'],
        conditional_all_mixed_trace_orders_zero=True,
        physical_trace_request=req,
        huge_independent_trace_samples_needed_if_global_input_contract_proved=False,
        complete_older_raw_norm_and_operator_request_graph_retained=True,
        actual_input_contract_or_same_band_gaps_verified=False)

def remaining_gate():
    return dict(target='for EVERY retained target Xi_j and requested I: |D^I Xi_j| <= C_jI epsilon^alpha_j zeta ya^-p_jI yb^-q_jI',
        source='full lifted regular source f_j/P_j with retained Gaussian quotient and mixed bounds, no additive tail sent to inverse',
        root='same original y_sigma >= c_sigma zeta/s0 and its weighted jets',
        floors='same-band full normal floor, root/cone gaps and target freezing, including ell^6>=8 Cfreeze/gamma0 when that frozen-direction argument is used',
        mean='whole shifted/ordinary primitive moment contracts and finite reserved hull; moving rho P and Jz+c_rho P retained',
        quantitative='original coefficient/profile/curl/kernel/domain constants at all nested raw requests; no new first band',
        actual_canonical_suprema_evaluated=False)

def certify_actual_inputs():
    raise NotImplementedError('Actual weighted mixed target/source jets and unchanged-band normal/root/freezing gates remain unverified; finite smooth-closure grammar is conditional')

def certificate():
    return dict(parent_commit='12603481241ee73dcfb3df6b773f8c3687b0dab9',
        result='conditional all-order mixed trace closure of the retained finite construction, reducing the seam gate to original weighted whole-path input contracts',
        sharp_weight_threshold='target zeta^theta gives signed zeta^(theta-1/2); theta>1/2 suffices, smooth compact target alone does not; fixed root comparison fails at/below equality',
        support=support_recipe(0,physical_power=0),remaining=remaining_gate(),
        pressure='full n.K.t-nprime.t+n.f; same K0/K1, rho P, Pi1/Pi2 and all six covariance entries; Q removes only auxiliary-independent assigned stress/compensator',
        caution='global finite recipe closure is proved under the displayed input contract; norm containers and protocol names do not certify that contract',
        full_regular_order='2-3kappa',actual_all_scale_inputs_absolute_decay_infinite_field_or_blowup_proved=False)
