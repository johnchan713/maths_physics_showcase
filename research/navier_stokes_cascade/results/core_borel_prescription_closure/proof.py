"""Exact algebra of the retained all-order prescription; no scales selected."""
import ast
from fractions import Fraction as Q
from pathlib import Path
import importlib.util

HERE=Path(__file__).resolve().parent


def load(name,path):
    s=importlib.util.spec_from_file_location(name,path);m=importlib.util.module_from_spec(s);s.loader.exec_module(m);return m


tangent=load('prescription_tangent',HERE.parent/'core_tangential_normal_closure/budgets.py')
exact=tangent.exact


def poly(*terms):
    """Sparse polynomial, (coefficient, variable-name tuple) terms."""
    out={}
    for c,variables in terms:
        key=tuple(sorted(variables));out[key]=out.get(key,Q(0))+Q(c)
    return {k:v for k,v in out.items() if v}


def add(a,b):return poly(*[(v,k) for row in (a,b) for k,v in row.items()])


def certificates():
    # g=n*h>0, ell=log(2)>0, v=log(C)>=0, b=log(a_n).
    normalized=poly((1,('g','b')),(-1,('g',)),(-1,('v',)),(-1,('n','ell')))
    separation=poly((1,('b',)),(-1,('previous',)),(-1,('ell',)),(-1,()))
    boundary=poly((1,('b',)),(-1,('domain',)),(-1,('ell',)),(-1,()))
    # Multiply the log-power requirement by gamma^2>0, clearing divisions.
    physical=poly((1,('gamma','gamma','b')),(-1,('gamma','gamma')),
                  (-2,('gamma','v')),(-2,('gamma','n','ell')),(-1,('P','P')))
    return {
        'normalized':dict(premise=normalized,remainder=poly((1,('g',))),
                          target=poly((1,('g','b')),(-1,('v',)),(-1,('n','ell')))),
        'separation':dict(premise=separation,remainder=poly((1,())),
                          target=poly((1,('b',)),(-1,('previous',)),(-1,('ell',)))),
        'boundary':dict(premise=boundary,remainder=poly((1,())),
                        target=poly((1,('b',)),(-1,('domain',)),(-1,('ell',)))),
        'physical_log_power':dict(premise=physical,remainder=poly((1,('gamma','gamma'))),
                                 target=poly((1,('gamma','gamma','b')),(-2,('gamma','v')),
                                             (-2,('gamma','n','ell')),(-1,('P','P')))),
    }


def verify_certificate(c):
    # All variable values are NONNEGATIVE; g and gamma are strictly positive.
    return add(c['premise'],c['remainder'])==c['target'] and bool(c['remainder']) and all(v>=0 for v in c['remainder'].values())


def source_rule_checks(source):
    """Identify the historical helper's expression, NOT its rounding behavior.

    The continuum result concerns the real formula in the frozen README.
    The helper evaluates manufactured data in mpmath and is not an interval
    implementation for an arbitrary actual infinite coefficient sequence.
    """
    tree=ast.parse(source)
    fn=next(x for x in tree.body if isinstance(x,ast.FunctionDef) and x.name=='cutoff_log_scale')
    expressions={ast.dump(x,include_attributes=False) for x in ast.walk(fn)}
    expected={
        'all_normalized_rows_through_n_required':'len(log_constants) != n+1',
        'fixed_gain_is_n_times_h':'mp.mpf(h.numerator)*n/h.denominator',
        'original_max_includes_each_row_separation_boundary_and_unit_padding':
            '1+max(previous+mp.log(2), boundary_log+mp.log(2), max((v+n*mp.log(2))/gain for v in values))',
    }
    return {name:ast.dump(ast.parse(expr,mode='eval').body,include_attributes=False) in expressions for name,expr in expected.items()}


def real_stage_residuals(g,n,ell,v,previous,domain,b):
    """Verify supplied rational stand-ins for real log variables; never choose b.

    ell stands for log(2) only in the theorem, not in these algebra controls.
    """
    g,ell,v,previous,domain,b=map(exact,(g,ell,v,previous,domain,b))
    if type(n) is not int or n<1 or g<=0 or ell<=0 or min(v,previous,domain,b)<0:raise ValueError('Valid positive gain/index and nonnegative real-log stand-ins required')
    slack=dict(normalized=b-1-(v+n*ell)/g,separation=b-1-previous-ell,boundary=b-1-domain-ell)
    return dict(prescription_slacks=slack,all_supplied_constraints_pass=all(s>=0 for s in slack.values()),
                normalized_log_margin=g*b-v-n*ell,scales_selected_or_changed=False,
                infinite_schedule_verified_by_finite_inputs=False)


def floor_inputs_from_rule(**kwargs):
    """Unfold the retained REAL prescription, relative to coefficient regularity.

    No choice of chi or a_n is made. Any fixed smooth chi and any sequence
    admitted by the original all-order rule obey the same low-order estimate.
    """
    if not all(verify_certificate(c) for c in certificates().values()):raise ValueError('Universal algebra certificate failed')
    return dict(floor_inputs=tangent.comparison(**kwargs),
                independent_diagonal_assumption_needed_for_rule_defined_construction=False,
                diagnostic_cutoff_identity_needed_for_this_floor=False,
                all_order_diagonal_proved_relative_to_original_real_rule_and_full_coefficient_regularity=True,
                actual_independently_fixed_numeric_sequence_identified=False,
                actual_band_or_full_residual_verified=False,
                scope='all choices admitted by the frozen real prescription with their OWN cutoff and coefficient seminorms; not a new scale selection')
