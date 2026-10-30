"""Guarded finite-head to ORIGINAL frame constants; no selector mutation."""
from fractions import Fraction as Q
from pathlib import Path
import importlib.util

HERE=Path(__file__).resolve().parent

def load(name,path):
    s=importlib.util.spec_from_file_location(name,path)
    m=importlib.util.module_from_spec(s);s.loader.exec_module(m);return m

head=load('q22_background',HERE.parent/'core_background_finite_head_budget/budgets.py')
frame=load('q22_frame',HERE.parent/'frame_column_foundation/frame_budget.py')

def required(background,Dbox):
    Dbox=head.nonnegative(Dbox)
    lead=background.leading_bound(2);b=background.bound(1)
    return {'Db':3*Dbox*lead+b['error_FG'],'Mz':b['full_FG']}

def join(background,original):
    p=frame.inputs(**original)
    if (background.Rmin,background.Rmax)!=(p['Rmin'],p['Rmax']):
        raise ValueError('The same closed radial chart is required')
    need=required(background,p['Dbox'])
    gates={k:need[k]<=p[k] for k in ('Db','Mz')}
    return {'required':need,'original':p,'dominance':gates,
            'original_constants_unchanged':True,
            'original_value_gates_still_required':True}

def require_reuse(result):
    p=frame.inputs(**result['original'])
    if not all(result['required'][k]<=p[k] for k in ('Db','Mz')):
        raise ValueError('Original uniform primitive bounds not justified by this recipe')
    b=frame.budgets(p);H=b['CE']+b['Cd']
    return {'CE':b['CE'],'Cd':b['Cd'],'H':H,
            'Kframe':(1+b['S0'])*max(Q(1),p['cmax'])/min(Q(1),p['cmin']),
            'power_of_3':-(-H.numerator//H.denominator),
            'scale_degree':0,'original_value_gates_still_required':True}

def reserve_log_lower(ell,h,reserve,M):
    vals=tuple(head.exact(x) for x in (ell,h,reserve,M));ell,h,reserve,M=vals
    if ell<1 or h<=0 or reserve<0 or M<0:raise ValueError('Invalid fixed reserve inputs')
    return ell*ell-(h*reserve+2*M)*ell

def record():
    return {'status':'compactness-quantifier-corrected-with-guarded-uniform-frame-recipe',
            'required_Db':'3 Dbox Jlead_2 + Jerr_1','required_Mz':'Jfull_1',
            'uniform_frame_condition':'log Kframe + CE + Cd <= A + B log S',
            'ODE_normalized_homogeneous':'exp(v)',
            'ODE_normalized_zero_entry':'exp(v)-1',
            'finite_band_constants_imply_polynomial_control':False,
            'fixed_epsilon_reserve_repairs_model':False,
            'original_selector_mutated':False,
            'actual_original_constant_dominance_proved':False,
            'actual_extended_source_bounds_proved':False,
            'unconditional_residual_or_blowup':False}
