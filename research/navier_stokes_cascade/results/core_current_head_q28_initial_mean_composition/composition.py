"""Two actual finite mean passes; positive norm recipes, not field evaluation."""
from pathlib import Path
import importlib.util
from fractions import Fraction as Q
HERE=Path(__file__).resolve().parent
spec=importlib.util.spec_from_file_location('q28_source',HERE.parent/'core_signed_prefix_source_budget/budgets.py')
source=importlib.util.module_from_spec(spec);spec.loader.exec_module(source)
P=source.P
KEYS={'fixed','wave','linear','KT1','KA1','time_inverse','slow'}

def compose(primary,signed,trace=None):
    """Initial means are derived, never supplied independently or reset later.

    Each dictionary has its own correctly padded operator archive. The signed
    wave is already formed from the PRE-signed target; no fixed point is solved.
    """
    if set(primary)!=KEYS or set(signed)!=KEYS:raise ValueError('Exact complete provider sets required; no initial callback')
    for data in (primary,signed):
        if not all(callable(v) for v in data.values()):raise TypeError('Finite-order provider callbacks required')
    initial=source.PrefixBudget(**primary,initial=lambda n:{k:P() for k in source.STATE},trace=trace)
    final=source.PrefixBudget(**signed,initial=initial.final,trace=trace)
    return initial,final

def poles(primary_wave,signed_wave):
    if any(type(p) is not int or p<0 for p in (primary_wave,signed_wave)):
        raise ValueError('Nonnegative integer poles required')
    p=primary_wave;r=max(2*p,signed_wave)
    return {'initial_input':p,'initial_mean':2*p,'signed_input':r,'final_mean':2*r,
            'source':3*r,'pressure_archive':sorted({2*p,4*p,2*r,4*r}),
            'fast_archive':sorted({2*p,2*r})}

def record():
    return {'result':'conditional-primary-to-signed-mean-composition',
        'initial_mean_callback_removed':True,'initial_old_means_zero':True,
        'signed_old_means_reset':False,'fresh_pressure_and_all_scalars':True,
        'primary_wave_order':'j+27','primary_root_order':'j+28 before target expansion',
        'signed_wave_order':'j+15','signed_root_order':'j+16 before target expansion',
        'source_pole':'3*max(2*p_primary,p_signed)',
        'fixed_map_and_background_premises_discharged':False,
        'complete_residual_or_blowup':False}
