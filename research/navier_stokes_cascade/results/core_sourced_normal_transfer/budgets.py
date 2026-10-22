"""Original normal floor in the unchanged sourced inverse; no band selection."""
from functools import lru_cache
from pathlib import Path
import importlib.util

HERE=Path(__file__).resolve().parent


def load(name,path):
    s=importlib.util.spec_from_file_location(name,path);m=importlib.util.module_from_spec(s);s.loader.exec_module(m);return m


normal=load('sourced_transfer_normal',HERE.parent/'core_normal_length_closure/budgets.py')
phase=load('sourced_transfer_phase',HERE.parent/'core_phase_coefficient_budget/budgets.py')
inverse=phase.inverse
P,order,exact=inverse.P,inverse.order,inverse.exact


def forced_from_original(p,selector_inputs,ell,*,F,G,btilde,path,source,tangent,length,harmonic=1):
    """Conditional bounds on the SAME full path and frozen label.

    F,G,btilde,path,source and the tangent propagator remain supplied upper
    bounds, not field enclosures certified by this function. Recompute the
    original normal certificate; do not accept a caller-written gate dict.
    """
    p=normal.frame_inputs(p)
    certificate=normal.verify_existing(p,selector_inputs,ell)
    if not certificate['conditional_normal_floor_implication_passes']:
        raise ValueError('Supplied original band and original length premises do not establish the normal floor')
    d=normal.frame.primitive(p)
    # Original L>=1 and original carrier triangle inequalities. No extra
    # comparison with the optional phase-beta or sqrt(lambda_max) recipe.
    b=phase.PhaseBudget(F=F,G=G,btilde=btilde,path=path,
        Rmin=p['Rmin'],Rmax=p['Rmax'],Lmin=1,u=p['u'],lambda_max=p['lam_max'],
        carriers=dict(theta=d['Pmax'],axial=d['Zmax'],radial=p['bmax']*p['u']/2))
    f=inverse.ForcedBudget(normal=b.normal,moving_normal=b.moving_normal,matrix=b.matrix,
        source=source,profile_rate=b.profile_rate,tangent=tangent,length=length,
        inverse_normal_square=certificate['inverse_normal_square'],harmonic=harmonic)
    return b,f,certificate


def degree_bound(kind,n):
    """Degree in the scalar reciprocal floor q, with all other inputs fixed."""
    order(n)
    if kind in ('square','pressure_row','damping'):return 0
    if kind in ('reciprocal','lift','undamped','normal_projector','generator'):return n+1
    if kind in ('projection','projected_source'):return n+1 if n else 0
    if kind=='ambient':return 2
    if kind=='velocity':return 4*n+2
    if kind in ('pressure','potential'):return 4*n+3
    raise ValueError('Known inverse coefficient kind required')


def transfer_factor(kind,n,q_old,q_new):
    q_old,q_new=map(exact,(q_old,q_new))
    if min(q_old,q_new)<=0:raise ValueError('Positive reciprocal bounds required')
    return max(1,q_new/q_old)**degree_bound(kind,n)


class FormalFloorBudget(inverse.ForcedBudget):
    """DIAGNOSTIC ONLY: use the existing polynomial variable as formal q.

    All non-floor input polynomials must first be evaluated at one S. The
    production constructor still requires an exact rational reciprocal.
    Only the reciprocal base case needs a symbolic lift; every positive
    recurrence at positive order is inherited verbatim.
    """
    def __init__(self,base,S):
        cv=lambda name:lambda n:P(base.data(name,n).evaluate(S))
        super().__init__(normal=cv('normal'),moving_normal=cv('moving_normal'),
            matrix=cv('matrix'),source=cv('source'),profile_rate=cv('profile_rate'),
            tangent=P(base.tangent.evaluate(S)),length=P(base.length.evaluate(S)),
            inverse_normal_square=1,harmonic=base.harmonic)
        self.inverse_normal_square=P({1:1})

    @lru_cache(None)
    def reciprocal(self,n):
        order(n)
        return P({1:1}) if n==0 else super().reciprocal(n)
