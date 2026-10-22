"""Exact denominator transfer and full-pressure controls; manufactured inputs."""
from fractions import Fraction as Q
from functools import lru_cache
from budgets import HERE,load,normal,phase,inverse,P,forced_from_original,degree_bound,transfer_factor,FormalFloorBudget


def fixture():
    p=normal.frame.inputs(bmin=Q(1,100),Dbox=6,Mz=1)
    c=normal.normal_length_certificate(p);sp=normal.selector.fixture();L=c['length_rule_floor']
    sp['length']=Q(-(-L.numerator//L.denominator))
    ell=normal.selector.select_tail(sp)['ell_tail'] # manufactured fixture ONLY
    kwargs=dict(F=lambda n:P(n+2),G=lambda n:P(n+3),btilde=lambda n:P(n+1),
        path=lambda n:P({1:1}) if n==0 else P(int(n==1)),source=lambda n:P(n+1),
        tangent=P(3),length=P({1:1}))
    return p,sp,ell,kwargs


def bridge_review():
    p,sp,ell,kw=fixture();b,f,c=forced_from_original(p,sp,ell,**kw)
    beta=p['lam_min']/(4*(1+p['u']**2)**2)
    d=normal.frame.primitive(p)
    checks=dict(
        same_original_normal_bound_enters_inverse=f.inverse_normal_square==4/p['bmin']**2==c['inverse_normal_square'],
        stronger_optional_phase_beta_bound_is_not_inferred=f.inverse_normal_square>4/beta,
        original_carrier_bounds_without_artificial_cross_gate=(b.p,b.pz,b.x0)==(d['Pmax'],d['Zmax'],p['bmax']*p['u']/2),
        original_shell_tilt_growth_and_L_ge_1_retained=(b.Rmin,b.Rmax,b.u,b.lam,b.Lmin)==(p['Rmin'],p['Rmax'],p['u'],p['lam_max'],1),
        background_and_path_callbacks_are_same_objects=all(getattr(b,k) is kw[k] for k in ('F','G','btilde','path')),
        unprojected_source_and_tangent_length_are_same_objects=f.source is kw['source'] and f.tangent is kw['tangent'] and f.length is kw['length'],
        no_actual_input_or_full_propagator_promotion=not c['actual_input_enclosures_or_arbitrary_historical_band_certified'] and not c['full_propagator_or_residual_certified'],
    )
    for n in range(4):
        checks[f'order{n}_unchanged_full_phase_callbacks_reach_inverse']=all(
            getattr(f,k)(n).record()==getattr(b,k)(n).record() for k in ('normal','moving_normal','matrix','profile_rate'))
    return dict(checks=checks,old_optional_q=str(4/beta),original_q=str(f.inverse_normal_square),existing_ell=str(ell))


def polynomial_review():
    checks={};rows=[]
    for S in (1,3):
        base=inverse.fixture();formal=FormalFloorBudget(base,S)
        for kind in ('reciprocal','lift','undamped','normal_projector','generator','projection','projected_source','velocity','pressure','potential'):
            for n in range(4):
                poly=getattr(formal,kind)(n);D=degree_bound(kind,n)
                checks[f'S{S}_{kind}_{n}_positive_formal_q_degree']=all(v>=0 for v in poly.c.values()) and max(poly.c,default=0)<=D
                # Nondegenerate manufactured callbacks attain the announced
                # highest degree for outputs, detecting accidental undercounts.
                if kind in ('velocity','pressure','potential'):
                    checks[f'S{S}_{kind}_{n}_top_degree_is_present']=poly.c.get(D,0)>0
                for q in (Q(1,3),Q(2),Q(8)):
                    f=inverse.fixture();f.inverse_normal_square=q
                    checks[f'S{S}_{kind}_{n}_q{q}_formal_matches_unchanged_inverse']=sum(v*q**j for j,v in poly.c.items())==getattr(f,kind)(n).evaluate(S)
                rows.append(dict(S=S,kind=kind,order=n,degree=max(poly.c,default=0),bound=D))
    for q0,q1 in ((Q(2),Q(8)),(Q(8),Q(2)),(Q(1,3),Q(2))):
        a=inverse.fixture();a.inverse_normal_square=q0;b=inverse.fixture();b.inverse_normal_square=q1
        for n in range(4):
            for kind in ('velocity','pressure','potential'):
                pa,pb=getattr(a,kind)(n),getattr(b,kind)(n)
                factor=transfer_factor(kind,n,q0,q1)
                checks[f'{kind}_{n}_q{q0}_to{q1}_coefficientwise_scale_transfer']=all(pb.c.get(j,0)<=factor*pa.c.get(j,0) for j in set(pa.c)|set(pb.c))
                checks[f'{kind}_{n}_q{q0}_to{q1}_scale_powers_unchanged']=set(pa.c)==set(pb.c)
    return dict(checks=checks,rows=rows,scope='formal positive-polynomial controls; analytical induction covers every fixed order and all S, no actual field measurement')


def pressure_review():
    checks={}
    # Constant n=e3, K=nprime=0, f=g*e3. Pi*f=0 gives t=0,
    # but exact ODE pressure is i*g/(km), not zero. Check the majorant
    # retains the unprojected source term even with velocity set to zero.
    f=inverse.ForcedBudget(normal=lambda n:P(int(n==0)),moving_normal=lambda n:P(),
        matrix=lambda n:P(),source=lambda n:P(7*int(n==0)),profile_rate=lambda n:P(),
        tangent=P(1),length=P(1),inverse_normal_square=1,harmonic=2)
    f.velocity=lambda n:P()
    checks['pure_normal_source_nonzero_pressure_with_zero_velocity']=f.pressure(0).evaluate(1)==Q(7,2)
    a=inverse.fixture();b=inverse.fixture();b.moving_normal=lambda n:P()
    checks['moving_normal_survives_pressure_and_ambient']=a.pressure(1).evaluate(1)>b.pressure(1).evaluate(1) and a.ambient().evaluate(1)>b.ambient().evaluate(1)
    cutoff=lambda n:P(n+1)
    for kind in ('pressure','potential'):
        checks[kind+'_same_cutoff_derivative_products_survive']=a.cut(kind,cutoff,2).evaluate(1)>cutoff(0).evaluate(1)*getattr(a,kind)(2).evaluate(1)
    return dict(checks=checks,scope='exact algebra and majorants; full Pi1/Pi2 pressure memory and fresh axial means remain required')


@lru_cache(None)
def review():
    parts=dict(original_floor_bridge=bridge_review(),denominator_polynomials=polynomial_review(),complete_pressure=pressure_review())
    return dict(checks={name+':'+k:v for name,r in parts.items() for k,v in r['checks'].items()},parts=parts)
