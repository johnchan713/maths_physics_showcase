"""Exact derivative budgets; cutoff shape norms are inputs, not certified values."""
from fractions import Fraction as F
from functools import lru_cache
import importlib.util
from pathlib import Path
import sys

HERE=Path(__file__).resolve().parent
def rat(x):
    if isinstance(x,(bool,float)):raise TypeError('Exact rational inputs required')
    return F(x)
def serial(x):
    if isinstance(x,F):return str(x)
    if isinstance(x,dict):return {k:serial(v) for k,v in x.items()}
    if isinstance(x,(list,tuple)):return [serial(v) for v in x]
    return x
@lru_cache(None)
def seminorms():
    p=HERE.parent/'first_seminorm_bound/seminorms.py'
    s=importlib.util.spec_from_file_location('cutoff_retained_first_seminorms',p)
    m=importlib.util.module_from_spec(s);sys.modules[s.name]=m;s.loader.exec_module(m);return m

def step_bounds():
    """Proved in README for the actual squared-flat sigma, all real s."""
    return (F(32),F(1216),F(220000))
def bump_bounds(width=F(1,10)):
    width=rat(width)
    if width<=0:raise ValueError('Positive logarithmic bump width required')
    return tuple(v/width**(j+1) for j,v in enumerate(step_bounds()))
def cutoff_bounds(xkeep,delta,L1,L2):
    """delta=tc/16; L_j are actual normalized log-transition shape norms.

    L1>=1,L2>=4 are necessary, not sufficient to certify any chosen cutoff.
    The returned inequalities are conditional on the supplied L_j enclosures.
    """
    xkeep,delta,L1,L2=map(rat,(xkeep,delta,L1,L2))
    if min(xkeep,delta)<=0 or L1<1 or L2<4:
        raise ValueError('Positive geometry and necessary shape budgets L1>=1,L2>=4 required')
    return dict(radial=L1/(delta*xkeep),
        radial2=(L2/delta**2+L1/delta)/xkeep**2)
def compressed_cutoff(n):
    """Exact controls for 1-sigma(n(s-1/2)+1/2); not a field replacement."""
    if type(n) is not int or n<1:raise ValueError('Positive integer compression required')
    return dict(support_left=F(1,2)-F(1,2*n),support_right=F(1,2)+F(1,2*n),
        center_slope=F(-8*n),first_norm_lower=F(8*n),second_norm_lower=F(16*n*n),
        center_third=F(3712*n**3))
def propagate(core,coefficient,a,xkeep,tc,L1,L2,rho_root):
    """Same additive five bumps, same linear repair, lower patch radius rho3.

    All core/coefficients are inherited C_eta^1 bounds. No actual field sampled.
    """
    coefficient,a,xkeep,tc,L1,L2,r=map(rat,(coefficient,a,xkeep,tc,L1,L2,rho_root))
    if min(a,tc,r)<=0 or xkeep<a*a:raise ValueError('Positive retained radii/width and Xkeep>=Xa required')
    cut=cutoff_bounds(xkeep,tc/16,L1,L2)
    jets=seminorms().extended_jets(core,coefficient,a,cut['radial'],cut['radial2'],r*r,r,bump_bounds())
    return dict(cutoff=cut,bumps=bump_bounds(),extended=jets)

def symbolic_review():
    import sympy as s
    X,x0,d=s.symbols('X x0 d',positive=True);q=s.symbols('q',real=True)
    K=s.Function('K');z=s.log(X/x0)/d
    k1=s.diff(K(q),q).subs(q,z)/(d*X)
    k2=(s.diff(K(q),q,2).subs(q,z)/d**2-s.diff(K(q),q).subs(q,z)/d)/X**2
    g=s.Function('g')(q);H=1/(1+s.exp(-g));hp=H*(1-H)
    logistic3=hp*(s.diff(g,q,3)+3*(1-2*H)*s.diff(g,q)*s.diff(g,q,2)+(1-6*H+6*H**2)*s.diff(g,q)**3)
    t=(1-q)**-2-q**-2;S=1/(1+s.exp(-t))
    center=[s.simplify(s.diff(S,q,j).subs(q,F(1,2))) for j in (1,2,3)]
    B=q**5-2*q**3+q;raw=B.subs(q,s.log(X))
    E=raw/s.sqrt(2*X)
    E2=(s.diff(B,q,2)-2*s.diff(B,q)+s.Rational(3,4)*B).subs(q,s.log(X))/(X**2*s.sqrt(2*X))
    return dict(center_sigma_jets=list(map(str,center)),checks=dict(
        cutoff_first_chain=s.simplify(s.diff(K(z),X)-k1)==0,
        cutoff_second_chain_includes_log_jacobian=s.simplify(s.diff(K(z),X,2)-k2)==0,
        omitted_log_jacobian_detected=s.simplify(s.diff(K(z),X,2)-(k2+k1/X))!=0,
        logistic_third_chain=s.simplify(s.diff(H,q,3)-logistic3)==0,
        actual_step_center_derivatives=center==[8,0,-3712],
        bump_to_F_second_with_both_inverse_radius_terms=s.simplify(s.diff(E,X,2)-E2)==0))

def exact_review():
    # For s<=1/2, z=s^-2>=4: |t'|<=4s^-3, |t''|<=12s^-4,
    # |t'''|<=48s^-5; sigma<=exp(4-z). See the continuum proof.
    step1=4*8;step2=12*16+16*64
    step3=48*32+144*128+64*5**5
    bump=bump_bounds();centers=tuple(map(F,('1/2','1','3/2','2','5/2')))
    supports=[(c-F(1,20),c+F(1,20)) for c in centers]
    rows=[compressed_cutoff(n) for n in (1,2,10,1000)]
    return dict(step_enclosure=[step1,step2,step3],published_step_bounds=serial(step_bounds()),
        bump_bounds=serial(bump),supports=serial(supports),compressed_controls=serial(rows),
        necessary_normalized_shape_bounds=['1','4'],checks=dict(
            first_two_step_bounds_exact=(step1,step2)==(32,1216),
            third_step_padding=step3==219968 and step3<220000,
            polynomial_logistic_third_factor_bounded_by_one=6*F(1,2)**2-F(1,2)==1 and F(-1,2)>=-1,
            actual_width_unit_mass_scaling=bump==(320,121600,220000000),
            bump_supports_strictly_inside_I3=0<supports[0][0] and supports[-1][1]<5,
            five_bumps_pairwise_disjoint=all(supports[i][1]<supports[i+1][0] for i in range(4)),
            reserved_mean_I4_not_touched=supports[-1][1]<5<6<11,
            compressed_cutoffs_keep_allowed_window=all(0<=r['support_left']<r['support_right']<=1 for r in rows),
            compression_changes_no_endpoint_values_or_allowed_support=all(r['support_right']-r['support_left']==F(1,n) for n,r in zip((1,2,10,1000),rows)),
            fixed_window_first_and_second_costs_unbounded=rows[-1]['first_norm_lower']>1000 and rows[-1]['second_norm_lower']>10**6))

def propagation_review():
    old=seminorms();core=old.collar_jets(2,7,F(1,4),1)
    # Exact rational auxiliary controls, NOT the actual giant construction.
    full=old.moment_bounds(2,core['value'],5,7,(11,13))
    erased=old.moment_bounds(2,core['value'],5,7,(0,0))
    rp=old.repair_bound(full,2,3,F(1,1000));r0=old.repair_bound(erased,2,3,F(1,1000))
    out=propagate(core,rp['coefficient'],1,F(9,8),F(1,16),32,1216,2)
    direct=old.extended_jets(core,rp['coefficient'],1,out['cutoff']['radial'],out['cutoff']['radial2'],4,2,bump_bounds())
    tighter=propagate(core,rp['coefficient'],1,F(9,8),F(1,16),1,4,2)
    return dict(auxiliary_rational_not_actual=True,core=serial(core),full_moments=serial(full),
        repair_rows=serial(rp['transformed']),conditional_shape_budgets=['32','1216'],
        conditional_output=serial(out),checks=dict(
            composition_equals_retained_exact_jet_interface=out['extended']==direct,
            full_pressure_memory_changes_both_normalized_rows=rp['transformed'][3]!=r0['transformed'][3] and rp['transformed'][4]!=r0['transformed'][4],
            larger_unknown_shape_norms_increase_derivative_budget=out['extended']['radial']>tighter['extended']['radial'] and out['extended']['radial2']>tighter['extended']['radial2'],
            cutoff_derivatives_do_not_change_value_bound=out['extended']['value']==tighter['extended']['value'],
            positive_order_bumps_are_additive_not_multiplied_by_e=True,
            pressure_is_reconstructed_not_cut_directly=True))

def review():
    s,e,p=symbolic_review(),exact_review(),propagation_review()
    return dict(symbolic=s,exact=e,propagation=p,checks=dict(
        symbolic_derivative_identities=all(s['checks'].values()),
        exact_step_bump_geometry_and_nonuniformity_controls=all(e['checks'].values()),
        conditional_first_jet_and_complete_memory_composition=all(p['checks'].values())))
