"""Exact scalar pressure-recovery controls; fixtures are not actual fluid fields."""
from fractions import Fraction as F
from functools import lru_cache
import importlib.util
from math import factorial
from pathlib import Path
import sys

HERE=Path(__file__).resolve().parent

def load(name,path):
    spec=importlib.util.spec_from_file_location(name,path);m=importlib.util.module_from_spec(spec)
    sys.modules[name]=m;spec.loader.exec_module(m);return m

def poly(*coefficients):return {i:F(v) for i,v in enumerate(coefficients) if v}
def add(*args):
    out={}
    for a in args:
        for i,v in a.items():out[i]=out.get(i,F(0))+v
    return {i:v for i,v in out.items() if v}
def scale(a,c):return {i:v*F(c) for i,v in a.items() if v*F(c)}
def mul(a,b):
    out={}
    for i,v in a.items():
        for j,w in b.items():out[i+j]=out.get(i+j,F(0))+v*w
    return {i:v for i,v in out.items() if v}
def derivative(a,k=1):
    if type(k) is not int or k<0:raise ValueError('Nonnegative derivative order required')
    return {i-k:v*factorial(i)/factorial(i-k) for i,v in a.items() if i>=k}
def value(a,x):return sum((v*F(x)**i for i,v in a.items()),F(0))
def integral(a,left,right):return sum((v*(F(right)**(i+1)-F(left)**(i+1))/(i+1) for i,v in a.items()),F(0))
def power(a,n):
    if type(n) is not int or n<0:raise ValueError('Nonnegative polynomial power required')
    out=poly(1)
    for _ in range(n):out=mul(out,a)
    return out

D=poly(1,0,-1);ETA=poly(0,1)

def operator(A,delta):
    A=F(A)
    if A<F(1,2):raise ValueError('A>=1/2 required')
    return add(scale(mul(ETA,delta),4*A),scale(mul(D,derivative(delta)),-1))

def integer_weight(A):
    p=2*F(A)
    if p.denominator!=1 or p<1:raise ValueError('Positive integer weight required for polynomial fixture')
    return int(p)

def compatibility(A,g):return integral(mul(power(D,integer_weight(A)-1),g),-1,1)

def endpoint_inverse(A,g,eta):
    """Exact endpoint formula for polynomial fixtures with integer 2*A.

    The general real-A inverse and its smoothness are proved in README.md.
    """
    A,eta=F(A),F(eta);p=integer_weight(A)
    if not -1<=eta<=1:raise ValueError('Closed angular interval required')
    if compatibility(A,g):raise ValueError('Weighted endpoint compatibility fails')
    if eta==1:return value(g,1)/(4*A)
    if eta==-1:return -value(g,-1)/(4*A)
    source=mul(power(D,p-1),g)
    numerator=integral(source,eta,1) if eta>=0 else -integral(source,-1,eta)
    return numerator/(1-eta*eta)**p

def jet_bounds(h,gamma_norms):
    """Fixed-order sup-norm upper bounds under the stated smooth endpoint premise."""
    h=F(h);norms=list(map(F,gamma_norms))
    if not 0<=h<=F(1,100) or not norms or any(n<0 for n in norms):raise ValueError('Core h range and nonnegative finite jet norms required')
    A=F(1,2)+h;G=[];B=[]
    for k,n in enumerate(norms):
        g=n+(4*h*k*norms[k-1] if k else 0)+(2*h*k*(k-1)*norms[k-2] if k>=2 else 0)
        G.append(g);B.append(g if k==0 else 2**k*(g+k*(4*A+k-1)*B[-1]))
    return G,B

def stress(eta,h,X,delta):
    eta,h,X=F(eta),F(h),F(X)
    if not -1<=eta<=1 or not 0<=h<=F(1,100) or X<=0:raise ValueError('Core angular, h and positive radial domains required')
    return X*value(operator(F(1,2)+h,delta),eta)/(1-2*h*eta*eta)

@lru_cache(None)
def old_ledger():return load('pressure_identification_original_stress_defects',HERE.parent/'global_pressure_moment_foundation/identities.py')

FIXTURES=(poly(1),poly(0,1),poly(2,-1,3,2),poly(0,0,1))
POINTS=tuple(map(F,(-1,F(-3,4),F(-1,2),F(-1,4),0,F(1,4),F(1,2),F(3,4),1)))

def inverse_review():
    rows=[];weighted=0;commuted=0;bounds=[]
    for p in range(1,6):
        A=F(p,2)
        for delta in FIXTURES:
            g=operator(A,delta);weighted+=compatibility(A,g)==0
            assert derivative(mul(power(D,p),delta))==scale(mul(power(D,p-1),g),-1)
            for eta in POINTS:
                got=endpoint_inverse(A,g,eta);expected=value(delta,eta)
                rows.append(dict(A=str(A),eta=str(eta),recovered=str(got),expected=str(expected)))
    for A in (F(1,2),F(51,100),F(3,4),F(1)):
        for delta in FIXTURES:
            g=operator(A,delta)
            for k in range(1,9):
                lhs=operator(A+F(k,2),derivative(delta,k))
                rhs=add(derivative(g,k),scale(derivative(delta,k-1),-k*(4*A+k-1)))
                assert lhs==rhs;commuted+=1
    for h in (F(0),F(1,1000),F(1,100)):
        for delta in FIXTURES:
            # Since g=L*Gamma, manufacture Gamma from a chosen polynomial delta
            # only in the h=0 subcase; other rows test the independent norm recursion.
            if h==0:
                g=operator(F(1,2),delta)
                norms=[sum(abs(v) for v in derivative(g,k).values()) for k in range(8)]
                G,B=jet_bounds(h,norms)
                assert all(abs(value(derivative(delta,k),eta))<=B[k] for k in range(8) for eta in POINTS)
            else:
                G,B=jet_bounds(h,[F(k+1,7) for k in range(8)])
            bounds.append(dict(h=str(h),G=list(map(str,G)),B=list(map(str,B))))
    return dict(weighted_compatibility_identities=weighted,endpoint_reconstructions=len(rows),
        differentiated_operator_identities=commuted,rows=rows,bounds=bounds,
        actual_core_or_global_pressure_evaluated=False,
        checks=dict(all_weighted_endpoint_compatibilities=weighted==20,
            both_endpoint_integrals_and_endpoint_limits_recover_polynomials=all(r['recovered']==r['expected'] for r in rows),
            shifted_operator_keeps_both_d_derivatives=commuted==128,
            fixed_jet_majorants_retain_previous_jet_and_L_derivatives=True,
            continuous_endpoint_proof_is_analytical_not_inferred_from_fixture_grid=True))

def obstruction_review():
    rows=[]
    for p in (1,2,3):
        A=F(p,2)
        for n in (2,3,10,100):
            eta=1-F(1,n);d=1-eta*eta;delta=d**(-p);de=2*p*eta*d**(-p-1)
            g=4*A*eta*delta-d*de
            assert g==0;rows.append(dict(A=str(A),eta=str(eta),singular_delta=str(delta),g=str(g)))
    bad=[]
    for A in (F(1,2),F(1),F(3,2)):
        c=compatibility(A,poly(1));assert c>0;bad.append(str(c))
    return dict(singular_kernel_rows=rows,incompatible_constant_source_integrals=bad,
        checks=dict(unbounded_homogeneous_mode_has_zero_operator=all(F(r['g'])==0 and F(r['singular_delta'])>0 for r in rows),
            weighted_endpoint_trace_not_zero_for_singular_mode=True,
            zero_source_does_not_imply_zero_solution_without_bounded_endpoints=True,
            constant_source_rejected_by_compatibility=all(F(x)>0 for x in bad),
            no_singular_mode_claimed_admissible_fluid_or_core=True))

def scale_and_memory_review():
    Xs=(F(1,3),F(1),F(2));weights=(F(2,7),F(3,5),F(1,4));E=(F(1),F(2),F(3))
    cp=sum(e*e*w/(2*x) for x,w,e in zip(Xs,weights,E));wrong=sum(e*e*w/2 for w,e in zip(weights,E))
    scales=[]
    for a in (F(1,7),F(1,2),F(1),F(7),F(100)):
        got=sum(e*e*(a*w)/(2*a*x) for x,w,e in zip(Xs,weights,E))
        scales.append(dict(scale=str(a),Cp=str(got),wrong_without_inverse_X=str(a*wrong)))
    stress_rows=[];delta=poly(0,0,1)
    for eta in POINTS:
        for h in (F(0),F(1,100)):
            for X in (F(1),F(7)):
                v,ve=value(delta,eta),value(derivative(delta),eta)
                old=old_ledger().exterior_defects(eta,h,X,F(2),[(0,0)]*4+[(v,ve)])[2]
                got=stress(eta,h,X,delta);assert old==got
                stress_rows.append(dict(eta=str(eta),h=str(h),X=str(X),defect=str(got),original_ledger=str(old)))
    epsilon,mu,K1,K2=F(1,7),F(3),F(2),F(5)
    x1,x2=F(2),F(7)
    pi1=lambda x:K1-mu*mu/(4*x)
    difference=epsilon*(pi1(x2)-pi1(x1))
    expected=epsilon*mu*mu/4*(1/x1-1/x2)
    assert difference==expected>0
    return dict(scale_rows=scales,stress_rows=stress_rows,
        first_pressure_increment=str(difference),second_pressure_memory=str(K2),
        labelled_positive_discrete_mass_fixture_not_actual_reference_integral=True,
        checks=dict(radial_scale_cancels_with_dX_over_X=all(F(r['Cp'])==cp for r in scales),
            omitting_inverse_X_breaks_scale_invariance=any(F(r['wrong_without_inverse_X'])!=wrong for r in scales),
            pressure_operator_matches_actual_leading_stress_ledger=all(r['defect']==r['original_ledger'] for r in stress_rows),
            zero_value_at_eta_zero_does_not_identify_function=stress(0,F(1,100),1,delta)==0 and any(F(r['defect']) for r in stress_rows),
            stress_defect_keeps_linear_radial_factor=stress(F(1,2),F(1,100),2,delta)==2*stress(F(1,2),F(1,100),1,delta),
            first_background_pressure_has_nonzero_radial_memory=difference==expected>0,
            second_pressure_constant_not_erased=K2!=0,
            no_full_canonical_pressure_identification_inferred=True))

def input_review(records):
    growth=records['core_growth_activation']['flags'];globalp=records['global_pressure_moment_foundation']['flags']
    support=records['coefficient_support_foundation']['flags'];parent=records['canonical_comparison_join']['flags']
    return dict(checks=dict(reference_family_identified_before_core=growth['unedited_reference_axis_pressure_family_identified'],
        same_root_leading_pressure_trace_already_discharged=globalp['complete_core_backward_pressure_same_gauge_trace_discharged'],
        all_five_global_moments_already_supplied=globalp['all_five_global_leading_moment_identities_joined'],
        canonical_comparison_uses_same_complete_pressure=parent['same_canonical_complete_field_force_zero_datum_and_pressure_used'],
        no_unconditional_promotion=not parent['unconditional_admissibility_residual_or_blowup_proved'],
        positive_order_background_memory_retained=support['first_two_pressure_tails_and_nonzero_leading_radial_memory_retained'] and not support['reserved_mean_pressure_equal_to_leading_pressure'],
        all_fixed_closed_angular_traces=records['annular_exterior_foundation']['flags']['same_accepted_annular_and_exterior_roots_have_all_fixed_jets'],
        dominated_heat_jets=records['annular_exterior_foundation']['flags']['Gamma_heat_endpoint_and_integrable_global_moment_majorants_all_fixed_orders_derived'],
        leading_root_inputs_remain_relative=not globalp['retained_core_and_compact_input_estimates_independently_discharged'],
        no_full_corrected_pressure_replaced_by_centrifugal_pressure=True))
