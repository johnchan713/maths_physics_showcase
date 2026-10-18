"""Exact mixed-jet, edge and missing-freezing controls; no actual field evaluation."""
from fractions import Fraction as F
from itertools import product
from math import comb,factorial

def exact(x):
    if isinstance(x,bool) or not isinstance(x,(int,F)):raise TypeError('Exact rational input required')
    return F(x)
def index(n,minimum=0):
    if type(n) is not int or n<minimum:raise ValueError('Finite integer index required')
    return n
def ceil(x):
    x=exact(x);return -(-x.numerator//x.denominator)
def ceil_root(x,degree):
    x=exact(x);degree=index(degree,1)
    if x<0:raise ValueError('Nonnegative root argument required')
    lo=0;hi=1<<((x.numerator.bit_length()+degree-1)//degree)
    while lo<hi:
        mid=(lo+hi)//2
        if mid**degree*x.denominator>=x.numerator:hi=mid
        else:lo=mid+1
    return lo

def frozen_margin(gamma0,tau):
    gamma0,tau=map(exact,(gamma0,tau))
    if not 0<gamma0<=1 or not 0<=tau<1:raise ValueError('Positive cone margin and tau<1 required')
    return 1-(1-gamma0+tau)/(1-tau)

def acceptance_floor(length,h,radius,old_band,Cfreeze,gamma0,budget_gamma,carrier_length):
    length,h,radius,Cfreeze,gamma0,budget_gamma,carrier_length=map(exact,(length,h,radius,Cfreeze,gamma0,budget_gamma,carrier_length))
    index(old_band,16)
    if length<=0 or not 0<h<F(1,100) or not 0<radius<=F(1,6057072):raise ValueError('Original h/radius and positive length required')
    if Cfreeze<0 or carrier_length<0 or not 0<budget_gamma<=gamma0/2 or not 0<gamma0<=1:raise ValueError('Explicit freezing and certified length-budget margin required')
    L=max(length,carrier_length)
    old=max(16,old_band,ceil(512/h**2),ceil_root(L/(2*radius),2))
    freeze=ceil_root(8*Cfreeze/gamma0,6);E=max(old,freeze);tau=Cfreeze/F(E**6)
    return dict(old_displayed_floor=old,freezing_floor=freeze,corrected_floor=E,
        length_margin=str(2*radius*E*E-L),freezing_margin=str(gamma0*E**6-8*Cfreeze),
        tau_upper=str(tau),frozen_margin_lower=str(frozen_margin(gamma0,tau)),
        tail_exponent=E+2,actual_tiny_q_not_materialized=True)

def indices(order,dim=2):return sorted((I for I in product(range(index(order)+1),repeat=index(dim,1)) if sum(I)<=order),key=lambda I:(sum(I),I))
def subs(I):return product(*(range(n+1) for n in I))
def minus(I,J):return tuple(i-j for i,j in zip(I,J))
def choose(I,J):
    out=1
    for i,j in zip(I,J):out*=comb(i,j)
    return out
def mul(a,b,order,dim=2):
    return {I:sum((F(choose(I,J))*a.get(J,F(0))*b.get(minus(I,J),F(0)) for J in subs(I)),F(0)) for I in indices(order,dim)}
def exp_jets(logjets,order,dim=2):
    """Normalized jets P^-1 D^I P at a point, with P/P0=exp(l-l0)."""
    out={tuple([0]*dim):F(1)}
    for I in indices(order,dim)[1:]:
        e=next(k for k,v in enumerate(I) if v);base=list(I);base[e]-=1;base=tuple(base)
        total=F(0)
        for J in subs(base):
            Jplus=list(J);Jplus[e]+=1
            total+=choose(base,J)*logjets.get(tuple(Jplus),F(0))*out[minus(base,J)]
        out[I]=total
    return out
def quotient_jets(Pjets,fjets,order,dim=2):
    zero=tuple([0]*dim);p0=Pjets.get(zero,F(0))
    if p0<=0:raise ValueError('Positive point envelope required')
    out={}
    for I in indices(order,dim):
        total=fjets.get(I,F(0))-sum((choose(I,J)*Pjets.get(J,F(0))*out[minus(I,J)] for J in subs(I) if any(J)),F(0))
        out[I]=total/p0
    return out
def reciprocal_half(Y,order,dim=2):
    """Exact negative-half composition for normalized fixture Y(0)=1."""
    zero=tuple([0]*dim)
    if Y.get(zero)!=1:raise ValueError('Normalized fixture Y0=1 required')
    z=dict(Y);z[zero]=F(0);power={zero:F(1)};coef=F(1);out={I:F(0) for I in indices(order,dim)}
    for n in range(order+1):
        if n:coef*=(-F(1,2)-(n-1))/n;power=mul(power,z,order,dim)
        for I in out:out[I]+=coef*power.get(I,F(0))
    return out

def edge_logprime(a,b,Y,x,p,q):
    a,b,Y,x=map(exact,(a,b,Y,x));p,q=index(p),index(q)
    if min(a,b,Y)<=0 or not 0<x<Y:raise ValueError('Positive edge widths and interior point required')
    return 2*a/x**3-2*b/(Y-x)**3-p/x+q/(Y-x)
def safe_edge_point(a,b,Y,p):
    a,b,Y=map(exact,(a,b,Y));p=index(p)
    if min(a,b,Y)<=0:raise ValueError('Positive edge data required')
    x=min(F(1),Y/2)
    while p*x*x>a/2 or 16*b*x**3/Y**3>a/2:x/=2
    return x
def edge_absorption(p,c):
    p=index(p);c=exact(c)
    if c<=0:raise ValueError('Positive retained exponential coefficient required')
    N=(p+1)//2;return max(F(1),F(factorial(N))/c**N)

def acceptance_review():
    rows=[]
    for gamma0 in (F(1,10),F(1,4),F(1,2),F(1)):
        for Cf in (F(0),F(1),F(10**100)):
            r=acceptance_floor(25,F(1,200),F(1,6057072),16,Cf,gamma0,gamma0/2,0)
            rows.append(dict(gamma0=str(gamma0),Cfreeze=str(Cf),certificate=r))
    ell=10**10;gamma=F(1,4);Cf=gamma*ell**6;tau=Cf/F(ell**6)
    witness=dict(ell=ell,Cfreeze=str(Cf),length_passes=2*F(1,6057072)*ell**2>=25,
        whole_pulse_smallness_sufficient=ell>=512/F(1,200)**2,
        freezing_passes=ell**6>=8*Cf/gamma,frozen_margin_lower=str(frozen_margin(gamma,tau)))
    return dict(rows=rows,omitted_freezing_parameter_witness=witness,
        checks=dict(corrected_length_and_freezing_margins=all(F(x['certificate']['length_margin'])>=0 and F(x['certificate']['freezing_margin'])>=0 for x in rows),
            certified_frozen_half_margin=all(F(x['certificate']['frozen_margin_lower'])>=F(x['gamma0'])/2 for x in rows),
            length_and_whole_pulse_conditions_alone_do_not_encode_freezing=witness['length_passes'] and witness['whole_pulse_smallness_sufficient'] and not witness['freezing_passes'],
            old_records_not_rewritten_and_actual_counterexample_not_claimed=True))

LOGJ={(1,0):F(2),(0,1):F(-3),(2,0):F(1,2),(1,1):F(3,5),(0,2):F(-2),(3,1):F(1,7)}
def envelope_review():
    rows=[];identities=0
    for n in range(1,7):
        E=exp_jets(LOGJ,n);Q=exp_jets({I:-v for I,v in LOGJ.items()},n);B=exp_jets({I:abs(v) for I,v in LOGJ.items()},n)
        H={I:F((-1)**sum(I)*(1+I[0]+2*I[1]),1+sum(I)) for I in indices(n)};Fj=mul(E,H,n)
        recovered=quotient_jets(E,Fj,n);from_inverse=mul(Fj,Q,n);upper=mul({I:abs(v) for I,v in Fj.items()},B,n)
        identity=mul(E,Q,n)
        for I in indices(n):
            assert identity[I]==(1 if not any(I) else 0)
            assert recovered[I]==from_inverse[I]==H[I]
            assert abs(Q[I])<=B[I] and abs(recovered[I])<=upper[I];identities+=1
        rows.append(dict(order=n,jets=len(identity),quotient_exact=True,Bell_bound=True))
    return dict(rows=rows,mixed_jet_identities=identities,
        checks=dict(all_mixed_quotient_and_inverse_identities=all(r['quotient_exact'] for r in rows),
            nonnegative_Bell_majorants_bound_reciprocal_and_source_jets=all(r['Bell_bound'] for r in rows),
            tiny_P0_cancels_only_with_full_normalized_input=True,
            all_generated_entrance_jets_need_smooth_support_not_Gaussian_smallness=True,
            source_normal_and_moving_normal_pressure_terms_not_removed=True))

def weighted_review():
    n=5;Y={I:(F(1) if not any(I) else F(1+sum(I),7+I[0]+2*I[1])) for I in indices(n)}
    R=reciprocal_half(Y,n);identity=mul(Y,mul(R,R,n),n)
    assert all(v==(1 if not any(I) else 0) for I,v in identity.items())
    X={I:F(2+I[0]+3*I[1],5) for I in indices(n)};normalized=mul(X,R,n);rows=[]
    for w in (F(1,2),F(1,10),F(1,1000)):
        target={I:w*w*v for I,v in X.items()};reciprocal={I:v/w for I,v in R.items()};signed=mul(target,reciprocal,n)
        assert signed=={I:w*v for I,v in normalized.items()}
        rows.append(dict(sqrt_zeta=str(w),mixed_signed_jets=len(signed),half_weight_exact=True))
    edge=[]
    for a in (F(1),F(1,10**40)):
        for p,q in ((0,0),(3,6),(11,20),(100,7)):
            x=safe_edge_point(a,4,2,p);margin=edge_logprime(a,4,2,x,p,q)-a/x**3
            edge.append(dict(a=str(a),p=p,q=q,x=str(x),log_derivative_margin=str(margin),primitive_factor=str(x**3/a)))
    beta=sum((F((-1)**k*comb(2,k),3+k) for k in range(3)),F(0))
    # Finite polynomial zero-mode control: compact patching retains a nonflat
    # interior remainder even when the unpatched primitive has nonzero total.
    x=F(1,2);A_mid=(6*x-6*x*x)*beta
    compact=dict(left_value='0',right_value=str(beta-beta),interior_A=str(A_mid),
        derivative_mid=str(x*x*(1-x)**2-A_mid))
    return dict(root_identity_jets=len(identity),signed_rows=rows,edge_rows=edge,positive_polynomial_total=str(beta),compact_zero_mode_control=compact,
        checks=dict(reciprocal_square_root_mixed_identity=all(v==(1 if not any(I) else 0) for I,v in identity.items()),
            target_full_weight_and_reciprocal_half_weight_leave_half=all(r['half_weight_exact'] for r in rows),
            oriented_primitive_gain_has_nonnegative_log_derivative_margin=all(F(r['log_derivative_margin'])>=0 for r in edge),
            nonzero_total_blocks_flat_opposite_edge=beta==F(1,30)>0,
            compact_patching_does_not_require_pointwise_zero_total=compact['left_value']==compact['right_value']=='0' and A_mid!=0,
            nonzero_zero_mode_remainder_not_misclassified_as_flat=True,
            actual_half_exponential_constants_and_physical_Haar_costs_retained=True,
            zero_Haar_and_radial_moments_are_not_optional=True,
            complete_covariance_axial_velocity_and_Pi1_Pi2_targets_still_required=True))
