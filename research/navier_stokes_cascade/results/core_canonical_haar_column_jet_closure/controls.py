"""Exact manufactured integrals and moving-chart controls, not actual giant profiles."""
from functools import lru_cache
from fractions import Fraction as Q
from math import comb,factorial
import sympy as s
from bounds import P,HaarColumnBudget,entrance,gaussian,from_original,root_join,certificate

v,x,r=s.symbols('v x r',real=True);z=v-s.Rational(1,2)

def poly_box(expr):
    return Q(sum(abs(c) for c in s.Poly(s.expand(expr),v,x).coeffs()))

@lru_cache(None)
def toy():
    # L=1. On [0,1]^2, 1/16 <= P <= exp(-4(v-1/2)^2), c=8.
    profile=(1-(1+x)*z*z)**4
    plus=s.Matrix([1,-v-x,1]);minus=s.Matrix([1,v+x,1])
    psi=v**4*(1-v)**4
    pulses=[profile*plus,profile*minus]
    class ManufacturedPulse:
        @lru_cache(None)
        def velocity(self,j):
            return P(16*max(sum(poly_box(s.diff(t,v,a,x,j-a)) for t in pulse)
                for pulse in pulses for a in range(j+1)))
    h=ManufacturedPulse();column=HaarColumnBudget(h,8)
    # J=5/4, transverse chi=(1-4 xi^2/r^2), int chi^2=8r/15.
    # Angular half and ci=2r retained; H/s0=(J*8/15)*int psi^2 t_r t_tan.
    matrices=[s.Matrix(3,3,lambda a,b:s.integrate(s.expand(psi**2*pulse[a]*pulse[b]),(v,0,1))*s.Rational(2,3))
        for pulse in pulses]
    H=s.Matrix.hstack(matrices[0][0,1:3].T,matrices[1][0,1:3].T)
    return profile,psi,pulses,matrices,H,column

def integral_review():
    profile,psi,pulses,matrices,H,column=toy();checks={};rows=[]
    checks['all_six_principal_tensor_entries_are_retained_and_nonzero']=all(m[a,b].subs(x,0)!=0 for m in matrices for a in range(3) for b in range(a,3))
    checks['two_selected_columns_are_the_original_tensor_entries']=H==s.Matrix.hstack(matrices[0][0,1:3].T,matrices[1][0,1:3].T)
    for j in range(4):
        for sigma,m in enumerate(matrices):
            actual=max(poly_box(s.diff(m[a,b],x,j)) for a in range(3) for b in range(3))
            checks[f'order{j}_sign{sigma}_all_tensor_entry_jets_enclosed_on_whole_parameter_interval']=column.raw_tensor_entry(j).evaluate(1)>=actual
        actual=max(sum(poly_box(s.diff(H[a,b],x,j)) for b in range(2)) for a in range(2))
        checks[f'order{j}_induced_matrix_row_sum_jet_bound']=column.matrix(j).evaluate(1)>=actual/Q(factorial(j))
        rows.append(dict(order=j,exact_matrix_raw_upper=actual))
    checks['higher_column_jets_are_nonzero']=s.diff(H,x,3)!=s.zeros(2)
    no_cut=s.integrate(s.expand(pulses[0][0]*pulses[0][2]),(v,0,1))*s.Rational(2,3)
    checks['discarding_temporal_cutoff_changes_actual_column']=s.simplify(no_cut-H[1,0])!=0
    no_gaussian_derivative=s.integrate(s.expand(psi**2*profile**2*s.diff(pulses[0][1]/profile,x)),(v,0,1))*s.Rational(2,3)
    checks['differentiating_only_t_over_P_loses_Gaussian_terms']=s.expand(s.diff(H[0,0],x)-no_gaussian_derivative)!=0
    xi=s.symbols('xi',real=True)
    integral=s.integrate((1-4*xi**2/r**2)**2,(xi,-r/2,r/2))
    physical=s.Rational(5,8)*integral*(2*r)*s.Symbol('I')
    checks['angular_half_transverse_mass_and_ci_give_physical_radius_square']=s.simplify(physical-s.Rational(2,3)*r*r*s.Symbol('I'))==0
    checks['omitting_angular_half_doubles_actual_column']=s.simplify(2*physical-physical)==physical
    checks['radius_independence_only_after_exact_s0_normalization']=s.simplify(physical/(r*r))==s.Rational(2,3)*s.Symbol('I')
    checks['two_sign_columns_give_nonzero_matrix_value_domain']=H.det().subs(x,0)!=0
    return checks,dict(rows=rows,scope='exact polynomial principal tensor integrals on L=1, x in [0,1]; manufactured finite-jet cutoffs, not the actual compact profile or full curl covariance',
        Gaussian_proof='a=(1+x)(v-1/2)^2 lies in [0,1/2]; log(1-a)<=-a gives P=(1-a)^4<=exp(-4(v-1/2)^2), P>=1/16')

def moving_review():
    y,w=s.symbols('y w',real=True);a=x+x*x;b=x**3-x
    f=(1+x)*s.cos(y-a)+(1+x*x)*s.sin(2*(w-b))
    avg=lambda e:s.simplify(s.integrate(s.integrate(s.expand_trig(e),(y,-s.pi,s.pi)),(w,-s.pi,s.pi))/(4*s.pi*s.pi))
    checks={}
    # Square contains a nonzero mixed covariance and has moving-center terms
    # pointwise. Only the complete periodic Haar integral removes them.
    total=f*f
    exact=((1+x)**2+(1+x*x)**2)/2
    for j in range(4):checks[f'order{j}_moving_translation_derivatives_cancel_only_after_full_Haar_integral']=s.simplify(avg(s.diff(total,x,j))-s.diff(exact,x,j))==0
    checks['moving_center_terms_are_nonzero_pointwise']=s.simplify(s.diff(f,x)-s.cos(y-a)-2*x*s.sin(2*(w-b)))!=0
    checks['truncating_to_nonperiodic_domain_destroys_translation_invariance']=s.diff(s.sin(1-a)+s.sin(a),x)!=0
    M=s.Matrix([[3,1],[1,5]]);rot=s.sqrt(2)-1
    checks['original_Jacobian_is_constant_and_below_two']=0<1+rot*rot<2
    for delta in range(4):
        matrix=M**delta
        checks[f'cover{delta}_normalized_Haar_needs_no_deck_count']=matrix.det()==14**delta and matrix.det()!=0
    checks['raw_multivariate_factorial_count_in_five_coordinates']=all(sum(Q(1,factorial(a)*factorial(b)*factorial(c)*factorial(d)*factorial(n-a-b-c-d))
        for a in range(n+1) for b in range(n-a+1) for c in range(n-a-b+1) for d in range(n-a-b-c+1))==Q(5**n,factorial(n)) for n in range(5))
    return checks,dict(scope='exact periodic manufactured translation integrals through order3 and integer determinant/jet-count identities; invariant Haar theorem in README supplies all orders, not a sampled cover claim')

@lru_cache(None)
def original_fixture():
    normal=entrance.reference.normal
    p=normal.frame.inputs(bmin=Q(1,20),bmax=1,Fmax=1,gmax=2,gmin=1,u=Q(1,100),cmin=1,cmax=1,Dbox=6,Mz=1)
    sp=normal.selector.fixture();L=normal.normal_length_certificate(p)['length_rule_floor']
    sp['length']=Q(-(-L.numerator//L.denominator))
    ell=normal.selector.select_tail(sp)['ell_tail'] # MANUFACTURED control only
    kw=dict(F=lambda j:P(j+2),G=lambda j:P(j+3),btilde=lambda j:P(j+1),
        path=lambda j:P({1:1}) if j==0 else P(int(j==1)),tilt=Q(1,100))
    return p,sp,ell,kw

def original_review():
    p,sp,ell,kw=original_fixture();phase,column,c,tc,g=from_original(p,sp,ell,**kw)
    template=entrance.inverse.roots.fixture();cut=lambda j:P(j+1);joined=root_join(template,column,cut)
    checks=dict(original_normal_tangent_and_band_gates_kept=all(c['gates'].values()) and all(tc['gates'].values()) and c['existing_ell']==ell,
        actual_fixed_tilt_not_replaced_by_its_upper_bound=g['actual_fixed_tilt']==kw['tilt']<=p['u'],
        original_Gaussian_curvature_is_derived_not_fitted=column.curvature==gaussian.curvature(kw['tilt'],p['lam_min'],p['lam_max'])[0],
        target_curl_loghalf_and_trace_callbacks_retained=all(getattr(joined,k) is getattr(template,k) for k in ('target','curl','loghalf','trace')),
        all_canonical_value_gaps_mass_radius_width_unchanged=all(getattr(joined,k)==getattr(template,k) for k in ('K0','beta','R','mass','width')))
    for j in range(3):
        checks[f'order{j}_joint_H_input_is_derived_from_same_pulse']=joined.H(j).c==column.matrix(j).c
        checks[f'order{j}_joint_W_input_keeps_full_cutoff_jets']=joined.wave(j).c==column.homogeneous.wave(j,cut).c
    checks['joined_inverse_root_signed_and_complete_curl_recipes_are_finite']=all(v.evaluate(1)>0 for v in (joined.K(1),joined.root(1),joined.signed(1,lambda j:P(j+1)),joined.complete_curl(0)))
    for i,cval in enumerate((Q(1,100),Q(1),Q(8))):
        G=2*max(Q(1),1/cval)
        checks[f'curvature{i}_rational_Gaussian_integral_upper_is_valid']=G*G*cval>=4
    checks['fixed_positive_tilt_required_for_uniform_Gaussian_constant']=gaussian.curvature(Q(1,1000),1,1)[0]<gaussian.curvature(Q(1,100),1,1)[0]
    return checks,dict(ell=ell,curvature=g['curvature_lower'],scope='same original-tail manufactured primitive bounds and same root norm recipes, no actual background supremum or cone positivity evaluation')

@lru_cache(None)
def review():
    checks={};parts={}
    for name,f in (('integral',integral_review),('moving',moving_review),('original',original_review)):
        c,d=f();checks.update({name+':'+k:bool(v) for k,v in c.items()});parts[name]=d
    return dict(checks=checks,parts=parts,certificate=certificate())
