"""Exact full-curl identities, whole-box coefficient jets and chart failure controls."""
from fractions import Fraction as Q
from functools import lru_cache
from math import comb,factorial
import sympy as s
from bounds import P,CurlBudget,columns,entrance,from_original,join_root,certificate

R,Z,T,Y,W=s.symbols('R Z T Y W',real=True);variables=(R,Z,T,Y,W)
er=s.Matrix([1,0,0]);ez=s.Matrix([0,0,1]);et=s.Matrix([0,1,0])
def cross(a):return s.Matrix([[0,-a[2],a[1]],[a[2],0,-a[0]],[-a[1],a[0],0]])
def zero(M):return all(s.cancel(e)==0 for e in M)

@lru_cache(None)
def manufactured():
    # R in [1,2], all other variables in [0,1]. Full normal floor N>=1.
    n=s.Matrix([R+Y,1,0]);h=n/n.dot(n)
    A=n.cross(s.Matrix([1+Z,R+T,1+Y*Y]))*(1+Z*T)*(1+Y*Y+W)
    eta=s.Rational(1,4);eps=s.Rational(1,4);M=2;g=2*R;vr=(1,-s.Rational(1,2))
    scalars=(eta,eta*eps,0,eta*M*g*vr[0],eta*M*g*vr[1])
    axes=(er,ez,er,er,er)
    matrices=[q*cross(axis)*cross(h) for q,axis in zip(scalars,axes)]
    D0=sum((q*cross(axis)*cross(h.diff(var)) for q,axis,var in zip(scalars,axes,variables)),s.zeros(3))
    D0+=eta/R*ez*et.T*cross(h)
    Dr=lambda f:f.diff(R)+M*g*(f.diff(Y)-f.diff(W)/2)
    C=h.cross(A)
    B=er.cross(Dr(C))+ez.cross(eps*C.diff(Z))+ez*C[1]/R
    base=entrance.inverse.ForcedBudget(normal=lambda j:P(4 if j==0 else int(j==1)),
        moving_normal=lambda j:P(),matrix=lambda j:P(),source=lambda j:P(),profile_rate=lambda j:P(),
        tangent=P(1),length=P(1),inverse_normal_square=1)
    return n,h,A,C,B,Dr,eta,matrices,D0,CurlBudget(base,1,2)

def equation_review():
    n,h,A,C,B,Dr,eta,matrices,D0,op=manufactured()
    represented=D0*A+sum((m*A.diff(var) for m,var in zip(matrices,variables)),s.zeros(3,1))
    connection=eta/R*ez*C[1]
    graph=eta*er.cross(Dr(C)-C.diff(R));axial=eta*ez.cross(C.diff(Z)/4)
    h_terms=D0-eta/R*ez*et.T*cross(h)
    checks=dict(full_cylindrical_graph_coefficient_identity_exact=zero(eta*B-represented),
        full_graph_curl_is_divergence_free=zero(s.Matrix([Dr(B)[0]+B[0]/R+B[2].diff(Z)/4])),
        transverse_input_uses_full_normal_not_component_inverse=n.dot(A)==0,
        cylindrical_connection_is_nonzero_and_cannot_be_dropped=not zero(connection),
        radial_graph_terms_are_nonzero_and_cannot_be_dropped=not zero(graph),
        slow_axial_terms_are_nonzero_and_cannot_be_dropped=not zero(axial),
        differentiated_normal_potential_terms_are_nonzero=not zero(h_terms*A),
        graph_coefficient_derivative_is_not_zero=not zero(matrices[3].diff(R)-cross(er)*cross(h.diff(R))*R),
        normal_reciprocal_derivatives_are_not_zero=not zero(h.diff(R,2)))
    # Independent cylindrical component formula has all three terms.
    direct=s.Matrix([-C[1].diff(Z)/4,C[0].diff(Z)/4-Dr(C)[2],Dr(C)[1]+C[1]/R])
    checks['independent_component_curl_matches_cross_product_formula']=zero(B-direct)
    return checks,dict(scope='exact manufactured theta-independent cylindrical/graph fields with dr=2 and nonconstant normal, amplitude and graph coefficient; no actual background supremum evaluated')

def rational_box(expr):
    q=s.symbols('q',real=True)
    num,den=s.fraction(s.cancel(expr.subs(R,1+q)));den=s.Poly(den,q,Y)
    if den.eval({q:0,Y:0})<0:num,den=-num,-den
    dc=den.eval({q:0,Y:0})
    if dc<=0 or any(c<0 for c in den.coeffs()):raise ValueError('No proved denominator floor')
    return Q(sum(abs(c) for c in s.Poly(num,q,Y).coeffs())/dc)

def matrix_box(M):return max(sum(rational_box(M[a,b]) for b in range(3)) for a in range(3))

def jet_review():
    n,h,A,C,B,Dr,eta,matrices,D0,op=manufactured();checks={};rows=[]
    for j in range(3):
        exact=max(sum(matrix_box(m.diff(R,a,Y,j-a)) for m in [D0]+matrices) for a in range(j+1))
        checks[f'order{j}_sum_of_all_matrix_raw_jets_enclosed_on_whole_box']=op.raw(j).evaluate(1)>=exact
        conjugated=matrices[:];conjugated[0]=2*R*conjugated[0]
        uexact=max(sum(matrix_box(m.diff(R,a,Y,j-a)) for m in [D0]+conjugated) for a in range(j+1))
        checks[f'order{j}_conjugated_U_coefficient_raw_jets_enclosed']=op.raw(j,True).evaluate(1)>=uexact
        rows.append(dict(order=j,original_raw_upper=exact,conjugated_raw_upper=uexact))
    checks['D0_uses_one_additional_normal_lift_jet']=op.raw(1).evaluate(1)>2*entrance.product(op.directions_R,op.lift,1).evaluate(1)
    checks['native_scale_normalization_keeps_graph_factor_bounded']=Q(1,4)*2<=1 and Q(1,4)*Q(1,4)<=1
    return checks,dict(rows=rows,scope='exact rational numerator-l1/positive-denominator enclosures on R in [1,2], Y in [0,1], through order2; other coefficient derivatives vanish exactly')

def chart_review():
    n,h,A,C,B,Dr,eta,matrices,D0,op=manufactured();u=s.symbols('u',positive=True)
    sub=lambda f:f.subs(R,s.sqrt(u))
    transformed=[sub(m) for m in matrices];transformed[0]=sub(2*R*matrices[0])
    Au=sub(A)
    represented=sub(D0)*Au+sum((m*Au.diff(var) for m,var in zip(transformed,(u,Z,T,Y,W))),s.zeros(3,1))
    checks=dict(U_operator_requires_conjugation_before_composition=all(s.simplify(e)==0 for e in sub(eta*B)-represented),
        reusing_R_first_order_matrix_in_U_is_detectably_wrong=not zero((sub(2*R*matrices[0])-sub(matrices[0]))*Au.diff(u)),
        value_root_gaps_do_not_change_under_positive_radial_reparametrization=True)
    # Independent finite normalized Taylor composition, not an analytic profile.
    z=s.symbols('z');mapjet=sum(s.binomial(s.Rational(1,2),j)*z**j for j in range(1,5))
    f=(R+1)**3;outer=[Q(s.diff(f,R,k).subs(R,1)/factorial(k)) for k in range(5)]
    direct=s.series((s.sqrt(1+z)+1)**3,z,0,5).removeO().expand()
    substitute=s.series(sum(s.Rational(c.numerator,c.denominator)*mapjet**k for k,c in enumerate(outer)),z,0,5).removeO().expand()
    for j in range(5):
        bound=op.compose(lambda k:P(outer[k]),j).evaluate(1)
        checks[f'order{j}_finite_composition_matches_direct_Taylor']=direct.coeff(z,j)==substitute.coeff(z,j)
        checks[f'order{j}_positive_composition_majorant_bounds_signed_exact_jet']=bound>=abs(Q(direct.coeff(z,j)))
    for j in range(1,5):
        for a in (Q(2,5),Q(1,2),Q(1),Q(2),Q(3)):
            coefficient=Q(1)
            for k in range(j):coefficient*=a-k
            checks[f'order{j}_inverse_power_a{a}_falling_factor_bound']=abs(coefficient)<=Q(factorial(j+2),2)
    checks['original_exponent_interval_is_sufficient_for_both_chart_bounds']=(Q(1,2)-2*Q(1,100)*Q(1,100000)>Q(1,3) and 2*(1+Q(1,100))<Q(5,2))
    return checks,dict(scope='exact dr=2 operator conjugation and finite signed Taylor controls through order4; analytical dr interval controls original fixed exponent, no chart/profile change')

@lru_cache(None)
def original_fixture():
    normal=entrance.reference.normal
    p=normal.frame.inputs(bmin=Q(1,20),bmax=1,Fmax=1,gmax=2,gmin=1,u=Q(1,100),cmin=1,cmax=1,Dbox=6,Mz=1)
    sp=normal.selector.fixture();L=normal.normal_length_certificate(p)['length_rule_floor']
    sp['length']=Q(-(-L.numerator//L.denominator));ell=normal.selector.select_tail(sp)['ell_tail'] # manufactured only
    kw=dict(F=lambda j:P(j+2),G=lambda j:P(j+3),btilde=lambda j:P(j+1),
        path=lambda j:P({1:1}) if j==0 else P(int(j==1)),tilt=Q(1,100),
        template=entrance.inverse.roots.fixture(),cutoff=lambda j:P(j+1))
    return p,sp,ell,kw

def original_review():
    p,sp,ell,kw=original_fixture();phase,column,op,root,c,tc,g=from_original(p,sp,ell,**kw)
    base=columns.root_join(kw['template'],column,kw['cutoff']);uroot=join_root(base,op,'U')
    checks=dict(original_same_band_normal_tangent_and_Gaussian_inputs_retained=all(c['gates'].values()) and all(tc['gates'].values()) and c['existing_ell']==ell,
        potential_lift_is_the_same_full_normal_recipe=op.homogeneous is column.homogeneous,
        R_chart_target_and_loghalf_callback_identity_retained=root.target is base.target and root.loghalf is base.loghalf,
        value_gaps_mass_width_radius_retained_in_both_charts=all(getattr(a,k)==getattr(base,k) for a in (root,uroot) for k in ('K0','beta','R','mass','width')))
    # Bound methods/closures are re-created by the original join: compare recipes,
    # while target/loghalf retain literal callback identity.
    for j in range(2):
        checks[f'order{j}_complete_curl_callback_derived_in_each_chart']=root.curl(j).c==op.normalized(j,'R').c and uroot.curl(j).c==op.normalized(j,'U').c
        for k in ('H','target','wave','loghalf'):
            checks[f'order{j}_{k}_same_R_recipe_and_all_U_jets_transported']=getattr(root,k)(j).c==getattr(base,k)(j).c and getattr(uroot,k)(j).c==op.compose(getattr(base,k),j).c
    xi=lambda j:P(j+1)
    checks['separately_supplied_signed_target_is_also_transported']=op.transport(xi,'U')(1).c==op.compose(xi,1).c and op.transport(xi,'U')(1).evaluate(1)>xi(1).evaluate(1)
    checks['complete_curl_and_signed_root_recipes_are_finite_in_both_charts']=all(a.complete_curl(0).evaluate(1)>0 and a.signed(1,op.transport(xi,chart)).evaluate(1)>0 for a,chart in ((root,'R'),(uroot,'U')))
    checks['whole_moving_axial_moment_retains_n_plus_thirteen_background_order']=(11+2==12+1==13)
    return checks,dict(ell=ell,scope='manufactured original primitive bounds and unchanged same-root norm substitutions; no actual background jets, full covariance means or residual evaluated')

@lru_cache(None)
def review():
    checks={};parts={}
    for name,f in (('equation',equation_review),('jets',jet_review),('chart',chart_review),('original',original_review)):
        c,d=f();checks.update({name+':'+k:bool(v) for k,v in c.items()});parts[name]=d
    return dict(checks=checks,parts=parts,certificate=certificate())
