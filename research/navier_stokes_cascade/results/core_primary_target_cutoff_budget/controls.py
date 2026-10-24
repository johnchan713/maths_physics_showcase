"""Finite exact edge/chart/cutoff controls; not giant-profile seminorm evaluations."""
from fractions import Fraction as Q
from functools import lru_cache
from math import comb,factorial
import sympy as s
from bounds import P,curl,ChartBudget,TargetBudget,CutoffBudget,from_original,certificate,load,RADIUS

tame=load('target_cutoff_tame','core_fixed_profile_tame/proof_checks.py')
x,e,d=s.symbols('x e d',real=True)

def l1(expr):return Q(sum(abs(c) for c in s.Poly(s.expand(expr),d,e).coeffs()))

@lru_cache(None)
def profile_fixture():
    F=s.Matrix([1+x+e*e,(2-x)**6*e*(1+x*x)])
    def norm(j):
        return sum((sum(abs(Q(c))*2**key[0] for key,c in s.Poly(s.expand(s.diff(f,x,a,e,j-a).subs(x,2-d)),d,e).terms())
            /Q(factorial(a)*factorial(j-a)) for f in F for a in range(j+1)),Q(0))
    return F,norm,TargetBudget(ChartBudget(1,Q(1,8)),norm)

def edge_review():
    F,M,b=profile_fixture();checks={};rows=[]
    quotient=F/(2-x)**3
    for j in range(4):
        actual=sum(l1(s.cancel((d**(j+3)*s.diff(f,x,a,e,j-a)).subs(x,2-d)))/Q(factorial(a)*factorial(j-a))
            for f in quotient for a in range(j+1))
        checks[f'order{j}_all_mixed_outer_collar_quotient_jets_enclosed']=b.quotient(j)>=actual
        rows.append(dict(order=j,weighted_exact_polynomial_upper=actual,profile_recipe=b.quotient(j)))
    a,bb=s.symbols('a b',positive=True)
    zeta=s.exp(-a/x**2-bb/(2-x)**2);T=zeta*quotient
    checks['outer_factor_reconstructs_same_T0_and_smooth_F']=s.simplify(T/zeta*(2-x)**3-F)==s.zeros(2,1)
    checks['inner_cancellation_leaves_smooth_opposite_weight_factor']=all(s.simplify(s.expand_power_exp(v))==0 for v in T/s.exp(-a/x**2)-s.exp(-bb/(2-x)**2)*F/(2-x)**3)
    checks['quotient_itself_keeps_cubic_outer_pole']=s.limit((2-x)**3*quotient[0],x,2)==3+e*e
    checks['replacing_inner_flat_constant_changes_exponential_ratio']=s.simplify(s.exp(-a/x**2)/s.exp(-1/x**2))==s.exp((1-a)/x**2)
    checks['sharper_pole_fits_original_root_allowance']=all(j+3<=3*j+3 for j in range(8))
    return checks,dict(rows=rows,scope='exact manufactured factored profile and every mixed jet through order3 on the outer collar; analytical edge factorization supplies actual fixed-profile finiteness, not these numbers')

def jet_expansion(jet):
    out=tame.Jet(jet.n);power=tame.Jet(jet.n,{(0,0):1})
    for k in range(1,jet.n+1):power=power*jet;out=out+Q((-1)**(k+1),k)*power
    return out

def integer_power(jet,n):
    p=tame.Jet(jet.n,{(0,0):1})
    for _ in range(n):p=p*jet
    return p

def degree_norm(jet,n):return sum(abs(c) for key,c in jet.c.items() if sum(key)==n)

def chart_review():
    _,_,budget=profile_fixture();chart=budget.chart;checks={}
    for case,(h,eta) in enumerate(((Q(1,200),Q(0)),(Q(1,200),Q(1)),(Q(1,100),Q(-1)))):
        inc=tame.geometry.implicit_jet(h,eta,4);u=jet_expansion(inc);p=tame.Jet(4,{(1,0):1})
        eta_jet=(eta+p)*inc.power_one(-(Q(1,2)-h));pref=inc.power_one(-(1+h))
        xx=1-u;dd=u # yb=1+u at the same x=1 center
        FF0=1+xx+eta_jet*eta_jet
        FF1=integer_power(1+u,6)*eta_jet*(1+xx*xx)
        targets=[pref*f*dd.power_one(-3) for f in (FF0,FF1)]
        loghalf=-chart.width**2/2*(-u).power_one(-2)-2*u.power_one(-2)
        for j in range(1,5):
            checks[f'case{case}_order{j}_implicit_log_coordinate_majorant']=degree_norm(u,j)<=chart.log_s(j)
            checks[f'case{case}_order{j}_eta_and_scale_factor_majorants']=degree_norm(eta_jet,j)<=2*chart.exponential('eta',j)+chart.exponential('eta',j-1) and degree_norm(pref,j)<=chart.exponential('scale',j)
            checks[f'case{case}_order{j}_original_band_target_composition_bound']=sum(degree_norm(f,j) for f in targets)<=budget.target(j).evaluate(1)
            checks[f'case{case}_order{j}_loghalf_composition_bound']=degree_norm(loghalf,j)<=budget.loghalf(j).evaluate(1)
        affine=tame.Jet(4,{key:c for key,c in u.c.items() if sum(key)==1})
        wrong=-chart.width**2/2*(-affine).power_one(-2)-2*affine.power_one(-2)
        checks[f'case{case}_nonlinear_chart_terms_cannot_be_omitted']=wrong.c!=loghalf.c
    checks['loghalf_order_zero_is_unused_not_a_finite_bound']=budget.loghalf(0).evaluate(1)==0
    changed=TargetBudget(ChartBudget(1,Q(1,4)),lambda j:1)
    checks['original_width_square_is_retained_in_loghalf']=changed.loghalf(1).c!=budget.loghalf(1).c
    return checks,dict(scope='exact fourth-order implicit two-slow-variable Taylor jets at three manufactured chart centers, including eta endpoints; no actual giant field sampled')

def fixed_profile(j):return (Q(3),Q(3),Q(1))[j] if j<3 else Q(0)

def cutoff_review():
    chart=ChartBudget(1,Q(1,8));c=CutoffBudget(chart,dyadic=fixed_profile,mesh=fixed_profile,transverse=fixed_profile,temporal=fixed_profile,radius=Q(1,5))
    H1,H2,rr,zz,tt,S,L,ci=s.symbols('H1 H2 rr zz tt S L ci',positive=True)
    f=lambda z:1+z+z*z
    transverse=s.Rational(9,10)*H1+s.Rational(3,10)*H2
    time=-s.Rational(3,10)*H1+s.Rational(9,10)*H2
    aux=s.expand(f(5*transverse)*f(s.Rational(5,2)*time))
    mesh=s.expand(f(S**3*rr)*f(S**3*zz)*f(S**3*tt));checks={}
    for j in range(4):
        av=sum(abs(Q(v)) for key,v in s.Poly(aux,H1,H2).terms() if sum(key)==j)
        mv=sum(abs(Q(v.subs(S,2))) for key,v in s.Poly(mesh,rr,zz,tt).terms() if sum(key)==j)
        checks[f'order{j}_dual_auxiliary_covector_cutoff_jets_enclosed']=av<=c.auxiliary_jets(j).evaluate(1)
        checks[f'order{j}_original_S_cubed_mesh_jets_enclosed']=mv<=c.mesh_jets(j).evaluate(2)
        checks[f'order{j}_complete_cutoff_polynomial_degree_at_most_3j']=all(k<=3*j for k in c.normalized(j).c)
    checks['relative_pulse_cutoff_cancels_fast_time_length_exactly']=s.simplify((time/(ci*L)).subs(ci,2*s.Rational(1,5)/L)-s.Rational(5,2)*time)==0
    checks['slow_partition_is_included_once_in_wave_norm']=c.normalized(1).evaluate(2)>c.normalized(1,False).evaluate(2)
    checks['raw_cutoff_is_sufficient_for_every_factorial_input_jet']=all(c.raw(j).evaluate(2)>=c.normalized(j).evaluate(2) for j in range(4))
    checks['original_fixed_geometry_radius_not_reselected']=RADIUS==Q(1,6057072)
    return checks,dict(scope='independent exact polynomial cutoff/mesh Taylor coefficients, not the actual compact profiles; finite fixed cutoff seminorms supplied by original smooth partition construction')

@lru_cache(None)
def original_fixture():
    normal=curl.entrance.reference.normal
    p=normal.frame.inputs(bmin=Q(1,20),bmax=1,Fmax=1,gmax=2,gmin=1,u=Q(1,100),cmin=1,cmax=1,Dbox=6,Mz=1)
    sp=normal.selector.fixture();L=normal.normal_length_certificate(p)['length_rule_floor'];sp['length']=Q(-(-L.numerator//L.denominator))
    ell=normal.selector.select_tail(sp)['ell_tail'] # manufactured only
    kw=dict(F=lambda j:P(j+2),G=lambda j:P(j+3),btilde=lambda j:P(j+1),path=lambda j:P({1:1}) if j==0 else P(int(j==1)),
        tilt=Q(1,100),template=curl.entrance.inverse.roots.fixture(),profile=lambda j:Q(j+1),
        dyadic=fixed_profile,mesh=fixed_profile,transverse=fixed_profile,temporal=fixed_profile)
    return p,sp,ell,kw

def original_review():
    p,sp,ell,kw=original_fixture();phase,column,op,root,c,tc,g,target,cut=from_original(p,sp,ell,**kw)
    uroot=curl.join_root(root,op,'U');checks=dict(
        original_tail_and_normal_tangent_gates_retained=all(c['gates'].values()) and all(tc['gates'].values()) and c['existing_ell']==ell,
        original_width_and_radius_retained=target.chart.width==root.width==kw['template'].width and cut.radius==RADIUS,
        all_original_root_value_gaps_and_support_width_unchanged=all(getattr(root,k)==getattr(kw['template'],k) for k in ('K0','beta','R','mass','width')))
    for j in range(2):
        checks[f'order{j}_target_and_loghalf_slots_derived']=root.target(j).c==target.target(j).c and root.loghalf(j).c==target.loghalf(j).c
        checks[f'order{j}_H_wave_and_complete_curl_slots_derived']=root.H(j).c==column.matrix(j).c and root.wave(j).c==column.homogeneous.wave(j,cut.raw).c and root.curl(j).c==op.normalized(j).c
        checks[f'order{j}_U_chart_transports_target_loghalf_and_cutoff_wave']=all(getattr(uroot,k)(j).c==op.compose(getattr(root,k),j).c for k in ('target','loghalf','wave'))
    checks['joined_amplitude_curl_and_signed_recipes_are_finite']=all(v.evaluate(1)>0 for v in (root.amplitude(1),root.complete_curl(0),root.signed(1,lambda j:P(j+1))))
    return checks,dict(ell=ell,scope='manufactured fixed-profile seminorms and original primitive bounds, same norm recipes only; no actual profile seminorms, root positivity, full covariance means or residual evaluated')

@lru_cache(None)
def review():
    checks={};parts={}
    for name,f in (('edge',edge_review),('chart',chart_review),('cutoff',cutoff_review),('original',original_review)):
        c,d=f();checks.update({name+':'+k:bool(v) for k,v in c.items()});parts[name]=d
    return dict(checks=checks,parts=parts,certificate=certificate())
