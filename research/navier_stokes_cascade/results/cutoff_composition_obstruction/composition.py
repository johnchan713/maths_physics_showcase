"""Exact local composition obstruction; manufactured quadrature diagnostics only."""
from functools import lru_cache
import mpmath as mp
import sympy as sp

def validate_parameters(delta,n):
    if not isinstance(n,int) or isinstance(n,bool) or n<1:
        raise ValueError('n must be an integer >=1')
    if isinstance(delta,bool) or not 0<mp.mpf(delta)<=1:raise ValueError('0<delta<=1 required')

@lru_cache(None)
def symbolic_review():
    d,n=sp.symbols('d n',positive=True)
    u=sp.symbols('u',real=True)
    sigma=1/(1+sp.exp(1/u**2-1/(1-u)**2))
    step0=sp.simplify(sigma.subs(u,sp.Rational(1,2)))
    step1=sp.simplify(sp.diff(sigma,u).subs(u,sp.Rational(1,2)))
    x=2+3*d/16;r=5*d/64;V=x+r;b=sp.sqrt(x-1)
    s=sp.sqrt(5*d)/8;s1=-sp.sqrt(d)*(2*sp.sqrt(5)*n+1/(8*sp.sqrt(5)))
    v1=d/8;t=-b;t1=-d/(16*b);V1=v1+2*s*s1
    bracket=t1*s+t*s1-t*s*V1/V
    positive=b*sp.sqrt(d)*(2*sp.sqrt(5)*n+1/(8*sp.sqrt(5)))*(x-r)/V+7*d**2*s/(1024*b*V)
    # Independent differentiation through the density inverse at theta=pi.
    ts,ss,vs,theta=sp.symbols('ts ss vs theta',real=True)
    vv=vs+ss**2
    extra=2*ts*sp.sqrt(2)*ss*(1-sp.cos(theta))-ss**2*sp.sin(2*theta)/2
    Phi=theta/(2*sp.pi)+extra/(2*sp.pi*vv)
    P=extra/(4*sp.pi*vv)
    subst={ts:t,ss:s,vs:x,theta:sp.pi}
    Phieta=sum(sp.diff(Phi,z).subs(subst)*j for z,j in ((ts,t1),(ss,s1),(vs,v1)))
    wpi=x/(2*sp.pi*V)
    moving=sp.simplify(Phieta/(4*sp.pi*wpi))
    frozen=sum(sp.diff(P,z).subs(subst)*j for z,j in ((ts,t1),(ss,s1),(vs,v1)))
    # Raw chain identity via independent normalized Taylor multiplication.
    e=sp.symbols('e');hj=sp.symbols('h0:5');vj=sp.symbols('v1:5');zj=sp.symbols('z0:5')
    h=sum(hj[j]*e**j/sp.factorial(j) for j in range(5))
    dv=sum(vj[j-1]*e**j/sp.factorial(j) for j in range(1,5))
    chain=sp.expand(h*sum(zj[j]*dv**j/sp.factorial(j) for j in range(5))).coeff(e,4)*24
    h0,h1,h2,h3,h4=hj;v_1,v2,v3,v4=vj;z0,z1,z2,z3,z4=zj
    claimed=h0*z4*v_1**4+z3*(6*h0*v_1**2*v2+4*h1*v_1**3)+z2*(h0*(3*v2**2+4*v_1*v3)+12*h1*v_1*v2+6*h2*v_1**2)+z1*(h0*v4+4*h1*v3+6*h2*v2+4*h3*v_1)+z0*h4
    # The scalar identities use vs=1+ts^2 exactly, not independent arbitrary data.
    density=(1+(ts+sp.sqrt(2)*ss*sp.sin(theta))**2)/(2*sp.pi*vv)
    relation={vs:1+ts**2}
    ga=-(vv/(1+(ts+sp.sqrt(2)*ss*sp.sin(theta))**2)-1)/2
    checks=dict(sigma_center_half=step0==sp.Rational(1,2),sigma_center_slope_eight=step1==8,
        exact_s_center=s**2==r,
        exact_s_first=sp.simplify((-8*n*sp.sqrt(5*d/16)-d/(32*sp.sqrt(5*d/16)))-s1)==0,
        exact_positive_bracket=sp.simplify(bracket-positive)==0,
        exact_moving_inverse_derivative=sp.simplify(moving-sp.sqrt(2)*bracket/(sp.pi*x))==0,
        frozen_angle_wrong_by_nonunit_factor=sp.simplify(moving-frozen*V/x)==0 and sp.simplify(V-x)!=0,
        exact_lift_density=sp.trigsimp((sp.diff(Phi,theta)-density).subs(relation))==0,
        exact_mean_one_lift=sp.simplify(Phi.subs(theta,2*sp.pi))==1,
        exact_anchored_primitive=sp.trigsimp((sp.diff(P,theta)-ga*density).subs(relation))==0,
        primitive_inverse_identity=sp.simplify(P-(Phi-theta/(2*sp.pi))/2)==0,
        exact_fourth_chain_all_terms=sp.expand(chain-claimed)==0,
        highest_shape_coefficient=sp.expand(chain).coeff(z4)==h0*v_1**4,
        missing_cross_term_rejected=sp.expand(chain-(claimed-12*z2*h1*v_1*v2))!=0)
    return dict(checks={k:bool(v) for k,v in checks.items()},raw_fourth_chain=str(sp.expand(chain)),
        first_target='-sqrt(delta)*(2sqrt(5)*n+1/(8sqrt(5)))',
        primitive_lower_bound='sqrt(10)*n*sqrt(delta)/(6pi)',
        global_pressure_compatibility_not_asserted=True)

def inequality_review():
    # Exact endpoint reductions; these are bounds for ALL 0<delta<=1,
    # not numerical sampling of the actual tiny deltaL.
    Q=sp.Rational
    checks=dict(native_vs_norm_upper_below_three=2+Q(5,16)+Q(1,8)<3,
        native_bs_value_upper_below_two=1+Q(5,16)<4,
        native_bs_positive_jet_sum_below_one_seventh=sum(Q(1,8)**j for j in range(1,5))<Q(1,7),
        native_bs_whole_norm_below_three=2+Q(1,7)<3,
        x_upper_below_three=2+Q(3,16)<3,
        b_lower_at_least_one=Q(3,16)>0, # b^2-1=(3/16)*delta
        twice_x_minus_r_exceeds_V=2-Q(3,64)>0,
        positive_bracket_remainder_coefficient=Q(7,1024)>0)
    return dict(checks={k:bool(v) for k,v in checks.items()},
        analytical_delta_domain='all 0<delta<=1; endpoints bound affine polynomials',
        actual_G_H_N_or_deltaL_materialized=False)

def diagnostic(delta,n):
    validate_parameters(delta,n)
    with mp.workdps(80):
        d=mp.mpf(delta);nn=mp.mpf(n)
        def state(eta):
            vs=2+3*d/16+d*eta/8
            z=mp.mpf('0.5')+nn*eta
            step=mp.mpf(0) if z<=0 else mp.mpf(1) if z>=1 else 1/(1+mp.exp(1/z**2-1/(1-z)**2))
            s=(1-step)*mp.sqrt(2+d/2-vs)
            return vs,-mp.sqrt(vs-1),s,vs+s*s
        def integrands(eta):
            vs,t,s,v=state(eta)
            def density(theta):return (1+(t+mp.sqrt(2)*s*mp.sin(theta))**2)/(2*mp.pi*v)
            def source(theta):return -(v/(1+(t+mp.sqrt(2)*s*mp.sin(theta))**2)-1)/2
            return density,source
        w0,_=integrands(mp.mpf(0));phi=mp.quad(w0,[0,mp.pi])
        def independently_solved_primitive(eta):
            w,g=integrands(eta)
            theta=mp.findroot(lambda th:mp.quad(w,[0,th])-phi,(mp.pi-mp.mpf('.01'),mp.pi+mp.mpf('.01')))
            return mp.quad(lambda th:g(th)*w(th),[0,theta])
        h=mp.mpf('1e-14')/nn
        measured=(independently_solved_primitive(h)-independently_solved_primitive(-h))/(2*h)
        x,t,s,V=state(mp.mpf(0));s1=-mp.sqrt(d)*(2*mp.sqrt(5)*nn+1/(8*mp.sqrt(5)))
        bracket=-d*s/(16*mp.sqrt(x-1))+t*s1-t*s*(d/8+2*s*s1)/V
        expected=mp.sqrt(2)*bracket/(mp.pi*x);frozen=expected*x/V
        return dict(delta=str(delta),n=n,digits=80,relative_error_limit='1e-20',
            checks=dict(independent_quadrature_inverse_derivative=bool(abs(measured/expected-1)<mp.mpf('1e-20')),
                frozen_phase_rejected=bool(abs(measured/frozen-1)>mp.mpf('1e-5')),
                primitive_oscillation_lower_bound=bool(expected>mp.sqrt(10)*nn*mp.sqrt(d)/(3*mp.pi)),
                positive_phase_density=bool(w0(mp.pi)>0),fixed_external_phase_in_period=bool(0<phi<1)))

@lru_cache(None)
def review():
    symbolic=symbolic_review();inequalities=inequality_review()
    rows=[diagnostic(d,n) for d in ('0.0625','0.00390625') for n in (1,4,16)]
    checks=dict(symbolic['checks']);checks.update(inequalities['checks'])
    checks.update({f'quadrature_delta{r["delta"]}_n{r["n"]}':all(r['checks'].values()) for r in rows})
    return dict(symbolic=symbolic,inequalities=inequalities,numerical=dict(rows=rows,manufactured_floating_not_interval_or_actual_profile=True),checks=checks)
