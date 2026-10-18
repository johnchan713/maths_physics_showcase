"""Explicit activated-collar envelopes; exact controls, not actual-core samples."""
from fractions import Fraction as F
from itertools import product
import importlib.util
from pathlib import Path
import sys

def rat(x):
    if isinstance(x,(bool,float)):raise TypeError('Exact rational inputs required')
    return F(x)

def serial(x):
    if isinstance(x,F):return str(x)
    if isinstance(x,dict):return {k:serial(v) for k,v in x.items()}
    if isinstance(x,(tuple,list)):return [serial(v) for v in x]
    return x

def bounds(S,rho,Lambda,t,ell):
    """Return factors times powers of two, without materializing exponentials.

    K=radial_factor*2**exponent; H=jet_factor*2**exponent;
    Ccol=operator_factor*2**(2*exponent). The actual input parameters
    remain symbolic in the proof; calls here are exact rational fixtures.
    """
    S,rho,Lambda,t,ell=map(rat,(S,rho,Lambda,t,ell))
    if S<2 or not 0<rho<=1 or Lambda<=0 or not 0<t<=F(1,100) or not 0<ell<t/20:
        raise ValueError('Need S>=2, 0<rho<=1, Lambda>0, 0<t<=.01, 0<ell<t/20')
    a=4/Lambda;x=F(41,10)/Lambda;b=1+x
    B=2*(256*S+9+S/Lambda)*(1+Lambda)**2
    j1=25*S;j2=1600*S/t+1000*S*S
    q=50*ell*S;n=(q.numerator+q.denominator-1)//q.denominator
    radial=B*(1+ell*x)*(1+(x+j1)/a+
        (2*x+x*x+2*x*j1+j2+j1*j1+j1+64*x/t)/(a*a))
    delta=rho/(800*S);loss=1+1/delta+2/delta**2
    H=radial*loss;T=1+b+x
    return dict(Xa=a,X_upper=x,xi_upper=b,reference_jet=B,
        log_shear=50*S,log_shear_y=2000*S*S,J_y=j1,J_yy=j2,
        exponential_power=n,radial_factor=radial,delta=delta,
        angular_loss=loss,jet_factor=H,operator_factor=4096*T**3*H**2,
        operator_power=2*n,sigma=2*delta,Delta=delta/4,d=delta/4)

def primitive_operator_bounds(x,b,H):
    """Triangle bounds on the exact six-variable rows, viscosity fixed at one."""
    x,b,H=map(rat,(x,b,H))
    if min(x,b)<0 or H<1:raise ValueError('Nonnegative radii and H>=1 required')
    W=512*(1+x)*H**2;Z2=128*(1+x)**2*H
    return dict(v=12*H,w0=W,axial_second=Z2,
        A0=128*(1+x+b)*H,A1=32*(1+x)*H,
        source=(b+3*x)*W+2*Z2,
        common=4096*(1+x+b)**3*H**2)

def symbolic_review():
    import sympy as s
    X,b,H=s.symbols('X b H',nonnegative=True);v=s.symbols('v',nonnegative=True)
    common=4096*(1+X+b)**3*H**2
    targets=dict(A0=128*(1+X+b)*H,A1=32*(1+X)*H,
        source=(b+3*X)*512*(1+X)*H**2+256*(1+X)**2*H)
    checks={name:all(c>=0 for c in s.Poly(s.expand((common-target).subs(H,1+v)),X,b,v).coeffs()) for name,target in targets.items()}
    # The inherited cutoff is exp(-1/s^2), NOT exp(-1/s).
    z=s.symbols('z',positive=True);J=s.Function('J')(z);f=s.Function('f')(z)
    checks['activated_swirl_second_product_rule']=s.simplify(s.diff(f*s.exp(J),z,2)-s.exp(J)*(s.diff(f,z,2)+2*s.diff(f,z)*s.diff(J,z)+f*(s.diff(J,z,2)+s.diff(J,z)**2)))==0
    Xp=s.symbols('Xp',positive=True)
    # Independent polynomial chain rule, not a rearranged derivative identity.
    y=s.symbols('y');p=y**4-2*y**2+3*y
    checks['independent_log_chain_rule']=s.simplify(s.diff(p.subs(y,s.log(Xp)),Xp,2)-(s.diff(p,y,2)-s.diff(p,y)).subs(y,s.log(Xp))/Xp**2)==0
    checks['dropping_lower_log_derivative_changes_result']=s.simplify(s.diff(p.subs(y,s.log(Xp)),Xp,2)-s.diff(p,y,2).subs(y,s.log(Xp))/Xp**2)!=0
    # All constants below are rational upper bounds on |eta|<=101/100.
    e=F(101,100);d=1+e*e;L=F(97,100);lp=F(1,25)*e
    numerator=2*e+e+d
    numerator_eta=2+2*e+1+e+2*e+d
    checks['regular_v_radial_bound']=numerator/L<12
    checks['regular_v_eta_bound']=numerator_eta/L+lp*numerator/L**2<12
    checks['reciprocal_L_no_pole']=1-F(1,50)*e*e>L
    checks['native_log_shear_bound']=2*F(41,10)*5<50
    checks['native_log_shear_y_bound']=(2*F(41,10)*5+2*F(41,10)**2*5)/2+2*F(41,10)**2*25<2000
    checks['flat_squared_cutoff_derivative_64']=2*(8+8)<=64
    checks['maximum_x3_exp_minus_x2_on_x_ge_2_at_2']=3-2*2**2<0
    def positive(poly):return all(c>=0 for c in s.Poly(s.expand(poly.subs(H,1+v)),X,b,v).coeffs())
    Q=lambda z:s.Rational(z.numerator,z.denominator)
    ee,dd,LL=map(Q,(e,d,L))
    raw_w=12*H*(1+X+ee/2)/LL+144*(s.Rational(1,2)+X)*H**2+12*H**2*(2*ee*(1+X)+dd)/LL+24*(2+X)*H
    checks['all_radial_source_terms_fit_W']=positive(512*(1+X)*H**2-raw_w)
    axial_raw=3*X*(3*X+12)*H+3*(6*X+18)*H+6*(3*X+9)*H
    checks['all_shifted_axial_product_terms_fit_Z']=positive(128*(1+X)**2*H-axial_raw)
    transport=b*((s.Rational(1,2)+ee*H)/LL+6*H)
    jf=(1+X)*H;ju=X*H;zg=(3*X+9)*H
    row4=2*(2*(2*ee*H+1)/LL+12*H)+2*(4*ee*jf/LL+zg)+8*ee*jf/LL+2*transport
    row5=8*ee*X*H/LL+2*(2*(2*ee*H+1)/LL+zg+4*ee*ju/LL)+8*ee*ju/LL+8*ee/LL+2*transport
    a1row4=6*H/LL+4*dd*jf/LL
    a1row5=2*(3*H+dd*ju)/LL+2*dd*ju/LL+2*dd/LL
    checks['A0_row_five_with_all_terms']=positive((73+24*X+16*b)*H-row4)
    checks['A0_row_six_with_pressure']=positive((40+33*X+16*b)*H-row5)
    checks['A1_row_five_with_angular_couplings']=positive((17+10*X)*H-a1row4)
    checks['A1_row_six_with_pressure_coupling']=positive((12+10*X)*H-a1row5)
    return dict(checks=checks,polynomial_domination_uses_H_ge_one=True,
        cutoff='exp(-1/s^2)/(exp(-1/s^2)+exp(-1/(1-s)^2))',
        regularity='real radial C-infinity and eta holomorphic, not radial analytic')

def legacy_review():
    """Compare envelope with independent retained system on polynomial fixtures."""
    import numpy as np
    root=Path(__file__).resolve().parent.parent;old=sys.path[:]
    sys.path[:0]=[str(root/'first_correction_construction'),str(root/'physical_residual_budget')]
    try:
        def load(n,p):
            spec=importlib.util.spec_from_file_location(n,p);m=importlib.util.module_from_spec(spec);sys.modules[n]=m;spec.loader.exec_module(m);return m
        system=load('collar_envelope_actual_system',root/'first_correction_construction/system.py')
        fixtures=load('collar_envelope_manufactured_profiles',root/'physical_residual_budget/review.py')
        base,_=fixtures.fixtures();rows=[]
        for xi,eta,h in product((0.,.3,1.2),(-1.,0.,1.),(.005,.01)):
            X=xi*xi;H=1.
            for p in (base.F,base.U,base.average_U):
                for r,k in product(range(3),repeat=2):
                    q=p
                    for _ in range(r):q=q.derivative(0)
                    for _ in range(k):q=q.derivative(1)
                    H=max(H,abs(q(X,eta)))
            A0,A1,f=system.matrices(base,xi,eta,h,1.)
            B=primitive_operator_bounds(F(str(X)),F(str(xi)),F(str(H)))
            good=[np.max(np.sum(np.abs(A0),axis=1))<=float(B['A0']),
                np.max(np.sum(np.abs(A1),axis=1))<=float(B['A1']),np.max(np.abs(f))<=float(B['source'])]
            rows.append(dict(xi=str(xi),eta=str(eta),h=str(h),checks=dict(zip(('A0','A1','full_source'),map(bool,good)))))
        # Nonzero source pressure entry and shifted axial terms are retained.
        A0,A1,f=system.matrices(base,.8,.4,.005,1.)
        return dict(manufactured_not_actual_core=True,rows=rows,
            checks=dict(all_18_matrix_source_controls=all(all(r['checks'].values()) for r in rows),
                nonzero_full_pressure_source=bool(abs(f[3])>1e-5),
                nonzero_shifted_axial_source=bool(abs(f[4])>1e-5),
                pressure_is_a_separate_system_variable=bool(A0[5,3]!=0 and A1[5,3]!=0)))
    finally:sys.path[:]=old

def review():
    r=bounds(2,F(1,8),4,F(1,100),F(1,10000))
    small=bounds(2,F(1,8),4,F(1,1000),F(1,100000))
    x,b,H=F(7,3),F(5,2),F(11,3);q=primitive_operator_bounds(x,b,H)
    sym=symbolic_review();legacy=legacy_review()
    checks=dict(symbolic_all=all(sym['checks'].values()),legacy_all=all(legacy['checks'].values()),
        width_cost_visible=small['radial_factor']>r['radial_factor'],
        shrunk_parameter_strips=r['Delta']+r['d']<r['sigma'],
        xi_enclosure=r['xi_upper']**2>=r['X_upper'],
        activated_endpoint_enclosure=F(4)*F(100,99)<F(41,10),
        exponent_not_materialized=r['operator_power']==2*r['exponential_power'],
        common_dominates_each_row=all(q['common']>=q[k] for k in ('A0','A1','source')))
    return serial(dict(exact_fixture=r,width_refinement=small,
        symbolic=sym,retained_system_diagnostic=legacy,checks=checks))
