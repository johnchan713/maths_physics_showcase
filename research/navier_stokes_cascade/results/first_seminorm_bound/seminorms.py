"""Explicit finite first-coefficient bounds; rational controls are not actual profiles."""
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

def collar_jets(C,M,d,a):
    """M=(z+beta)exp(beta) on a remaining eta strip of width >d.

    a is sqrt(Xa), NOT the larger analytic-collar endpoint in z.
    C is the actual collar operator envelope, not the protected-core C^16.
    Bounds are C_eta^1 factorial norms. No radial analyticity is assumed.
    """
    C,M,d,a=map(rat,(C,M,d,a))
    if C<1 or min(M,d,a)<=0 or d>1:raise ValueError('Positive collar bounds and 0<d<=1 required')
    J=M*(1+1/d);Jeta=M*(1/d+2/d**2)
    H=C*(1+1/d)*(J+Jeta+1)
    return dict(value=J,radial=J/(2*a),radial2=H/(4*a*a)+J/a**3,
                angular_derivative_norm=Jeta,system_rhs_norm=H)

def moment_bounds(xcut,seed,f0,u0,omega_integrals):
    """Tentatively cut F1,U1; full w0=Omega0/X integrates through Xv."""
    xcut,seed,f0,u0=map(rat,(xcut,seed,f0,u0));w,v=map(rat,omega_integrals)
    if xcut<=0 or min(seed,f0,u0,w,v)<0:raise ValueError('Positive cutoff radius and nonnegative norms required')
    return (xcut*seed,xcut*xcut*seed,2*f0*xcut*seed+w/2,
            xcut*xcut*(u0+f0)*seed,2*u0*xcut*seed+f0*xcut*xcut*seed+v/2)

def repair_bound(moments,rho_root,e_star,lam):
    """Use sqrt(2)>1 to enclose K=1/(sqrt(2)rho^(3/2)).

    rho is specified as r^2, so rho_root=r=sqrt(rho) is the exact input.
    Full C_eta^1 reciprocal cost is 4/e_star, not 1/e or 2/e_star.
    """
    r,e_star,lam=map(rat,(rho_root,e_star,lam));mu=tuple(map(rat,moments))
    if len(mu)!=5 or min(mu)<0 or min(r,e_star,lam)<=0 or lam>F(1,100):
        raise ValueError('Five nonnegative moments and actual positive patch scales required')
    R=r*r;K=1/r**3;w=4/e_star;m,i,cp,j,s=mu
    z=(m/R,(m/R+K*w*j)/lam,K*i,w*s/R,w*cp)
    return dict(transformed=z,coefficient=800*max(z))

def extended_jets(core,coefficient,a,k1,k2,p,rp,bumps):
    """Xa=a^2; p lower patch radius; 0<rp<=sqrt(2p); Bj are y-jet bounds."""
    coefficient,a,k1,k2,p,rp=map(rat,(coefficient,a,k1,k2,p,rp));B0,B1,B2=map(rat,bumps)
    if min(a,p,rp)<=0 or rp*rp>2*p or min(coefficient,k1,k2,B0,B1,B2)<0:
        raise ValueError('Nonnegative jet bounds and valid positive radius enclosures required')
    J,Jx,Jxx=(core[n] for n in ('value','radial','radial2'))
    return dict(value=J+coefficient*B0*max(F(1),1/rp),
        radial=Jx+k1*J+coefficient*max(B1/p,(B1+B0/2)/(p*rp)),
        radial2=Jxx+2*k1*Jx+k2*J+coefficient*max((B2+B1)/p**2,(B2+2*B1+F(3,4)*B0)/(p*p*rp)))

def T(k,f,fx,fe,X,eta,h):
    D=F(1,2)-h;L=1-2*h*eta*eta
    return (-k*f+X*fx+D*eta*fe)/L

def Z(k,f,fx,fe,X,eta,h):
    return (2*k*eta*f-2*eta*X*fx+(1-eta*eta)*fe)/(1-2*h*eta*eta)

def first_sources(h,nu,X,eta,leading,first):
    """Exact order-one coefficients, with rtheta1=R*theta and rz1=axial.

    F/U tuples are (value, X, XX, eta); P is (value, X, eta).
    CF,CU are full shifted axial second derivatives of the leading fields.
    """
    h,nu,X,eta=map(rat,(h,nu,X,eta));A=F(1,2)+h;lam=2*h;b=-A-F(1,2);c=-A
    f,fx,fxx,fe=leading['F'];u,ux,uxx,ue=leading['U'];v=leading['V']
    g,gx,gxx,ge=first['F'];w,wx,wxx,we=first['U'];v1=first['V'];p,px,pe=first['P']
    theta=(T(b+lam,g,gx,ge,X,eta,h)-2*nu*(X*gxx+2*gx)
        +v*(gx+g/X)+v1*(fx+f/X)+u*Z(b+lam,g,gx,ge,X,eta,h)
        +w*Z(b,f,fx,fe,X,eta,h)-nu*leading['CF'])
    axial=(T(c+lam,w,wx,we,X,eta,h)-2*nu*(X*wxx+wx)
        +Z(-2*A+lam,p,px,pe,X,eta,h)+v*wx+v1*ux
        +u*Z(c+lam,w,wx,we,X,eta,h)+w*Z(c,u,ux,ue,X,eta,h)-nu*leading['CU'])
    return theta,axial

def source_bounds(h,nu,a,c,jets,f0,u0,f0x,u0x,v0,omega,CF,CU):
    """Leading f0,u0 are C_eta^1 norms on [0,Xc]; other data per README."""
    h,nu,a,c,f0,u0,f0x,u0x,v0,omega,CF,CU=map(rat,(h,nu,a,c,f0,u0,f0x,u0x,v0,omega,CF,CU))
    if not 0<h<F(1,100) or min(nu,a,c)<=0 or c<a or min(f0,u0,f0x,u0x,v0,omega,CF,CU)<0:
        raise ValueError('Admissible h, viscosity, radii and leading bounds required')
    S,Sx,Sxx=(rat(jets[n]) for n in ('value','radial','radial2'))
    if min(S,Sx,Sxx)<0:raise ValueError('Nonnegative first jet bounds required')
    A=F(1,2)+h;D=F(1,2)-h;lam=2*h;b=-A-F(1,2);k=-A;L=1-2*h;X=c*c
    time=lambda t:((abs(t)+D)*S+X*Sx)/L
    axial=lambda t:((2*abs(t)+1)*S+2*X*Sx)/L
    ZF=((2*abs(b)+1)*f0+2*X*f0x)/L
    ZU=((2*abs(k)+1)*u0+2*X*u0x)/L
    V1=X*(3+2*A)*S/L
    PX=2*f0*S+omega/2;P=X*PX;ZP=5*X*PX/L
    theta=time(b+lam)+2*nu*(X*Sxx+2*Sx)+v0*(Sx+S/(a*a))+V1*(f0x+f0/(a*a))+u0*axial(b+lam)+S*ZF+nu*CF
    z=time(k+lam)+2*nu*(X*Sxx+Sx)+ZP+v0*Sx+V1*u0x+u0*axial(k+lam)+S*ZU+nu*CU
    return dict(theta=theta,axial=z,radial_velocity=V1,pressure=P,pressure_X=PX,pressure_axial=ZP)

def compact_gate(a,c,jets,sources,stress_floor):
    """T1=0 through Xinner>Xa; floor applies ONLY on [Xinner,Xc]."""
    a,c,tau=map(rat,(a,c,stress_floor))
    if min(a,c,tau)<=0 or c<a:raise ValueError('Positive ordered radii and active stress floor required')
    E=2*c*jets['value'];EX=2*c*jets['radial']+jets['value']/a;UX=jets['radial']
    W=(c*c*sources['theta']/2+c*sources['axial'])/tau
    CE,CEX,CUX,CW=E+F(1,2),EX+F(1,2),UX+F(1,2),W+F(1,2)
    Cval=max(F(1),5*CE/a,5*(2*c*(CEX+CUX)+CE/a),CW)
    return dict(A_E=E,A_EX=EX,A_UX=UX,A_W=W,Cval=Cval)

def legacy_source_review():
    """Independent full-series coefficient assembly; floating manufactured diagnostics."""
    root=Path(__file__).resolve().parent.parent/'physical_residual_budget';old=sys.path[:]
    sys.path.insert(0,str(root))
    try:
        def load(name,path):
            spec=importlib.util.spec_from_file_location(name,path);m=importlib.util.module_from_spec(spec)
            sys.modules[name]=m;spec.loader.exec_module(m);return m
        r=load('first_seminorm_physical_residual',root/'residual.py')
        fixtures=load('first_seminorm_physical_fixtures',root/'review.py')
        base,correction=fixtures.fixtures();worst=0.;cases=0
        for h,nu,X,eta in product((.005,.009),(.4,1.,1.7),(.2,.7,1.4),(-.6,0.,.5)):
            J0,J1=(r.profile_jets(p,X,eta,h,lam) for p,lam in ((base,0.),(correction,2*h)))
            cv=lambda x:F(str(float(x)))
            get=lambda j,key:tuple(cv(x) for x in (j[key].v,j[key].g[0],j[key].H[0,0],j[key].g[1]))
            A=.5+h
            leading=dict(F=get(J0,'F'),U=get(J0,'U'),V=cv(J0['V'].v),
                CF=cv(r.axial_second(-A-.5,J0['F'],X,eta,h)),CU=cv(r.axial_second(-A,J0['U'],X,eta,h)))
            first=dict(F=get(J1,'F'),U=get(J1,'U'),V=cv(J1['V'].v),
                P=tuple(cv(x) for x in (J1['Pi'].v,J1['Pi'].g[0],J1['Pi'].g[1])))
            actual=first_sources(cv(h),cv(nu),cv(X),cv(eta),leading,first)
            expected=r.series_coefficients([base,correction],X,eta,h,nu)
            for got,key in zip(actual,('theta','axial')):
                target=expected[key][1];worst=max(worst,abs(float(got)-target)/max(1.,abs(target)))
            cases+=1
        return dict(manufactured=True,case_count=cases,tested_relative_error_upper='2e-12',
            checks=dict(full_first_residual_matches_independent_series=bool(worst<2e-12)))
    finally:sys.path[:]=old

def review():
    # Symbolic derivative identities: no actual profile or gigantic parameter is sampled.
    import sympy as s
    X,xi=s.symbols('X xi',positive=True);P=s.Function('P')(xi);Q=s.Function('Q')(xi)
    HP,HQ=s.symbols('HP HQ');E=s.Function('E')(X)
    fxx=s.diff(P/(2*xi),xi)/(2*xi);uxx=s.diff(Q/(2*xi),xi)/(2*xi)
    radial=dict(F_regular_system=s.simplify(fxx.subs(s.diff(P,xi),HP-3*P/xi)-(HP/(4*xi**2)-P/xi**3))==0,
        U_regular_system=s.simplify(uxx.subs(s.diff(Q,xi),HQ-Q/xi)-(HQ/(4*xi**2)-Q/(2*xi**3)))==0,
        E_to_F_first=s.simplify(s.diff(E/s.sqrt(2*X),X)-(s.diff(E,X)-E/(2*X))/s.sqrt(2*X))==0,
        E_to_F_second=s.simplify(s.diff(E/s.sqrt(2*X),X,2)-(s.diff(E,X,2)-s.diff(E,X)/X+3*E/(4*X**2))/s.sqrt(2*X))==0)
    # Verify the logarithmic chain rule using independent polynomial b(y), not a preset True.
    y=s.symbols('y');bp=y**5-2*y**3+3*y
    radial['logarithmic_bump_second']=s.simplify(s.diff(bp.subs(y,s.log(X)),X,2)-(s.diff(bp,y,2)-s.diff(bp,y)).subs(y,s.log(X))/X**2)==0
    # Physical radial measure, including the extra R in rtheta1.
    R=s.symbols('R',positive=True);q=s.symbols('q');v=q**3-2*q+1
    It=s.integrate(q*v,(q,0,X));Iz=s.integrate(v,(q,0,X))
    theta=-It/X;axial=-Iz/s.sqrt(2*X)
    jacobian=dict(theta_measure=s.simplify(s.diff((R**2*theta).subs(X,R**2/2),R)+R**3*v.subs(q,R**2/2))==0,
        axial_measure=s.simplify(s.diff((R*axial).subs(X,R**2/2),R)+R*v.subs(q,R**2/2))==0,
        missing_theta_X_weight_detected=s.simplify(theta+Iz/X)!=0,
        missing_axial_sqrt_detected=s.simplify(axial+Iz/(2*X))!=0)
    # Coherent upper bound for the diagnostic Ccol=2,a_col=2,Delta=1/2:
    # e<3 gives beta<384 and (8+beta)exp(beta)<2^9*2^768<2^800.
    # The radius/strip data are illustrative, NOT the selected construction.
    picard=dict(Ccol='2',a_col='2',Delta='1/2',d='1/4',sigma='1',
        M_upper='2^800',beta_upper='384',
        checks=dict(remaining_strip=F(1,2)+F(1,4)<1,
                    exact_elementary_exponential_enclosure=8+384<2**9 and 9+768<800,
                    cutoff_inside_analytic_radius=2<2**2))
    core=collar_jets(2,2**800,F(1,4),1)
    moments=moment_bounds(2,core['value'],5,7,(11,13));rep=repair_bound(moments,2,3,F(1,1000))
    ext=extended_jets(core,rep['coefficient'],1,2,5,3,2,(7,11,17))
    sources=source_bounds(F(1,200),1,1,4,ext,5,7,11,13,17,19,23,29)
    gate=compact_gate(1,4,ext,sources,F(1,10))
    mem0=moment_bounds(2,core['value'],5,7,(0,0));rep0=repair_bound(mem0,2,3,F(1,1000))
    # Changes to pressure memory need not change the max-row amplitude bound;
    # compare the actual normalized row bounds, not just their common maximum.
    memory=dict(full_m3=str(moments[2]),erased_m3=str(mem0[2]),full_m5=str(moments[4]),erased_m5=str(mem0[4]),
        different_third_and_fifth_rows=rep['transformed'][4]!=rep0['transformed'][4] and rep['transformed'][3]!=rep0['transformed'][3])
    L=dict(F=(F(2),F(1),F(0),F(3)),U=(F(1),F(-2),F(1),F(2)),V=F(-3),CF=F(5),CU=F(-7))
    rows=[]
    for eta,sign in product((F(-1),F(-1,3),F(0),F(2,5),F(1)),(-1,1)):
        J=dict(F=tuple(F(sign)*v for v in (1,2,3,4)),U=tuple(F(sign)*v for v in (2,3,4,5)),
            V=F(sign)*7,P=tuple(F(sign)*v for v in (3,4,5)))
        exact=first_sources(F(1,200),1,F(2),eta,L,J)
        B=source_bounds(F(1,200),1,1,2,dict(value=7,radial=4,radial2=4),5,3,1,2,3,100,5,7)
        noP=dict(J,P=(F(0),)*3);omitted=first_sources(F(1,200),1,F(2),eta,L,noP)
        rows.append(dict(eta=str(eta),sign=sign,exact=list(map(str,exact)),
            checks=dict(theta_triangle=abs(exact[0])<=B['theta'],axial_triangle=abs(exact[1])<=B['axial'],
                        pressure_does_not_enter_theta=exact[0]==omitted[0],pressure_enters_axial=exact[1]!=omitted[1])))
    legacy=legacy_source_review()
    return dict(exact_radial_and_E_to_F_identities=radial,exact_stress_Jacobians=jacobian,independent_series_diagnostic=legacy,
        illustrative_Picard_enclosure=picard,
        rational_propagation_example=serial(dict(core=core,moments=moments,repair=rep,extended=ext,sources=sources,gate=gate)),
        full_pressure_memory_control=memory,source_triangle_controls=rows,
        checks=dict(radial_and_cutoff_derivative_algebra=all(radial.values()),stress_integral_measures=all(jacobian.values()),
            illustrative_Picard_upper_enclosure=all(picard['checks'].values()),
            all_rational_source_triangle_controls=all(all(r['checks'].values()) for r in rows),
            complete_first_source_matches_legacy_series=all(legacy['checks'].values()),
            future_pressure_source_changes_normalized_repair_rows=memory['different_third_and_fifth_rows'],
            all_enormous_but_finite_example_constants_positive=gate['Cval']>gate['A_W']>0,
            endpoint_stress_floor_not_taken_at_Xa=F(9,8)>1))
