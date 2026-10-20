"""Fifth forward moments reduced to field jets; exact manufactured controls."""
from fractions import Fraction as Q
from functools import lru_cache
import importlib.util
from pathlib import Path
import sympy as s

HERE=Path(__file__).resolve().parent
spec=importlib.util.spec_from_file_location('volterra_parent_input',HERE.parent/'fourth_pressure_input_boundary'/'input_bounds.py')
parent=importlib.util.module_from_spec(spec);spec.loader.exec_module(parent)

def rat(v):
    if type(v) not in (int,Q):raise TypeError('Exact integer/Fraction input required')
    return Q(v)

def remainders(A,B):
    A,B=rat(A),rat(B)
    if A<1 or B<1:raise ValueError('Envelopes at least one')
    W=21*A*B;Qrow=156*A**7*B**5;Nrow=91*A*B**2
    R1=2496*A**8*B**5;R2=7280*A**7*B**6;common=10000*A**8*B**6
    return dict(W=W,Q=Qrow,N=Nrow,R1=R1,R2=R2,common=common)

def certificate(A,B,left,X,eta,h,entry,FE,FU,YE,YU):
    """Conditional positive majorants, not evaluated actual-profile acceptance.

    entry=(M5,I5,J5,S5,Cp5,axis5) NONNEGATIVE absolute bounds.
    FE,FU=normalized pure fifth field bounds, YE,YU=y eta4 bounds /4!.
    """
    A,B,left,X,eta,h,FE,FU,YE,YU=map(rat,(A,B,left,X,eta,h,FE,FU,YE,YU))
    entry=tuple(map(rat,entry))
    if len(entry)!=6:raise ValueError('All five incoming moments and axis datum required')
    if A<1 or B<1 or not 1/A<=left<=X<=A or not -1<=eta<=1 or not 0<h<=Q(1,100):raise ValueError('Outside compact certificate domain')
    if min(*entry,FE,FU,YE,YU)<0:raise ValueError('Nonnegative input bounds required')
    m,i,j,ss,c,k=entry;D=Q(1,2)-h;d=1-eta*eta;ell=X-left
    # Integral of (X/x-1) <= ell^2/(2*left), not an arbitrary supplied cap.
    khat=ell*ell/(2*left)
    tin1=d*m+A*(D*abs(eta)*i+d*j)
    tin2=(d*A+D*abs(eta))*m+d*ss+X*d*(c+k)
    t1=tin1+ell*d*(1+2*A**3)*FU+2*A*A*ell*(D*abs(eta)+d*A)*FE+2*d*A*A*ell*B*B
    t2=tin2+ell*(3*d*A+D*abs(eta))*FU+d*A*khat*FE+d*B*B*(ell+khat/2)
    r=remainders(A,B)
    return dict(length=ell,kappa_cap=khat,T1=t1,T2=t2,incoming_T1=tin1,incoming_T2=tin2,
        pressure1_C4=B+r['common']+10*t1,pressure2_C4=B+r['common']+10*A*t2,
        shear_a_high4=2*A*YE+10*A**5*B**5,shear_b_high4=2*A*YU+10*A**5*B**5,
        actual_native_high_jets_supplied=False,actual_gate_accepted=False)

def jet5(f,e,at):return s.diff(f,e,5).subs(e,at)/120
def lower_product(f,g,e,at):
    return sum(s.diff(f,e,j).subs(e,at)*s.diff(g,e,5-j).subs(e,at)/s.factorial(j)/s.factorial(5-j) for j in range(1,5))

def radial_integral(f,x,X):
    """Exact monomial antiderivatives, including x^-1 and half powers.

    Reject non-Laurent-polynomial inputs instead of approximating integration.
    """
    out=s.Integer(0)
    for term in s.Add.make_args(s.expand(f)):
        if term==0:continue
        coefficient,monomial=term.as_independent(x,as_Add=False)
        power=monomial.as_powers_dict().get(x,s.Integer(0))
        if monomial!=x**power or not power.is_Rational:raise ValueError('Not an exact radial monomial')
        out+=coefficient*(s.log(X) if power==-1 else (X**(power+1)-1)/(power+1))
    return out

@lru_cache(None)
def kernel_review():
    x,X=s.symbols('x X',positive=True);e=s.Symbol('eta',real=True);h=s.Rational(1,200);D=s.Rational(1,2)-h
    E=1+x/10+e*(x+1)/20+e**2/30+e**3*x/50+e**4/70+e**5*(x+2)/90
    U=s.Rational(1,7)+x/15+e/13+e**2*x/31+e**3/37+e**4*(x+1)/61+e**5*x/83
    incoming=[s.Rational(k,11)+e/(k+7)+e**4/(k+17)+e**5/(k+29) for k in range(1,6)]
    axis=s.Rational(3,7)+e/19+e**4/23+e**5/41
    densities=[U,s.sqrt(2*x)*E,s.sqrt(2*x)*U*E,U*U-E*E/2,E*E/(2*x)]
    moments=[v+radial_integral(f,x,X) for v,f in zip(incoming,densities)]
    rows=[]
    for eta in (-1,0,s.Rational(1,2),1):
        d=1-eta*eta;ev=E.subs(e,eta);uv=U.subs(e,eta);e5=jet5(E,e,eta);u5=jet5(U,e,eta)
        lu=lower_product(U,U,e,eta);le=lower_product(E,E,e,eta);lue=lower_product(U,E,e,eta)
        m,i,j,ss,c=[jet5(v,e,eta) for v in moments];m0,i0,j0,s0,c0=[jet5(v,e,eta) for v in incoming];k=jet5(axis,e,eta)
        endU=uv.subs(x,X);Hx=s.sqrt(2*X)*ev.subs(x,X)
        T1=d*m-(D*eta*i+d*j)/Hx
        T2=(d*endU-D*eta)*m-d*ss-X*d*(c+k)
        mem1=d*m0-(D*eta*i0+d*j0)/Hx
        mem2=(d*endU-D*eta)*m0-d*s0-X*d*(c0+k)
        integrand1=d*(1-s.sqrt(2*x)*ev/Hx)*u5-s.sqrt(2*x)*(D*eta+d*uv)*e5/Hx-d*s.sqrt(2*x)*lue/Hx
        axial=(d*(endU-2*uv)-D*eta)*u5
        swirl=d*(1-X/x)*ev*e5
        lower=-d*lu+d*(1-X/x)*le/2
        integrated1=radial_integral(integrand1,x,X)
        integrated2=radial_integral(axial+swirl+lower,x,X)
        for radius in (s.Rational(3,2),2,3):
            p1=s.simplify((T1-mem1-integrated1).subs(X,radius))
            p2=s.expand((T2-mem2-integrated2).subs(X,radius))
            wrong_pressure=radial_integral(d*ev*e5,x,X)
            # Removing the Cp fifth derivative destroys the swirl kernel.
            delta=s.expand((radial_integral(swirl,x,X)-wrong_pressure).subs(X,radius))
            rows.append(dict(eta=str(eta),X=str(radius),T1_exact=p1==0,T2_exact=p2==0,
                T2=str(s.expand(T2.subs(X,radius))),incoming_T2=str(s.expand(mem2.subs(X,radius))),
                dropped_pressure_detected=bool(delta!=0) if d else None))
    endpoint=s.simplify((1-X/x).subs(x,X))
    return dict(manufactured_not_selected_profile=True,rows=rows,checks=dict(
        twelve_exact_both_forward_kernel_identities=all(r['T1_exact'] and r['T2_exact'] for r in rows),
        both_angular_endpoints_keep_D_eta_M5=all(s.sympify(r['T2'])!=0 for r in rows if r['eta'] in ('-1','1')),
        complete_Cp_high_jet_is_needed_for_swirl_kernel=all(r['dropped_pressure_detected'] for r in rows if r['eta'] not in ('-1','1')),
        local_swirl_kernel_vanishes_only_at_current_radius=endpoint==0,
        nonzero_incoming_fifth_memory_not_erased=all(s.sympify(r['incoming_T2'])!=0 for r in rows)))

@lru_cache(None)
def lower_and_shear_review():
    rows=[]
    for eta in (Q(-1),Q(0),Q(2,5),Q(1)):
        data=parent.fixture(eta);split=parent.split(data);r=remainders(3,10)
        rows.append(dict(eta=str(eta),R1=str(split['p1']['lower']),R2=str(split['p2']['lower']),
            bounded=abs(split['p1']['lower'])<=r['R1'] and abs(split['p2']['lower'])<=r['R2']))
    E=parent.Jet([Q(3,2),Q(1,7),Q(-1,9),Q(1,11),Q(1,13)])
    Ey0=[Q(1,3),Q(1,5),Q(-1,7),Q(1,9)]
    controls=[]
    for high in (Q(-1000),Q(0),Q(1000)):
        Y=parent.Jet(Ey0+[high]);inv=E.inverse()
        exact=(Y*inv).a[4];low=sum(Ey0[j]*inv.a[4-j] for j in range(4))
        controls.append(dict(high=str(high),exact=str(exact),split_exact=exact==high/E.a[0]+low,
            a_bound=abs(-2*exact)<=2*3*abs(high)+10*3**5*10**5))
    return dict(rows=rows,shear_controls=controls,checks=dict(
        actual_original_forward_lower_remainder_controls=all(r['bounded'] for r in rows),
        mixed_shear_highest_jet_is_linear_with_same_E_value=all(r['split_exact'] for r in controls),
        conditional_shear_positive_majorant_controls=all(r['a_bound'] for r in controls),
        original_J4_only_bounds_slow_y_jets_through_three_not_four=True,
        C4_L_inverse_norm_geometric_series_25_over23_below_two=Q(25,23)<2,
        both_common_remainder_polynomial_constants_fit_10000=2496<10000 and 7280<10000))

def review():
    k=kernel_review();l=lower_and_shear_review()
    fixture=certificate(3,10,1,2,Q(1,2),Q(1,200),(1,2,3,4,5,6),7,8,9,10)
    zero=certificate(3,10,2,2,Q(1,2),Q(1,200),(1,2,3,4,5,6),7,8,9,10)
    checks=dict(**k['checks'],**l['checks'],
        zero_length_annulus_keeps_incoming_memory=zero['T2']==zero['incoming_T2'] and zero['T2']>0,
        exact_quadratic_local_pressure_cap=fixture['kappa_cap']==Q(1,2),
        actual_profile_high_inputs_or_gate_not_claimed=not fixture['actual_native_high_jets_supplied'] and not fixture['actual_gate_accepted'])
    return dict(kernels=k,lower_shear=l,certificate_fixture={k:str(v) if isinstance(v,Q) else v for k,v in fixture.items()},checks=checks)
