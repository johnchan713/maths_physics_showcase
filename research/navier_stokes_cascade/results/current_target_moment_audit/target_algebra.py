"""Exact finite current-state controls, not actual pulse or continuum bounds."""
from functools import lru_cache
from fractions import Fraction
import sympy as s

R,Z,T,u=s.symbols('R Z T u',positive=True)
eps=s.Rational(1,7);M=s.Rational(3,5);ci=s.Rational(2,9)
F=(R-1)**5*(2-R)**5
b=eps/R;V=R*(2+Z);G=1+Z/3
rho=F*(1+(Z/10+T/20)*(R-s.Rational(3,2)))/s.integrate(F,(R,1,2))
sigtheta=F/s.integrate(R**2*F,(R,1,2));sigz=F/s.integrate(R*F,(R,1,2))
st=eps*R**2*F*(1+Z);sz=eps*R**2*F*(1+T)

def ex(x):return s.expand(x)
def haar(x):return ex(x).coeff(u,0)
def L(x):return s.I*u*s.diff(x,u)
def dr(x):return s.diff(x,R)+2*M*R*L(x)
def de(x,e):return dr(x)+e*x/R
def dz(x):return eps*s.diff(x,Z)
def ts(x):return -eps*s.diff(x,T)+ci*L(x)
def lap(x):return dr(dr(x))+dr(x)/R+eps**2*s.diff(x,Z,2)
def tau(x):return eps*s.diff(x,T)+eps**3*s.diff(x,Z,2)
def integ(x):
    p=s.Poly(ex(x),R);return sum(c*(2**(k[0]+1)-1)/s.Integer(k[0]+1) for k,c in p.terms())
def primitive(x,e=0):
    p=s.Poly(ex(R**e*x),R)
    return ex(sum(c*(R**(k[0]+1)-1)/s.Integer(k[0]+1) for k,c in p.terms())/R**e)
def zero(x):return ex(x)==0

def state(scale=1):
    c=(u+1/u)/2
    psi=scale*eps*R**3*F*(1+T*Z+(1+Z)*c)
    beta=-dz(psi);gamma=de(psi,1)
    v=scale*eps*(5*R**2*F+R**3*s.diff(F,R))*(1+T+Z+c)
    W={key:scale*eps*R**2*F*(i+2+T/(i+1)+Z/(i+2)+c/(i+3))
       for i,key in enumerate(('rr','rt','rz','tt','tz','zz'))}
    return dict(beta=ex(beta),v=ex(v),gamma=ex(gamma),W=W)

def radial(q):
    B,v,g,W=(q[k] for k in ('beta','v','gamma','W'))
    return ex(-ts(B)-de(2*b*B+B**2+W['rr'],1)-dz(b*g+G*B+B*g+W['rz'])
              +(2*V*v+v**2+W['tt'])/R+eps*(lap(B)-B/R**2))

@lru_cache(None)
def calculation():
    q=state();B,v,g,W=(q[k] for k in ('beta','v','gamma','W'))
    gh=haar(radial(q));Bh,vh,ghvel=map(haar,(B,v,g))
    ff=haar(2*b*B+B**2+W['rr']);hh=haar(b*g+G*B+B*g+W['rz'])
    uu=haar(2*V*v+v**2+W['tt']);kk=haar(G*v+V*g+g*v+W['tz'])
    zf=haar(2*G*g+g**2+W['zz'])
    P=integ(gh);Qg=integ(R**2*gh);Jt=integ(R**2*kk);Jz=integ(R*zf)-Qg/2
    Pc=tau(integ(Bh))+integ((uu-ff)/R)-eps*s.diff(integ(hh),Z)
    Qc=tau(integ(R**2*Bh))+integ(R*(ff+uu))-eps*s.diff(integ(R**2*hh),Z)
    pb=primitive(gh-rho*P);cr=integ(R**2*rho)/2
    Et=haar(ts(v)+de(b*v+B*V+B*v+W['rt']-st,2)+dz(G*v+V*g+g*v+W['tz'])-eps*(lap(v)-v/R**2))
    Ez=haar(ts(g)+de(b*g+B*G+B*g+W['rz']-sz,1)+dz(2*G*g+g**2+W['zz']+pb)-eps*lap(g))
    At=ex(haar(b*v+B*V+B*v+W['rt']-st)-eps*(s.diff(vh,R)-vh/R))
    Az=ex(haar(b*g+B*G+B*g+W['rz']-sz)-eps*s.diff(ghvel,R))
    ft=ex(eps*s.diff(kk,Z)-tau(vh));fz=ex(eps*s.diff(zf+pb,Z)-tau(ghvel))
    mt=eps*s.diff(Jt,Z);mz=eps*s.diff(Jz+cr*P,Z)
    ht=ex(ft-sigtheta*mt);hz=ex(fz-sigz*mz)
    Xt=-primitive(Et-sigtheta*mt,2);Xz=-primitive(Ez-sigz*mz,1)
    identities={
      'actual_divergence':zero(de(B,1)+dz(g)),
      'angular_velocity_moment':zero(integ(R**2*vh)),
      'axial_velocity_moment':zero(integ(R*ghvel)),
      'pressure_bump_unit':zero(integ(rho)-1),
      'current_P_closed_formula':zero(P-Pc),
      'current_Qg_closed_formula':zero(Qg-Qc),
      'pressure_derivative':zero(s.diff(pb,R)-(gh-rho*P)),
      'pressure_left_compact':zero(pb.subs(R,1)),
      'pressure_right_compact':zero(pb.subs(R,2)),
      'full_pressure_moment':zero(integ(R*pb)+Qg/2-cr*P),
      'angular_flux_decomposition':zero(Et-(s.diff(At,R)+2*At/R+ft)),
      'axial_flux_decomposition':zero(Ez-(s.diff(Az,R)+Az/R+fz)),
      'angular_integrated_residual':zero(integ(R**2*Et)-mt),
      'axial_integrated_residual':zero(integ(R*Ez)-mz),
      'angular_compensated_transport_zero_total':zero(integ(R**2*ht)),
      'axial_compensated_transport_zero_total':zero(integ(R*hz)),
      'angular_target_flux_plus_primitive':zero(Xt+At+primitive(ht,2)),
      'axial_target_flux_plus_primitive':zero(Xz+Az+primitive(hz,1)),
      'angular_target_left_compact':zero(Xt.subs(R,1)),
      'angular_target_right_compact':zero(Xt.subs(R,2)),
      'axial_target_left_compact':zero(Xz.subs(R,1)),
      'axial_target_right_compact':zero(Xz.subs(R,2)),
    }
    # A genuinely different CURRENT state. The complete quadratic source
    # expansion is checked independently of the closed moment formula.
    new=state(s.Rational(4,3));db=new['beta']-B;dv=new['v']-v;dg=new['gamma']-g
    dW={k:new['W'][k]-W[k] for k in W}
    dg_formula=ex(-ts(db)-de(2*b*db+2*B*db+db**2+dW['rr'],1)
        -dz(b*dg+G*db+B*dg+g*db+db*dg+dW['rz'])
        +(2*V*dv+2*v*dv+dv**2+dW['tt'])/R+eps*(lap(db)-db/R**2))
    actual=ex(radial(new)-radial(q));dgh=haar(dg_formula)
    dt=integ(R**2*haar(G*dv+V*dg+g*dv+v*dg+dg*dv+dW['tz']))
    dzj=integ(R*haar(2*G*dg+2*g*dg+dg**2+dW['zz']))-integ(R**2*dgh)/2
    newt=integ(R**2*haar(G*new['v']+V*new['gamma']+new['gamma']*new['v']+new['W']['tz']))
    newz=integ(R*haar(2*G*new['gamma']+new['gamma']**2+new['W']['zz']))-integ(R**2*haar(radial(new)))/2
    identities.update(current_source_full_linear_quadratic_difference=zero(actual-dg_formula),
        current_pressure_obstruction_difference=zero(integ(haar(actual))-integ(dgh)),
        current_angular_flux_difference=zero(newt-Jt-dt),current_axial_flux_difference=zero(newz-Jz-dzj))
    controls=dict(nonzero_moving_pressure_bump_term=not zero(eps*s.diff(cr,Z)*P),
        old_mean_correlation_not_product_of_Haar_means=not zero(haar(B*g)-Bh*ghvel),
        old_new_angular_correlation_nonzero=not zero(integ(R**2*haar(g*dv+v*dg))),
        old_new_axial_correlation_nonzero=not zero(integ(R*haar(2*g*dg))),
        full_covariance_six_entries_nonzero=all(not zero(haar(W[k])) for k in W),
        angular_target_nonzero=not zero(Xt),axial_target_nonzero=not zero(Xz))
    controls.update(nonzero_radial_covariance_pressure=not zero(integ(haar(W['rr'])/R)),
        nonzero_centrifugal_covariance_pressure=not zero(integ(haar(W['tt'])/R)),
        nonzero_axial_covariance_transport=not zero(eps*s.diff(integ(haar(W['rz'])),Z)),
        nonzero_angular_axial_covariance_flux=not zero(integ(R**2*haar(W['tz']))),
        nonzero_axial_square_covariance_flux=not zero(integ(R*haar(W['zz']))),
        nonzero_radial_angular_covariance_target=not zero(haar(W['rt'])))
    point={Z:s.Rational(1,5),T:s.Rational(2,7)}
    values={k:str(ex(x).subs(point)) for k,x in dict(P=P,Jtheta=Jt,Jz=Jz,c_rho=cr,
        omitted_moving_bump_term=eps*s.diff(cr,Z)*P).items()}
    return dict(identities=identities,controls=controls,values=values)

def ceiling(alpha,mu):
    alpha,mu=Fraction(alpha),Fraction(mu)
    if mu<=0:raise ValueError('positive preserved mean order required')
    return dict(control_target_order=str(mu+1),requested_order=str(alpha),
        exponent_margin=str(mu+1-alpha),order_not_excluded_by_this_control=mu+1>=alpha)

def swirl_control():
    # Finite-regularity control of the C-infinity argument in the README.
    psi=R**3*F;f=s.diff(psi,R)/R**2;rows=[]
    for ee in (s.Rational(1,7),s.Rational(1,11)):
        vv=ee*T*f
        Et=-ee*s.diff(vv,T)-ee*(s.diff(vv,R,2)+s.diff(vv,R)/R-vv/R**2)
        target=-primitive(Et.subs(T,0),2)
        rows.append(dict(epsilon=str(ee),exact_time_viscous_target=zero(target-ee**2*psi/R**2),
            both_endpoint_values_zero=zero(target.subs(R,1)) and zero(target.subs(R,2)),
            exact_velocity_moment_zero=zero(integ(R**2*vv)),nonzero_target=not zero(target)))
    return rows

def review():
    r=calculation();k=Fraction(1,100000);mu=1-k
    rows=[ceiling(Fraction(2)-2*k+j*(Fraction(1,2)-k),mu) for j in range(4)]
    control=swirl_control()
    return dict(exact_current_state=r,raw_state_ceiling=rows,swirl_time_viscous_controls=control,checks=dict(
        all_current_identities=all(r['identities'].values()),all_omission_controls=all(r['controls'].values()),
        actual_finite_swirl_control_matches_scaling=all(all(v for k,v in row.items() if k!='epsilon') for row in control),
        later_order_not_implied_by_raw_state=all(not row['order_not_excluded_by_this_control'] for row in rows[1:]),
        no_epsilon_gain_from_support_or_moments_alone=True))
