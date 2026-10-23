"""Leading directional enclosures from the unchanged regional construction.

Rational evaluations are controls, not numerical enclosures of the selected
giant profile. K below is the explicit collar expression in README (1).
"""
from fractions import Fraction as Q
from pathlib import Path
import importlib.util

HERE=Path(__file__).resolve().parent
def load(name,path):
    s=importlib.util.spec_from_file_location(name,HERE.parent/path)
    m=importlib.util.module_from_spec(s);s.loader.exec_module(m);return m
parent=load('leading_direction_parent','core_native_axial_transport_floor/bounds.py')
frame=load('leading_direction_frame','leading_frame_foundation/frame_floors.py')
exact=parent.exact

def positive(*xs):
    xs=tuple(map(exact,xs))
    if min(xs)<=0:raise ValueError('Positive primitive bounds required')
    return xs

def collar_recipe(core_S,rho,Lambda,t,ell_axis):
    S,rho,Lam,t,w=positive(core_S,rho,Lambda,t,ell_axis)
    if S<2 or rho>1 or t>Q(1,100) or w>=t/20:
        raise ValueError('Original core/collar range required')
    x=Q(41,10)/Lam;Xa=4/Lam;delta=rho/(800*S)
    B=2*(256*S+9+S/Lam)*(1+Lam)**2
    j1=25*S;j2=1600*S/t+1000*S*S
    pref=B*(1+w*x)*(1+(x+j1)/Xa+
        (2*x+x*x+2*x*j1+j2+j1*j1+j1+64*x/t)/Xa**2)
    return dict(x=x,Xa=Xa,delta=delta,K_prefactor=pref,K_exp_argument=25*w*S,
        K_expression='K_prefactor * exp(K_exp_argument)',
        actual_parameters_evaluated=False)

def compact_majorants(Ahi,H,N):
    a,H,N=positive(Ahi,H,N)
    if H<16 or a>=H or not H**32<N<2*H**32:
        raise ValueError('Retained hierarchy Ahi<H, H>=16, H^32<N<2H^32 required')
    repair=H**20/N
    # Joint J1 norm of P(y,eta,Ny)/N is <=H. The value is <=H/N<1.
    E=3*a*(1+H)+repair;U=a+H+repair
    return dict(E_direct=E,U_direct=U,E_budget=H**4,U_budget=H**2,
        repair_upper=repair,checks=dict(first_E_below_H4=E<H**4,
            first_U_below_H2=U<H**2,leading_repair_below_one=repair<1),
        actual_profile_measured=False)

def leading_enclosures(h,x,delta,K,H):
    h,x,d,K,H=positive(h,x,delta,K,H)
    if h>=Q(1,100) or H<16:raise ValueError('Original h and retained H required')
    A=Q(1,2)+h;b=1+x
    JE=b*K*(Q(3,2)+x+1/d);JU=K*(1+x+1/d)
    DEaxis=b*K*parent.directional_envelope(1+h+x,1/d)
    DUaxis=K*parent.directional_envelope(A+x,1/d)
    return dict(E_direction=max(DEaxis,2*H**4),U_direction=max(DUaxis,2*H**2),
        axis=dict(E_pointwise_jet_sum=JE,U_pointwise_jet_sum=JU,
            E_direction=DEaxis,U_direction=DUaxis),
        later=dict(E_pointwise_jet_sum=H**4,U_pointwise_jet_sum=H**2),
        scope='global positive native X, |eta|<=1, retaining all actual regional joins',
        K_required='same explicit collar K at ell_axis=tc/2; not a desired phase derivative',
        actual_parameters_evaluated=False)

def reference_enclosures(fmin,fmax,alpha,abar,bbar,nu,u_squared):
    vals=tuple(map(exact,(fmin,fmax,alpha,abar,bbar,nu,u_squared)))
    if vals[-1]<0:raise ValueError('Nonnegative original squared tilt required')
    r=frame.squared_floors(*vals[:-1]);u2=vals[-1]
    # 1/8 <= s^(-1-h) <= 8 on the original enlarged chart.
    # min(1,t)<=sqrt(t)<=max(1,t), t>0, avoids algebraic roots.
    return dict(gmin=min(Q(1),r['gmin2']/64),
        lambda_min=min(Q(1),r['lambda_min2']/64),
        lambda_max=max(Q(1),64*r['lambda_max2']),tilt_upper=1+u2,
        native_squared_floors=r,representatives='same original closed active annulus',
        actual_parameters_evaluated=False,original_tilt_reselected=False)

def native_first_jets(x,delta,K,H,Xlo):
    x,d,K,H,lo=positive(x,delta,K,H,Xlo)
    if H<16:raise ValueError('Retained H required')
    b=1+x
    return dict(E=dict(value=max(b*K,H**4),
        radial=max(b*K*(1+1/(2*lo)),H**4/lo),angular=max(b*K/d,H**4)),
        U=dict(value=max(K,H**2),radial=max(K,H**2/lo),angular=max(K/d,H**2)),
        original_positive_hull_Xlo=lo,actual_parameters_evaluated=False)

def comparison(enclosure,h,*,Rmin,Rmax,tilt_upper,lambda_max,gmin,Llow,Lhigh,Dbox):
    # Reuse the original reference/carrier/geometry validation, then replace
    # its independent ordinary-jet estimate with the DERIVED regional bound.
    zero=dict(value=0,radial=0,angular=0)
    r=parent.comparison(h,zero,zero,Rmin=Rmin,Rmax=Rmax,tilt_upper=tilt_upper,
        lambda_max=lambda_max,gmin=gmin,Llow=Llow,Lhigh=Lhigh,Dbox=Dbox)
    h=exact(h);lo=exact(Rmin);Lhi=exact(Lhigh)
    DE,DU=positive(enclosure['E_direction'],enclosure['U_direction'])
    c=r['carriers'];Hz=4/(1-2*h)*(c['theta']*DE/lo+c['axial']*DU)
    r['HZ_leading']=Hz;r['native_direction_E']=DE;r['native_direction_U']=DU
    r['coefficients']['A']=r['groups']['tilt']+Lhi*Hz
    r['leading_bound_source']='regional original-field/loop/repair/heat bounds, not fine H_Z callback'
    return r

def selector_margin(result,lambda_min,Lreq,r0,Cslow,tau):
    lam,L,r,Cslow,tau=positive(lambda_min,Lreq,r0,Cslow,tau)
    if lam>result['lambda_max'] or r>Q(1,6057072) or Cslow<1 or tau>1:
        raise ValueError('Original growth, radius and selector ranges required')
    Sfloor=max(Q(1),L/(2*r));c=result['coefficients']
    error=c['A']/Sfloor+c['B']/Sfloor**2+c['C']/Sfloor**3
    beta=lam/(4*(1+result['tilt_upper']**2)**2)
    return dict(existing_S_lower=Sfloor,retained_power_error=error,beta=beta,
        normal_margin=4*error**2<=beta,slow_upper_from_existing_selector=tau/Cslow,
        slow_gate_from_existing_selector=tau/Cslow<=1,inverse_normal_square=4/beta,
        supplied_Lreq=L,existing_selector_required=True,new_band_or_Lreq_selected=False,
        actual_same_band_floor_verified=False)

def certificate():
    return dict(parent_commit='102f6fed681ccffc474887df75c99f40d16ed259',
        compact='first native logarithmic jets E<H^4,U<H^2 with no N loss; second radial jets still see N',
        axis='K=explicit collar expression at tc/2; E=sqrt(2X)F, E_y=sqrt(2X)(F/2+XF_X)',
        whole_hull='take maxima of pointwise directional bounds, not max regional sums as a global sum of suprema',
        diagonal='same background_jet_foundation (10)-(11) includes m=1 for EVERY n>=1 and full cutoff-curl tuple',
        selector='original core_primary_tail_recheck (2)-(3) supplies slow gate; no ell>=512/h^2 identification needed on that already defined tail',
        obstruction='4[A/Sfloor+B/Sfloor^2+C/Sfloor^3]^2<=beta, Sfloor=max(1,Lreq/(2r0)), with SAME Lreq',
        actual_parameter_substitution_complete_residual_or_blowup_verified=False)

def certify_same_band():
    raise NotImplementedError('Actual same-Lreq normal-margin inequality, target freezing/root/cone and source bounds remain to be verified')
