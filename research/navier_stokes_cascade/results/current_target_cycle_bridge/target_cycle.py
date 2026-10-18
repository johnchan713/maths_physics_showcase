"""Exact manufactured current updates and rational all-half-line target budgets."""
from fractions import Fraction as F
from functools import lru_cache
import importlib.util
from pathlib import Path
import sys

@lru_cache(None)
def algebra():
    p=Path(__file__).resolve().parent.parent/'current_target_moment_audit/target_algebra.py'
    spec=importlib.util.spec_from_file_location('cycle_current_equation_algebra',p)
    m=importlib.util.module_from_spec(spec);sys.modules[spec.name]=m;spec.loader.exec_module(m);return m

def certificate(k=F(1,100000)):
    k=F(k)
    if not 0<k<F(1,100):raise ValueError('Use the reviewed positive small kappa range')
    C0=2-2*k;gain=F(1,2)-k
    orders=dict(signed_old_error=(1,F(1,2)-k),signed_full_curl=(1,F(1,2)-k),
        signed_square=(2,-1),signed_axial_pressure=(1,1),
        fast_old_new_means=(1,2-2*k),fast_square=(2,1-2*k),
        slow_linear_pressure_viscosity=(1,1),slow_old_new_means=(1,2-k),slow_square=(2,1))
    rows={n:dict(order_slope=str(a),order_offset=str(b),margin_slope=str(a-1),
        boundary_margin=str((a-1)*C0+b-gain)) for n,(a,b) in orders.items()}
    mean_viscosity_margin=F(1,2)-2*k
    return dict(kappa=str(k),C0=str(C0),gain=str(gain),rows=rows,full_mean_viscosity_margin=str(mean_viscosity_margin),
        checks=dict(all_half_line_slopes_nonnegative=all(F(r['margin_slope'])>=0 for r in rows.values()),
            all_boundary_margins_nonnegative=all(F(r['boundary_margin'])>=0 for r in rows.values()),
            both_generic_crosses_tie=rows['signed_old_error']['boundary_margin']==rows['signed_full_curl']['boundary_margin']=='0',
            fast_full_mean_viscosity_strictly_above_next=mean_viscosity_margin>0))

def prefix():
    k=F(1,100000);mu=rho=1-k;alpha=2-3*k
    real=min(alpha+1-2*k,2*alpha-1,alpha+1)
    realfast=min(mu+(alpha-k)+1,2*(alpha-k)+1);realslow=alpha+1
    B=1-k;C=B+F(1,2)
    generic=min(rho+C-F(1,2),C+F(1,2)-k,2*C-1,C+1)
    fast=min(mu+(C-k)+1,2*(C-k)+1);slow=C+1
    values=dict(initial_target=alpha,postreal_target=min(real,realfast,realslow),
        real_fast_target_change=realfast,real_slow_target_change=realslow,
        first_sourced_covariance=C,first_generic_target=generic,first_fast_target_change=fast,
        first_slow_target_change=slow,final_target_seed=min(generic,fast,slow))
    return dict(orders={n:str(v) for n,v in values.items()},checks=dict(
        real_prefix_current_target_stronger_than_first_sourced=values['postreal_target']>C,
        first_sourced_square_no_worse_than_cross=2*B>=C,
        current_generic_and_means_give_exact_seed=values['final_target_seed']==2-2*k,
        real_quadrature_shortcut_only_in_real_prefix=True))

@lru_cache(None)
def current_review():
    a=algebra();s=a.s;R,Z,T,u=a.R,a.Z,a.T,a.u;ep=a.eps
    ex,haar,I,integ,zero=a.ex,a.haar,a.primitive,a.integ,a.zero
    sig={}
    for e in (1,2):
        mass=integ(R**e*a.F);center=integ(R**(e+1)*a.F)/mass
        sig[e]=a.F*(1+(Z/13+T/17)*(R-center))/mass
    def Q(x,e):return ex(x-sig[e]*integ(R**e*x))
    def D(x,e):return ex(s.diff(x,R)+e*x/R)
    def engine(q,stress):
        beta,v,g,W=(q[n] for n in ('beta','v','gamma','W'))
        gh=haar(a.radial(q));P=integ(gh);pb=I(ex(gh-a.rho*P));cr=integ(R**2*a.rho)/2
        Jt=integ(R**2*haar(a.G*v+a.V*g+g*v+W['tz']))
        Jz=integ(R*haar(2*a.G*g+g**2+W['zz']))-integ(R**2*gh)/2
        Et=haar(a.ts(v)+a.de(a.b*v+beta*a.V+beta*v+W['rt']-stress[2],2)
            +a.dz(a.G*v+a.V*g+g*v+W['tz'])-ep*(a.lap(v)-v/R**2))
        Ez=haar(a.ts(g)+a.de(a.b*g+beta*a.G+beta*g+W['rz']-stress[1],1)
            +a.dz(2*a.G*g+g**2+W['zz']+pb)-ep*a.lap(g))
        return dict(gr=gh,P=P,p=pb,Jt=Jt,Jz=Jz,cr=cr,E={2:Et,1:Ez},X={2:-I(Q(Et,2),2),1:-I(Q(Ez,1),1)})
    q=a.state();start=engine(q,{1:a.sz,2:a.st})
    desired={2:ep**2*R**2*a.F*(1+Z+T),1:ep**2*R**2*a.F*(2+Z-T)}
    # Manufacture a FIXED stress, not an actual positive pulse covariance.
    stress={1:ex(a.sz+desired[1]-start['X'][1]),2:ex(a.st+desired[2]-start['X'][2])}
    old=engine(q,stress)
    error={2:ep**3*R**2*a.F*(Z+2*T),1:ep**3*R**2*a.F*(Z-T)}
    dw=dict(rr=ep**2*R**2*a.F*(1+Z**2+T),tt=ep**2*R**2*a.F*(3+Z*T),
        tz=ep**2*R**2*a.F*(2+Z+T),zz=ep**2*R**2*a.F*(5+Z**2+T),
        rt=desired[2]+error[2],rz=desired[1]+error[1])
    wave=dict(q,W={k:q['W'][k]+dw[k] for k in q['W']});postwave=engine(wave,stress)
    dg=ex(-D(dw['rr'],1)-a.dz(dw['rz'])+dw['tt']/R)
    dp=I(ex(dg-a.rho*integ(dg)))
    changeH={2:dw['tz'],1:dw['zz']+dp}
    predicted={e:-error[e]-I(Q(a.dz(changeH[e]),e),e) for e in (1,2)}
    ids=dict(engineered_current_target_angular=zero(old['X'][2]-desired[2]),
        engineered_current_target_axial=zero(old['X'][1]-desired[1]),
        complete_covariance_source_change=zero(postwave['gr']-old['gr']-dg),
        complete_covariance_pressure_change=zero(postwave['p']-old['p']-dp),
        actual_signed_angular_target_cancellation=zero(postwave['X'][2]-predicted[2]),
        actual_signed_axial_target_cancellation=zero(postwave['X'][1]-predicted[1]))
    for e in (1,2):
        H=changeH[e]
        ids[f'moving_compensator_unit_{e}']=zero(integ(R**e*sig[e])-1)
        ids[f'moving_compensator_commutator_{e}']=zero(Q(a.dz(H),e)-a.dz(Q(H,e))-ep*s.diff(sig[e],Z)*integ(R**e*H))
    controls=dict(axial_pressure_target_response_nonzero=not zero(I(Q(a.dz(dp),1),1)),
        moving_angular_compensator_commutator_nonzero=not zero(ep*s.diff(sig[2],Z)*integ(R**2*changeH[2])),
        moving_axial_compensator_commutator_nonzero=not zero(ep*s.diff(sig[1],Z)*integ(R*changeH[1])),
        moving_pressure_moment_product_nonzero=not zero(ep*s.diff(old['cr'],Z)*integ(dg)),
        signed_square_error_not_discarded=all(not zero(v) for v in error.values()))
    def curl_increment(scale,fast):
        c=(u+1/u)/2 if fast else 1
        psi=scale*ep**2*R**3*a.F*(1+T*Z+Z)*c
        return dict(beta=ex(-a.dz(psi)),gamma=ex(a.de(psi,1)),
            v=ex(scale*ep**2*(5*R**2*a.F+R**3*s.diff(a.F,R))*(1+T+Z)*c))
    def add_mean(q,d):return dict(beta=q['beta']+d['beta'],v=q['v']+d['v'],gamma=q['gamma']+d['gamma'],W=q['W'])
    fast=curl_increment(1,True);qfast=add_mean(wave,fast);postfast=engine(qfast,stress)
    beta,v,g=(wave[k] for k in ('beta','v','gamma'));db,dv,dga=(fast[k] for k in ('beta','v','gamma'))
    nonlinear_g=ex(-D(haar(2*beta*db+db**2),1)-a.dz(haar(beta*dga+g*db+db*dga))+haar(2*v*dv+dv**2)/R)
    pfast=I(ex(nonlinear_g-a.rho*integ(nonlinear_g)))
    dF={2:haar(beta*dv+db*v+db*dv),1:haar(beta*dga+db*g+db*dga)}
    dH={2:haar(g*dv+dga*v+dga*dv),1:haar(2*g*dga+dga**2)+pfast}
    ids['fast_complete_source_has_only_nonlinear_Haar_terms']=zero(postfast['gr']-postwave['gr']-nonlinear_g)
    ids['fast_pressure_recomputed_from_nonlinear_source']=zero(postfast['p']-postwave['p']-pfast)
    for e in (1,2):
        ids[f'fast_current_residual_difference_{e}']=zero(postfast['E'][e]-postwave['E'][e]-D(dF[e],e)-a.dz(dH[e]))
        ids[f'fast_current_target_difference_{e}']=zero(postfast['X'][e]-postwave['X'][e]+dF[e]+I(Q(a.dz(dH[e]),e),e))
    for name in ('beta','v','gamma'):ids['fast_exact_zero_Haar_'+name]=zero(haar(fast[name]))
    ids['fast_linear_time_Haar_zero']=zero(haar(a.ts(dv)))
    ids['fast_linear_viscosity_Haar_zero']=zero(haar(a.lap(dv)-dv/R**2))
    controls['fast_old_new_correlation_nonzero']=not zero(haar(2*v*dv))
    controls['fast_defect_actually_changed']=not zero(postfast['P']-postwave['P'])
    slow=curl_increment(ep,False);qslow=add_mean(qfast,slow);postslow=engine(qslow,stress)
    beta,v,g=(qfast[k] for k in ('beta','v','gamma'));db,dv,dga=(slow[k] for k in ('beta','v','gamma'))
    fullg=ex(-a.ts(db)-a.de(2*a.b*db+2*beta*db+db**2,1)
        -a.dz(a.b*dga+a.G*db+beta*dga+g*db+db*dga)
        +(2*a.V*dv+2*v*dv+dv**2)/R+ep*(a.lap(db)-db/R**2))
    dpfull=I(ex(haar(fullg)-a.rho*integ(haar(fullg))))
    et=haar(a.ts(dv)+a.de(a.b*dv+db*a.V+beta*dv+db*v+db*dv,2)
        +a.dz(a.G*dv+a.V*dga+g*dv+dga*v+dga*dv)-ep*(a.lap(dv)-dv/R**2))
    ez=haar(a.ts(dga)+a.de(a.b*dga+db*a.G+beta*dga+db*g+db*dga,1)
        +a.dz(2*a.G*dga+2*g*dga+dga**2+dpfull)-ep*a.lap(dga))
    ids['slow_full_current_source_difference']=zero(postslow['gr']-postfast['gr']-haar(fullg))
    ids['slow_full_current_pressure_difference']=zero(postslow['p']-postfast['p']-dpfull)
    for e,dE in ((2,et),(1,ez)):
        ids[f'slow_full_residual_difference_{e}']=zero(postslow['E'][e]-postfast['E'][e]-dE)
        ids[f'slow_current_target_difference_{e}']=zero(postslow['X'][e]-postfast['X'][e]+I(Q(dE,e),e))
    for name,d in (('fast',fast),('slow',slow)):
        ids[name+'_divergence']=zero(a.de(d['beta'],1)+a.dz(d['gamma']))
        ids[name+'_angular_velocity_moment']=zero(integ(R**2*haar(d['v'])))
        ids[name+'_axial_velocity_moment']=zero(integ(R*haar(d['gamma'])))
    for name,r in (('old',old),('wave',postwave),('fast',postfast),('slow',postslow)):
        ids[name+'_current_angular_moment_identity']=zero(integ(R**2*r['E'][2])-a.dz(r['Jt']))
        ids[name+'_current_axial_pressure_moment_identity']=zero(integ(R*r['E'][1])-a.dz(r['Jz']+r['cr']*r['P']))
    controls['slow_old_nonzero_Haar_means_retained']=not zero(integ(R**2*haar(g*dv+v*dga)))
    controls['slow_changed_pressure_nonzero']=not zero(dpfull)
    # Angular frequency and auxiliary Haar are deliberately DIFFERENT symbols.
    ang=s.symbols('angular_mode',positive=True)
    for m in (1,2,3):
        f=(1+R+Z+u)*ang**m+(2+T+1/u)*ang**(-m)
        ids[f'linear_cutoff_tail_zero_angular_mean_mode_{m}']=ex(s.diff(f,R)+s.diff(f,Z)+s.diff(f,T)+ang*s.diff(f,ang)).coeff(ang,0)==0
    return dict(identities=ids,controls=controls,
        covariance_fixture='manufactured complete tensor; fixed stress engineered once; not a positive canonical pulse realization',
        fast_fixture='exact zero-Haar curl, including its full axial/radial components',
        slow_fixture='moment-preserving slow curl, not a numerical current five-map solve')

def review():
    c=certificate();p=prefix();r=current_review()
    return dict(current_equation_cancellations=r,all_half_line_certificate=c,canonical_prefix_routing=p,
        checks=dict(all_exact_current_update_identities=all(r['identities'].values()),
            all_nonzero_omission_controls=all(r['controls'].values()),
            whole_half_line_target_budget=all(c['checks'].values()),accepted_initial_target_routes_to_seed=all(p['checks'].values())))
