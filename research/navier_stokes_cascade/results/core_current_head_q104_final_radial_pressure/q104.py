"""Complete final means/covariance and full shifted radial pressure bounds."""
from pathlib import Path
from fractions import Fraction as Q
import importlib.util
p=Path(__file__).resolve().parent.parent/'core_current_head_q103_current_haar_postslow/q103.py'
s=importlib.util.spec_from_file_location('q104_postslow',p);post=importlib.util.module_from_spec(s);s.loader.exec_module(post)
full=post.full;EP=post.EP;cv=post.cv
rows=post.slow.rows;pressure=rows.changes.pressure;covariance=pressure.covariance
edge=full.refresh.maps.update.signed.maps.edge_absorption

def total_covariance(n,bundle):
    full.order(n)
    if len(bundle['old'])<=n or len(bundle['new'])<=n:raise ValueError('Complete old and sourced raw wave jets required')
    U=[cv(a)+cv(b) for a,b in zip(bundle['old'][:n+1],bundle['new'][:n+1])]
    return pressure.through(covariance.prod(U,U),n)

def field_sum(parts,L,pole):
    out={k:EP() for k in ('x','u','z')}
    for name,source_pole,fields in parts:
        if set(fields)!=set(out):raise ValueError('All final mean components required in '+name)
        B=rows.embedding(L,source_pole,pole)
        for k in out:out[k]+=243*B*cv(fields[k])
    return out

def radial_envelope(n,mean,W,a,fixed,absorption):
    full.order(n);a=Q(a)
    if a<=0 or set(mean)!={'x','u','z'} or set(fixed)!={'b','V','G','c','graph'}:
        raise ValueError('Positive whole hull and complete unweighted fixed multipliers required')
    x,u,z=(cv(absorption)*cv(mean[k]) for k in ('x','u','z'))
    b,V,G,c,g=(cv(fixed[k]) for k in ('b','V','G','c','graph'));g=1+2*g
    W=cv(W);d=n+1;dd=(n+1)*(n+2);r=sum((a**(-j-1) for j in range(n+1)),Q(0))
    S=2*b*x+x*x+W;H=b*z+G*x+x*z+W;C=2*V*u+u*u+W
    time=d*(1+4*c)*x
    visc=(dd*(g*g+1)+r*d*g+r*r)*x
    reduced=time+r*S+d*H+r*C+dd*x
    local=S+(d*g+r)*x
    direct=time+(d*g+r)*S+d*H+r*C+visc
    return dict(S=S,H=H,C=C,local_pressure=local,reduced_source=reduced,
                direct_source=direct,radial_viscosity=visc)

class FinalRadial:
    def __init__(self,postslow,*,maximum,pole,c0,c1,shape,time_inverse):
        if not isinstance(postslow,post.CurrentPostslow):raise TypeError('Q103 current postslow reset required')
        full.order(maximum);full.order(pole)
        if postslow.fresh.rows.maximum<maximum+9:raise ValueError('Final radial budget needs current row reserve maximum+9')
        self.postslow=postslow;self.state=postslow.state
        self.current=postslow.current;self.maximum,self.pole=maximum,pole
        d=self.current.fast.rows.pressure.d
        self.operator=full.compact(self.current.current,d,c0,c1,shape,time_inverse,maximum,pole,0)
    def common_pole(self,bundle):
        p=max(self.current.initial.budget.pole,self.postslow.fresh.rows.common_pole(bundle))
        if p>self.pole:raise ValueError('Final radial source exceeds new compact proof archive')
        return p
    def fields(self,n,bundle):
        full.order(n)
        if n>self.maximum+2:raise ValueError('Outside final field reserve')
        if bundle['source_state'] is not self.state or bundle['means_updated'] or bundle['stage']!='post_first_sourced_wave_fixed_means':
            raise ValueError('Same original sourced wave bundle required for final mean construction')
        if n+7>self.postslow.fresh.rows.maximum or bundle['order']<n+9:
            raise ValueError('Complete final slow field requires rows n+7 and covariance n+9')
        p=self.common_pole(bundle);c=self.current;b=c.initial.budget
        f0=b._fast(n);s0=c.initial.slow_fields(n,True)
        f1=c.fast(n);s1=c.slow(n,True)
        f2=self.postslow.fresh.rows.fast(n,bundle);s2=self.postslow.fresh.slow(n,bundle,True)
        if f1['order']<n or f1['pole']!=c.fast.rows.pole or f2['order']<n or f2['pole']!=self.postslow.fresh.rows.common_pole(bundle) or f2['source_state'] is not self.state:
            raise ValueError('Actual old and new complete fast domains required')
        if s2['order']<n or s2['source_state'] is not self.state or s2['weighted'] is not True:
            raise ValueError('Current weighted complete fresh slow field required')
        select=lambda v:{k:cv(v[k]) for k in ('x','u','z')}
        parts=[('initial_fast',b.pole,select(f0)),
               ('initial_slow',0,{k:cv(s0[j]) for k,j in (('x','X'),('u','U'),('z','Z'))}),
               ('first_real_fast',f1['pole'],select(f1)),('first_real_slow',0,select(s1)),
               ('sourced_fast',f2['pole'],select(f2)),('sourced_slow',0,select(s2['fields']))]
        return dict(fields=field_sum(parts,c.current.target.Lhi,p),parts=parts,pole=p,order=n)
    def bound(self,n,bundle):
        full.order(n)
        if n>self.maximum:raise ValueError('Outside final radial output reserve')
        self.postslow.check(n,bundle)
        scalar=self.postslow.scalars(n,bundle)['scalars'];N=n+2
        mean=self.fields(N,bundle);p=mean['pole'];cp=self.current.current.profile
        W=rows.embedding(self.current.current.target.Lhi,bundle['covariance_pole'],p)*total_covariance(N,bundle)
        b=self.current.initial.budget;data=b._data(N,b.pole)
        fixed={k:243*cv(data[k]) for k in ('b','V','G','c','graph')}
        E=1+cv(edge(p,self.state.target.width))
        out=radial_envelope(n,mean['fields'],W,cp.Rmin,fixed,E)
        op=self.operator(n,p)
        rho=243*rows.embedding(self.current.current.target.Lhi,0,p)*cv(cp.weighted_rho(n))
        comp=rho*scalar['P'];source=out['reduced_source']+comp
        out.update(pressure=out['local_pressure']+op['KT']*source,
                   A0=op['KA']*source,compensator=comp,
                   radial_residual=comp+op['KA']*source,mean=mean,total_wave_covariance=W,
                   source_state=self.state,pole=p,order=n,norm_stage='final_sourced_postslow',
                   pressure_and_direct_radial_order='1-kappa',compensator_order='5/2-kappa',
                   whole_compact_flux_trace_contract_required=True,
                   actual_native_contracts_certified=False,complete_residual_improvement=False)
        return out

def record():
    return dict(status='complete-final-mean-and-full-radial-pressure-envelope',
                final_mean_contributions=6,complete_old_plus_sourced_covariance=True,
                local_viscous_and_covariance_pressure_retained=True,
                full_source_pressure_and_A0_coarse_order='1-kappa',
                compensator_order='5/2-kappa',whole_compact_trace_contract_required=True,
                actual_native_contracts_certified=False,complete_residual_or_blowup=False)
