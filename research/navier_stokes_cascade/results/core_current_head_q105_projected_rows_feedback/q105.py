"""Final projected mean rows and mean-wave feedback, with source stages kept."""
from pathlib import Path
from fractions import Fraction as Q
import importlib.util
p=Path(__file__).resolve().parent.parent/'core_current_head_q104_final_radial_pressure/q104.py'
s=importlib.util.spec_from_file_location('q105_radial',p);radial=importlib.util.module_from_spec(s);s.loader.exec_module(radial)
full=radial.full;EP=radial.EP;cv=radial.cv

def projected_join(first,second,tail):
    return dict(Qtheta=2*(first['Etheta']+second['Etheta']),
                Qz=2*(cv(tail)+first['Ez']+second['Ez']))

class ProjectedRows:
    def __init__(self,final_radial):
        if not isinstance(final_radial,radial.FinalRadial):raise TypeError('Q104 complete final radial factory required')
        if final_radial.maximum<1:raise ValueError('Reserve the next full pressure output')
        self.radial=final_radial;self.state=final_radial.state
        self.reset=final_radial.postslow;self.fresh=self.reset.fresh
        self.current=final_radial.current
    def pre_fields(self,N,bundle,p):
        c=self.current;b=c.initial.budget
        f0=b._fast(N);s0=c.initial.slow_fields(N,True)
        f1=c.fast(N);s1=c.slow(N,True)
        if f1['order']<N or f1['pole']!=c.fast.rows.pole:raise ValueError('Original pre-sourced-fast mean domain changed')
        take=lambda f:{k:cv(f[k]) for k in ('x','u','z')}
        parts=[('initial_fast',b.pole,take(f0)),('initial_slow',0,{k:cv(s0[j]) for k,j in (('x','X'),('u','U'),('z','Z'))}),
               ('first_real_fast',f1['pole'],take(f1)),('first_real_slow',0,take(s1))]
        return radial.field_sum(parts,c.current.target.Lhi,p)
    def bound(self,n,bundle):
        full.order(n)
        if n+1>self.radial.maximum:raise ValueError('Projected rows need next pressure output')
        self.reset.check(n+1,bundle);N=n+3
        p=self.radial.common_pole(bundle);c=self.current;b=c.initial.budget;cp=c.current.profile
        L=c.current.target.Lhi;edge=radial.edge
        E=1+cv(edge(p,self.state.target.width))
        data=b._data(N,b.pole)
        fixed={k:243*cv(data[k]) for k in ('b','V','G','c','graph')}
        op=self.radial.operator(n+1,p)
        rho=243*radial.rows.embedding(L,0,p)*cv(cp.weighted_rho(N))
        # This source is BEFORE the new fast/slow means: only the first four means.
        oldw=self.pre_fields(N,bundle,p)
        W=radial.rows.embedding(L,bundle['covariance_pole'],p)*radial.total_covariance(N,bundle)
        pre=radial.radial_envelope(n+1,oldw,W,cp.Rmin,fixed,E)
        preP=self.fresh.changed.bound(n+1,bundle)['changed']['P']
        primitive=op['KT']*(pre['reduced_source']+rho*preP)
        Z=E*oldw['z']
        H1=2*fixed['G']*Z+Z*Z+W+pre['S']+primitive
        axial=self.fresh.rows.axial(n,p)
        if axial['order']<n or axial['pole']!=p:raise ValueError('Exact pre-fast axial proof domain required')
        tail=2*(n+1)*axial['KA']*H1
        # Complete mean differences, with full shifted pressure at the SAME pole.
        old={k:cv(v) for k,v in self.fresh.means.fields(N).items()}
        fast=self.fresh.rows.fast(N,bundle);slow=self.fresh.slow(N,bundle,True)
        if fast['source_state'] is not self.state or fast['order']<N or fast['pole']!=self.fresh.rows.common_pole(bundle):
            raise ValueError('Actual complete fresh fast field required')
        if slow['source_state'] is not self.state or slow['order']<N or slow['weighted'] is not True:
            raise ValueError('Actual weighted fresh slow field required')
        F={k:243*E*radial.rows.embedding(L,fast['pole'],p)*cv(fast[k]) for k in ('x','u','z')}
        S={k:243*E*radial.rows.embedding(L,0,p)*cv(slow['fields'][k]) for k in ('x','u','z')}
        FU={k:243*cv(edge(fast['pole'],self.state.target.width))*cv(fast[k]) for k in ('x','u','z')}
        f=dict(fixed,graph=1+2*fixed['graph'],r=sum((cp.Rmin**(-j-1) for j in range(N+1)),Q(0)),
               Kp=op['KT'],rho=rho,Kmom=(cp.Rmax-cp.Rmin)*(1/cp.Rmin+1+cp.Rmax+cp.Rmax**2)*cv(edge(p,self.state.target.width)))
        first=full.refresh.post.fast.difference(n,old,F,f)
        second=full.refresh.post.fast.difference(n,{k:old[k]+FU[k] for k in old},S,f)
        out=projected_join(first,second,tail)
        out.update(fast_difference=first,slow_difference=second,pre_fast_H1=H1,pre_fast_A1=tail,
                   source_state=self.state,pole=p,order=n,tail_order='2-kappa',
                   angular_order='5/2-4*kappa',axial_order='min(2-kappa,5/2-4*kappa)',
                   norm_stage='final_sourced_postslow',pointwise_Haar_rows_supplied=False,actual_native_contracts_certified=False,
                   complete_residual_improvement=False)
        return out
    def feedback(self,n,bundle):
        full.order(n);self.reset.check(n+1,bundle);N=n+3
        if bundle['covariance_pole']!=2*bundle['wave_pole']:raise ValueError('Common total-wave pole required')
        fast=self.fresh.rows.fast(N,bundle);slow=self.fresh.slow(N,bundle,False)
        p=self.fresh.rows.common_pole(bundle)
        if fast['source_state'] is not self.state or fast['order']<N or fast['pole']!=p or slow['source_state'] is not self.state or slow['order']<N or slow['weighted'] is not False:
            raise ValueError('Complete current mean increments required')
        inc={k:243*(cv(radial.edge(p,self.state.target.width))*cv(fast[k])+cv(slow['fields'][k])) for k in ('x','u','z')}
        if len(bundle['old'])<=N or len(bundle['new'])<=N:raise ValueError('Complete total-wave jets through n+3 required')
        U=radial.pressure.through([cv(a)+cv(b) for a,b in zip(bundle['old'][:N+1],bundle['new'][:N+1])],N)
        wave=self.fresh.changed.changed.response.inputs.source.waves
        normal=wave.through(N)
        if normal['wave_state'] is not self.state or normal['order']<N:raise ValueError('Original same-state phase normal required')
        b=self.current.initial.budget;a=b._data(N,b.pole)
        meanwave=radial.covariance.response.inputs.source.waves.linear.meanwave
        coef=meanwave.envelope(n,sum(inc.values(),EP()),U,243*cv(a['graph']),cv(normal['normal']),243*cv(a['r']),harmonic=2)
        return dict(coefficient=coef,source_state=self.state,wave_pole=bundle['wave_pole'],order=n,
                    epsilon_order='3/2-2*kappa',complete_regular_source_bound_supplied=False,
                    actual_native_contracts_certified=False)

def record():
    return dict(status='final-projected-rows-and-mean-wave-feedback',
                angular_order='5/2-4*kappa',axial_order='min(2-kappa,5/2-4*kappa)',
                actual_pre_fast_potential_A1_order='2-kappa',mean_wave_feedback_order='3/2-2*kappa',
                tail_source='old_final_means_plus_total_sourced_wave_before_fresh_fast',
                both_axial_viscosities_and_full_pressure_memory_retained=True,
                pointwise_Haar_rows_or_full_regular_residual_supplied=False,
                actual_native_contracts_certified=False,unconditional_pulse_residual_or_blowup=False)
