"""Incoming final rows plus sourced changes; complete conditional fast map."""
from pathlib import Path
from fractions import Fraction as Q
import importlib.util
p=Path(__file__).resolve().parent.parent/'core_current_head_q100_sourced_full_row_changes/q100.py'
s=importlib.util.spec_from_file_location('q101_changes',p);changes=importlib.util.module_from_spec(s);s.loader.exec_module(changes)
full=changes.full;EP=changes.EP;cv=changes.cv

def embedding(L,source,destination):
    full.order(source);full.order(destination);L=Q(L)
    if L<=0 or destination<source:raise ValueError('Positive original hull and upward pole embedding required')
    return max(Q(1),L)**(destination-source)

def join(fast,slow,tail,delta,tail_embedding=1,delta_embedding=1):
    oldt=2*(fast['Etheta']+slow['Etheta'])
    oldz=2*(cv(tail_embedding)*cv(tail)+fast['Ez']+slow['Ez'])
    return dict(incoming_theta=oldt,incoming_z=oldz,
                Etheta=oldt+cv(delta_embedding)*delta['Qtheta'],
                Ez=oldz+cv(delta_embedding)*delta['Qz'])

class FinalRows:
    def __init__(self,sourced_changes,*,maximum,pole,c0,c1,shape,time_inverse):
        if not isinstance(sourced_changes,changes.SourcedRowChanges):raise TypeError('Q100 same-state sourced row changes required')
        full.order(maximum);full.order(pole)
        self.changes=sourced_changes;self.state=sourced_changes.state
        self.means=sourced_changes.scalars.means;self.current=self.means.current
        if self.means.postslow.max_output<maximum+3:raise ValueError('Incoming final differences need factory reserve maximum+3')
        self.old_rows=self.current.fast.rows
        self.maximum,self.pole=maximum,pole
        shape=tuple(shape)
        d=self.old_rows.pressure.d;maps=self.current.current
        self.pressure=full.compact(maps,d,c0,c1,shape,time_inverse,maximum+1,pole,0)
        self.axial=full.compact(maps,d,c0,c1,shape,time_inverse,maximum+1,pole,1)
    def common_pole(self,bundle):
        p=max(bundle['covariance_pole'],self.old_rows.pole)
        if p>self.pole:raise ValueError('Total row exceeds declared new operator archive')
        return p
    def bound(self,n,bundle):
        full.order(n)
        if self.means.state is not self.state or self.old_rows.current is not self.current.current:raise ValueError('Incoming final pressure state changed')
        if n>self.maximum:raise ValueError('Outside total changed-row reserve')
        p=self.common_pole(bundle)
        delta=self.changes.bound(n,bundle) # includes exact state/stage/pressure guard
        N=n+3;c=self.current;b=c.initial.budget;cp=c.current.profile
        L=c.current.target.Lhi;edge=full.refresh.maps.update.signed.maps.edge_absorption
        # Raw normalized field coefficients, not arbitrary outgoing-row callbacks.
        old={k:243*cv(v) for k,v in c.old_full(N).items()}
        fw=c.fast(N);full.order(fw['order']);full.order(fw['pole'])
        if fw['order']<N or fw['pole']!=self.old_rows.pole:raise ValueError('Original complete fast domain changed')
        fu={k:243*cv(v) for k,v in c.fast_unweighted(N).items()}
        sw={k:243*cv(v) for k,v in c.slow(N,True).items()}
        E=1+cv(edge(p,self.state.target.width))
        ef=embedding(L,fw['pole'],p);es=embedding(L,0,p)
        F={k:243*E*ef*cv(fw[k]) for k in ('x','u','z')}
        S={k:E*es*sw[k] for k in ('x','u','z')}
        fixed={k:243*cv(v) for k,v in b._data(N,b.pole).items()}
        fixed['graph']=1+2*fixed['graph']
        a,z=cp.Rmin,cp.Rmax
        fixed.update(Kp=self.pressure(n+1,p)['KT'],rho=243*es*cv(cp.weighted_rho(N)),
                     Kmom=(z-a)*(1/a+1+z+z*z)*cv(edge(p,self.state.target.width)))
        first=full.refresh.post.fast.difference(n,old,F,fixed)
        after={k:old[k]+fu[k] for k in old}
        second=full.refresh.post.fast.difference(n,after,S,fixed)
        # The source here is POSTSIGNED, not the primary initial_A1 callback.
        before=self.old_rows.rows(n)
        tail=self.current.fast.compact(n,self.old_rows.pole)['KA']*243*before['Ez']
        out=join(first,second,tail,delta,embedding(L,self.old_rows.pole,p),
                 embedding(L,bundle['covariance_pole'],p))
        out.update(pole=p,order=n,fast_difference=first,slow_difference=second,
                   post_signed_A1=tail,source_state=self.state,row_order='3/2-2*kappa',
                   actual_native_contracts_certified=False,means_updated=False)
        return out
    def fast(self,n,bundle):
        full.order(n)
        if n+5>self.maximum:raise ValueError('Fast inverse needs row reserve n+5')
        p=self.common_pole(bundle)
        out=full.fast.fast_recipe(n,lambda j:self.bound(j,bundle),self.axial,p)
        out.update(source_state=self.state,velocity_order='3/2-2*kappa',
                   radial_order='5/2-2*kappa',actual_native_contracts_certified=False,
                   fresh_slow_or_postfast_scalars_supplied=False,means_updated=False)
        return out

def record():
    return dict(status='incoming-final-row-and-sourced-fast-map-join',
                row_and_fast_order='3/2-2*kappa',fast_radial_order='5/2-2*kappa',
                same_pole_product_absorption_explicit=True,post_signed_A1_retained=True,
                old_operators_evaluated_at_their_exact_pole=True,
                original_time_inverse_S_and_axial_one_plus_KA_retained=True,
                actual_native_contracts_certified=False,fresh_slow_or_scalar_reset_supplied=False,
                unconditional_pulse_residual_or_blowup=False)
