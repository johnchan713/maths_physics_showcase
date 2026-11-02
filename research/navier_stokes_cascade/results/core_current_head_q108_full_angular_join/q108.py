"""Full current angular row from old cancellation, sourced Haar, fresh means."""
from pathlib import Path
from fractions import Fraction as Q
import importlib.util
HERE=Path(__file__).resolve().parent

def load(name,path):
    s=importlib.util.spec_from_file_location(name,HERE.parent/path);m=importlib.util.module_from_spec(s);s.loader.exec_module(m);return m
projected=load('q108_projected','core_current_head_q105_projected_rows_feedback/q105.py')
angular=load('q108_angular','core_current_head_q87_final_angular_adapter/q87.py')
haar=load('q108_haar','core_current_head_q107_sourced_haar_flux/q107.py')
radial=projected.radial;full=projected.full;EP=projected.EP;cv=projected.cv
native=radial.covariance.response.inputs.source.waves.linear.meanwave.means

def join(old,sourced,fast,slow,source_edge=1,final_edge=1):
    return cv(old)+cv(source_edge)*cv(sourced)+cv(final_edge)*(cv(fast)+cv(slow))

class FullAngular:
    def __init__(self,rows):
        if not isinstance(rows,projected.ProjectedRows):raise TypeError('Q105 current projected rows required')
        self.rows=rows;self.means=rows.fresh.means;self.state=rows.state
        self.check()
    def check(self):
        if native.FinalMeans(self.means.postslow).state is not self.state or self.means is not self.rows.fresh.means or self.rows.state is not self.state:
            raise ValueError('Same original first-real mean and sourced state required')
    def old(self,n):
        """Q87 coefficient on this EXACT Q86 instance, before the sourced wave."""
        full.order(n);self.check();f=self.means
        if n+3>f.postslow.max_output:raise ValueError('Old angular mean reserve n+3 required')
        c=f.current;N=n+3;b=c.initial.budget
        fixed={k:cv(v) for k,v in b._data(N,b.pole).items()}
        O=c.old_full(N);F=c.fast_unweighted(N);S=c.slow(N,False)
        first=angular.angular_difference(n,O,F,fixed)
        second=angular.angular_difference(n,{k:O[k]+F[k] for k in O},S,fixed)
        h=self.state.target.source.haar;edge=radial.edge;w=self.state.target.width
        initial=cv(edge(2*h.pole,w))*cv(h.rows(n)['Etheta'])
        signed=native.covariance.defect(self.state,n+1,target_rate=Q(2))
        if signed['order']<n+1:raise ValueError('Old signed covariance order n+1 required')
        delta=cv(edge(signed['pole'],w))*cv(signed['signed'])
        signed_haar=((n+1)*(fixed['graph']+1)+2*fixed['r'])*delta
        H=initial+signed_haar
        return dict(coefficient=243*(H+first+second),initial=initial,signed_haar=signed_haar,
                    first_difference=first,second_difference=second,order=n,source_state=self.state,
                    norm='unweighted-factorial-through-R',epsilon_order='2-2*kappa')
    def bound(self,n,bundle):
        full.order(n);self.check();self.rows.reset.check(n+1,bundle)
        old=self.old(n);pm=self.rows.fresh.rows.changes.pressure
        if pm.state is not self.state:raise ValueError('Current changed pressure state required')
        P=pm.bound(n+1,bundle)
        delta=haar.envelope(n,bundle['change'],pm.maps.profile.Rmin,P['pressure'])
        E=cv(radial.edge(bundle['covariance_pole'],self.state.target.width))
        return dict(pre_fast_haar=join(old['coefficient'],delta['theta'],0,0,E),old=old,
                    sourced_haar=delta,source_edge=E,order=n,source_state=self.state,
                    norm='unweighted-factorial-through-R',epsilon_order='3/2-kappa')
    def final(self,n,bundle):
        pre=self.bound(n,bundle);out=self.rows.bound(n,bundle)
        if out['source_state'] is not self.state or out['order']<n:raise ValueError('Same current full mean differences required')
        E=cv(radial.edge(out['pole'],self.state.target.width))
        coefficient=join(pre['pre_fast_haar'],0,out['fast_difference']['Etheta'],out['slow_difference']['Etheta'],1,E)
        return dict(coefficient=coefficient,pre_fast=pre,final_edge=E,source_state=self.state,order=n,
                    norm='unweighted-factorial-through-R',epsilon_order='3/2-kappa',
                    full_angular_row=True,full_axial_row_supplied=False,
                    actual_native_contracts_certified=False,full_cone_or_residual_certified=False)

def rates(kappa=Q(1,100000)):
    k=Q(kappa)
    if not 0<k<Q(1,4):raise ValueError('Original kappa range required')
    vals=dict(old=2-2*k,sourced=Q(3,2)-k,fast=Q(5,2)-4*k,slow=Q(5,2)-3*k)
    return dict(vals,full_angular=min(vals.values()),gap_to_old_target=Q(1,2)-k)

def record():
    return dict(status='current-full-angular-row-from-incoming-cancellation-and-sourced-haar',
                full_angular_order='3/2-kappa',old_angular_order='2-2*kappa',
                old_program='Q87 on identical Q86 state',sourced_change_added_once=True,
                each_weight_unweighted_at_its_own_pole=True,full_axial_row_supplied=False,
                actual_native_contracts_certified=False,full_cone_pulse_residual_or_blowup=False)
