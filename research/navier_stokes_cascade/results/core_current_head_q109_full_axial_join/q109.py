"""Full axial mean row; Haar removes old compact remainder, not its velocity."""
from pathlib import Path
from fractions import Fraction as Q
from math import factorial
import importlib.util
p=Path(__file__).resolve().parent.parent/'core_current_head_q108_full_angular_join/q108.py'
s=importlib.util.spec_from_file_location('q109_angular',p);angular=importlib.util.module_from_spec(s);s.loader.exec_module(angular)
radial=angular.radial;full=angular.full;EP=angular.EP;cv=angular.cv;haar=angular.haar

class FullAxial:
    def __init__(self,angular_rows):
        if not isinstance(angular_rows,angular.FullAngular):raise TypeError('Q108 identical-state angular/mean graph required')
        self.angular=angular_rows;self.rows=angular_rows.rows;self.state=angular_rows.state
    def check(self,n,bundle):
        full.order(n);self.angular.check()
        if self.angular.state is not self.state or self.rows is not self.angular.rows:raise ValueError('One current axial/mean state required')
        self.rows.reset.check(n+1,bundle)
    def old(self,n,bundle):
        self.check(n,bundle);fresh=self.rows.fresh.rows;c=self.rows.current
        olddiff=fresh.bound(n,bundle)
        if olddiff['source_state'] is not self.state or olddiff['order']<n:raise ValueError('Same old full difference program required')
        p=olddiff['pole'];L=self.state.target.Lhi;cp=c.current.profile
        signed=[]
        for j in range(n+3):
            d=angular.native.covariance.defect(self.state,j,target_rate=Q(2))
            if d['order']<j:raise ValueError('Old signed covariance jet reserve required')
            B=radial.rows.embedding(L,d['pole'],p)
            signed.append(factorial(j)*B*cv(d['signed']))
        op=fresh.pressure(n+1,p)
        if op['order']<n+1 or op['pole']!=p:raise ValueError('Exact common-pole full pressure archive required')
        edge=cv(radial.edge(p,self.state.target.width))
        rho=243*radial.rows.embedding(L,0,p)*cv(cp.weighted_rho(n+1))
        pm=haar.pressure.envelope(n+1,signed,cp.Rmin,cp.Rmax,edge,rho,op['KT'],op['KA'])
        delta=haar.envelope(n,signed,cp.Rmin,pm['pressure'])
        h=self.state.target.source.haar
        initial=cv(radial.edge(2*h.pole,self.state.target.width))*cv(h.rows(n)['Ez'])
        coefficient=243*(initial+edge*(cv(delta['axial'])+cv(olddiff['fast_difference']['Ez'])+cv(olddiff['slow_difference']['Ez'])))
        return dict(coefficient=coefficient,initial=initial,signed_haar=delta,signed_pressure=pm,
                    old_differences=olddiff,pole=p,order=n,source_state=self.state,
                    signed_covariance_request=n+2,pressure_request=n+1,
                    old_A1_haar=0,old_A1_velocity_retained=True,epsilon_order='2-2*kappa')
    def bound(self,n,bundle):
        old=self.old(n,bundle);pm=self.rows.fresh.rows.changes.pressure
        if pm.state is not self.state:raise ValueError('Current sourced pressure state required')
        newpressure=pm.bound(n+1,bundle)
        delta=haar.envelope(n,bundle['change'],pm.maps.profile.Rmin,newpressure['pressure'])
        edge=cv(radial.edge(bundle['covariance_pole'],self.state.target.width))
        return dict(pre_fast_haar=old['coefficient']+edge*cv(delta['axial']),old=old,
                    sourced_haar=delta,source_state=self.state,order=n,epsilon_order='3/2-kappa')
    def final(self,n,bundle):
        pre=self.bound(n,bundle);out=self.rows.bound(n,bundle)
        if out['source_state'] is not self.state or out['order']<n:raise ValueError('Same new full differences and A1 source required')
        edge=cv(radial.edge(out['pole'],self.state.target.width))
        coefficient=pre['pre_fast_haar']+edge*(cv(out['pre_fast_A1'])+cv(out['fast_difference']['Ez'])+cv(out['slow_difference']['Ez']))
        return dict(coefficient=coefficient,pre_fast=pre,source_state=self.state,order=n,
                    norm='unweighted-factorial-through-R',epsilon_order='3/2-kappa',
                    full_axial_mean_row=True,radial_cone_closed=False,
                    actual_native_contracts_certified=False,full_residual_certified=False)

def rates(kappa=Q(1,100000)):
    k=Q(kappa)
    if not 0<k<Q(1,4):raise ValueError('Original kappa range required')
    v=dict(old=2-2*k,sourced=Q(3,2)-k,current_A1=2-k,fast=Q(5,2)-4*k,slow=Q(5,2)-3*k)
    return dict(v,full_axial=min(v.values()),target_gap=Q(1,2)-k)

def record():
    return dict(status='current-full-axial-mean-row-with-complete-haar-pressure',
                old_haar_order='2-2*kappa',full_axial_order='3/2-kappa',
                old_A1_haar_zero_but_velocity_retained=True,
                initial_signed_and_sourced_pressure_memory_retained=True,
                radial_cone_or_stronger_target_closed=False,
                actual_native_contracts_certified=False,unconditional_pulse_residual_or_blowup=False)
