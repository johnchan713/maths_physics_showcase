"""Current radial split and a deferred exact A0 request, not tail certification."""
from pathlib import Path
from fractions import Fraction as Q
import importlib.util
HERE=Path(__file__).resolve().parent

def load(name,path):
    s=importlib.util.spec_from_file_location(name,HERE.parent/path);m=importlib.util.module_from_spec(s);s.loader.exec_module(m);return m
axial=load('q110_axial','core_current_head_q109_full_axial_join/q109.py')
tail=load('q110_tail','core_final_mean_row_closure/budgets.py')
radial=axial.radial;full=axial.full;cv=axial.cv;EP=axial.EP

def request(n,*,kappa=Q(1,100000),target=None):
    full.order(n);k=Q(kappa)
    if not 0<k<Q(1,4):raise ValueError('Original kappa range required')
    B=1-k;C=Q(3,2)-k;t=C if target is None else Q(target)
    m=max(1,(t-B)//k+1);J=n+2*m+3
    return dict(output_order=n,integration_order=m,raw_source_order=J,
                final_radial_archive=J,mean_fields=J+2,rows=J+9,covariance=J+11,
                mean_factory=J+14,target_jets=J+28,profile_jets=J+41,cutoff_shape=J+11,
                source_rate=B,target_rate=t,tail_rate=B+m*k,strict_margin=B+m*k-t,
                compensator_rate=Q(5,2)-k,compensator_sufficient=Q(5,2)-k>=t,
                full_graph_common_pole_required=True,actual_tail_bound_certified=False)

class RadialAudit:
    def __init__(self,axial_rows):
        if not isinstance(axial_rows,axial.FullAxial):raise TypeError('Q109 identical-state mean graph required')
        self.axial=axial_rows;self.radial=axial_rows.rows.radial;self.state=axial_rows.state
    def check(self):
        self.axial.angular.check()
        if self.axial.state is not self.state or self.radial is not self.axial.rows.radial:raise ValueError('Same completed radial/mean state required')
    def coarse(self,n,bundle):
        self.check();out=self.radial.bound(n,bundle)
        if out['source_state'] is not self.state or out['norm_stage']!='final_sourced_postslow':raise ValueError('Final six-mean radial stage required')
        return out
    def source(self,J,bundle):
        """Q104 CURRENT FINAL r0 through J, all six means and TOTAL covariance."""
        full.order(J);self.check()
        if J>self.radial.maximum:raise ValueError('Rebuild the same-state radial archive through raw source order')
        self.radial.postslow.check(J,bundle);N=J+2
        mean=self.radial.fields(N,bundle);p=mean['pole'];c=self.radial.current;cp=c.current.profile
        W=radial.rows.embedding(c.current.target.Lhi,bundle['covariance_pole'],p)*radial.total_covariance(N,bundle)
        b=c.initial.budget;data=b._data(N,b.pole)
        fixed={k:243*cv(data[k]) for k in ('b','V','G','c','graph')}
        E=1+cv(radial.edge(p,self.state.target.width))
        out=radial.radial_envelope(J,mean['fields'],W,cp.Rmin,fixed,E)
        return dict(coefficient=out['reduced_source'],order=J,pole=p,source_state=self.state,
                    mean=mean,total_wave_covariance=W,epsilon_order='1-kappa',norm='weighted-factorial-through-R',
                    pressure_memory_removed_only_under_A0=True,actual_native_contracts_certified=False)
    def recipe(self,n,bundle,*,geometry,kappa=Q(1,100000),target=None):
        """Requires full high-jet metadata before issuing even a deferred recipe."""
        self.check();r=request(n,kappa=kappa,target=target);J=r['raw_source_order']
        if self.radial.maximum<J:raise ValueError('Low-output radial archive cannot cover the A0 raw source')
        self.radial.postslow.check(J,bundle)
        if any(key not in bundle or len(bundle[key])<=r['covariance'] for key in ('old','new','change')):
            raise ValueError('Complete old/new/changed wave jets through the full covariance request required')
        p=self.radial.common_pole(bundle)
        if Q(geometry['width'])!=self.state.target.width:raise ValueError('Original frozen width required')
        out=tail.oscillatory_recipe(n,r['integration_order'],e=0,
            source='Qcirc CURRENT FINAL reduced r0 from RadialAudit.source',
            source_rate=r['source_rate'],target_rate=r['target_rate'],kappa=kappa,
            geometry=geometry,pole=p)
        out.update(request=r,source_method='RadialAudit.source(raw_source_order, rebuilt_bundle)',
                   source_state=self.state,whole_graph_pole_proof_still_required=True,
                   original_source_and_output_chart_transfers_required=True,chart_factors_materialized=False,
                   materialized_actual_input_bound=False,radial_cone_closed=False)
        return out

def record():
    r=request(0)
    return dict(status='current-radial-split-and-full-graph-A0-request-audit',
                default_integrations=r['integration_order'],default_raw_source_order=r['raw_source_order'],
                covariance_request=r['covariance'],profile_request=r['profile_jets'],
                scalar_sufficient_for_current_mean_target=True,
                low_output_archive_not_a_high_jet_proof=True,
                actual_tail_or_native_inputs_certified=False,full_cone_pulse_residual_or_blowup=False)
