"""Current post-signed source, full homogeneous pressure and physical tails."""
from fractions import Fraction as Q
from math import factorial
from pathlib import Path
import importlib.util

HERE=Path(__file__).resolve().parent

def load(name,path):
    s=importlib.util.spec_from_file_location(name,path)
    m=importlib.util.module_from_spec(s);s.loader.exec_module(m);return m

mean=load('post_residual_mean',HERE.parent/'core_post_signed_mean_closure/budgets.py')
mult=load('post_residual_multiharmonic',HERE.parent/'core_multiharmonic_mean_residual_budget/budgets.py')
physical=load('post_residual_physical',HERE.parent/'core_physical_tail_join/bounds.py')
signed,prefix=mean.signed,mult.prefix
P,exact,order=mult.P,mult.exact,mult.order


def rates(kappa=Q(1,100000)):
    r=mean.rates(kappa);k=exact(kappa);B,C,D,C0=(r[x] for x in ('B','C','D','C0'))
    groups=dict(inherited_current_source=D,signed_complete_linear=B+Q(1,2)-k,
        both_old_signed_wave_crosses=Q(1,2)+B-k,signed_square=2*B-k,
        old_mean_signed_wave=B+B-Q(1,2),new_fast_feedback=D,new_slow_feedback=C)
    return dict(B=B,C=C,D=D,C0=C0,groups=groups,nonzero_source=min(groups.values()),
        regular_mean=C0,whole_regular=min(C0,*groups.values()),
        seed_relation=C0-Q(1,2),postslow_scalar=C+1,
        absolute_decay_or_actual_inputs_certified=False)


def homogeneous_pressure(T,*,normal,K,c,normal_time,inverse_normal_square):
    """Full f=0 principal pressure, with 1/|omega| factored as epsilon^1/2.

    Pi = i/(omega |n|^2) (n.K-c Nt n).t. Inputs dominate through n,
    normal_time through n needs the original normal through n+1.
    Scalar signed coefficients are auxiliary-time independent; otherwise
    principal homogeneous cancellation must include c(Nt alpha)t.
    """
    return inverse_normal_square*(normal*K+c*normal_time)*T


def regular_source_envelope(n,old_source,data):
    """Full signed-wave difference plus the NEW two actual mean increments.

    All entries through n+3. T/R are UNCUT RootSignedEnvelope L/R;
    Pi is derived above, never a fine-pressure callback. newU below is
    the actual complete cut velocity, not a residual tail. No old source
    cancellation is available at this homogeneous update.
    """
    order(n)
    keys=('b','G','graph','c','r')+prefix.LINEAR+('T','R','oldU','oldM','fastM','slowM','normal_time','inverse_normal_square')
    d=prefix.require(data,keys)
    pi=homogeneous_pressure(d['T'],normal=d['normal'],K=d['K'],c=d['c'],
        normal_time=d['normal_time'],inverse_normal_square=d['inverse_normal_square'])
    v=d['cutoff']*(d['T']+d['R'])
    supplied={k:d[k] for k in ('b','G','graph','c','r')+prefix.LINEAR+('T','R','oldU','oldM')}
    supplied.update(Pi=pi,newU=v,dM=d['fastM']+d['slowM'])
    out=mult.increment_recipe(n,supplied,(-2,-1,1,2))
    out.update(old_regular_source=old_source,signed_pressure=pi,complete_cut_signed=v,
        F=old_source+out['F'],raw_source_factorial=factorial(n),
        pressure_source_normal_term='zero because f=0; moving-normal and K pressure retained')
    return out


def physical_request(m,N,h,kappa=Q(1,100000)):
    """Trace current compact sources at PHYSICAL, not normalized, margins."""
    a0=physical.compact_orders(0,m,N,h,kappa);a1=physical.compact_orders(1,m,N,h,kappa)
    trace=mean.dependency_trace(m,a1['integration_order'],a0['integration_order'])
    # The nonzero signed linear expression reads uncut curl through m+3,
    # hence its principal/target through m+4. The compact graph dominates.
    trace['maxima']['signed_linear_target_input']=m+4
    trace['maxima']['next_signed_target_input']=max(trace['maxima']['next_signed_target_input'],m+4)
    return dict(A0=a0,A1=a1,trace=trace)


def whole_domain(trace,previous_wave_pole,original_pole):
    p,q0=map(exact,(original_pole,previous_wave_pole))
    if min(p,q0)<0:raise ValueError('Nonnegative finite poles required')
    maxima=trace['maxima'];N=max(maxima['new_covariance'],maxima['next_signed_target_input']-1)
    groups=dict(retained_wave=q0,old_complete_mean_domain=6*p,
                original_curl=6*N+9,new_signed_curl=6*N+12*p+6)
    q=max(groups.values())
    return dict(q=str(q),groups={k:str(v) for k,v in groups.items()},
        maximum_covariance_input=N,maximum_current_target_input=maxima['next_signed_target_input'],
        post_signed_rows=str(2*q),new_fast_fields=str(2*q),final_radial_source=str(4*q),
        full_nonzero_regular_source=str(3*q),old_nonzero_source=str(9*p),
        required_pressure_domains=[str(2*q),str(4*q)],oriented_target_primitive_pole=str(2*q),
        original_p_covers_ALL_nested_old_jets_required=True,actual_input_bounds_evaluated=False)


class CurrentResidualJoin:
    """One finite SAME-prefix join; all giant arithmetic remains deferred."""
    def __init__(self,mean_recipes,h):
        self.mean=mean_recipes;self.h=physical.scales(h,mean_recipes.kappa)['h']
        self.orders=rates(mean_recipes.kappa)

    def domains(self,m,N):
        request=physical_request(m,N,self.h,self.mean.kappa)
        domain=whole_domain(request['trace'],self.mean.wave_pole,self.mean.original_pole)
        geometry=self.mean.geometry;q=Q(domain['q'])
        operators=dict(pressure_and_fast=dict(method='fresh.shift_recipe',parameter=str(q/3),geometry=geometry,materialized=False),
            oriented_target=dict(method='edge.primitive_recipe',pole=str(2*q),geometry=geometry,materialized=False,
                upper_rule='K_Ie=sum(branch coefficient * 3^ceil(branch exponent)); retain TargetBudget F/K radial and derivative factors'),
            archived_estimate_shifts='K_old: sum of all finite shifts already required by the SAME archived lower recipes',
            common_argument='S+K_old+K_mean+K_Ie; neither archived nor oriented constants are inferred from a higher-pole pressure constant',
            new_band_or_actual_geometry_changed=False)
        return dict(request=request,domain=domain,operators=operators)

    def regular_recipe(self,n):
        order(n)
        return dict(method='regular_source_envelope',order=n,input_order=n+3,
            old_source='core_fast_row_gain_join: mult.increment_recipe with SAME derived fast D and scalar slow C, not an outgoing fine-source input',
            signed_principal_curl='generic RootSignedEnvelope.through(n+3) on SAME original root and paired current TargetBudget',
            full_signed_pressure='homogeneous_pressure from SAME t,n,K,c,Nt n and inverse |n|^2 through n+3; normal jet n+4',
            old_means='previous FinalMeanRowClosure full current means',
            new_means='PostSignedMeanRecipes.fast_recipe(n+3) and fresh scalar_recipe(n+3), including actual A1 velocity',
            retained_angular_inputs=[-2,-1,1,2],inherited_nonzero_outputs=list(range(-4,0))+list(range(1,5)),
            signed_cross_outputs=[-3,-2,-1,0,1,2,3],source_order=str(self.orders['D']),
            new_tail='homogeneous temporal commutator only; f=0, so no new (1-psi) source term',
            no_actual_input_or_materialized_source_certificate=True)

    def compact_recipes(self,m,N):
        d=self.domains(m,N);q=Q(d['domain']['q']);out={}
        for e,key in ((0,'A0'),(1,'A1')):
            a=d['request'][key]
            source='Qcirc NEW FINAL full gr' if e==0 else 'Qcirc post-generic-signed Ez'
            r=mean.final.oscillatory_recipe(m,a['integration_order'],e=e,source=source,
                source_rate=a['source_rate'],target_rate=a['normalized_target'],kappa=self.mean.kappa,
                geometry=self.mean.geometry,pole=(4 if e==0 else 2)*q)
            out[key]=dict(orders={k:str(v) for k,v in a.items()},raw_kernel=r,
                exact_tail='-A0 Qcirc gr_final' if e==0 else '+A1 Qcirc EpostSigned_z',
                source_method='prefix.full_mean on SAME old final + new fast + new slow' if e==0 else 'PostSignedMeanRecipes.post_signed_rows',
                raw_source_order=a['raw_source_order'],whole_older_DAG_preserved=True,
                polynomial_absorption='sum_j c_j (16+K_old+K_mean+K_Ie)^j (2j)! / delta^(2j)',
                evaluated_actual_high_jet_bound=False)
        return out

    def certificate(self,m,N):
        d=self.domains(m,N)
        return dict(regular=self.regular_recipe(m),domains=d,compact=self.compact_recipes(m,N),
            mean_regular_order=str(self.orders['C0']),nonzero_regular_order=str(self.orders['D']),
            full_identity='Rphysical=Rregular_current+(-A0 Qcirc gr_final,0,+A1 Qcirc EpostSigned_z)+TGaussian_old+Tpsi_signed+Tbase_exterior',
            Gaussian='original additive tail ledger plus full UNCUT signed complete curl times psi_t; never inverted',
            base_exterior='retained exact flat remainder after bare stress divergence assigned once',
            physical=physical.regular_transfer(m,self.h,self.mean.kappa),
            all_retained_actual_weighted_support_background_inputs_required=True,
            complete_regular_grammar_mean_tail_join_is_conditional=True,
            complete_actual_residual_improvement_absolute_decay_or_blowup_certified=False)

    def certify_actual_inputs(self):
        raise NotImplementedError('Actual weighted high jets/support/background/base tail and later higher regular orders/controlled diagonal remain unproved')
