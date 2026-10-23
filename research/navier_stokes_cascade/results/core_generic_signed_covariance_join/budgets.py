"""Complete generic signed covariance, new pressure domains and target gain."""
from fractions import Fraction as Q
from functools import lru_cache
from pathlib import Path
import importlib.util

HERE=Path(__file__).resolve().parent


def load(name,path):
    s=importlib.util.spec_from_file_location(name,path)
    m=importlib.util.module_from_spec(s);s.loader.exec_module(m);return m


preserved=load('generic_signed_preserved',HERE.parent/'core_preserved_haar_scalar_join/budgets.py')
final,scalar,prefix=preserved.final,preserved.scalar,preserved.prefix
edge=load('generic_signed_edge',HERE.parent/'core_signed_target_edge_budget/budgets.py')
fresh=final.haar.fresh
P,order,exact=prefix.P,prefix.order,prefix.exact


def convert(p):return P(p.c)


def rates(kappa=Q(1,100000)):
    k=exact(kappa);a=preserved.rates(k)
    B,C,D=a['B'],a['C'],a['D'];C0=2*B
    terms=dict(old_error_cross=2*B,complete_curl_cross=Q(1,2)+D,signed_square=2*B)
    return dict(B=B,C=C,D=D,C0=C0,principal_signed=B,curl_signed=D,
        principal_covariance=C,complete_covariance=C,error_terms=terms,
        covariance_error=min(terms.values()),Haar_pressure_change=C,
        pointwise_pressure_change=D,scalar_change=C,
        axial_covariance_pressure_response=C+1,post_signed_target=C0,
        post_signed_Haar_tangential=C0,post_signed_radial_Haar=C,
        complete_regular_improvement_or_actual_input_validity=False)


def covariance_budget(T,E,U,L,R):
    """Through-order harmonic-l1 majorants, same weighted coefficient chart.

    T bounds ORIGINAL principal W0, E bounds U-W0 from its actual
    summands, U bounds current complete wave. Neither E nor T is inferred
    from the size of U or a difference of unrelated norm bounds.
    """
    lam=2*T*L
    err=2*E*L+2*U*R+(L+R)*(L+R)
    return dict(principal=lam,error=err,complete=lam+err)


def pressure_target_budget(n,D,Dnext,Dnext2,error,*,r,rho,Kp,Kmom,primitive,bump):
    """Full six-entry covariance -> SAME reconstructed Haar pressure.

    All inputs are through-order norms, Dnext/Dnext2 at n+1/n+2;
    fixed/operator inputs dominate all calls through n+1. Output separates
    the C-order pressure from its (C+1)-order axial target response.
    No graph loss survives Haar; the full pointwise pressure is only D.
    """
    order(n)
    gr=2*(n+1)*Dnext+2*r*D
    Pdef=Kmom*gr
    # Pressure through n+1 requires the covariance through n+2. This
    # conservative envelope also covers gr at n and all fixed products.
    grnext=2*(n+2)*Dnext2+2*r*Dnext
    pm=Kp*(grnext+rho*Kmom*grnext)
    projection=1+bump*Kmom
    axial_theta=primitive*projection*(n+1)*Dnext
    axial_z=primitive*projection*(n+1)*(Dnext+pm)
    target_theta=error+axial_theta;target_z=error+axial_z
    moments=scalar.wave_change(n,D,Dnext,Kmom)
    return dict(gr=gr,P=Pdef,pm_next=pm,axial_theta=axial_theta,axial_z=axial_z,
        target_theta=target_theta,target_z=target_z,
        Jtheta=moments['Jtheta'],Jz=moments['Jz'],moment_P=moments['P'])


class PairedTargetBudget:
    """Keep e=2 and e=1 targets separate, then bound their vector by a sum."""
    def __init__(self,theta,axial):self.theta=theta;self.axial=axial

    def V(self,n):return edge.P(self.theta.V(n).c)+edge.P(self.axial.V(n).c)

    def edge(self,N):
        a,b=self.theta.edge(N),self.axial.edge(N)
        return dict(d=max(a['d'],b['d']),theta=a,axial=b)


class RootSignedEnvelope:
    """Use the ORIGINAL RootBudget.signed and identical full curl operator.

    target combines the TWO original e=2/e=1 TargetBudget recipes with
    PairedTargetBudget, including their distinct moments. transfer(n) is a retained
    common-chart/label norm conversion, including all weight/chart jets.
    These lower inequalities are not certified by evaluating polynomials.
    """
    def __init__(self,root,target,transfer):
        if not isinstance(target,PairedTargetBudget):raise ValueError('Both e=2/e=1 current target recipes are required')
        self.root=root;self.target=target;self.transfer=transfer
        self.RootP=type(root.data('wave',0))

    @lru_cache(None)
    def principal_jet(self,n):
        order(n)
        return convert(self.root.signed(n,lambda j:self.RootP(self.target.V(j).c)))

    @lru_cache(None)
    def curl_jet(self,n):
        order(n)
        return sum((convert(self.root.data('curl',j))*(self.principal_jet(n-j)+(n-j+1)*self.principal_jet(n-j+1)) for j in range(n+1)),P())

    def through(self,n):
        order(n);w=self.transfer(n)
        return dict(L=w*sum((self.principal_jet(j) for j in range(n+1)),P()),
                    R=w*sum((self.curl_jet(j) for j in range(n+1)),P()))


def pressure_domains(n,p,target):
    """Finite proof-domain expansion for the NEXT signed pressure call.

    At target output n, full Haar pressure through n+1 uses covariance
    through n+2, signed curl through n+2 and target input through n+3.
    q is a sufficient COMMON wave pole for those finite inputs. Original
    physical cutoffs/roots/width/band do not change when estimates shift.
    """
    order(n);p=exact(p)
    if p<0:raise ValueError('Nonnegative original pole required')
    N=n+2;te=target.edge(n+3);d=te['d']
    groups=dict(current_complete_wave=3*p,original_primary_curl=6*N+9,
                new_signed_curl=6*N+d+6)
    q=max(groups.values());required=2*q
    old_domains=[6*p,12*p]
    return dict(output_order=n,covariance_input_order=N,current_target_input_order=n+3,
        target_quotient_allowance=str(d),wave_pole_groups={k:str(v) for k,v in groups.items()},
        common_wave_pole=str(q),covariance_and_pressure_source_pole=str(required),
        retained_old_pressure_domains=[str(v) for v in old_domains],
        old_domain_proof_covers_new_call=required in old_domains,
        replacement_estimate_parameter=str(q/6),
        new_operator_domains=[str(q),str(2*q)],
        original_raw_pressure_profiles_retained=True,materialized_actual_input_bound=False)


def signed_pressure_shift(n,p,target,geometry):
    d=pressure_domains(n,p,target)
    recipe=fresh.shift_recipe(Q(d['replacement_estimate_parameter']),Q(geometry['L']),
        Q(geometry['width']),Q(geometry['gap']))
    primitive=edge.primitive_recipe(Q(d['covariance_and_pressure_source_pole']),Q(geometry['L']),Q(geometry['width']))
    return dict(domains=d,shift_recipe=recipe,oriented_target_primitive=dict(recipe=primitive,
        upper_rule='K_Ie=sum(branch coefficient * 3^ceil(branch exponent)); original TargetBudget.F/K derivative and radial factors retained'),
        expression='reevaluate every nonnegative input/response polynomial at S+K_new+K_Ie; no S>=K gate',
        fields_or_accepted_domain_changed=False)


class CurrentSignedRecipe:
    """Current final rows and complete moments, not a new fine-row callback.

    This deferred recipe references the SAME latest actual mean prefix.
    It does not materialize the giant original A0/A1 source jets, certify
    original weighted inputs, or implement the subsequent fresh means.
    """
    def __init__(self,final_mean,haar_join):
        if final_mean.means is not haar_join.means:raise ValueError('One and the same current mean prefix required')
        self.final=final_mean;self.haar=haar_join

    def target_inputs(self,n):
        order(n)
        return dict(stage='after-latest-fast-and-slow-before-next-generic-signed',
            current_row_recipe=self.final.tangential_recipe(n),
            scalar_recipe=self.haar.certificate(n+1),
            moment_orders=dict(theta=str(self.haar.orders['C']+2),axial=str(self.haar.orders['C']+2)),
            theta_moment='(n+1) PreservedHaarScalarJoin.postslow(n+1).Jtheta',
            axial_moment='PreservedHaarScalarJoin.axial_residual(n,c_rho_next)',
            target_routing='TWO TargetBudget recipes e=2/e=1 on current weighted rows and distinct whole moments; PairedTargetBudget sums their V bounds and takes the larger edge allowance',
            original_root_routing='RootSignedEnvelope on SAME positive original roots; no updated-root inverse',
            actual_weighted_input_or_support_validity=False)

    def certify_actual_inputs(self):
        raise NotImplementedError('Current weighted target/root/common-chart inputs and next pressure source domains require actual input validation')
