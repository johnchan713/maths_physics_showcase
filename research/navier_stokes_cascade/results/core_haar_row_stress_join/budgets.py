"""Complete Haar rows and postslow radial moments; explicit principal-stress gap."""
from fractions import Fraction as Q
from functools import lru_cache
from pathlib import Path
import importlib.util

HERE=Path(__file__).resolve().parent


def load(name,path):
    s=importlib.util.spec_from_file_location(name,path);m=importlib.util.module_from_spec(s);s.loader.exec_module(m);return m


fast=load('haar_stress_fast',HERE.parent/'core_fast_row_gain_join/budgets.py')
scalar,fresh,prefix=fast.scalar,fast.fresh,fast.prefix
P,order,exact=prefix.P,prefix.order,prefix.exact


def haar_rows(n,state,fixed,W,leading_flux,fine_flux):
    """Current full Haar equations; means B, W order1, fine flux order C.

    Leading_flux is the unverified order-one stress contribution. The
    other row terms have order >=B+1; neither c Nt nor graph losses
    survive normalized Haar averaging. All correlations remain.
    """
    order(n)
    if set(state)!=set(prefix.STATE):raise ValueError('All complete mean components required')
    x,u,z=(state[k] for k in prefix.STATE)
    b,V,G,r,rho,kp,km=(fixed[k] for k in ('b','V','G','r','rho','Kp','Kmom'))
    d=Q(n+4);vis=2*d*d+r*d+r*r
    S=2*b*x+x*x+W;H=b*z+G*x+x*z+W;T=2*V*u+u*u+W
    gr=d*x+(d+r)*S+d*H+r*T+vis*x
    pp=km*gr;pm=kp*(gr+rho*pp)
    rt=b*u+V*x+x*u;rz=b*z+G*x+x*z
    at=G*u+V*z+z*u+W;az=2*G*z+z*z+W+pm
    non_t=d*u+(d+2*r)*rt+d*at+vis*u
    non_z=d*z+(d+r)*rz+d*az+vis*z
    return dict(gr=gr,P=pp,pm=pm,
                mean_and_axial_theta=non_t,mean_and_axial_z=non_z,
                fine_theta=non_t+(d+2*r)*fine_flux,
                fine_z=non_z+(d+r)*fine_flux,
                leading_theta=(d+2*r)*leading_flux,leading_z=(d+r)*leading_flux,
                input_order=n+3,pressure_source_order=n+1)


def flux_remainder(real_curl,first_signed,sourced,higher_stress):
    """Original real curl square, both signed crosses/squares, sourced change.

    These are already finite coefficient majorants for the COMPLETE
    differences, not differences of two unrelated norm bounds.
    The separate principal stress identity is not asserted by this sum.
    """
    return real_curl+first_signed+sourced+higher_stress


def rates(kappa=Q(1,100000)):
    r=fast.rates(kappa);k=exact(kappa);B,C,D=r['B'],r['C'],r['D']
    beta=2-2*k;alpha=2-3*k
    terms=dict(real_primary_curl_square=beta,first_signed_covariance=alpha,
               sourced_covariance=C,higher_background_stress=Q(3))
    return dict(B=B,C=C,D=D,principal_flux_required_gain=D-1,
                fine_flux=min(terms.values()),flux_remainder_terms=terms,
                nonflux_Haar_row=min(B+1,2*B+1,Q(2)),
                same_principal_identity_Haar_row=C,
                postslow_scalar=min(B+C,2*C,C+2),radial_Haar_compensator=B+C,
                actual_leading_stress_identity_certified=False,
                full_Haar_from_coarse_or_scalar_bounds_alone=False,
                complete_residual_or_blowup_certified=False)


class HaarStressJoin:
    """One-way current-state and radial-moment join; no stress identity flag."""
    def __init__(self,means,kappa=Q(1,100000)):
        self.means=means;self.scalars=scalar.ScalarGainBudget(means)
        self.orders=rates(kappa);self.embedding=self.scalars.embedding

    @lru_cache(None)
    def postslow_scalars(self,n):
        order(n);N=n+2;s=self.scalars;m=self.means
        # Old postfast means still have coarse B; SAME actual slow step
        # has the parent's sharper C bound. Exact linear scalar rows cancel.
        old=s.embed(m.postfast(N));inc=s.embed(s.slow_increment(N))
        fixed=m.data('fixed',N)
        return scalar.slow_reset(n,old,inc,fixed['b'],fixed['G'],s.moment(N))

    @lru_cache(None)
    def radial_Haar_compensator(self,n):
        order(n);a=self.postslow_scalars(n)
        return fresh.monomial_upper(self.means.data('fixed',n)['rho']*a['P'])

    def whole_axial_moment(self,n,c_rho):
        a=self.postslow_scalars(n)
        return a['Jz']+c_rho*a['P']

    def current_row_envelope(self,n,*,principal_defect,fine_flux):
        """Two-rate envelope; a supplied norm never certifies exact matching."""
        order(n);N=n+3;m=self.means;E=self.embedding
        state={k:E*v for k,v in m.final(N).items()}
        W=m.data('wave',N)['U'];W=E*W*W
        out=haar_rows(n,state,m.data('fixed',N),W,E*principal_defect,E*fine_flux)
        return dict(out,leading_order=Q(1),fine_order=self.orders['C'],
                    leading_stress_is_unverified=True)

    def certify_actual_principal_matching(self,*args,**kwargs):
        raise NotImplementedError('Prove the principal covariance/stress identity for the SAME actual primary/root family; coarse norms or frame closeness do not certify it')
