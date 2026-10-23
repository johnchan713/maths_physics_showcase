"""Finite complete-scalar refinement; original weighted inputs stay conditional."""
from fractions import Fraction as Q
from functools import lru_cache
from pathlib import Path
import importlib.util

HERE=Path(__file__).resolve().parent


def load(name,path):
    s=importlib.util.spec_from_file_location(name,path)
    m=importlib.util.module_from_spec(s);s.loader.exec_module(m);return m


final=load('preserved_haar_final',HERE.parent/'core_final_mean_row_closure/budgets.py')
haar,scalar,prefix=final.haar,final.haar.scalar,final.prefix
seed_module=load('preserved_haar_initial',HERE.parent/'core_initial_scalar_seed_budget/proof_checks.py')
P,order,exact=prefix.P,prefix.order,prefix.exact
SCALARS=('P','Jtheta','Jz')


def convert(p):return P(p.c)


def rates(kappa=Q(1,100000)):
    k=exact(kappa)
    if not 0<k<Q(1,4):raise ValueError('Retained 0<kappa<1/4 required')
    B=1-k;C=B+Q(1,2);D=C-k;alpha=2-3*k
    return dict(B=B,C=C,D=D,first_signed=alpha,initial_scalars=Q(2),
                first_postfast=min(Q(2),alpha,2*B),first_slow=alpha,
                current_postfast=C,current_slow=C,old_Haar_tangent=Q(1),
                old_Haar_radial=Q(2),postslow=min(C+1,2*C,C+2),
                coarse_postslow=B+C,gain_over_coarse=k,
                whole_axial_residual=C+2,radial_compensator=C+1,
                full_regular_order=D,absolute_physical_decay_certified=False)


def response(n,a,an,k,kn):
    """SAME fixed five-map, on fresh complete scalar majorants."""
    order(n)
    return dict(u=k['H']*(k['Kp']*a['P']+k['Kz']*a['Jz']),
                z=k['H']*k['Kg']*a['Jtheta'],
                x=(n+1)*kn['H']*kn['Kpsi']*an['Jtheta'])


def signed_covariance(primary,signed):
    """Full l1 tensor difference; no real-quadrature mixed-term deletion."""
    return 2*primary*signed+signed*signed


def haar_slow_reset(n,old_haar,inc,b,M):
    """Old tangential order 1/radial 2, new C/C+1, reserved G=0.

    This is only for an auxiliary-independent SLOW increment. It is not
    legal for a zero-Haar fast increment or for a total pressure field.
    """
    return scalar.slow_reset(n,old_haar,inc,b,b*0,M)


class PreservedHaarScalarJoin:
    """Derive the missing Haar envelope from LOWER original input recipes.

    seed is the ORIGINAL canonical SeedBudget; primary(n) bounds the
    original complete CUT primary at order 1/2. A bound for an updated
    total wave is not a bound for that original summand. first_signed(n) bounds the
    complete cut original-root real signed velocity at alpha-1/2, not an
    incoming fine scalar or mean callback. transfer(n) is a bound for the
    original raw/slow-chart to current normalized weighted closed-hull
    norm conversion, including ALL chart/weight derivatives. It must be
    independently justified on the same actual hull. These retained
    lower input bounds are NOT verified by constructing this object.

    Supplied seed/transfer data cannot be silently inferred from the
    generic PrefixBudget.initial callback. Diagnostics use manufactured
    norm polynomials only. No physical fields or high jets are evaluated.
    """
    def __init__(self,means,seed,primary,first_signed,transfer,kappa=Q(1,100000)):
        self.means=means;self.before=means.before;self.seed=seed
        self.primary=primary;self.first_signed=first_signed;self.transfer=transfer
        self.orders=rates(kappa);self.current=scalar.ScalarGainBudget(means)
        self.embedding=self.current.embedding

    def seed_fields(self,n):
        f=self.seed.fields(n);w=self.transfer(n)
        return {key:self.embedding*w*convert(f[other])
                for key,other in (('x','X'),('u','U'),('z','Z'))}

    @lru_cache(None)
    def initial_defects(self,n):
        order(n);w=self.transfer(n);a=self.seed.output(n)
        return {key:self.embedding*w*convert(a[key]) for key in SCALARS}

    @lru_cache(None)
    def first_wave_change(self,n):
        order(n)
        def cov(j):
            # Both are original lower coefficient recipes. A bound for
            # the updated wave alone cannot bound its original summand.
            U=self.primary(j);V=self.first_signed(j)
            return self.embedding*signed_covariance(U,V)
        return scalar.wave_change(n,cov(n),cov(n+1),self.current.moment(n+1))

    @lru_cache(None)
    def first_fast_change(self,n):
        order(n);N=n+1
        return scalar.fast_change(n,self.current.embed(self.before.data('initial',N)),
            self.current.embed(self.before.fast(N)),self.current.moment(N))

    @lru_cache(None)
    def first_postfast(self,n):
        order(n);parts=(self.initial_defects(n),self.first_wave_change(n),self.first_fast_change(n))
        return {key:sum((a[key] for a in parts),P()) for key in SCALARS}

    @lru_cache(None)
    def first_slow(self,n):
        order(n)
        return response(n,self.first_postfast(n),self.first_postfast(n+1),
            self.before.data('slow',n),self.before.data('slow',n+1))

    @lru_cache(None)
    def old_Haar(self,n):
        """BEFORE the latest slow step: seed Haar + first signed slow.

        All intervening complete fast increments have zero Haar in every
        component, including A1; signed/sourced waves change no mean
        velocities. Both displayed slow orders exceed one.
        """
        order(n);a,b=self.seed_fields(n),self.first_slow(n)
        return {key:a[key]+b[key] for key in prefix.STATE}

    @lru_cache(None)
    def postslow(self,n):
        order(n);N=n+2;fixed=self.means.data('fixed',N)
        return haar_slow_reset(n,self.old_Haar(N),
            self.current.embed(self.current.slow_increment(N)),fixed['b'],self.current.moment(N))

    def whole_axial_moment(self,n,c_rho):
        """c_rho bounds the WHOLE through-n moving product multiplier."""
        a=self.postslow(n);return a['Jz']+c_rho*a['P']

    def axial_residual(self,n,c_rho_next):
        order(n)
        return (n+1)*self.whole_axial_moment(n+1,c_rho_next)

    def radial_compensator(self,n):
        return self.means.data('fixed',n)['rho']*self.postslow(n)['P']

    def certificate(self,n):
        order(n)
        return dict(order=n,postslow_scalar_order=str(self.orders['postslow']),
            axial_residual_moment_order=str(self.orders['whole_axial_residual']),
            radial_Haar_compensator_order=str(self.orders['radial_compensator']),
            initial_recipe='original SeedBudget.fields/output with same-hull weighted norm transfer',
            first_slow_recipe='response(initial S2 + complete first signed covariance change + full-old-mean fast change)',
            old_Haar_recipe='seed_fields + first_slow; never the full current postfast norm relabelled',
            sharp_current_slow_recipe='ScalarGainBudget.slow_increment at n+2',
            missing_actual_input_inequalities=[
                '||initial Haar fields||_n,weighted <= epsilon^(1,2) transfer_n * SeedBudget.fields(n)',
                '||initial current defects||_n,weighted <= epsilon^2 transfer_n * SeedBudget.output(n)',
                '||original complete cut primary||_n,weighted <= epsilon^(1/2) primary(n)',
                '||complete cut first original-root signed velocity||_n,weighted <= epsilon^(alpha-1/2) first_signed(n)',
                'same actual G=0 reserved potential hull, whole-line compact fast operators and complete zero-Haar identities'],
            deferred_lower_inputs_independently_verified=False,
            actual_current_fine_norms_certified=False,complete_residual_improvement_or_blowup=False)

    def certify_actual_inputs(self):
        raise NotImplementedError('Original weighted seed/complete signed coefficient/support inequalities require actual lower input validation')
