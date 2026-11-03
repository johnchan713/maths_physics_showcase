"""Current frozen-phase/rate/linear leaves from the original full-background recipe."""
from pathlib import Path
from fractions import Fraction as Q
from functools import lru_cache
from math import factorial
import importlib.util
HERE=Path(__file__).resolve().parent
def load(name,path):
    s=importlib.util.spec_from_file_location(name,HERE.parent/path);m=importlib.util.module_from_spec(s);s.loader.exec_module(m);return m
bound=load('q118_response','core_current_head_q117_current_covariant_response/q117.py')
background=load('q118_background','core_background_finite_head_budget/budgets.py')
charts=load('q118_charts','core_current_head_q111_closed_chart_majorant/q111.py')
native=bound.native;EP,cv=bound.EP,bound.cv
ExpSum=type(EP(1).c[0])

def request(j):
    native.full.order(j)
    return dict(output=j,source_ceiling=j+1,phase_common_ceiling=j+3,
                background_ceiling=j+4,Borel_Euler_cutoff_ceiling=j+5,
                source_factory_reserve=j+3,actual_inputs_certified=False)

class BackgroundPhase:
    def __init__(self,postslow,budget,*,maximum,path,Lmin,U,lambda_max,gmin):
        if not isinstance(budget,background.BackgroundBudget):raise TypeError('Original full finite-head background recipe required')
        native.full.order(maximum)
        self.postslow=postslow
        self.state=native.inputs.source.waves.linear.meanwave.means.FinalMeans(postslow).state
        self.budget=budget;self._maximum=maximum
        chart=self.state.primary.chart
        if (budget.Rmin,budget.Rmax)!=(chart.Rmin,chart.Rmax):
            raise ValueError('Same original current closed R hull required')
        self.hull=charts.hull(chart.Rmin,chart.Rmax,postslow.current.fast.rows.pressure.d)
        self.recipe=budget.phase_inputs(path,Lmin,U,lambda_max,gmin)
        base=16*(1+self.hull['Rmax'])/min(Q(1),self.hull['Umin'])
        # j! times the Q111 forward transfer, dominated at one FIXED ceiling.
        self.forward_height=Q(486)+2*maximum*maximum+maximum*base
        self.forward=EP(ExpSum([(1,self.forward_height)]))
    @property
    def maximum(self):return self._maximum
    def check(self,n):
        native.full.order(n)
        if n>self.maximum:raise ValueError('Rebuild the same background phase at the full fixed ceiling')
        if native.inputs.source.waves.linear.meanwave.means.FinalMeans(self.postslow).state is not self.state:
            raise ValueError('Current first-real state changed')
        c=self.state.primary.chart
        if (self.budget.Rmin,self.budget.Rmax)!=(c.Rmin,c.Rmax):raise ValueError('Original background hull changed')
    @lru_cache(None)
    def common(self,kind):
        # The old PhaseBudget has RAW R jets. First form its factorial R norm,
        # then pay R->U and raw conversion, before Q93 later returns to R.
        self.check(self.maximum)
        raw=getattr(self.recipe,kind)
        norm=sum((Q(5**j,factorial(j))*cv(raw(j)) for j in range(self.maximum+1)),EP())
        return self.forward*norm
    def jet(self,kind,n):self.check(n);return self.common(kind)
    def normal(self,n):return self.jet('normal',n)
    def matrix(self,n):return self.jet('matrix',n)
    def moving_normal(self,n):return self.jet('moving_normal',n)
    def profile_rate(self,n):return self.jet('profile_rate',n)
    def r_through(self,kind,n):
        self.check(n)
        norm=sum((Q(5**j,factorial(j))*self.jet(kind,j) for j in range(n+1)),EP())
        factor=2*cv(self.state.primary.chart.transfer(n))*16*188**n*factorial(n)
        return factor*native.inputs.source.waves.rescale(norm)
    def linear_table(self,R,psi):
        self.check(R+2)
        if len(psi)!=R+1:raise ValueError('Original cutoff/psi bounds at every source order required')
        return [dict(Eik=self.r_through('eikonal',j+2),bR=self.r_through('bR',j+2),
                     Bz=self.r_through('Bz',j+2),psi=cv(psi[j])) for j in range(R+1)]
    def waves(self,*,regular_wave,inverse_normal_square):
        self.check(0)
        return native.inputs.source.waves.CurrentWaves(self.postslow,regular_wave=regular_wave,
                                                     phase=self,inverse_normal_square=inverse_normal_square)

class BackgroundResponse(bound.CovariantResponse):
    def __init__(self,inputs):
        super().__init__(inputs)
        if not isinstance(inputs.source.waves.phase,BackgroundPhase):
            raise TypeError('Current same-hull full-background phase binding required')
        self.background_phase=inputs.source.waves.phase
        if self.background_phase.state is not self.state:raise ValueError('Identical background/response state required')
    def build(self,j,linear_table,log_table,rate_table,cutoff,*,psi_table,**kw):
        native.full.order(j);R=j+1
        if linear_table is not None or rate_table is not None:
            raise ValueError('Linear phase and reference-rate tables are derived, not independent inputs')
        if self.inputs.source.waves.phase is not self.background_phase or self.background_phase.state is not self.state:
            raise ValueError('Original bound phase object and source state required')
        self.background_phase.check(request(j)['phase_common_ceiling'])
        if type(kw.get('entrance_flat_through')) is not int or kw['entrance_flat_through']<j:
            raise ValueError('Original mixed entrance compatibility through output order required')
        bound.damping_lower(self.inputs.source.waves.q,kw.get('damping_floor'))
        bound.mixed.first.harmonic(kw.get('harmonic',1))
        linear=self.background_phase.linear_table(R,psi_table)
        rate=[self.phase('profile_rate',k) for k in range(R+1)]
        out=super().build(j,linear,log_table,rate,cutoff,**kw)
        out.update(background_phase_request=request(j),derived_rate_jets=rate,
                   phase_provider_inequalities_certified=False)
        return out

def record():
    return dict(status='current-full-background-phase-rate-and-linear-binding',
                raw_R_to_common_U_to_R_transfers_retained=True,
                forward_factorial_height='486+2*N^2+N*B; fixed before all jets',
                normal_matrix_moving_normal_and_varying_profile_rate_derived=True,
                Eik_bR_Bz_derived_from_full_background=True,
                complete_output_j_background_order='j+4; Borel Euler cutoff j+5',
                actual_native_background_or_source_inputs_certified=False,
                full_cone_pulse_residual_or_blowup=False)
