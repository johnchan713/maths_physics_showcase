"""Closed slow-coordinate geometry factors; actual tuple/schedule inputs remain."""
from pathlib import Path
from fractions import Fraction as Q
from math import factorial
from functools import lru_cache
import importlib.util
p=Path(__file__).resolve().parent.parent/'core_current_head_q118_background_phase_binding/q118.py'
s=importlib.util.spec_from_file_location('q119_phase',p);phase=importlib.util.module_from_spec(s);s.loader.exec_module(phase)
base=phase.background

def order(m):
    if type(m) is not int or m<0:raise ValueError('Nonnegative integer derivative order required')

def chart_recipe(m,a,b):
    order(m);a,b=Q(a),Q(b)
    if not 0<a<=b:raise ValueError('Original positive closed normalized R hull required')
    tau=min(Q(1),a)/(1024*(1+b)**2)
    return dict(order=m,tau=tau,factorial_argument=m,geometric_base=1/tau,
                prefactors=dict(F=48/a,G=Q(24),btilde=Q(12)),
                whole_tuple_inputs_required=True,actual_profile_inputs_certified=False)

def materialize(recipe,maximum=40):
    order(maximum);m=recipe['factorial_argument'];order(m)
    if m>maximum:raise ValueError('Keep large derivative factorials in the analytical recipe')
    common=factorial(m)*recipe['geometric_base']**m
    return dict(kind='closed-real-slow-chart-majorant',tau=recipe['tau'],
                **{k:c*common for k,c in recipe['prefactors'].items()})

class ClosedBackground(base.BackgroundBudget):
    """Same complete tuple/head/tail inputs; preserve old low-order value budgets."""
    def chart_recipe(self,m):return chart_recipe(m,self.Rmin,self.Rmax)
    def leading_bound(self,m):
        order(m)
        if m<=2:return super().leading_bound(m)
        if self.trace is not None:self.trace.add(('leading',m))
        chart=materialize(self.chart_recipe(m),maximum=m)
        return (chart['F']+chart['G'])*base.nonnegative(self.leading(m))
    @lru_cache(None)
    def bound(self,m):
        order(m)
        if m<=2:return super().bound(m)
        if self.trace is not None:self.trace.add(('leading',m));self.trace.add(('cutoff',m+1))
        L=base.nonnegative(self.leading(m));K=base.nonnegative(self.cutoff(m+1))
        head={}
        for n in range(1,max(2,m)):
            if self.trace is not None:self.trace.add(('profile',n,m))
            head[n]=base.increment_constant(self.h,n,m,self.profiles(n,m),K)
        out=base.finite_head(m,head);chart=materialize(self.chart_recipe(m),maximum=m)
        total=L+out['correction'];fg=chart['F']+chart['G']
        return dict(F=chart['F']*total,G=chart['G']*total,btilde=chart['btilde']*total,
                    leading_FG=fg*L,full_FG=fg*total,error_FG=4*fg*out['correction'],
                    head=out,chart=chart,cutoff_input_order=m+1)

def record():
    return dict(status='closed-real-slow-chart-background-majorants',
                tau='min(1,Rmin)/(1024*(1+Rmax)^2)',
                prefactors='F:48/Rmin; G:24; btilde:12',
                raw_through_m_factor='m!*tau^-m',
                original_low_orders_zero_one_two_preserved=True,
                original_full_tuple_Borel_head_tail_and_schedule_retained=True,
                analytic_coordinates_only_not_analytic_profiles=True,
                actual_native_profile_or_source_inequalities_certified=False,
                full_cone_pulse_residual_or_blowup=False)
