"""Same-weight common-atlas covariance substitution into the original mean rows."""
from fractions import Fraction as Q
from functools import lru_cache
from math import factorial
from pathlib import Path
import importlib.util
HERE=Path(__file__).resolve().parent
def load(name,path):
    sp=importlib.util.spec_from_file_location(name,HERE.parent/path);m=importlib.util.module_from_spec(sp);sp.loader.exec_module(m);return m
primary=load('weighted_join_primary','core_primary_target_cutoff_budget/bounds.py')
rows=load('weighted_join_rows','core_initial_haar_row_budget/budgets.py')
P,order=primary.P,primary.order

def rescale(p):return P({k:v*Q(25,16)**k for k,v in p.c.items()})

class WeightedPrimary:
    """Use the complete original U-chart root recipe, before edge unweighting.

    Actual weight equality, matched coarsest normalization and closed-support
    separation are retained geometry hypotheses, not new field choices.
    """
    def __init__(self,root):
        if not isinstance(root,primary.curl.entrance.inverse.roots.RootBudget):raise TypeError('Complete same-family U-chart primary root recipe required')
        self.root=root

    @lru_cache(None)
    def coefficients(self,n):
        order(n)
        a=16*188**n*factorial(n)*rescale(sum((self.root.amplitude(j) for j in range(n+1)),P()))
        b=16*188**n*factorial(n)*rescale(sum((self.root.complete_curl(j) for j in range(n+1)),P()))
        return a,b

    def covariance(self,n):
        order(n);a,b=self.coefficients(n)
        return dict(W=Q(2**n,2)*(a*a+b*b),Dcurl=Q(2**n,2)*b*b,
            pole=12*n+18,epsilon_order_W=Q(1),epsilon_order_Dcurl=2-Q(2,100000))

    def initial_rows(self,n,data,higher_stress,row_transfer):
        """Replace W,D only, preserving every mean/pressure/base input.

        row_transfer bounds the SAME-WEIGHT U-to-R tensor jet transfer through
        n+3. It is required, not inferred to equal one. higher_stress bounds
        the original weighted Sigma-Sigma0 in the ROW chart at order3,
        with pole <=12*(n+3)+18; its existing Borel/profile validity remains.
        """
        order(n)
        if set(data)!=set(rows.KEYS)-{'W','D'}:raise ValueError('All remaining base, mean, axial and full pressure inputs required')
        if not isinstance(higher_stress,P) or not isinstance(row_transfer,P) or not all(isinstance(v,P) for v in data.values()):raise TypeError('Positive polynomial lower-input and same-weight chart-transfer recipes required')
        c=self.covariance(n+3)
        convert=lambda v:rows.P(dict(v.c))
        inputs={k:convert(v) for k,v in data.items()};inputs.update(W=convert(row_transfer*c['W']),D=convert(row_transfer*c['Dcurl']+higher_stress))
        return dict(rows=rows.row_budget(n,inputs),covariance=c,
            input_order=n+3,input_pole=c['pole'],remaining_inputs=tuple(sorted(data)),
            complete_mean_or_physical_residual_certified=False)

def from_original(p,selector_inputs,ell,**kw):
    if 'chart' in kw:raise ValueError('This bridge explicitly constructs the U-chart recipe once')
    bundle=primary.from_original(p,selector_inputs,ell,chart='U',**kw)
    return bundle,WeightedPrimary(bundle[3])

def review():
    import sympy as s
    q,Qb,r,z,D=s.symbols('q Q r z D',positive=True)
    X=(r/s.sqrt(Qb))**2/(2*q/Qb);eta=(z/Qb**D)/(q/Qb)**D
    t=s.symbols('t',real=True);A=s.Matrix([1+t,2-t,1+t*t]);B=s.Matrix([t,1+t,2-t])
    theta=s.symbols('theta',real=True);wave=A*s.cos(theta)-B*s.sin(theta)
    tensor=s.simplify(s.integrate(wave*wave.T,(theta,0,2*s.pi))/(2*s.pi))
    checks=dict(common_X_weight_coordinate_invariant=s.simplify(X-r*r/(2*q))==0,
        common_eta_weight_coordinate_invariant=s.powsimp(eta*q**D/z,force=True)==1,
        all_six_real_covariances_retained=all(s.expand(v)==0 for v in tensor-(A*A.T+B*B.T)/2),
        all_six_curl_square_entries_nonzero=all((B*B.T)[i,j]!=0 for i in range(3) for j in range(i,3)),
        mixed_real_quadrature_averages_zero=s.integrate(s.sin(theta)*s.cos(theta),(theta,0,2*s.pi))==0,
        generic_complex_sine_cross_must_remain=s.integrate(s.sin(theta)**2,(theta,0,2*s.pi))==s.pi,
        covariance_pole_covers_higher_stress=all(3*n+3<=12*n+18 for n in range(8)),
        full_curl_square_retains_fine_order=2-Q(2,100000)>1,
        normalized_cover_has_no_copy_count=all(sum(Q(1,14**d) for _ in range(14**d))==1 for d in range(4)))
    rr=s.symbols('R',positive=True)
    checks['unit_reverse_chart_transfer_is_not_valid_in_general']=s.diff(rr**2,rr).subs(rr,1)==2>1
    for n in range(5):
        checks[f'order{n}_Leibniz_binomial_sum']=sum(s.binomial(n,i) for i in range(n+1))==2**n
        checks[f'order{n}_square_weight_pole_is_twelve_n_plus_eighteen']=2*(6*n+9)==12*n+18
    root=primary.curl.entrance.inverse.roots.fixture();budget=WeightedPrimary(root)
    class Adapter:
        # Exact coefficient copy between independently loaded diagnostic classes.
        def amplitude(self,j):return rows.P(dict(root.amplitude(j).c))
        def complete_curl(self,j):return rows.P(dict(root.complete_curl(j).c))
    for n in range(3):
        old=rows.weighted_covariance(n,Adapter());new=budget.covariance(n)
        checks[f'order{n}_original_weighted_covariance_recipe_reproduced']=all(old[k].c==new[k].c for k in ('W','Dcurl'))
    data={k:P(i+1) for i,k in enumerate(rows.KEYS) if k not in ('W','D')}
    result=budget.initial_rows(0,data,P(3),P(1));cc=budget.covariance(3)
    direct=rows.row_budget(0,{k:rows.P(dict(v.c)) for k,v in dict(data,W=cc['W'],D=cc['Dcurl']+P(3)).items()})
    checks['complete_mean_row_substitution_matches_original_formula']=all(result['rows'][k].c==direct[k].c for k in ('gr','pm','Etheta','Ez'))
    changed=dict(data);changed['PiP']=P(0)
    checks['whole_pressure_moment_input_cannot_be_deleted']=budget.initial_rows(0,changed,P(3),P(1))['rows']['Ez'].c!=result['rows']['Ez'].c
    checks['nontrivial_row_chart_transfer_changes_the_row_bound']=budget.initial_rows(0,data,P(3),P(2))['rows']['Etheta'].c!=result['rows']['Etheta'].c
    return dict(checks={k:bool(v) for k,v in checks.items()},scope='exact manufactured tensors, atlas identities and lower-input substitutions, not actual continuum suprema',
        parent_commit='b919fe2154652bfc4df6586ef22918a61c147c07',
        remaining='derive same-weight U-to-R row transfer (not unit by default); weighted fast/slow/full pressure and higher-stress inputs, same hull/value gates and full residual remain conditional')
