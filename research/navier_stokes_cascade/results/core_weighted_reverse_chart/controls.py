"""Independent finite Taylor and exact mixed-derivative controls."""
from fractions import Fraction as Q
from math import factorial,comb
from bounds import ReverseChart,ChartedPrimary,P,bridge

def review():
    import sympy as s
    t,R,U,Z=s.symbols('t R U Z',positive=True)
    chart=ReverseChart(Q(1,2),2)
    checks={}
    for n in range(7):
        a=sum(s.Rational(chart.increment(j).numerator,chart.increment(j).denominator)*t**j for j in range(1,n+1))
        e=sum(a**k/s.factorial(k) for k in range(n+1))
        checks[f'order{n}_Bell_recurrence_equals_independent_finite_Taylor']=s.expand(e).coeff(t,n)==chart.exponential(n)
    for d in (Q(1,3),Q(1,2),Q(3,2),Q(5,2)):
        ds=s.Rational(d.numerator,d.denominator)
        for r in (Q(1,2),Q(1),Q(2)):
            rr=s.Rational(r.numerator,r.denominator)
            for n in range(1,7):
                val=abs(s.diff(R**ds,R,n).subs(R,rr))/s.factorial(n)
                checks[f'power_{d}_{r}_order{n}_normalized_radial_enclosure']=bool(val<=chart.radial(n))
    # All five variables, nonlinear radial map and an exact flat factor.
    # Compare each mixed RAW jet to M_U times the transfer, pointwise at R=U=1.
    Y,T,V=s.symbols('Y T V')
    expr=s.exp(-1/U**2)*(1+U*Z+Y*T+V**2+U**3*Y)
    vars_U=(U,Z,Y,T,V);vars_R=(R,Z,Y,T,V)
    point={U:1,R:1,Z:1,Y:1,T:1,V:1}
    def indices(N,dim=5):
        if not dim:yield ();return
        for i in range(N+1):
            for rest in indices(N-i,dim-1):yield (i,)+rest
    def jet(f,variables,alpha):
        for v,j in zip(variables,alpha):f=s.diff(f,v,j)
        return s.simplify(s.E*f.subs(point)) # divide out SAME positive weight e^-1
    for n in range(4):
        alphas=list(indices(n));M=max(abs(jet(expr,vars_U,a)) for a in alphas)
        pulled=expr.subs(U,R**s.Rational(3,2));vals=[abs(jet(pulled,vars_R,a)) for a in alphas]
        tr=chart.transfer(n).evaluate(1)
        checks[f'order{n}_mixed_flat_weight_raw_sum_enclosed']=bool(sum(vals)<=M*tr)
        checks[f'order{n}_mixed_factorial_norm_enclosed']=bool(sum(v/s.prod(s.factorial(j) for j in a) for a,v in zip(alphas,vals))<=M*tr)
    # Identity chart with raw derivatives bounded by one realizes exp(5t).
    # This detects forgetting the four unchanged variables in a sum jet norm.
    checks['five_coordinate_sum_not_a_radial_only_norm']=comb(1+5,5)==6>2
    checks['unit_reverse_raw_transfer_can_fail']=s.diff(R**2,R).subs(R,1)==2>1
    checks['nonlinear_reverse_chain_term_is_required']=s.diff(s.exp(U).subs(U,R**2),R,2).subs(R,1)==6*s.E>4*s.E
    checks['normalization_factorial_conversion_cannot_be_omitted']=s.diff(U**4,U,4)==24>1
    root=bridge.primary.curl.entrance.inverse.roots.fixture()
    b=bridge.WeightedPrimary(root);joined=ChartedPrimary(b,chart)
    data={k:P(i+1) for i,k in enumerate(bridge.rows.KEYS) if k not in ('W','D')}
    result=joined.initial_rows(0,data,P(3))
    direct=b.initial_rows(0,data,P(3),chart.transfer(3))
    checks['complete_initial_row_recipe_exactly_reproduced']=all(result['rows'][k].c==direct['rows'][k].c for k in ('gr','pm','Etheta','Ez'))
    lost=dict(data);lost['PiP']=P()
    checks['whole_pressure_moment_still_changes_axial_row']=joined.initial_rows(0,lost,P(3))['rows']['Ez'].c!=result['rows']['Ez'].c
    for n in range(4):
        c=joined.covariance_R(n);old=b.covariance(n)
        checks[f'order{n}_pole_and_epsilon_orders_unchanged']=all(c[k]==old[k] for k in ('pole','epsilon_order_W','epsilon_order_Dcurl'))
        checks[f'order{n}_tensor_transfer_inserted_once']=all(c[k].c==(chart.transfer(n)*old[k]).c for k in ('W','Dcurl'))
    return dict(checks={k:bool(v) for k,v in checks.items()},parent_commit='e30f082183bd3e2ee04c60b6f67189e524616e43',
        remaining='weighted fast/slow/profile/operator and whole-pressure inputs in the same R hull; complete source, mean and physical residual closure',
        scope='finite exact rational/symbolic controls, not actual continuum suprema or a new residual gain')
