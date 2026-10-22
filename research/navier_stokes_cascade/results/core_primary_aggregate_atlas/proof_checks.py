"""Exact atlas, edge, common-normalization and primary-aggregate controls."""
from fractions import Fraction as Q
from math import ceil,factorial
from pathlib import Path
import importlib.util

HERE=Path(__file__).resolve().parent


def load(name,path):
    spec=importlib.util.spec_from_file_location(name,path)
    module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module)
    return module


seed=load("aggregate_initial_seed",HERE.parent/"core_initial_scalar_seed_budget/proof_checks.py")
atlas=load("aggregate_actual_atlas",HERE.parent/"auxiliary_geometry_audit/geometry.py")
P=seed.P;order=seed.order;exact=seed.coupled.cycle.exact


def edge_constant(b,width):
    b,width=exact(b),exact(width)
    if b<0 or not 0<width<=1:raise ValueError("Nonnegative edge loss and fixed width in (0,1] required")
    n=ceil(b/2);alpha=width*width/2
    return max(Q(1),Q(factorial(n))/alpha**n)


def rescale(poly,ratio):
    ratio=exact(ratio)
    if ratio<=0:raise ValueError("Positive slow-scale comparison required")
    if not isinstance(poly,P):raise TypeError("Positive scale polynomial required")
    return P({power:value*ratio**power for power,value in poly.c.items()})


def aggregate(n,piA,piB,bA,bB,width,arbitrary=False):
    order(n)
    if type(arbitrary) is not bool:raise TypeError("Reference convention must be explicit")
    a=8*188**n*edge_constant(bA,width)*rescale(piA,Q(25,16))
    b=8*188**n*edge_constant(bB,width)*rescale(piB,Q(25,16))
    ratio=2 if arbitrary else 1
    return dict(a=a,b=b,Omega=Q(2**n,2)*(a*a+ratio*b*b),
                carrier_ratio=ratio,time_floor=1296 if arbitrary else 6)


def atlas_review():
    parts={"covering_gap":atlas.covering_gap_review(),"coloring":atlas.coloring_review(),
           "separation":atlas.separation_review(),"covering":atlas.covering_review(),"counting":atlas.counting_review()}
    checks={name+":"+key:value for name,row in parts.items() for key,value in row["checks"].items()}
    for delta in range(4):
        J=atlas.matrix_power(delta);norm=atlas.norm_infinity(J)
        checks[f"delta{delta}_symmetric_covering_has_same_column_norm"]=J[0][1]==J[1][0] and max(sum(J[i][j] for i in range(2)) for j in range(2))==norm<=188
        degree=14**delta
        checks[f"delta{delta}_normalized_deck_average_has_no_degree_factor"]=sum((Q(1,degree) for _ in range(degree)),Q(0))==1
        if delta:checks[f"delta{delta}_unnormalized_copy_sum_is_different"]=degree!=1
        # Expand the two constant-coefficient pullback derivative operators.
        for p,q in ((0,0),(1,0),(1,1),(2,1),(0,4)):
            coefficients={(0,0):Q(1)}
            for column,count in ((0,p),(1,q)):
                for _ in range(count):
                    new={}
                    for (a,b),v in coefficients.items():
                        new[a+1,b]=new.get((a+1,b),0)+v*J[0][column]
                        new[a,b+1]=new.get((a,b+1),0)+v*J[1][column]
                    coefficients=new
            checks[f"delta{delta}_p{p}_q{q}_raw_jet_pullback_bound"]=sum(abs(v) for v in coefficients.values())<=188**(p+q)
    return dict(checks=checks,scope="replayed exact actual atlas certificates and finite pullback expansions; all-order chain rule proved in note")


def edge_review():
    checks={};rows=[]
    for b in (Q(0),Q(1),Q(3,2),Q(2),Q(5),Q(8)):
        n=ceil(b/2)
        for width in (Q(1),Q(1,8),Q(1,1000)):
            alpha=width*width/2;E=edge_constant(b,width)
            checks[f"b{b}_t{width}_edge_series_coefficient_controls_x_power"]=E>=Q(factorial(n))/alpha**n and n>=b/2
            checks[f"b{b}_t{width}_interior_delta_one_is_covered"]=E>=1
            for delta in (Q(1,2),Q(1,7)):
                x=delta**-2
                checks[f"b{b}_t{width}_d{delta}_exponential_series_term"]=E*(alpha*x)**n/factorial(n)>=x**n
        rows.append(dict(edge_loss=str(b),series_order=n,width_eighth_bound=str(edge_constant(b,Q(1,8)))))
    checks["width_dependence_cannot_be_declared_uniform"]=edge_constant(8,Q(1,1000))>edge_constant(8,Q(1,8))
    return dict(checks=checks,rows=rows,scope="rational exponential-series sufficient constants, not numerical evaluation of the selected width or profile")


def normalization_review():
    checks={}
    for h in (Q(1,200),Q(1,2**40)):
        kappa=Q(1,100000);A=Q(1,2)+h
        pA=A-h/2;pB=pA+h*kappa
        checks[f"h{h}_coefficient_conversion_exponents_below_three_quarters"]=0<pA<pB<Q(3,4)
        checks[f"h{h}_radial_U_multiplier_bound"]=1+h<=Q(101,100)<Q(5,4)
        checks[f"h{h}_arbitrary_reference_carrier_ratio_below_two"]=h<Q(1,4)
    checks["actual_band_window_slow_scale_ratio"]=(1+Q(4,16))**2==Q(25,16)
    checks["raw_slow_chart_multipliers_below_covering_constant"]=max(32,4,16)<188
    checks["amplitude_normalization_bound_eight"]=16**3==8**4
    checks["radial_U_normalization_bound_thirty_two"]=16**5==32**4
    checks["coarsest_index_is_common_minimum"]=Q(1,2)-Q(2,16)>0
    # At h=1/200, ell_i=40000, the native carrier is exactly 2^100.
    h=Q(1,200);ell=40000;ref=ell+1
    checks["finer_reference_literal_carrier_bound_fails"]=(h*ell/2).denominator==1 and h*ref>h*ell
    checks["native_rounding_counterexample_is_exact"]=2**int(h*ell/2)==2**100
    # c_common = Tg^-Delta c_ref; use rational Tg=5 to test the scale algebra.
    for delta in range(4):
        S=Q(37);Tg=Q(5);c_ref=Q(3,4)/S;c_common=c_ref/Tg**delta
        checks[f"delta{delta}_arbitrary_time_floor_constant_suffices"]=1/c_common<Tg**4*S
        if delta:checks[f"delta{delta}_native_time_floor_cannot_be_reused"]=1/c_common>Tg*S
    checks["rational_time_floor_majorants_match_two_conventions"]=6**4==1296
    return dict(checks=checks,scope="exact scaling exponents and representative floor arithmetic; no selected band or carrier evaluation")


def aggregate_fixture(arbitrary=False,Svalue=None,trace=None):
    source=seed.fixture_budget()
    def omega(n):
        if trace is not None:trace.add(n)
        row=aggregate(n,P({0:n+1,1:1}),P({0:2,1:n+1}),2,4,Q(1,8),arbitrary)
        return row["Omega"] if Svalue is None else row["Omega"].evaluate(Svalue)
    base=seed.fixture_budget(Svalue)
    base.omega=omega;base.tg=1296 if arbitrary else 6
    return base


def chain_review():
    checks={};rows=[]
    for arbitrary in (False,True):
        for n in (0,2):
            abstract=aggregate_fixture(arbitrary);out=abstract.output(n);moment=abstract.axial_moment(n,Q(5))
            for S in (1,4,16):
                concrete=aggregate_fixture(arbitrary,S)
                checks[f"arbitrary{arbitrary}_n{n}_S{S}_aggregate_to_seed_substitution"]=all(poly.evaluate(S)==concrete.output(n)[key] for key,poly in out.items())
                checks[f"arbitrary{arbitrary}_n{n}_S{S}_whole_pressure_moment_substitution"]=moment.evaluate(S)==concrete.axial_moment(n,Q(5))
            trace=set();aggregate_fixture(arbitrary,trace=trace).output(n)
            checks[f"arbitrary{arbitrary}_n{n}_scalar_jet_ledger_unchanged"]=max(trace)==n+10
            trace=set();aggregate_fixture(arbitrary,trace=trace).axial_moment(n,Q(5))
            checks[f"arbitrary{arbitrary}_n{n}_axial_moment_jet_ledger_unchanged"]=max(trace)==n+11
            rows.append(dict(arbitrary_reference=arbitrary,output_order=n,
                             ThetaP_S4_numerator_bits=out["P"].evaluate(4).numerator.bit_length(),
                             ThetaZ_S4_numerator_bits=out["Jz"].evaluate(4).numerator.bit_length()))
    for n in (0,1,3):
        a=P({0:2,1:1});b=P({1:3})
        c=aggregate(n,a,b,0,0,1);other=aggregate(n,a,b,0,0,1,True)
        for S in (1,4):
            nativeA=a.evaluate(Q(25,16)*S);nativeB=b.evaluate(Q(25,16)*S)
            checks[f"n{n}_S{S}_no_label_or_deck_factor_in_aggregate"]=c["a"].evaluate(S)==8*188**n*nativeA and c["b"].evaluate(S)==8*188**n*nativeB
            checks[f"n{n}_S{S}_arbitrary_carrier_ratio_retained"]=other["Omega"].evaluate(S)==Q(2**n,2)*(c["a"].evaluate(S)**2+2*c["b"].evaluate(S)**2)
    return dict(checks=checks,rows=rows,scope="positive polynomial substitution under retained native per-label bounds; selected primary constants remain unevaluated")


def review():
    parts={"actual_atlas":atlas_review(),"edge_unweighting":edge_review(),"matched_normalization":normalization_review(),
           "aggregate_seed_chain":chain_review(),"pressure_tails":seed.means.pressure_tail_review(),
           "complete_covariance":seed.means.covariance_review()}
    return dict(checks={p+":"+k:v for p,row in parts.items() for k,v in row["checks"].items()},parts=parts)
